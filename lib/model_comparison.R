library(optparse)
library(randomForest)
#library(ranger)
library(tidyr)
library(dplyr)
library(sf)
library(fs)
library(stringr)
library(rockchalk)
library(terra)
library(paws.storage)

get_height_column_names <- function(in_data){
  return(
    names(in_data)[grep('^RH_[0-9]{2}$', names(in_data))]
  )
}

rename_height_columns_to_match_pretrained_models <- function(in_data){
  return(
    in_data |>
      rename_with(~gsub('\\brh([0-9]{2})\\b', 'RH_\\1', .x), matches='^rh[0-9]{2}$') |>
      rename(RH_98=h_canopy)
  )
}

offset_RH_columns <- function(all_train_data, offset){
  RH_columns <- get_height_column_names(all_train_data)

  return(
    all_train_data |>
      mutate(across(all_of(RH_columns), ~. + offset))
  )
}

set_model_id_for_AGB_prediction <- function(in_data, offset){
  # TODO: uncomment correct model ids once tested against old results
  return(
    in_data |>
      mutate(model_id = case_when(
        segment_landcover %in% c(111, 113, 121, 123) ~ "m3", # needle leaf
        segment_landcover %in% c(112, 114, 122, 124) ~ "m1", # broad leaf
        TRUE ~ "m8"
      )
      )
  )
}

randomize_AGB_model <- function(model){
  # modify coeffients through sampling variance covariance matrix
  model_coeffs <- mvrnorm(n=50, mu=model$coefficients, Sigma=vcov(model))
  model$coefficients <- model_coeffs[1,]

  return(model)
}

stratified_random_sample <- function(df, max_n) {
  df |> group_by(segment_landcover) |> summarise(c=n()) |> filter(c < 1000)
  #group_size = max_n / 

}

GEDI2AT08AGB<-function(rds_models, df, randomize=FALSE, max_n=10000, sample=TRUE){
  if (sample && nrow(df) > max_n)
    df <- reduce_sample_size(df, max_n)

  df$AGB <- NA
  df$SE <- NA

  ids<-unique(df$model_id)
  n_models <- length(ids)

  for (i in ids){
    model_i <- rds_models[[i]]

    # Modify coeffients through sampling variance covariance matrix
    if(randomize)
      model_i <- randomize_AGB_model(model_i)

    # Predict AGB and SE
    df$AGB[df$model_id==i] <- predict(model_i, newdata=df[df$model_id==i,])
    df$SE[df$model_id==i] <- summary(model_i)$sigma^2

    df$AGB[df$AGB < 0] <- 0.0

    # Calculate Correction Factor C
    C <- mean(model_i$fitted.values^2)/mean(model_i$model$`sqrt(AGBD)`^2)

    # Bias correction in case there is a systematic over or under estimation in the model
    df$AGB[df$model_id==i] <- C*(df$AGB[df$model_id==i]^2)
  }

  # Apply slopemask, validmask and landcover masks
  bad_lc <- c(0, 60, 80, 200, 50, 70)
  df$AGB[df$slopemask == 0 |
                df$ValidMask == 0 |
                df$segment_landcover %in% bad_lc] <- 0.0
  return(df)
}

DOY_and_solar_filter <- function(tile_data, start_DOY, end_DOY, solar_elevation){
  filter <- which((tile_data$doy >= start_DOY) &
                    (tile_data$doy <= end_DOY) &
                    (tile_data$solar_elevation < solar_elevation)
                  )
  return(filter)
}

late_season_filter <- function(tile_data, minDOY, maxDOY,
                               default_maxDOY, min_icesat2_samples, max_sol_el){
  n_late <- 0
  for(late_months in 0:3) {
    if(n_late < min_icesat2_samples) {

      default_maxDOY <- default_maxDOY + 30 * late_months

      if(default_maxDOY < maxDOY){
        filter <- DOY_and_solar_filter(minDOY, default_maxDOY, max_sol_el)
        n_late <- length(filter)
      }
    }
  }
  return(list(filter=filter, default_maxDOY=default_maxDOY))
}

early_and_late_season_filter <- function(tile_data, minDOY,
                                         default_minDOY, default_maxDOY,
                                         min_icesat2_samples, max_sol_el){
  n_early <- 0
  for(early_months in 0:3){
    if(n_early < min_icesat2_samples){
      default_minDOY <- default_minDOY - 30 * early_months

      if(default_minDOY > minDOY){
        filter <- DOY_and_solar_filter(default_minDOY, default_maxDOY, max_sol_el)
        n_early <- length(filter)
      }
    }

  }
  return(list(filter=filter, default_minDOY=default_minDOY))
}

expand_training_around_season <- function(tile_data, minDOY, maxDOY,
                                          default_minDOY, default_maxDOY,
                                          max_sol_el, min_icesat2_samples){
  filter <- DOY_and_solar_filter(tile_data, minDOY, maxDOY, max_sol_el)
  if(length(filter) >= min_icesat2_samples){
    cat('Found nough data with max_solar_elevation:', max_sol_el, '\n')
    return(tile_data[filter,])
  }
  # next try expanding 1 month later in growing season, iteratively, up to 3 months
  filter <- late_season_filter(
    tile_data, minDOY, maxDOY, default_maxDOY, min_icesat2_samples, max_sol_el
  )
  if(length(filter$filter) >= min_icesat2_samples){
    cat('Found enough data when expanding into late season DOY:', filter$default_maxDOY, '\n')
    return(tile_data[filter$filter,])
  }

  # next try expanding 1 month earlier in growing season, iteratively, up to 3 months
  # Note that the upper window might be later in the growing season from the previous call
  current_default_maxDOY <- filter$default_maxDOY
  filter <- early_and_late_season_filter(
    tile_data, minDOY, default_minDOY, current_default_maxDOY, min_icesat2_samples, max_sol_el
  )
  if(length(filter$filter) >= min_icesat2_samples){
    cat(
      'Found enough data when expanding into early and late season DOY:[',
      filter$default_minDOY, ' ', default_minDOY,  ']\n'
    )
    return(tile_data[filter$filter,])
  }

  print('Search into late and early season did not return enough data')
  print('applying basic filter')
  tile_data <- tile_data[DOY_and_solar_filter(tile_data, default_minDOY, default_maxDOY, 0),]
  return(tile_data)
}

reduce_sample_size <- function(df, sample_size){
  return(df[sample(row.names(df), sample_size, replace=FALSE),])
}

remove_stale_columns <- function(df, column_names) {
  columns_to_remove <- intersect(names(df), column_names)
  df <- df[, !names(df) %in% columns_to_remove, drop=FALSE]

  return(df)
}

sample_broad_data_within_latitude <- function(broad_data, lat, threshold, samples_needed){
  broad_within_lat <- which(broad_data$lat > (lat-threshold) & broad_data$lat < (lat+threshold))
  broad_data <- broad_data[broad_within_lat,]
  return(broad_data[sample(row.names(broad_data), samples_needed, replace=TRUE), ])
}

expand_training_with_broad_data <- function(broad_data, tile_data, samples_needed){
  broad_data <- sample_broad_data_within_latitude(broad_data, min(tile_data$lat), 5, samples_needed)
  # TODO this check may no longer be needed, ask Paul
  if (!setequal(names(broad_data), names(tile_data))) {
    only_in_broad <- setdiff(names(broad_data), names(tile_data))
    only_in_local <- setdiff(names(tile_data), names(broad_data))
    diff <- union(only_in_broad, only_in_local)
    print('Warning!')
    print('Boreal wide training data and local data have non matching columns!')
    print('Will drop the following non matching columns and continue:')
    print(diff)
    common <- intersect(names(broad_data), names(tile_data))
    tile_data <- tile_data[common]
    broad_data <- broad_data[common]
    print(common)
  }

  return(rbind(tile_data, broad_data))
}

remove_height_outliers <- function(all_train_data){
  # remove height outliers based on more than 3SD from the landcover mean
  return(
  all_train_data |>
    group_by(segment_landcover) |>
    summarise(thresh=mean(h_canopy, na.rm=T) + 3 * sd(h_canopy, na.rm=T)) |>
    right_join(all_train_data, by='segment_landcover') |>
    filter(h_canopy <= thresh)
  )
}

set_output_file_names <- function(predict_var, tile_num, year){
  key <- if (predict_var == 'AGB') 'agb' else 'ht'
  out_fn_stem = paste(
    paste0('output/boreal_', key, '_', year), format(Sys.time(),"%Y%m%d%s"), str_pad(tile_num, 7, pad = "0"),
    sep="_"
  )

  fn_suffix <- c('.tif', '_summary.csv', '_train_data.csv', '_stats.Rds',
  '_model.Rds', '_comparison.csv', '_cv.csv', '_lvis_r2_rmse.csv')
  names <- c('map', 'summary', 'train', 'stats', 'model', 'comparison', 'cv', 'r2_rmse_lvis')

  output_file_names <- paste0(out_fn_stem, fn_suffix)
  names(output_file_names) <- names

  return(output_file_names)
}

write_single_model_summary <- function(model, df, target, out_fns){
  target <- if(target == 'AGB') df$AGB else df$RH_98
  local_model <- lm(model$predicted ~ target)
  saveRDS(model, file=out_fns['model'])

  rsq <- max(model$rsq, na.rm=T)
  cat('rsq_model: ', rsq, '\n')

  rsq_local <- summary(local_model)$r.squared
  cat('rsq_local: ', rsq_local, '\n')

  na_data <- which(is.na(local_model$predicted==TRUE))

  if(length(na_data) == 0)
    rmse_local <- sqrt(mean(local_model$residuals^2))

  cat('rmse_local: ', rmse_local, '\n')

  imp_vars <- model$importance
  out_accuracy <- list(RSQ=rsq_local, RMSE=rmse_local, importance=imp_vars)
  saveRDS(out_accuracy, file=out_fns['stats'])
}

read_and_filter_training_data <- function(ice2_30_atl08_path, expand_training, min_icesat2_samples, minDOY, maxDOY, max_sol_el){
  default_maxDOY <- 273
  default_minDOY <- 121
  tile_data <- read.csv(ice2_30_atl08_path)

  night_time_in_season <- DOY_and_solar_filter(tile_data, default_minDOY, default_maxDOY, 0)
  cat('train data size before any filtering:', nrow(tile_data), '\n')
  cat('length(night_time_in_season):', length(night_time_in_season), 'expand_training:', expand_training, ' min_n:', min_icesat2_samples, '\n')

  if (length(night_time_in_season) < min_icesat2_samples && expand_training){
    cat('running expansion:', length(night_time_in_season), '<', min_icesat2_samples, '\n')
    tile_data <- expand_training_around_season(tile_data, minDOY, maxDOY, default_minDOY, default_maxDOY, max_sol_el, min_icesat2_samples)
  }
  else{
    print('night time filter only')
    tile_data <- tile_data[night_time_in_season, ]
  }
  cat('training data size after filtering:', nrow(tile_data), '\n')
  tile_data <- remove_stale_columns(tile_data, c("binsize", "num_bins"))
  return(tile_data)
}

augment_training_data_with_broad_data <- function(tile_data, ice2_30_sample_path, local_train_perc, min_icesat2_samples){
  broad_data <- read.csv(ice2_30_sample_path)
  broad_data <- remove_stale_columns(broad_data, c("X__index_level_0__", "geometry"))

  # take proportion of broad data we want based on local_train_perc
  sample_local <- ceiling(nrow(tile_data) * local_train_perc / 100)
  cat('sample_local:', sample_local, '\n')

  if (sample_local < min_icesat2_samples){
    cat('reducing sample size to', sample_local, ' from ', nrow(tile_data), 'to complete with broad data \n')
    tile_data <- reduce_sample_size(tile_data, sample_local)
  }

  # sample from broad data to complete sample size
  # this will work if either there aren't enough local samples for n_min OR if there is forced broad sampling
  n_broad <- min_icesat2_samples - nrow(tile_data)
  if(n_broad > 1){
    tile_data <- expand_training_with_broad_data(broad_data, tile_data, n_broad)
    cat('training data size after augmenting with broad data:', nrow(tile_data), '\n')
  }
  return(tile_data)
}

reformat_training_data_for_AGB_modeling <- function(tile_data, offset){
  tile_data <- rename_height_columns_to_match_pretrained_models(tile_data)
  tile_data$h_canopy <- tile_data$RH_98
  tile_data <- offset_RH_columns(tile_data, offset)
  tile_data <- set_model_id_for_AGB_prediction(tile_data)
  return(tile_data)
}

prepare_training_data <- function(ice2_30_atl08_path, ice2_30_sample_path,
                                  expand_training, minDOY, maxDOY, max_sol_el,
                                  min_icesat2_samples, local_train_perc, offset, stack_vars){
  print('preparing training data ...')

  tile_data <- read_and_filter_training_data(
    ice2_30_atl08_path, expand_training,
    min_icesat2_samples, minDOY, maxDOY, max_sol_el
  )

  tile_data <- augment_training_data_with_broad_data(
    tile_data, ice2_30_sample_path, local_train_perc, min_icesat2_samples
  )
  print(stack_vars)
  needed_cols <- union(
    c('lat', 'lon', 'segment_landcover', 'h_canopy', 'rh25', 'rh50', 'rh60',
      'rh70', 'rh75', 'rh80', 'rh85', 'rh90', 'rh95'),
    stack_vars
  )

  tile_data <- tile_data |> select(all_of(needed_cols))

  tile_data <- reformat_training_data_for_AGB_modeling(tile_data, offset)

  tile_data <- remove_height_outliers(tile_data)
  cat('training data size after removing height outliers:', nrow(tile_data), '\n')

  tile_data <- tile_data |> filter(if_all(everything(), ~ !is.na(.x) & .x != -9999))
  cat('training data size after removing NAs:', nrow(tile_data), '\n')

  str(tile_data)
  cat('table for model training generated with ', nrow(tile_data), ' observations\n')

  if (nrow(tile_data) <= 1) {
    # TODO another option could be to drop SAR columns and continue
    stop('No traing data available, likley due to SAR being all -9999')
  }
  return(tile_data)
}

get_rds_models <- function(){
  rds_model_fns <- list.files(pattern='*.rds')
  rds_models <- lapply(rds_model_fns, readRDS)
  names(rds_models) <- paste0("m",1:length(rds_models))
  print(rds_models)
  return(rds_models)
}

fit_model <- function(model, model_config, train_df, pred_vars, predict_var){
  y_fit <- if (predict_var == 'Ht') train_df$h_canopy else train_df$AGB
  x_fit <- train_df[pred_vars]

  model_fit <- do.call(model, modifyList(model_config, list(y=y_fit, x=x_fit)))
  return(model_fit)
}

setup_kfold <- function(train_df, k) {
  train_df <- train_df[sample(nrow(train_df)), ]
  fold_size <- floor(nrow(train_df) / k)
  train_df$folds <- 0
  n <- 1

  for (fold in 1:k) {
    end <- min(n + fold_size - 1, nrow(train_df))
    train_df[n:end, 'folds'] <- fold
    n <- end + 1
  }
  return(train_df)
}

validate <- function(model, val_df, predict_var){
print(names(model))
  target_var <- if (predict_var == 'Ht') 'h_canopy' else 'AGB'
  predictions <- predict(model, val_df)
  val_df$residuals <- predictions - val_df[[target_var]]
  metrics <- val_df %>%
    group_by(segment_landcover) %>%
    summarise(
      RMSE = sqrt(mean(residuals^2, na.rm = TRUE)),
      R2 = 1 - (sum(residuals^2, na.rm = TRUE) / sum((.data[[target_var]] - mean(.data[[target_var]], na.rm = TRUE))^2, na.rm = TRUE))
    ) %>%
    ungroup()
  return(metrics)
}

download_rasters_from_s3 <- function(s3_paths){
  bucket <- 'ornl-cumulus-prod-protected'
  local_dir <- '/tmp'
  s3 <- s3()
  local_paths <- c()

  for (s3_path in s3_paths) {
    s3_parts <- strsplit(sub("s3://", "", s3_path), "/")[[1]]
    bucket <- s3_parts[1]
    key <- paste(s3_parts[-1], collapse = "/")
    local_path <- file.path(local_dir, basename(s3_path))
    cat("Bucket:", bucket, "Key:", key, "local_path:", local_path, "\n")
    response <- s3$get_object(Bucket = bucket, Key = key)
    writeBin(response$Body, local_path)
    cat("Downloaded:", s3_path, "->", local_path, "\n")
    local_paths <- c(local_paths, local_path)
  }

  return(local_paths)
}

create_lvis_mosaic <- function(stack, lvis_footprints_gpkg, year){
  lvis_footprints <- st_read(lvis_footprints_gpkg, layer=paste0('LVISF3_', year))
  # Transform stack extent to the CRS of lvis_footprints
  print('AOI Before transform:')
  print(ext(stack))
  raster_extent <- as.polygons(ext(stack), crs = crs(stack)) |> st_as_sf()
  print(raster_extent)
  raster_extent <- st_transform(raster_extent, st_crs(lvis_footprints))
  print('AOI after transform:')
  print(raster_extent)
  intersecting_footprints <- st_intersection(lvis_footprints, raster_extent)
  intersecting_footprints <- intersecting_footprints |>
    mutate(
      FILE = str_replace(FILE, "lvis_pt_cnt", "RH098_mean"),
      FILE = paste0("s3://ornl-cumulus-prod-protected/above/ABoVE_LVIS_VegetationStructure/data/", FILE)
    )
  lvis_files <- download_rasters_from_s3(intersecting_footprints$FILE)
  print(lvis_files)
  if (length(lvis_files) == 0) {
    stop("No LVIS rasters found in /tmp")
  }
  lvis_stack <- lapply(lvis_files, rast)
  lvis_mosaic <- do.call(merge, lvis_stack)

  # Ensure LVIS mosaic matches the CRS, resolution, and extent of the stack
  lvis_mosaic <- project(lvis_mosaic, stack)
  lvis_mosaic <- crop(lvis_mosaic, ext(stack))
  lvis_mosaic <- resample(lvis_mosaic, stack)

  output_mosaic <- "/tmp/lvis_mosaic.tif"
  writeRaster(lvis_mosaic, output_mosaic, overwrite = TRUE)
  cat("Mosaic saved to:", output_mosaic)

  return(lvis_mosaic)
}


predict_stack <- function(model, stack, majority_lc_classes){
  stack <- na.omit(stack)
  lc_classes_matrix <- cbind(majority_lc_classes, 1)  # Assign value 1 to valid classes
  landcover_mask <- classify(stack$esa_worldcover_v100_2020, lc_classes_matrix, others=NA)
  map <- predict(stack, model, na.rm=TRUE)
  map[!landcover_mask] <- NA
  # set slope and valid mask to zero
  # TODO maybe mask can skip over these pixels by default?
  map <- mask(map, stack$slopemask, maskvalues=0, updatevalue=0)
  map <- mask(map, stack$ValidMask, maskvalues=0, updatevalue=0)

  return(map)
}

subset_to_majority_lc_classes <- function(df){
  print('creating AGB traing data frame.')
  counts <- df |> group_by(segment_landcover) |> summarise(num_samples=n()) |> ungroup()
  n <- nrow(df)
  nunique <- n_distinct(df$segment_landcover)
  cat("n=",n, " nunique classes=", nunique, " inclusion pct=", n/nunique, "\n")
  df <- df |>
    group_by(segment_landcover) |>
    summarise(count=n()) |>
    right_join(df, by='segment_landcover') |>
    filter(count >= n/nunique)

  print(df |> group_by(segment_landcover) |> summarise(count=n()))
  return(df)
}

cross_validate <- function(df, pred_vars, pred_vars_nosar, folds, model, model_config, tile_num){
  df <- setup_kfold(df, 10)
  results <- data.frame()

  for (fold in 1:folds) {
    t1 <- Sys.time()
    cat('fold:', fold, '\n')

    train_df <- df[df['folds'] != fold, ]
    val_df <- df[df['folds'] == fold, ]

    print('fitting model with SAR')
    print(train_df |> group_by(segment_landcover) |> summarise(count=n()))
    
    model_sar <- fit_model(model, model_config, train_df, pred_vars, 'Ht')
    print('predicting with SAR')
    result <- validate(model_sar, val_df, 'Ht')
    result$predictors <- 'with_sar'
    result$fold <- fold
    results <- rbind(results, result)

    print('fitting model without SAR')
    model_nosar <- fit_model(model, model_config, train_df, pred_vars_nosar, 'Ht')
    print('predicting with SAR')
    result <- validate(model_nosar, val_df, 'Ht')
    result$predictors <- 'without_sar'
    result$fold <- fold
    results <- rbind(results, result)
    print(result)
    t2 <- Sys.time()
    cat('Fold runtime:', difftime(t2, t1, units="mins"), ' (m)\n')
  }

  counts <- df |> group_by(segment_landcover) |> summarise(num_samples=n()) |> ungroup()
  print(counts)
  cv_results <- results |>
    group_by(segment_landcover, predictors) |>
    summarise(
      R2=mean(R2, na.rm=TRUE),
      RMSE=mean(RMSE, na.rm=TRUE)
    ) |>
    ungroup() |>
    mutate(tile_num=as.integer(tile_num)) |>
    pivot_wider(names_from = predictors, values_from = c(R2, RMSE)) |>
    left_join(counts, by='segment_landcover')
  print(cv_results)
  return(cv_results)
}

evaluate_rmse_r2 <- function(predictions, truth, landcover, majority_lc_classes) {
  pred_values <- values(predictions, mat=FALSE)
  truth_values <- values(truth, mat=FALSE)
  lc_values <- values(landcover, mat=FALSE)
  print('-------------------')
  print(majority_lc_classes)
  print('-------------------')
  vals <- data.frame(
    predictions = pred_values,
    truth = truth_values,
    landcover = lc_values
  )
  print(vals |> group_by(landcover) |> summarise(count=n()))
  vals <- na.omit(vals)
  print(vals |> group_by(landcover) |> summarise(count=n()))

  # overal RMSE and R2 on all classes
  overall_metrics <- data.frame(
    landcover = -1,
    RMSE = sqrt(mean((vals$predictions - vals$truth)^2)),
    R2 = cor(vals$predictions, vals$truth, method = "pearson")^2
  )

  # Filter by majority land cover classes
  vals <- vals[vals$landcover %in% majority_lc_classes, ]
  print(vals |> group_by(landcover) |> summarise(count=n()))
  if (nrow(vals) == 0) {
    stop("No valid data for RMSE and R² calculation. All values may be NA.")
  }
  filtered_metrics <- data.frame(
    landcover = -2,
    RMSE = sqrt(mean((vals$predictions - vals$truth)^2)),
    R2 = cor(vals$predictions, vals$truth, method = "pearson")^2
  )

  print("Compute RMSE and R² by land cover classes")
  metrics <- vals %>%
    group_by(landcover) %>%
    summarise(
      RMSE = sqrt(mean((predictions - truth)^2)),
      R2 = cor(predictions, truth, method = "pearson")^2
    ) %>%
    ungroup()

  metrics <- bind_rows(metrics, overall_metrics, filtered_metrics)
  return(metrics)
}

fit_predict_evaluate <- function(model, model_config, training_df, pred_vars, predict_vars, with_without_sar){
  cat('fitting model :', with_without_sar)
  model <- fit_model(model, model_config, df, pred_vars, predict_var)
  cat('predicting :', with_without_sar)
  sar_map <- predict_stack(model, stack, majority_lc_classes)

  print('evaluating r2 and rmse with SAR')
  sar_lvis <- evaluate_rmse_r2(sar_map, lvis_mosaic, stack$esa_worldcover_v100_2020, majority_lc_classes_list)
  sar_lvis$predict_vars <- 'with_sar'

}

run_modeling_pipeline <-function(rds_models, all_train_data, model, model_config,
                                 max_samples, sample, pred_vars, pred_vars_nosar,
                                 predict_var, tile_num, stack, lvis_footprints_gpkg, year, folds=10){

  all_train_data_AGB <- GEDI2AT08AGB(rds_models, all_train_data, randomize=FALSE, max_samples, sample)
  df <- subset_to_majority_lc_classes(all_train_data_AGB)
  majority_lc_classes <- df |> distinct(segment_landcover)
  majority_lc_classes_list <- majority_lc_classes |> pull(segment_landcover))
  cv_results <- cross_validate(df, pred_vars, pred_vars_nosar, folds, model, model_config, tile_num)
  # print(cv_results)
  print(freq(stack$esa_worldcover_v100_2020))
  print('creating LVIS mosaic')
  lvis_mosaic <- create_lvis_mosaic(stack, lvis_footprints_gpkg, year)

  print('fitting model with sar')
  sar_model <- fit_model(model, model_config, df, pred_vars, predict_var)

  print('predicting with sar')
  sar_map <- predict_stack(sar_model, stack, majority_lc_classes)

  print('evaluating r2 and rmse with SAR')
  sar_lvis <- evaluate_rmse_r2(sar_map, lvis_mosaic, stack$esa_worldcover_v100_2020, majority_lc_classes_list)
  sar_lvis$predict_vars <- 'with_sar'

  print('fitting model without sar')
  nosar_model <- fit_model(model, model_config, df, pred_vars_nosar, predict_var)

  print('predicting without sar')
  nosar_map <- predict_stack(nosar_model, stack, majority_lc_classes)

  print('evaluating r2 and rmse without SAR')
  nosar_lvis <- evaluate_rmse_r2(nosar_map, lvis_mosaic, stack$esa_worldcover_v100_2020, majority_lc_classes_list)
  nosar_lvis$predict_vars <- 'without_sar'

  r2_rmse_lvis <- bind_rows(nosar_lvis, sar_lvis)
  print(r2_rmse_lvis)

  return(list(cv_results=cv_results, r2_rmse_lvis=r2_rmse_lvis))
}


resample_if_needed <- function(src, des){
  if (nrow(src) != nrow(des) || ncol(src) != ncol(des)){
    src <- resample(src, des, method = 'near')
    ext(src) <- ext(des)
  }
  return(src)
}

parse_pred_vars <- function(pred_vars){
  pred_vars <- unlist(strsplit(pred_vars, split = " "))
  sar_vars <- list(
    "vv_median_frozen", "vh_median_frozen",
    "vv_median_summer", "vh_median_summer",
    "vv_median_shoulder", "vh_median_shoulder"
  )
  pred_vars_nosar <- pred_vars[!pred_vars %in% sar_vars]

  return(list(pred_vars=pred_vars, pred_vars_nosar=pred_vars_nosar))
}

resample_if_needed <- function(src, des){
  if (nrow(src) != nrow(des) || ncol(src) != ncol(des)){
    src <- resample(src, des, method = 'near')
    ext(src) <- ext(des)
  }
  return(src)
}

prepare_raster <- function(path, subset_bands=NULL, extra_bands=NULL, dest_raster=NULL){
  raster <- rast(path)
  raster_bands <- names(raster)

  if (!is.null(subset_bands))
    raster_bands <- intersect(raster_bands, subset_bands)

  if (!is.null(extra_bands))
    raster_bands <- c(raster_bands, extra_bands)

  raster <- subset(raster, raster_bands)

  if (!is.null(dest_raster))
    raster <- resample_if_needed(raster, dest_raster)

  return(raster)
}

resample_reproject_and_mask <- function(topo_path, hls_path, lc_path, pred_vars, mask, sar_path){
  hls <- prepare_raster(hls_path, subset_bands=pred_vars, extra_bands='ValidMask')
  topo <- prepare_raster(topo_path, subset_bands=pred_vars, extra_bands='slopemask', dest_raster=hls)
  lc <- prepare_raster(lc_path, dest_raster=hls)
  sar <- prepare_raster(sar_path, subset_bands=pred_vars, dest_raster=hls)
  stack <- c(hls, sar, topo, lc)

  if(mask)
    stack <- mask_input_stack(stack)
  # stack <- subset(stack, names(stack) != 'esa_worldcover_v100_2020')
  return(stack)
}

mask_input_stack <- function(stack){
  MASK_LYR_NAMES = c('slopemask', 'ValidMask')
  MASK_LANDCOVER_NAMES = c(50, 60, 70, 80)

  print("Masking stack...")
  # Bricking the stack will make the masking faster (i think)
  # brick = rast(stack)
  for(LYR_NAME in MASK_LYR_NAMES){
    m <- terra::subset(stack, grep(LYR_NAME, names(stack), value = T))
    stack <- mask(stack, m == 0, maskvalue=TRUE)
  }

  for(LC_NAME in MASK_LANDCOVER_NAMES){
    n <- terra::subset(stack, grep('esa_worldcover_v100_2020', names(stack), value=LC_NAME))
    stack <- mask(stack, n == LC_NAME, maskvalue=TRUE)
  }

  return(stack)
}

mapBoreal<-function(atl08_path, broad_path, hls_path, topo_path, lc_path, year,
                    lvis_gpkg, sar_path, mask=TRUE, max_sol_el=0, offset=100, minDOY=1, maxDOY=365,
                    expand_training=TRUE, local_train_perc=100, min_samples=5000, max_samples=10000,
                    predict_var='AGB', pred_vars=c('elevation', 'slope', 'NDVI')){

  tile_num = tail(unlist(strsplit(path_ext_remove(atl08_path), "_")), n=1)
  cat("Modelling and mapping boreal AGB tile: ", tile_num, "\n")

  pred_vars_all <- parse_pred_vars(pred_vars)
  pred_vars <- pred_vars_all[['pred_vars']]
  pred_vars_nosar <- pred_vars_all[['pred_vars_nosar']]
  print(pred_vars)
  print(pred_vars_nosar)

  stack <- resample_reproject_and_mask(topo_path, hls_path, lc_path, pred_vars, mask, sar_path)
  print(names(stack))

  all_train_data <- prepare_training_data(
    atl08_path, broad_path, expand_training, minDOY,
    maxDOY, max_sol_el, min_samples, local_train_perc, offset, c(pred_vars, c('slopemask', 'ValidMask'))
  )

  fixed_modeling_pipeline_params <- list(
    rds_models=get_rds_models(), all_train_data=all_train_data,
    pred_vars=pred_vars, pred_vars_nosar=pred_vars_nosar, predict_var=predict_var, year=year,
    model=randomForest, sample=TRUE, tile_num=tile_num, stack=stack, lvis_footprints_gpkg=lvis_gpkg
  )

  results <- do.call(run_modeling_pipeline, modifyList(
    fixed_modeling_pipeline_params,
    list(max_samples=max_samples, model_config=list(ntree=100, mtry=6))
  ))

  output_fns <- set_output_file_names(predict_var, tile_num, year)
  write.csv(results[['cv_results']], output_fns[['cv']])
  write.csv(results[['r2_rmse_lvis']], output_fns[['r2_rmse_lvis']])
  print('Model comparison finished successfully!')
}

option_list <- list(
  make_option(
    c("-a", "--atl08_path"), type = "character",
    help = "Path to the atl08 training data"
  ),
  make_option(
    c("-b", "--broad_path"), type = "character",
    help = "Path to the boreal wide training data"
  ),
  make_option(
    c("-t", "--topo_path"), type = "character",
    help = "Path to the topo stack file"
  ),
  make_option(
    c("-h", "--hls_path"), type = "character",
    help = "Path to the HLS stack file"
  ),
  make_option(
    c("-l", "--lc_path"), type = "character",
    help = "Path to the land cover mask file"
  ),
  make_option(
    c("-s", "--sar_path"), type = "character", default = NULL,
    help = "Path to the land cover mask file"
  ),
  make_option(
    c("--lvis_gpkg"), type = "character", default = NULL,
    help = "Path to lvis"
  ),
  make_option(
    c("-y", "--year"), type = "character",
    help = "Year of the input HLS imagery"
  ),
  make_option(
    c("-m", "--mask"), type = "logical", default = TRUE,
    help = "Whether to mask imagery [default: %default]"
  ),
  make_option(
    c("--max_sol_el"), type = "numeric", default = 5,
    help = "Maximum solar elevation degree allowed in training data [default: %default]"
  ),
  make_option(
    c("--minDOY"), type = "integer", default = 130,
    help = "Minimum day of year allowed in training data [default: %default]"
  ),
  make_option(
    c("--maxDOY"), type = "integer", default = 250,
    help = "Maximum day of year allowed in training data [default: %default]"
  ),
  make_option(
    c("--min_samples"), type = "integer", default = 5000,
    help = "Minimum number of samples to avoid augmenting the training data with broad data [default: %default]"
  ),
  make_option(
    c("--max_samples"), type = "integer", default = 10000,
    help = "Maximum number of samples used for training [default: %default]"
  ),
  make_option(
    c("-e", "--expand_training"), type = "logical", default = TRUE,
    help = "Whether to expand training around the season [default: %default]"
  ),
  make_option(
    c("-p", "--local_train_perc"), type = "integer", default = 100,
    help = "Percent of atl08 data to be used in case it is augmented with broad data [default: %default]"
  ),
  make_option(
    c("--predict_var"), type = "character", default = "AGB",
    help = "Variable to predict, it can be either AGB or Ht [default: %default]"
  ),
  make_option(
    c("--pred_vars"), type = "character",
    default = paste(
      "Red Green elevation slope tsri tpi NIR SWIR SWIR2 NDVI",
      "SAVI MSAVI NDMI EVI NBR NBR2 TCB TCG TCW",
      "vv_median_frozen vh_median_frozen vv_median_summer",
      "vh_median_summer vv_median_shoulder vh_median_shoulder"
    ),
    help = paste(
      "List of predictor variables, must be a subset from the default options",
      "seperated by space, e.g, NDVI slope\n [default: %default]"
    )
  ),
  make_option(
    c("--help"), action = "store_true",
    help = "Show help message"
  )
)

opt_parser <- OptionParser(option_list = option_list, add_help_option = FALSE)
opt <- parse_args(opt_parser)

cat("Parsed arguments:\n")
print(opt)

if (!is.null(opt$help)) {
  print_help(opt_parser)
}
set.seed(123)
do.call(mapBoreal, opt)
