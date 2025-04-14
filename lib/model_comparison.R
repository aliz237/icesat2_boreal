library(optparse)
#library(randomForest)
library(ranger)
library(tidyr)
library(dplyr)
library(fs)
library(stringr)
library(rockchalk)
library(terra)

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

  fn_suffix <- c('.tif', '_summary.csv', '_train_data.csv', '_stats.Rds', '_model.Rds', '_comparison.csv')
  names <- c('map', 'summary', 'train', 'stats', 'model', 'comparison')

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
  rds_model_fns <- list.files(path='~/Downloads/bio_models_noground', pattern='*.rds', full.names=TRUE)
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

validate <- function(df, predict_var){
  y_true <- if (predict_var == 'Ht') df$h_canopy else df$AGB
  predict_var <- if (predict_var == 'Ht')  'h_canopy' else 'AGB'
  df$residuals_SAR <- df[['SAR']] - y_true
  df$residuals_no_SAR <- df[['NO_SAR']] - y_true
  df <- na.omit(df)
  result <- data.frame()
  for (var in c('segment_landcover', 'height_class', 'slope_class')){
    result <- rbind(
      result,
      df |>
        group_by(across(all_of(var))) |>
        summarise(
          num_samples = n(),
          RMSE_SAR = sqrt(mean(residuals_SAR^2, na.rm = TRUE)),
          RMSE_PCT_SAR=sqrt(mean(residuals_SAR^2, na.rm = TRUE)) / mean(.data[[predict_var]], na.rm=TRUE) * 100,
          MAE_SAR=mean(abs(residuals_SAR)),
          MAE_PCT_SAR=mean(abs(residuals_SAR))/mean(.data[[predict_var]], na.rm=TRUE) * 100,
          bias_SAR=mean(residuals_SAR),
          bias_PCT_SAR=mean(residuals_SAR)/mean(.data[[predict_var]], na.rm=TRUE) * 100,
          R2_SAR = 1 - (sum(residuals_SAR^2, na.rm = TRUE) / sum((.data[[predict_var]] - mean(.data[[predict_var]], na.rm = TRUE))^2, na.rm = TRUE)),
          RMSE_no_SAR = sqrt(mean(residuals_no_SAR^2, na.rm = TRUE)),
          RMSE_PCT_no_SAR=sqrt(mean(residuals_no_SAR^2, na.rm = TRUE)) / mean(.data[[predict_var]], na.rm=TRUE) * 100,
          MAE_no_SAR=mean(abs(residuals_no_SAR)),
          MAE_PCT_no_SAR=mean(abs(residuals_no_SAR))/mean(.data[[predict_var]], na.rm=TRUE) * 100,
          bias_no_SAR=mean(residuals_no_SAR),
          bias_PCT_no_SAR=mean(residuals_no_SAR)/mean(.data[[predict_var]], na.rm=TRUE) * 100,
          R2_no_SAR = 1 - (sum(residuals_no_SAR^2, na.rm = TRUE) / sum((.data[[predict_var]] - mean(.data[[predict_var]], na.rm = TRUE))^2, na.rm = TRUE))
        ) |>
        ungroup() |>
        mutate(group_type=var) |>
        rename(group_value=all_of(var))
    )
  }
  return(result)
}

add_height_and_slope_class <- function (df){

  # Add height class and counts using RH_98 column
  labels <- c("0-1 m", "1-2 m", "2-3 m", "3-4 m", "4-5 m", "5-10 m", "10-15 m", "15-20 m", "20-25 m", "25-30 m", "30> m")
  breaks <- c(0, 1, 2, 3, 4, 5, 10, 15, 20, 25, 30, Inf)
  df$height_class <- cut(df$h_canopy, breaks = breaks , include.lowest = TRUE, labels = labels, right = FALSE)
  df <- df |>
    group_by(height_class) |>
    summarise(n_height_class=n()) |>
    right_join(df, by='height_class') |>
    filter(n_height_class >= 1) |>
    select(-n_height_class)

  # Add Slope Class and counts
  labels <- c("0-5", "5-10", "10-15", "15-20", "20-25", "25-30", "30-35", "35-40", "40-45", "45-90")
  df$slope_class <- cut(df$slope, breaks = c(seq(0, 45, 5), 90), right = FALSE, include.lowest = TRUE, labels=labels)
  df <- df |>
    group_by(slope_class) |>
    summarise(n_slope_class=n()) |>
    right_join(df, by='slope_class') |>
    filter(n_slope_class >= 1) |>
    select(-n_slope_class)

  # Add Landcover class counts
  df <- df |>
    group_by(segment_landcover) |>
    summarise(n_lc_class=n()) |>
    right_join(df, by='segment_landcover') |>
    select(-n_lc_class)

  return(df)
}

run_modeling_pipeline <-function(rds_models, all_train_data, model, model_config,
                                 max_samples, sample, pred_vars, pred_vars_nosar, predict_var, tile_num, folds=10){

  print('creating AGB traing data frame.')
  #stop('ss')
  all_train_data_AGB <- GEDI2AT08AGB(rds_models, all_train_data, randomize=FALSE, max_samples, sample)
  all_train_data_AGB <- add_height_and_slope_class(all_train_data_AGB)
  df <- setup_kfold(all_train_data_AGB, folds)
  # run a 10-fold CV with and without SAR features
  results <- data.frame()
  df$SAR <- NA
  df$NO_SAR <- NA

  for (fold in 1:folds) {
    t1 <- Sys.time()
    cat('fold:', fold, '\n')

    train_df <- df[df['folds'] != fold, ]
    val_df <- df[df['folds'] == fold, ]

    print('fitting model with SAR')
    model_sar <- fit_model(model, model_config, train_df, pred_vars, predict_var)
    print('predicting with SAR')
    df[df['folds'] == fold, ]$SAR<- predict(model_sar,  df[df['folds'] == fold, ])$predictions
    print(head(df[df['folds'] == fold, c('folds', 'SAR', 'h_canopy', 'RH_98')]))
    print('fitting model without SAR')
    model_nosar <- fit_model(model, model_config, train_df, pred_vars_nosar, predict_var)
    print('predicting without SAR')
    df[df['folds'] == fold, ]$NO_SAR<- predict(model_nosar,  df[df['folds'] == fold, ])$predictions

    t2 <- Sys.time()
    cat('Fold runtime:', difftime(t2, t1, units="mins"), ' (m)\n')
  }
  

  final_results <- validate(df, predict_var)
  
  return(final_results)
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

mapBoreal<-function(atl08_path, broad_path, boreal_vector_path, year,
                    max_sol_el=0, offset=100, minDOY=1, maxDOY=365,
                    expand_training=TRUE,
                    local_train_perc=100, min_samples=5000, max_samples=10000,
                    predict_var='AGB', pred_vars=c('elevation', 'slope', 'NDVI')){

  tile_num = tail(unlist(strsplit(path_ext_remove(atl08_path), "_")), n=1)
  cat("Modelling and mapping boreal AGB tile: ", tile_num, "\n")

  pred_vars_all <- parse_pred_vars(pred_vars)
  pred_vars <- pred_vars_all[['pred_vars']]
  pred_vars_nosar <- pred_vars_all[['pred_vars_nosar']]
  print(pred_vars)
  print(pred_vars_nosar)

  all_train_data <- prepare_training_data(
    atl08_path, broad_path, expand_training, minDOY,
    maxDOY, max_sol_el, min_samples, local_train_perc, offset, c(pred_vars, c('slopemask', 'ValidMask'))
  )

  fixed_modeling_pipeline_params <- list(
    rds_models=get_rds_models(), all_train_data=all_train_data,
    pred_vars=pred_vars, pred_vars_nosar=pred_vars_nosar, predict_var=predict_var,
    model=ranger, sample=TRUE, tile_num=tile_num
  )

  results <- do.call(run_modeling_pipeline, modifyList(
    fixed_modeling_pipeline_params,
    list(max_samples=max_samples, model_config=list(num.trees=50, mtry=6))
  ))

  output_fns <- set_output_file_names(predict_var, tile_num, year)
  write.csv(results, output_fns[['comparison']])
  print('AGB successfully predicted!')
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
    c("-y", "--year"), type = "character",
    help = "Year of the input HLS imagery"
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
setwd('~/')
do.call(mapBoreal, opt)
