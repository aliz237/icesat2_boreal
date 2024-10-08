
#----------------------------------------------#
# loading packages
#----------------------------------------------#

#require(pacman)
#pacman::p_load(randomForest, raster,rgdal)


# 2.6 ICESat-2 biomass 30m ATL08

#### i) Description
#### This notebook use GEDI-derived AGB models for estimating sparse ICESat-2 AGB using 30 m ATL08 product.

#### ii) How the algorithm works?
#### GEDI-derived AGB models are loaded in R and applied over the ICESat-2 30m ATL08 data. 

#### iii) Inputs
####  - rds_models: list of GEDI-derived AGB model paths
####  - models_id: models id
####  - ice2_30_atl08: list containing the path to the ATL08 files
####  - offset: offset applied in the model

#### iii) Outputs
####  -ICEsat-2 AGB estimates are data.frame object in R

#----------------------------------------------#
############# functions ########################
#----------------------------------------------#
GEDI2AT08AGB<-function(rds_models,models_id,ice2_30_atl08_path, offset=100){
  
  # rds_models
  names(rds_models)<-models_id
  # read table
  xtable<-read.table(ice2_30_atl08_path, sep=",", head=T)
  xtable_i<-na.omit(as.data.frame(xtable))
  names(xtable_i)[1:11]<-c("lon","lat","RH_25","RH_50","RH_60","RH_70","RH_75","RH_80","RH_90","RH_95","RH_98")
  
  #
  xtable_sqrt<-xtable_i[3:11]+offset
  #names(xtable_sqrt)<-paste0("sqrt(",names(xtable_sqrt),")")
  #xtable_sqrt<-cbind(xtable_i,xtable_sqrt)
  #xtable_sqrt<-xtable_i
  #xtable_sqrt[,3:12]<-xtable_sqrt[,3:12]
  
  # get unique ids

  # apply models by id
  xtable_sqrt$AGB<-NA
  xtable_sqrt$SE<-NA
    
  #Assign model id based on landcover

    # 1 = DBT trees (boreal-wide), 2=Evergreen needleleaf trees (boreal-wide), 12 = boreal-wide all PFT
    
    # for the seg_landcov: {0: "water", 1: "evergreen needleleaf forest", 2: "evergreen broadleaf forest", \ 3: "deciduous needleleaf forest", 4: "deciduous broadleaf forest", \ 5: "mixed forest", 6: "closed shrublands", 7: "open shrublands", \ 8: "woody savannas", 9: "savannas", 10: "grasslands", 11: "permanent wetlands", \ 12: "croplands", 13: "urban-built", 14: "croplands-natural mosaic", \ 15: "permanent snow-ice", 16: "barren"})
    
    xtable_sqrt$model_id <- NA
    xtable_sqrt$model_id[xtable_sqrt$seg_landcov==0] <- 12
    xtable_sqrt$model_id[xtable_sqrt$seg_landcov==1] <- 4
    xtable_sqrt$model_id[xtable_sqrt$seg_landcov==2] <- 4
    xtable_sqrt$model_id[xtable_sqrt$seg_landcov==3] <- 1
    xtable_sqrt$model_id[xtable_sqrt$seg_landcov==4] <- 1
    xtable_sqrt$model_id[xtable_sqrt$seg_landcov==5] <- 12
    xtable_sqrt$model_id[xtable_sqrt$seg_landcov==6] <- 12
    xtable_sqrt$model_id[xtable_sqrt$seg_landcov==7] <- 12
    xtable_sqrt$model_id[xtable_sqrt$seg_landcov==8] <- 12
    xtable_sqrt$model_id[xtable_sqrt$seg_landcov==9] <- 12
    xtable_sqrt$model_id[xtable_sqrt$seg_landcov==10] <-12
    xtable_sqrt$model_id[xtable_sqrt$seg_landcov==11] <- 12
    xtable_sqrt$model_id[xtable_sqrt$seg_landcov==12] <- 12
    xtable_sqrt$model_id[xtable_sqrt$seg_landcov==13] <- 12
    xtable_sqrt$model_id[xtable_sqrt$seg_landcov==14] <- 12
    xtable_sqrt$model_id[xtable_sqrt$seg_landcov==15] <- 12
    xtable_sqrt$model_id[xtable_sqrt$seg_landcov==16] <- 12

  xtable_sqrt$model_id<-names(rds_models)[1]
  ids<-unique(xtable_sqrt$model_id)

    
  #i=ids[1]
  for ( i in ids){
    
    # subset data for model id
    model_i<-readRDS(rds_models[names(rds_models)==i])
    
    # SE
    xtable_sqrt$SE[xtable_sqrt$model_id==i] <- summary(model_i)$sigma
    
    # AGB prediction
    xtable_sqrt$AGB[xtable_sqrt$model_id==i]<-predict(model_i, newdata=xtable_sqrt[xtable_sqrt$model_id==i,])
    
    #adjust for offset in model fits (note, this was 20 for ages, and now it's 100; essentially we added to all the height metrics so they would never be negative)
    #xtable_sqrt$AGB[xtable_sqrt$model_id==i] <- xtable_sqrt$AGB[xtable_sqrt$model_id==i]+offset
    
    #define C
    C <- mean(model_i$fitted.values)/mean(model_i$model$`sqrt(AGBD)`)
    
    #we multiply by C in case there is a systematic over or under estimation in the model (bias correction)
    xtable_sqrt$AGB[xtable_sqrt$model_id==i]<-C*(xtable_sqrt$AGB[xtable_sqrt$model_id==i]^2)
    #print(head(xtable_sqrt$AGB[xtable_sqrt$model_id==i]))
  }
  xtable2<-cbind(xtable, xtable_sqrt$AGB, xtable_sqrt$SE)#[,c(names(xtable_i),"AGB","SE")]
    ncol <- ncol(xtable2)
  colnames(xtable2)[(ncol-1):ncol]<-c('AGB', 'SE')
  return(xtable2)
}



# 2.7 ICESat-2 biomass algorithm development

#### i) Description
#### This notebook uses the sparse ICESat-2 AGB estimates and remote sensing covariates data (ex. extracted from ALOS-2, Landsat 8 OLI, Sentinel 2A...etc) for calibrating and testing AGB models using parametric and non-parametric statistical modeling approaches, such as OLS, RF, k-NN, SVM and ANN. 

#### ii) How the algorithm works?
#### Users must to select the number of bootstrapping runs (ex. 100). In each run, the original dataset is divided into training (ex. 70%) and testing (ex. 30%) datastes for model calibration and validation. Users can select if they want to create the training and testing dataset using a random or stratified random sampling approach. R2, RMSE and MD are computed based on the training dataset. 

#### iii) Inputs
####  -rep: number of bootstrapping runs
####  -y: vector of the response variable in the dataset (ex. "g_agbm")
####  -x: dataframe containing the name of the covariates
####  -s_train: percentage used selected the training dataset (ex.70), 
####  -stack: stack of covar
####  -strat_random: If TRUE the original dataset will be splitted into training and testing using the stratified random approach (bins are defined by the quartiles). Otherwise, the random approach will be used. 

#### iii) Outputs
#### Maps of AGB (mean and sd)
#### Maps (mean and sd) of the prediction accuracy (RMSE, bias and R2 for the testig datasets)
#### Legend:
#### agb_mean = mean of the AGB predictions 
#### agb_sd = sd of the AGB predictions (n = rep)
#### armse_map_mean = mean of absolute RMSE values
#### armse_map_sd = sd of absolute RMSE values
#### rrmse_map_mean = mean of relative RMSE values
#### rrmse_map_sd = sd of relative RMSE values
#### abias_map_mean = mean of absolute RMSE values
#### abias_map_sd = sd of absolute RMSE values
#### rbias_map_mean = mean of relative of bias values
#### rbias_map_sd = sd of relative bias values
#### r2_map_mean = mean of r2 values
#### r2_map_sd = sd of r2 values 


#----------------------------------------------#
############# functions ########################
#----------------------------------------------#

# stats
StatModel <- function( obs, est){
  xy<-na.omit(cbind(obs,est))
  obs<-xy[,1]
  est<-xy[,2]
  rmse <- sqrt( sum(( est - obs )^2)/length(obs) ) # Root mean square error
  bias <- mean( est - obs ) # bias
  rmseR <- 100 * sqrt( sum(( est - obs )^2)/length(obs) ) / mean( obs )
  biasR <- 100 * mean( est - obs ) / mean( obs )
  r <- cor(est,obs)
  r2<-summary(lm(obs~est))$r.squared
  Stats<-data.frame( Stat=c("rmse","rmseR","bias","biasR","r","r2"),
                     Values=round(c(rmse,rmseR,bias,biasR,r,r2),2)) 
  return(Stats)
}

stratRandomSample<-function(agb=y,breaks, p){
  #require(data.table)
  n<-length(agb)
  ids<-1:n
  s<-round(n*p)
  agb[agb==0]<-0.0000000001
  ids_cut<-cut(agb,breaks=breaks, labels=F)
  df<-cbind(agb,ids,ids_cut)
  df<-data.table(df[!is.na(df[,1]),])
  number_sample<-ceiling(s/(length(breaks)-1))
  sel_all<-df[,.SD[sample(.N, min(number_sample,.N), replace = T)],by=ids_cut]
  return(ids_selected=sel_all$ids)
}

# modeling 
agbModelingMapping<-function(x=x,y=y,se=NULL,s_train=70, rep=10,stack=stack,strat_random=TRUE,output){
  
  pb <- txtProgressBar(min = 0, max = rep, style = 3)
  
  stats_df<-NULL
  n<-nrow(x)
  ids<-1:n
  map_pred<-NULL
  
  #p_load("VSURF")
  #varsel <- VSURF(y=y, x=x, ntree=1000)
  #varnames<-colnames(x)[varsel$varselect.interp]
  stack_df <- na.omit(as.data.frame(stack, xy=TRUE))
  stack_df$grid_id<-1:nrow(stack_df)
  
  #i=10
  i.s=0
  #j=1
  
  for (j in 1:rep){
    i.s<-i.s+1
    setTxtProgressBar(pb, i.s)
    set.seed(j)
    if ( strat_random==TRUE){
      trainRowNumbers<-stratRandomSample(agb=y,breaks=quantile(y, na.rm=T), p=s_train/100)
    } else {
      trainRowNumbers<-sort(sample(ids,round(n*s_train/100), T))
    }
    # Step 2: Create the training  dataset
    # select % of data to training and testing the models
    trainData.x <- x[trainRowNumbers,]
    trainData.y <- y[trainRowNumbers]
    
    
    # Step 3: Create the test dataset
    # select % of the data for validation
    testData.x <- x[!row.names(x) %in% trainRowNumbers,]
    testData.y <- y[!row.names(x) %in% trainRowNumbers]
    
    # rf modeling
    set.seed(j)
    
    fit.rf <- randomForest(y=trainData.y, x=trainData.x, ntree=1000)
    pred.rf<-predict(fit.rf, newdata=testData.x) # 
    #fit.rf <- randomForest(y=trainData.y, x=trainData.x[,varnames], ntree=1000)
    #pred.rf<-predict(fit.rf, newdata=testData.x[,varnames]) # 
    stats.rf<-cbind(method=rep("RF",6), rep=rep(j,6), StatModel(testData.y,pred.rf))
    
    # model validation stats
    stats_df<-rbind(stats_df,
                    stats.rf)
    row.names(stats_df)<-1:nrow(stats_df)
    
    # mapping
    #map_i<-cbind(stack_df[,1:2],agb=predict(fit.rf, newdata=stack_df[,varnames]), rep=i)
    map_i<-cbind(stack_df[,1:2],agb=predict(fit.rf, newdata=stack_df), PI=predict(fit.rf, newdata=stack_df, interval='predict', CI=predict(fit.rf, newdata=stack_df, interval='confidence'), rep=j, grid_id=stack_df$grid_id)
    #map_pred<-rbind(map_pred,map_i)
    
  }
  
  
  #mean_map<-tapply(map_pred$agb,map_pred$grid_id,mean)
  #sd_map<-tapply(map_pred$agb,map_pred$grid_id,sd)
  
  
  #agb_maps <- rasterFromXYZ(cbind(stack_df[,1:2],
  #                                agb_mean=mean_map,
  #                                agb_sd=sd_map,
  #                                armse_map_mean=mean(stats_df[stats_df[,3]=="rmse",4]),
  #                                armse_map_sd=sd(stats_df[stats_df[,3]=="rmse",4]),
  #                                
  #                                rrmse_map_mean=mean(stats_df[stats_df[,3]=="rmseR",4]),
  #                                rrmse_map_sd=sd(stats_df[stats_df[,3]=="rmseR",4]),
  #                                
  #                                abias_map_mean=mean(stats_df[stats_df[,3]=="bias",4]),
  #                                abias_map_sd=sd(stats_df[stats_df[,3]=="bias",4]),
  #                                
  #                                rbias_map_mean=mean(stats_df[stats_df[,3]=="bias",4]),
  #                                rbias_map_sd=sd(stats_df[stats_df[,3]=="bias",4]),
  #                                
  #                                r2_map_mean=mean(stats_df[stats_df[,3]=="r2",4]),
  #                                r2_map_sd=sd(stats_df[stats_df[,3]=="r2",4])))
  
    
    agb_maps <- rasterFromXYZ(cbind(stack_df[,1:2],
                                   agb_pred=agb,
                                   agb_CI_low=CI$lwr,
                                   agb_CI_high=CI$upr,
                                   agb_PI_low=PI$lwr,
                                   agb_PI_high=PI$upr))
  close(pb)
  return(agb_maps)
}

mapBoreal<-function(rds_models,
                    models_id,
                    ice2_30_atl08_path, 
                    offset=100,
                    s_train=70, 
                    rep=10,
                    stack,
                    strat_random=TRUE,
                    output){
  
  # apply GEDI models  
  xtable<-GEDI2AT08AGB(rds_models=rds_models,
                       models_id=models_id,
                       ice2_30_atl08_path=ice2_30_atl08_path, 
                       offset=offset)
  hist(xtable$AGB)
  # run 
    pred_vars <- c('elevation', 'slope', 'tsri', 'tpi', 'NDVI', 'SAVI', 'MSAVI', 'NDMI', 'EVI', 'NBR', 'NBR2', 'TCW')

maps<-agbModelingMapping(x=xtable[pred_vars],
                           y=xtable$AGB,
                           se=xtable$SE,
                           s_train=s_train,
                           rep=rep,
                           stack=stack,
                           strat_random=strat_random)

  writeRaster(maps,output,overwrite=T)
  return(maps)
}


