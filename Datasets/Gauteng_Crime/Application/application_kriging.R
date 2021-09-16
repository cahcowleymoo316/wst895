setwd("H:/Personal/Data Science/MSc/WST 895")
ds.name <- "Gauteng_Crime"
section <- "Application"
library(sp)
library(spdep)
library(gstat)
library(tidyverse)
library(scales) 
library(RColorBrewer)
library(Metrics)
library(rgeos)
library(caret)
# Load Functions
source("Datasets/esda_functions.R")


# Load Data
poly.sp <- readRDS(paste0("Datasets/", ds.name, "/", ds.name, "_Data.RDS"))
point.sp <- readRDS(paste0("Datasets/", ds.name, "/", ds.name, "_Point_Data.RDS"))

### KRIGING/GP MODEL ###

# Empirical variogram
par(mfrow=c(1,1))
vgm <- gstat::variogram(log(crime_vehicle_1000hh_use)~1, data = point.sp) 
plot(x=vgm$dist
     , y = vgm$gamma
     , pch = 20
     , main = "Empirical Semivariogram"
     , sub = "Vehicle Crime"
     , xlab = "h"
     , ylab = expression(paste(gamma,"(h)", sep = "")))

fit.sph <- gstat::fit.variogram(vgm, model=vgm(psill = 0.5, model = "Sph", range = 0.4, nugget = 0.1))
line.sph <- gstat::variogramLine(fit.sph,maxdist = max(vgm$dist))

fit.exp <- gstat::fit.variogram(vgm, model=vgm(psill = 0.5, model = "Exp", range = 0.4, nugget = 0.1))
line.exp <- gstat::variogramLine(fit.exp,maxdist = max(vgm$dist))


fit.gau <- gstat::fit.variogram(vgm, model=vgm(model = "Gau", range = 0.3))
line.gau <- gstat::variogramLine(fit.gau,maxdist = max(vgm$dist))


fit.Mat <- gstat::fit.variogram(vgm, model=vgm(psill = 0.5, "Mat", kappa = 0.4))
line.Mat <- gstat::variogramLine(fit.Mat,maxdist = max(vgm$dist))

all.fits <- data.frame(dist = line.exp$dist
                       , Exponential = line.exp$gamma
                       , Spherical = line.sph$gamma
                       , Gaussian = line.gau$gamma
                       , Matern = line.Mat$gamma)

png(paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name), "_kriging_0_vehicle_semivariograms.png"), width = 750, height = 500)
par(mfrow = c(1,1), mar = rep(4.4,4), bg = "Grey80", col = "black")
plot(x=vgm$dist
     , y = vgm$gamma
     , pch = 20
     , main = "Fitted Semivariogram Models"
     , xlab = "h"
     , ylab = expression(paste(gamma,"(h)", sep = "")))
lines(x  = all.fits$dist, y = all.fits$Exponential, col = "blue")
lines(x  = all.fits$dist, y = all.fits$Spherical, col = "red")
lines(x  = all.fits$dist, y = all.fits$Gaussian, col = "Grey20")
lines(x  = all.fits$dist, y = all.fits$Matern, col = "orange")
legend("bottomright"
       , title = "Semivariogram"
       , legend = c("Exponential", "Spherical", "Gaussian", "Matérn")
       , fill = c("blue", "red", "grey20", "orange")
       , inset = 0.0125)
dev.off()
# 100 times 10-fold cross validation
set.seed(12345)
crime_veh_kriged <- gstat::krige.cv(formula = log(crime_vehicle_1000hh_use) ~ 1
                              , locations = point.sp
                              , model = fit.sph
                              , nfold = 10
                              , nmax = 100)
summary(crime_veh_kriged)

# Residual analysis
# Translate residual to original units
crime_veh_kriged@data$eresidual <- exp(crime_veh_kriged$observed) - exp(crime_veh_kriged$var1.pred)

png(paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name), "_kriging_1_vehicle_crime_res_plots.png"), width = 750, height = 500)
par(mfrow = c(2,2), mar = rep(3,4), bg = "grey80")
hist(crime_veh_kriged$eresidual
     , xlab = "Residual"
     , main = "Histogram of Residuals")
text(x = -20
     ,y = 800
     , paste("MSE =\n", round(mean(crime_veh_kriged$eresidual^2), 4)))
# Mean Square Normalised Error
hist(crime_veh_kriged$zscore
     , xlab = "Normalised Residual"
     , main = "Histogram of \nNormalised Residuals")
text(x = -4
     ,y = 600
     , paste("MSE =", round(mean(crime_veh_kriged$zscore^2), 4)))
# Correlation observed and predicted, ideally close to 1
plot(x = exp(crime_veh_kriged$observed)
     , y = exp(crime_veh_kriged$var1.pred)
     #, col = crime_veh_kriged$fold
     , pch = 20
     , main = "Observed vs. Predicted"
     , xlab = "Observed"
     , ylab = "Predicted")
text(x = 47.50, y = 7.5
     , paste("Corr =", round(cor(exp(crime_veh_kriged$observed), exp(crime_veh_kriged$var1.pred)),4)))
# Correlation predicted and residual, ideally 0
plot(x = exp(crime_veh_kriged$var1.pred)
     , y = crime_veh_kriged$eresidual
     #, col = crime_veh_kriged$fold
     , pch = 20
     , main = "Predicted vs. Residual"
     , xlab = "Precicted"
     , ylab = "Residual")
text(x = 46, y = -25
     , paste("Corr =", round(cor(exp(crime_veh_kriged$var1.pred), exp(crime_veh_kriged$eresidual)),4)))
dev.off()

crime_veh_kriged_est_colours <- fn_mapcolours(var = exp(crime_veh_kriged@data$var1.pred), brewerpal = "YlGn", num_colours = 5, dig.lab = 2)

crime_veh_kriged_var_colours <- fn_mapcolours(var = exp(crime_veh_kriged@data$var1.var), brewerpal = "YlGn", num_colours = 5, dig.lab = 2)
png(paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name), "_kriging_5_estimate_variance.png"), width = 375, height = 250)
par(mfrow = c(1,1), mar = rep(3,4), bg = "Grey80")
fn_spdfplot(spdf = poly.sp
            , colours = crime_veh_kriged_var_colours
            , main = "(c), Kriging estimates")
dev.off()


# Fitted values
crime_veh_kriged_pred_bins <- data.frame(class = cut(exp(crime_veh_kriged$var1.pred), breaks = seq(from = 0, to = 60, by = 5))) %>% 
        dplyr::group_by(class) %>% 
        dplyr::summarise(n = n()) %>% 
        dplyr::mutate(midpoint = seq(2.5, 57.5, by = 5))

png(paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name), "_kriging_2_vehicle_crime_model_fit_plots.png"), width = 750, height = 500)
par(mfrow = c(2,2), mar = rep(3.8,4), bg = "Grey80")
hist(poly.sp$crime_vehicle_1000hh_use
     , xlab = "Vehicle Crime"
     , ylim = c(0, 850)
     , main = "(a), Histogram of \nActual vs. Fitted Values")
lines(x = crime_veh_kriged_pred_bins$midpoint
      , y = crime_veh_kriged_pred_bins$n
      , col = "red")
text(x = 50
     , y = 100
     , "Fitted"
     , col = "red")
plot(x = poly.sp@data$crime_res_1000hh_use
     , y = poly.sp@data$crime_vehicle_1000hh_use
     , pch = 20
     , ylab = "Vehicle Crime"
     , xlab = "Residential Crime"
     , main = "(b), Actual vs. Fitted")
points(x = poly.sp@data$crime_res_1000hh_use
       , y = exp(crime_veh_kriged@data$var1.pred)
       , col = scales::alpha("red", 0.3)
       , pch = 19)
fn_spdfplot(spdf = poly.sp
            , colours = crime_veh_kriged_est_colours
            , main = "(c), Kriging estimates")
plot(0,type='n',axes=FALSE,ann=FALSE)
legend("top"
       , title = "Estimate"
       , legend = levels(crime_veh_kriged_est_colours$class)
       , fill=crime_veh_kriged_est_colours$uniquecolour)
dev.off()


### ACCURACY MEASURES ###

# Calculate MSE
crime_veh_kriged_mse <- Metrics::mse(actual = poly.sp$crime_vehicle_1000hh_use
                              , predicted = exp(crime_veh_kriged$var1.pred))
crime_veh_kriged_mse
# Add plot of individual squared errors


# Calculate Global Moran's I - Actual
crime_veh_global_actual <- read.csv(paste0("Datasets/", ds.name, "/esda/esda_", tolower(ds.name) , "_3_global_moran_i.csv")) %>% 
        dplyr::filter(Variable == "crime_vehicle_1000hh_use")
# Calculate Global Moran's I - Predicted
neighbour <- readRDS(paste0("Datasets/", ds.name, "/esda/esda", "_", tolower(ds.name), "_neighbours.RDS"))
crime_veh_global_pred <- spdep::moran(x = exp(crime_veh_kriged$var1.pred)
                                , listw = neighbour$lw
                                , n = 1000
                                , S0 = spdep::Szero(neighbour$lw)
                                , zero.policy = T)
# Calculate difference
crime_veh_global_diff <- abs(crime_veh_global_actual$Global.Moran.s.I - crime_veh_global_pred$I)
crime_veh_global_diff_df <- crime_veh_global_actual %>% 
  dplyr::mutate(Predicted = crime_veh_global_pred$I) %>% 
  dplyr::rename(Actual = Global.Moran.s.I) %>% 
  reshape2::melt(id.vars = "Variable") %>% 
  dplyr::select(-Variable)
png(paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name), "_kriging_3_vehicle_crime_global_corr_plot.png"), width = 375, height = 250)
par(mfrow=c(1,1), mar = rep(4,4), bg = "grey80")
barplot(crime_veh_global_diff_df$value
        , names.arg = paste0(crime_veh_global_diff_df$variable, "\n(", round(crime_veh_global_diff_df$value,4),")")
        , col = c("black", "red"), 0.4
        , main = "Global Moran's I:\n Actual vs.Predicted")
dev.off()

# Calculate Local Moran's I - Actual
crime_veh_local_actual <- readRDS(paste0("Datasets/", ds.name, "/esda/esda_", tolower(ds.name) , "_vehicle_crime_local_obj.RDS"))
# Calculate Local Moran's I - Predicted
crime_veh_kriged$expvar1.pred <- exp(crime_veh_kriged$var1.pred)
crime_veh_local_pred <-  fn_localmoran(spdf = crime_veh_kriged
                                       , lw = neighbour$lw
                                       , var = "expvar1.pred"
                                       , zero.policy = F)



# Get  LISA clusters from predicted Moran's I

crime_veh_local_quadrant <- readRDS(paste0("Datasets/", ds.name, "/esda/esda_", tolower(ds.name) , "_vehicle_crime_local_quadrant.RDS"))
crime_veh_local_quadrant_pred <- fn_quadrant(spdf = crime_veh_local_pred
                                        , lw = neighbour$lw
                                        , var = "expvar1.pred")
table(crime_veh_local_quadrant_pred$quadrant)
quadrant_compare <- data.frame(Actual = crime_veh_local_quadrant$quadrant
                               , Predicted = crime_veh_local_quadrant_pred$quadrant)

quadrant_compare_matrix <- matrix(nrow = 5, ncol = 5)

crime_veh_confusion <- caret::confusionMatrix(as.factor(quadrant_compare$Actual)
                       , as.factor(quadrant_compare$Predicted)
                       , dnn = c("Predicted", "Actual"))
crime_veh_confusion_df <- data.frame(t(crime_veh_confusion$byClass)) %>% 
  dplyr::mutate(`Metric` = NA)
names(crime_veh_confusion_df) <- c("insignificant (I)","low-low (L-L)","high-high (H-H)", "Metric")
crime_veh_confusion_df$Metric <- rownames(crime_veh_confusion_df)
write.csv(crime_veh_confusion_df
          , paste0("Datasets/"
                   , ds.name
                   , "/"
                   , section
                   , "/"
                   , section
                   , "_"
                   , tolower(ds.name)
                   , "_kriging_4_vehicle_crime_confusion_metrics.csv")
          , row.names = F)

for(i in 0:4){
  for(j in 0:4){
    df <- quadrant_compare %>% 
      dplyr::filter(Actual == i & Predicted == j) %>% 
      dplyr::summarise(n= n())
    quadrant_compare_matrix[i+1,j+1] <- df$n
  }
}

quadrant_compare_df <- data.frame(quadrant_compare_matrix)
names(quadrant_compare_df) <- c("insignificant (I)","low-low (L-L)","low-high (L-H)","high-low (H-L)","high-high (H-H)")
rownames(quadrant_compare_df) <- c("insignificant (I)","low-low (L-L)","low-high (L-H)","high-low (H-L)","high-high (H-H)")
quadrant_compare_df_long <- quadrant_compare_df %>% 
  dplyr::mutate(Actual = rownames(quadrant_compare_df)) %>% 
  reshape2::melt(id.vars = "Actual", variable.name = "Predicted", value.name = "n") %>% 
  dplyr::mutate(x = ifelse(Actual == "insignificant (I)", 0
                           , ifelse(Actual == "low-low (L-L)", 0.25
                                    , ifelse(Actual == "low-high (L-H)", 0.5
                                             , ifelse(Actual == "high-low (H-L)", 0.75, 1))))) %>% 
  dplyr::mutate(y = ifelse(Predicted == "insignificant (I)", 1
                           , ifelse(Predicted == "low-low (L-L)", 0.75
                                    , ifelse(Predicted == "low-high (L-H)", 0.5
                                             , ifelse(Predicted == "high-low (H-L)", 0.25, 0)))))
# Get Bivariate Moran's I - Actual
crime_veh_res_bivariate_actual <- readRDS(paste0("Datasets/"
                                             , ds.name
                                             , "/esda/esda_"
                                             , tolower(ds.name)
                                             , "_crime_veh_res_bivmoran_obj.RDS"))
# Calculate Bivariate Moran's I - Predicted
crime_veh_res_bivariate_pred <-  fn_bivariatemoran(x = crime_veh_kriged$var1.pred
                                           , y = poly.sp@data$crime_res_1000hh_use
                                           , nb = neighbour$nb)
cor(crime_veh_res_bivariate_actual$local, crime_veh_res_bivariate_pred$local, use = "complete.obs")

# Calculate difference
crime_veh_bivariate_quadrant_actual <- readRDS(paste0("Datasets/"
               , ds.name
               , "/esda/esda_", tolower(ds.name) , "_vehicle_res_crime_bivariate_quadrant.RDS"))
crime_veh_res_bivariate_quadrant_pred <- fn_bivariatequadrant(x = crime_veh_kriged$var1.pred
                                                  , y = poly.sp@data$crime_res_1000hh_use
                                                  , sig = crime_veh_res_bivariate_pred$local_sig
                                                  , W = crime_veh_res_bivariate_pred$W)


quadrant_bivariate_compare <- data.frame(Actual = crime_veh_bivariate_quadrant_actual$quadrant
                               , Predicted = crime_veh_res_bivariate_quadrant_pred$quadrant)

quadrant_bivariate_compare_matrix <- matrix(nrow = 5, ncol = 5)

crime_veh_res_bivariate_confusion <- caret::confusionMatrix(reference = as.factor(quadrant_bivariate_compare$Actual)
                                                            , data = as.factor(quadrant_bivariate_compare$Predicted)
                                                            , dnn = c("Actual", "Predicted"))

crime_veh_res_bivariate_confusion_df <- data.frame(t(crime_veh_res_bivariate_confusion$byClass)) %>% 
  dplyr::mutate(`Metric` = rownames(crime_veh_res_bivariate_confusion_df))
names(crime_veh_res_bivariate_confusion_df) <- c("insignificant (I)","low-low (L-L)","low-high (L-H)","high-low (H-L)","high-high (H-H)", "Metric")
write.csv(crime_veh_res_bivariate_confusion_df
          , paste0("Datasets/"
                   , ds.name
                   , "/"
                   , section
                   , "/"
                   , section
                   , "_"
                   , tolower(ds.name)
                   , "_kriging_4_veh_res_bivariate_confusion_metrics.csv")
                   , row.names = F)

for(i in 0:4){
  for(j in 0:4){
    df <- quadrant_bivariate_compare %>% 
      dplyr::filter(Actual == i & Predicted == j) %>% 
      dplyr::summarise(n= n())
    quadrant_bivariate_compare_matrix[i+1,j+1] <- df$n
  }
}
quadrant_bivariate_compare_df <- data.frame(quadrant_bivariate_compare_matrix)
names(quadrant_bivariate_compare_df) <- c("insignificant (I)","low-low (L-L)","low-high (L-H)","high-low (H-L)","high-high (H-H)")
rownames(quadrant_bivariate_compare_df) <- c("insignificant (I)","low-low (L-L)","low-high (L-H)","high-low (H-L)","high-high (H-H)")
quadrant_bivariate_compare_df_long <- quadrant_bivariate_compare_df %>% 
  dplyr::mutate(Actual = rownames(quadrant_bivariate_compare_df)) %>% 
  reshape2::melt(id.vars = "Actual", variable.name = "Predicted", value.name = "n") %>% 
  dplyr::mutate(x = ifelse(Actual == "insignificant (I)", 0
                           , ifelse(Actual == "low-low (L-L)", 0.25
                                    , ifelse(Actual == "low-high (L-H)", 0.5
                                             , ifelse(Actual == "high-low (H-L)", 0.75, 1))))) %>% 
  dplyr::mutate(y = ifelse(Predicted == "insignificant (I)", 1
                           , ifelse(Predicted == "low-low (L-L)", 0.75
                                    , ifelse(Predicted == "low-high (L-H)", 0.5
                                             , ifelse(Predicted == "high-low (H-L)", 0.25, 0)))))



png(paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name), "_kriging_4_vehicle_crime_local_bivariate_corr_plot.png"), width = 750, height = 500)
par(mfrow=c(2,2), mar = rep(4,4), bg = "grey80")
plot(x = crime_veh_local_actual$Ii
     , y = crime_veh_local_pred$Ii
     , pch = 20
     , main = "(a), Local Moran's I:\n Actual vs. Predicted"
     , xlab = "Actual"
     , ylab = "Predicted")
text(x = 2
     , y = 11
     , paste("Corr =", round(cor(crime_veh_local_actual$Ii, crime_veh_local_pred$Ii, use = "complete.obs"),4)))
image(t(quadrant_compare_matrix[nrow(quadrant_compare_matrix):1,])
      , col = colorRampPalette(brewer.pal(8, "YlOrRd"))(10000)
      , xlab = "Actual"
      , xaxt = "n"
      , ylab = "Predicted"
      , yaxt = "n"
      , main = "(b), Local LISA Clusters: \nConfusion Matrix")
axis(1
     , at = unique(quadrant_compare_df_long$x)
     , labels = c("I","L-L","L-H","H-L","H-H")
     , side = 1
     , )
axis(2
     , at = unique(quadrant_compare_df_long$y)
     , labels = c("I","L-L","L-H","H-L","H-H")
     , )
text(x = quadrant_compare_df_long$x
     , y = quadrant_compare_df_long$y
     ,  quadrant_compare_df_long$n)

plot(x = crime_veh_res_bivariate_actual$local
     , y = crime_veh_res_bivariate_pred$local
     , pch = 20
     , main = "(c), Bivariate Moran's I:\n Actual vs. Predicted"
     , xlab = "Actual"
     , ylab = "Predicted")
text(x = 4.5
     , y = -1.5
     , paste("Corr =", round(cor(crime_veh_res_bivariate_actual$local, crime_veh_res_bivariate_pred$local, use = "complete.obs"),4)))

image(t(quadrant_bivariate_compare_matrix[nrow(quadrant_bivariate_compare_matrix):1,])
      , col = colorRampPalette(brewer.pal(8, "YlOrRd"))(10000)
      , xlab = "Actual"
      , xaxt = "n"
      , ylab = "Predicted"
      , yaxt = "n"
      , main = "(d), Bivariate LISA Clusters: \nConfusion Matrix")
axis(1
     , at = unique(quadrant_bivariate_compare_df_long$x)
     , labels = c("I","L-L","L-H","H-L","H-H")
     , side = 1
     , )
axis(2
     , at = unique(quadrant_bivariate_compare_df_long$y)
     , labels = c("I","L-L","L-H","H-L","H-H")
     , )
text(x = quadrant_bivariate_compare_df_long$x
     , y = quadrant_bivariate_compare_df_long$y
     ,  quadrant_bivariate_compare_df_long$n)

dev.off()