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
library(boot)
# Load Functions
source("Datasets/esda_functions.R")
# Load Data
poly.sp <- readRDS(paste0("Datasets/", ds.name, "/", ds.name, "_Data.RDS"))
point.sp <- readRDS(paste0("Datasets/", ds.name, "/", ds.name, "_Point_Data.RDS"))
### BASELINE MODEL ###
# Poisson regression
set.seed(103207)
crime_veh_poiss <- glm(crime_vehicle_1000hh_use ~ crime_res_1000hh_use, poly.sp@data, family = "poisson")
crime_veh_poiss_cv <- boot::cv.glm(data = poly.sp@data
                                   , glmfit = crime_veh_poiss
                                   , K = 10)
crime_veh_poiss_cv$delta



# Residual analysis
# Translate residual to original units
crime_veh_poiss$data$eresidual <- crime_veh_poiss$data$crime_vehicle_1000hh_use - crime_veh_poiss$fitted.values
crime_veh_poiss$data$zscore <- scale(crime_veh_poiss$residuals)
mean(crime_veh_poiss$data$eresidual^2)


png(paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name), "_poisson_1_vehicle_crime_res_plots.png"), width = 750, height = 500)
par(mfrow = c(2,2), mar = rep(4,4), bg = "grey80")
hist(crime_veh_poiss$data$eresidual
     , xlab = "Residual"
     , main = "(a), Histogram of Residuals")
text(x = -30
     ,y = 550
     , paste("MSE =\n", round(mean(crime_veh_poiss$data$eresidual^2), 4)))
# Mean Square Normalised Error
hist(crime_veh_poiss$data$zscore
     , xlab = "Normalised Residual"
     , main = "(b), Histogram of \nNormalised Residuals")
text(x = 3
     ,y = 300
     , paste("MSE =", round(mean(crime_veh_poiss$data$zscore^2), 4)))
# Correlation observed and predicted, ideally close to 1
plot(x = crime_veh_poiss$data$crime_vehicle_1000hh_use
     , y = crime_veh_poiss$fitted.values
     #, col = crime_veh_poiss$fold
     , pch = 20
     , main = "Observed vs. Predicted"
     , xlab = "Observed"
     , ylab = "Predicted")
text(x = 45, y = 60
     , paste("Corr =", round(cor(crime_veh_poiss$data$crime_vehicle_1000hh_use, crime_veh_poiss$fitted.values),4)))
# Correlation predicted and residual, ideally 0
plot(x = crime_veh_poiss$fitted.values
     , y = crime_veh_poiss$data$eresidual
     #, col = crime_veh_poiss$fold
     , pch = 20
     , main = "Predicted vs. Residual"
     , xlab = "Precicted"
     , ylab = "Residual")
text(x = 50, y = 35
     , paste("Corr =", round(cor(crime_veh_poiss$fitted.values, crime_veh_poiss$data$eresidual),4)))
dev.off()




crime_veh_poiss_est_colours <- fn_mapcolours(var = crime_veh_poiss$fitted.values, brewerpal = "YlOrRd", num_colours = 5, dig.lab = 2)

# Fitted values
crime_veh_poiss_pred_bins <- data.frame(class = cut(crime_veh_poiss$fitted.values, breaks = seq(from = 0, to = 70, by = 5))) %>% 
        dplyr::group_by(class) %>% 
        dplyr::summarise(n = n()) %>% 
        dplyr::mutate(midpoint = c(seq(7.5, 57.5, by = 5), 67.5))

png(paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name), "_poisson_2_vehicle_crime_model_fit_plots.png"), width = 750, height = 500)
par(mfrow = c(2,2), mar = rep(4,4), bg = "Grey80")
hist(poly.sp$crime_vehicle_1000hh_use
     , xlab = "Vehicle Crime"
     , ylim = c(0, 850)
     , main = "(a), Histogram of \nActual vs. Fitted Values")
lines(x = crime_veh_poiss_pred_bins$midpoint
      , y = crime_veh_poiss_pred_bins$n
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
       , y = crime_veh_poiss$fitted.values
       , col = scales::alpha("red", 0.3)
       , pch = 19)
fn_spdfplot(spdf = poly.sp
            , colours = crime_veh_poiss_est_colours
            , main = "(c), Poisson Regression \nestimates")
plot(0,type='n',axes=FALSE,ann=FALSE)
legend("top"
       , title = "Estimate"
       , legend = levels(crime_veh_poiss_est_colours$class)
       , fill=crime_veh_poiss_est_colours$uniquecolour)
dev.off()

# Correlation Plots

# Calculate Global Moran's I - Actual
crime_veh_global_actual <- read.csv(paste0("Datasets/", ds.name, "/esda/esda_", tolower(ds.name) , "_3_global_moran_i.csv")) %>% 
        dplyr::filter(Variable == "crime_vehicle_1000hh_use")
# Calculate Global Moran's I - Predicted
neighbour <- readRDS(paste0("Datasets/", ds.name, "/esda/esda", "_", tolower(ds.name), "_neighbours.RDS"))
crime_veh_global_pred <- spdep::moran(x = crime_veh_poiss$fitted.values
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
png(paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name), "_poisson_3_vehicle_crime_global_corr_plot.png"), width = 750, height = 500)
par(mfrow=c(1,1), mar = rep(4,4), bg = "grey80")
barplot(crime_veh_global_diff_df$value
        , names.arg = paste0(crime_veh_global_diff_df$variable, "\n(", round(crime_veh_global_diff_df$value,4),")")
        , col = c("black", "red"), 0.4
        , main = "Global Moran's I:\n Actual vs.Predicted")
dev.off()

# Calculate Local Moran's I - Actual
crime_veh_local_actual <- readRDS(paste0("Datasets/", ds.name, "/esda/esda_", tolower(ds.name) , "_vehicle_crime_local_obj.RDS"))
# Calculate Local Moran's I - Predicted
poly.sp@data$poiss_pred <- crime_veh_poiss$fitted.values
crime_veh_local_pred <-  fn_localmoran(spdf = poly.sp
                                       , lw = neighbour$lw
                                       , var = "poiss_pred"
                                       , zero.policy = F)



# Get  LISA clusters from predicted Moran's I
crime_veh_local_quadrant_pred <- fn_quadrant(spdf = crime_veh_local_pred
                                             , lw = neighbour$lw
                                             , var = "poiss_pred")
table(crime_veh_local_quadrant_pred$quadrant)
quadrant_compare <- data.frame(Actual = crime_veh_local_quadrant$quadrant
                               , Predicted = crime_veh_local_quadrant_pred$quadrant)

quadrant_compare_matrix <- matrix(nrow = 5, ncol = 5)

crime_veh_confusion <- caret::confusionMatrix(reference = as.factor(quadrant_compare$Actual)
                                              , data = as.factor(quadrant_compare$Predicted)
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
                   , "_poisson_4_vehicle_crime_confusion_metrics.csv")
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
crime_veh_res_bivariate_pred <-  fn_bivariatemoran(x = poly.sp@data$poiss_pred
                                                   , y = poly.sp@data$crime_res_1000hh_use
                                                   , nb = neighbour$nb)

# Calculate difference
crime_veh_bivariate_quadrant_actual <- readRDS(paste0("Datasets/"
                                                      , ds.name
                                                      , "/esda/esda_", tolower(ds.name) , "_vehicle_res_crime_bivariate_quadrant.RDS"))
crime_veh_res_bivariate_quadrant_pred <- fn_bivariatequadrant(x = poly.sp@data$poiss_pred
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
                   , "_poisson_4_veh_res_bivariate_confusion_metrics.csv")
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



png(paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name), "_poisson_4_vehicle_crime_local_bivariate_corr_plot.png"), width = 750, height = 500)
par(mfrow=c(2,2), mar = rep(4,4), bg = "grey80")
plot(x = crime_veh_local_actual$Ii
     , y = crime_veh_local_pred$Ii
     , pch = 20
     , main = "(a), Local Moran's I:\n Actual vs. Predicted"
     , xlab = "Actual"
     , ylab = "Predicted")
text(x = 7
     , y = 40
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
text(x = 0
     , y = 25
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
