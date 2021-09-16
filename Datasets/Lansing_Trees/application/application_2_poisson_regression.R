setwd("H:/Personal/Data Science/MSc/WST 895")
ds.name <- "Lansing_Trees"
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
### BASELINE MODEL ###
# Poisson regression
set.seed(103207)
maple_poiss <- glm(maple ~ hickory, poly.sp@data, family = "poisson")
#saveRDS(maple_poiss, paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name), "_poisson_maple_glm.RDS"))
maple_poiss_cv <- boot::cv.glm(data = poly.sp@data
                                   , glmfit = maple_poiss
                                   , K = 10)
maple_poiss_cv$delta



# Residual analysis
# Translate residual to original units
maple_poiss$data$eresidual <- maple_poiss$data$maple - maple_poiss$fitted.values
maple_poiss$data$zscore <- scale(maple_poiss$residuals)
mean(maple_poiss$data$eresidual^2)


png(paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name), "_poisson_1_maple_res_plots.png"), width = 750, height = 500)
par(mfrow = c(2,2), mar = rep(4,4), bg = "grey80")
hist(maple_poiss$data$eresidual
     , xlab = "Residual"
     , main = "(a), Histogram of Residuals")
text(x = 6
     ,y = 100
     , paste("MSE =\n", round(mean(maple_poiss$data$eresidual^2), 4)))
# Mean Square Normalised Error
hist(maple_poiss$data$zscore
     , xlab = "Normalised Residual"
     , main = "(b), Histogram of \nNormalised Residuals")
text(x = 3
     ,y = 80
     , paste("MSE =", round(mean(maple_poiss$data$zscore^2), 4)))
# Correlation observed and predicted, ideally close to 1
plot(x = maple_poiss$data$maple
     , y = maple_poiss$fitted.values
     #, col = maple_poiss$fold
     , pch = 20
     , main = "Observed vs. Predicted"
     , xlab = "Observed"
     , ylab = "Predicted")
text(x = 8, y = 0.5
     , paste("Corr =", round(cor(maple_poiss$data$maple, maple_poiss$fitted.values),4)))
# Correlation predicted and residual, ideally 0
plot(x = maple_poiss$fitted.values
     , y = maple_poiss$data$eresidual
     #, col = maple_poiss$fold
     , pch = 20
     , main = "Predicted vs. Residual"
     , xlab = "Precicted"
     , ylab = "Residual")
text(x = 0.75, y = 8
     , paste("Corr =\n", round(cor(maple_poiss$fitted.values, maple_poiss$data$eresidual),4)))
dev.off()




maple_poiss_est_colours <- fn_mapcolours(var = maple_poiss$fitted.values, brewerpal = "YlGn", num_colours = 5, dig.lab = 2)

# Fitted values
maple_poiss_pred_bins <- data.frame(class = cut(maple_poiss$fitted.values, breaks = seq(from = 0, to = 5, by = 1))) %>% 
        dplyr::group_by(class) %>% 
        dplyr::summarise(n = n()) %>% 
        dplyr::mutate(midpoint = seq(0.5, 3.5, by = 1))

png(paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name), "_poisson_2_maple_model_fit_plots.png"), width = 750, height = 500)
par(mfrow = c(2,2), mar = rep(4,4), bg = "Grey80")
hist(poly.sp$maple
     , xlab = "Maple"
     , ylim = c(0, 200)
     , main = "(a), Histogram of \nActual vs. Fitted Values")
lines(x = maple_poiss_pred_bins$midpoint
      , y = maple_poiss_pred_bins$n
      , col = "red")
text(x = 5.5
     , y = 50
     , "Fitted"
     , col = "red")
plot(x = poly.sp@data$hickory
     , y = poly.sp@data$maple
     , pch = 20
     , ylab = "Maple"
     , xlab = "Hickory"
     , main = "(b), Actual vs. Fitted")
points(x = poly.sp@data$hickory
       , y = maple_poiss$fitted.values
       , col = scales::alpha("red", 0.3)
       , pch = 19)
fn_spdfplot(spdf = poly.sp
            , colours = maple_poiss_est_colours
            , main = "(c), Poisson Regression \nestimates")
plot(0,type='n',axes=FALSE,ann=FALSE)
legend("top"
       , title = "Estimate"
       , legend = levels(maple_poiss_est_colours$class)
       , fill=maple_poiss_est_colours$uniquecolour)
dev.off()

# Correlation Plots

# Calculate Global Moran's I - Actual
maple_global_actual <- read.csv(paste0("Datasets/", ds.name, "/esda/esda_", tolower(ds.name) , "_3_global_moran_i.csv")) %>% 
        dplyr::filter(Variable == "maple")
# Calculate Global Moran's I - Predicted
neighbour <- readRDS(paste0("Datasets/", ds.name, "/esda/esda", "_", tolower(ds.name), "_neighbours.RDS"))
maple_global_pred <- spdep::moran.mc(x = maple_poiss$fitted.values
                                     , listw = neighbour$lw
                                     , nsim = 1000
                                     , zero.policy = T)
# Calculate difference
maple_global_diff <- abs(maple_global_actual$statistic - maple_global_pred$statistic)
maple_global_diff_df <- maple_global_actual %>% 
        dplyr::mutate(Predicted = maple_global_pred$statistic) %>% 
        dplyr::rename(Actual = statistic) %>% 
        reshape2::melt(id.vars = "Variable") %>% 
        dplyr::select(-Variable)
png(paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name), "_poisson_3_maple_global_corr_plot.png"), width = 375, height = 250)
par(mfrow=c(1,1), mar = rep(4,4), bg = "grey80")
barplot(maple_global_diff_df$value
        , names.arg = paste0(maple_global_diff_df$variable, "\n(", round(maple_global_diff_df$value,4),")")
        , col = c("black", "red"), 0.4
        , main = "Global Moran's I:\n Actual vs.Predicted")
dev.off()

# Calculate Local Moran's I - Actual
maple_local_actual <- readRDS(paste0("Datasets/", ds.name, "/esda/esda_", tolower(ds.name) , "_maple_local_obj.RDS"))
# Calculate Local Moran's I - Predicted
poly.sp@data$poiss_pred <- maple_poiss$fitted.values
maple_local_pred <-  fn_localmoran(spdf = poly.sp
                                       , lw = neighbour$lw
                                       , var = "poiss_pred"
                                       , zero.policy = F)



# Get  LISA clusters from predicted Moran's I
maple_local_quadrant_pred <- fn_quadrant(spdf = maple_local_pred
                                             , lw = neighbour$lw
                                             , var = "poiss_pred")
table(maple_local_quadrant_pred$quadrant)
quadrant_compare <- data.frame(Actual = maple_local_quadrant$quadrant
                               , Predicted = maple_local_quadrant_pred$quadrant)

quadrant_compare_matrix <- matrix(nrow = 5, ncol = 5)

maple_confusion <- caret::confusionMatrix(reference = as.factor(quadrant_compare$Actual)
                                              , data = as.factor(quadrant_compare$Predicted)
                                              , dnn = c("Predicted", "Actual"))
maple_confusion_df <- data.frame(t(maple_confusion$byClass)) %>% 
        dplyr::mutate(`Metric` = NA)
names(maple_confusion_df) <- c("insignificant (I)","low-low (L-L)","high-high (H-H)", "Metric")
maple_confusion_df$Metric <- rownames(maple_confusion_df)
write.csv(maple_confusion_df
          , paste0("Datasets/"
                   , ds.name
                   , "/"
                   , section
                   , "/"
                   , section
                   , "_"
                   , tolower(ds.name)
                   , "_poisson_4_maple_confusion_metrics.csv")
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


png(paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name), "_poisson_4_maple_local_corr_plot.png"), width = 750, height = 250)
par(mfrow=c(1,2), mar = rep(4,4), bg = "grey80")
plot(x = maple_local_actual$Ii
     , y = maple_local_pred$Ii
     , pch = 20
     , main = "(a), Local Moran's I:\n Actual vs. Predicted"
     , xlab = "Actual"
     , ylab = "Predicted")
text(x = 4
     , y = -1
     , paste("Corr =", round(cor(maple_local_actual$Ii, maple_local_pred$Ii, use = "complete.obs"),4)))
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
dev.off()
