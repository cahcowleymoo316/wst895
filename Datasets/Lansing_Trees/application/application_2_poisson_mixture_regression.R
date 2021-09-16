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
library(flexmix)
# Load Functions
source("Datasets/esda_functions.R")
# Load Data
poly.sp <- readRDS(paste0("Datasets/", ds.name, "/", ds.name, "_Data.RDS"))

# Choosing the right number of mixtures
mix.dat <- NULL
for(k in 2:10){
        maple_mix <- flexmix::stepFlexmix(maple ~ hickory
                                              , k = k
                                              , nrep = 5
                                              , data = poly.sp@data
                                              , model = flexmix::FLXMRglm(family = "poisson"))
        mix.dat <- rbind(mix.dat, data.frame(k = k, AIC = AIC(maple_mix)))
        print(mix.dat)
}

mix.dat$col = ifelse(mix.dat$AIC == min(mix.dat$AIC), "red", "black")
png(paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name), "_mixture_0_maple_choosing_k.png"), width = 750, height = 500)
par(mfrow = c(1,1), bg = "Grey80")
barplot(height = mix.dat$AIC
        , col = mix.dat$col
        , names.arg = paste0(mix.dat$k, "\n(", round(mix.dat$AIC,1), ")")
        , xlab = "Number of Mixtures (AIC)"
        , ylab = "AIC"
        , main = "Choosing k - Maple")
dev.off()
# Poisson mixture regression
set.seed(12345)
maple_mix <- flexmix::stepFlexmix(maple ~ hickory
                                      , k = 2
                                      , nrep = 5
                                      , data = poly.sp@data
                                      , model = flexmix::FLXMRglm(family = "poisson"))
# saveRDS(maple_mix, paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name), "_mixture_maple_mix.RDS"))
maple_mix@components
summary(maple_mix)
clusters(maple_mix)
maple_mix_pred <- fitted(maple_mix, newdata = poly.sp@data)
maple_mix_pred_which <- data.frame(maple = poly.sp@data$maple
                                   , clusters = as.numeric(clusters(maple_mix))
                                   , Comp.1 = maple_mix_pred[,1]
                                   , Comp.2 = maple_mix_pred[,2]
                                 ) %>% 
        dplyr::mutate(estimate = ifelse(clusters == 1, Comp.1, 2))

# Residual analysis
# Translate residual to original units
maple_mix_eresidual <- poly.sp@data$maple - maple_mix_pred_which$estimate
maple_mix_zscore <- scale(maple_mix_eresidual)
mean(maple_mix_eresidual^2)


png(paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name), "_mixture_1_maple_res_plots.png"), width = 750, height = 500)
par(mfrow = c(2,2), mar = rep(4,4), bg = "grey80")
hist(maple_mix_eresidual
     , xlab = "Residual"
     , main = "(a), Histogram of Residuals")
text(x = 3
     ,y = 125
     , paste("MSE =\n", round(mean(maple_mix_eresidual^2), 4)))
# Mean Square Normalised Error
hist(maple_mix_zscore
     , xlab = "Normalised Residual"
     , main = "(b), Histogram of \nNormalised Residuals")
text(x = 3
     ,y = 75
     , paste("MSE =\n", round(mean(maple_mix_zscore^2), 4)))
# Correlation observed and predicted, ideally close to 1
plot(x = poly.sp@data$maple
     , y = maple_mix_pred_which$estimate
     #, col = maple_mix$fold
     , pch = 20
     , main = "Observed vs. Predicted"
     , xlab = "Observed"
     , ylab = "Predicted")
text(x = 10, y = 2
     , paste("Corr =\n", round(cor(poly.sp@data$maple, maple_mix_pred_which$estimate),4)))
# Correlation predicted and residual, ideally 0
plot(x = maple_mix_pred_which$estimate
     , y = maple_mix_eresidual
     #, col = maple_mix$fold
     , pch = 20
     , main = "Predicted vs. Residual"
     , xlab = "Precicted"
     , ylab = "Residual")
text(x = 2, y = 6
     , paste("Corr =\n", round(cor(maple_mix_pred_which$estimate, maple_mix_eresidual),4)))
dev.off()



maple_mix_est_colours <- fn_mapcolours(var = maple_mix_pred_which$estimate, brewerpal = "YlGn", num_colours = 5, dig.lab = 2)

# Fitted values
maple_mix_pred_bins <- data.frame(class = cut(maple_mix_pred_which$estimate, breaks = seq(from = 0, to = 7, by = 0.5))) %>% 
        dplyr::group_by(class) %>% 
        dplyr::summarise(n = n()) %>% 
        dplyr::mutate(midpoint = c(seq(0.75, 3.25, by = 0.5), 4.25, 5.25, 6.75))

png(paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name), "_mixture_2_maple_model_fit_plots.png"), width = 750, height = 500)
par(mfrow = c(2,2), mar = rep(3.9,4), bg = "Grey80")
hist(poly.sp$maple
     , xlab = "Maple"
     , ylim = c(0, 200)
     , main = "(a), Histogram of \nActual vs. Fitted Values")
lines(x = maple_mix_pred_bins$midpoint
      , y = maple_mix_pred_bins$n
      , col = "red")
text(x = 8
     , y = 20
     , "Fitted"
     , col = "red")


plot(x = poly.sp@data$hickory
     , y = poly.sp@data$maple
     , col = clusters(maple_mix)
     , pch = 20
     , ylab = "Maple"
     , xlab = "Hickory"
     , main = "(b), Actual vs. Fitted")
points(x = poly.sp@data$hickory
       , y = maple_mix_pred[,1]
       , col = "black"
       , pch = 19)
points(x = poly.sp@data$hickory
       , y = maple_mix_pred[,2]
       , col = "red"
       , pch = 19)

legend("topright"
       , legend = c("1", "2")
       , fill = c("black", "red")
       , inset = 0.01
       , bty = "\n")
fn_spdfplot(spdf = poly.sp
            , colours = maple_mix_est_colours
            , main = "(c), Poisson Mixture\n Regression estimates")
plot(0,type='n',axes=FALSE,ann=FALSE)
legend("top"
       , title = "Estimate"
       , legend = levels(maple_mix_est_colours$class)
       , fill=maple_mix_est_colours$uniquecolour)
dev.off()

# Correlation Plots

# Calculate Global Moran's I - Actual
maple_global_actual <- read.csv(paste0("Datasets/", ds.name, "/esda/esda_", tolower(ds.name) , "_3_global_moran_i.csv")) %>% 
        dplyr::filter(Variable == "maple")
# Calculate Global Moran's I - Predicted
neighbour <- readRDS(paste0("Datasets/", ds.name, "/esda/esda", "_", tolower(ds.name), "_neighbours.RDS"))
maple_global_pred <-spdep::moran.mc(x = maple_mix_pred_which$estimate
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
png(paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name), "_mixture_3_maple_global_corr_plot.png"), width = 750, height = 500)
par(mfrow=c(1,1), mar = rep(4,4), bg = "grey80")
barplot(maple_global_diff_df$value
        , names.arg = paste0(maple_global_diff_df$variable, "\n(", round(maple_global_diff_df$value,4),")")
        , col = c("black", "red"), 0.4
        , main = "Global Moran's I:\n Actual vs.Predicted")
dev.off()

# Calculate Local Moran's I - Actual
maple_local_actual <- readRDS(paste0("Datasets/", ds.name, "/esda/esda_", tolower(ds.name) , "_maple_local_obj.RDS"))
# Calculate Local Moran's I - Predicted
poly.sp@data$mix_pred <- maple_mix_pred_which$estimate
maple_local_pred <-  fn_localmoran(spdf = poly.sp
                                       , lw = neighbour$lw
                                       , var = "mix_pred"
                                       , zero.policy = F)



# Get  LISA clusters from predicted Moran's I
maple_local_quadrant_pred <- fn_quadrant(spdf = maple_local_pred
                                             , lw = neighbour$lw
                                             , var = "mix_pred")
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
                   , "_mixture_4_maple_confusion_metrics.csv")
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

png(paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name), "_mixture_4_maple_local_corr_plot.png"), width = 750, height = 250)
par(mfrow=c(1,2), mar = rep(4,4), bg = "grey80")
plot(x = maple_local_actual$Ii
     , y = maple_local_pred$Ii
     , pch = 20
     , main = "(a), Local Moran's I:\n Actual vs. Predicted"
     , xlab = "Actual"
     , ylab = "Predicted")
text(x = 4
     , y = -0.5
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
