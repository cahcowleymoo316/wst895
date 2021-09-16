setwd("H:/Personal/Data Science/MSc/WST 895")
library(RODBC)
library(sp)
library(gstat)
library(tidyverse)
library(rgeos)
library(tmap)
library(spdep)
library(data.table)

# 1. Load data #
ds.name <- "Gauteng_Crime"

#subs.df0 can be made available upon request

subs.df.estimate <- subs.df0[(grepl("MULTIPOLYGON", subs.df0$sub_wkt)),]

subs.df <- subs.df0[!(grepl("MULTIPOLYGON", subs.df0$sub_wkt)),]

# Create SpatialPointsDataFrame for Model cross validation
# point.sp <- subs.df %>% 
#   dplyr::select(suburb_id, suburb, crime_res_1000hh_use, crime_vehicle_1000hh_use, x, y)
# coordinates(point.sp) = ~x+y
# saveRDS(point.sp, paste0("Datasets/", ds.name, "/", ds.name, "_Point_Data.RDS"))
# 2. Convert data frame to SpatialPolygonsDataFrame
poly.sp <- sp::SpatialPolygonsDataFrame(readWKT(subs.df$sub_wkt[1]), data=data.frame(OBJECTID=subs.df$suburb_id[1]
                                                                                     , suburb = subs.df$suburb[1]
                                                                                     , municipality = subs.df$municipality[1]
                                                                                     , adults = subs.df$adults[1]
                                                                                     , crime_vehicle_1000hh_use=subs.df$crime_vehicle_1000hh[1]
                                                                                     , crime_res_1000hh_use=subs.df$crime_res_1000hh[1]))

for (n in 2:length(subs.df$suburb_id)) {
  poly.sp <- rbind(poly.sp,
                   sp::SpatialPolygonsDataFrame(readWKT(subs.df$sub_wkt[n]), data=data.frame(OBJECTID=subs.df$suburb_id[n]
                                                                                                  , suburb = subs.df$suburb[n]
                                                                                                  , municipality = subs.df$municipality[n]
                                                                                                  , adults = subs.df$adults[n]
                                                                                                  , crime_vehicle_1000hh_use=subs.df$crime_vehicle_1000hh[n]
                                                                                                  , crime_res_1000hh_use=subs.df$crime_res_1000hh[n])))
}
saveRDS(poly.sp, file = paste0("Datasets/", ds.name, "/", ds.name, "_Data.RDS"))
section <- "esda"
# 3. Plot variables of interest:
# 3a. Create colours for the variables
crime_veh_colours <- fn_mapcolours(var =poly.sp@data$crime_vehicle_1000hh_use)
crime_res_colours <- fn_mapcolours(var = poly.sp@data$crime_res_1000hh_use, brewerpal = "Blues")
# 3b. Set plotting parameters
# Uncomment png and dev.off to save plot
png(paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name) , "_1_variable_plots.png"), width = 750, height = 500)
par(mfrow = c(2,1), mar = c(1,1,1,1), bg = "Grey80")
# 3c. Plots
fn_spdfplot(poly.sp
            , colours = crime_veh_colours
            , main = "(a), Vehicle crimes per 1000 households in Gauteng"
            , legend.pos = "topleft"
            , legend.title = "Vehicle crimes per 1000 households")
fn_spdfplot(poly.sp
            , colours = crime_res_colours
            , main = "(a), Residential crimes per 1000 households in Gauteng"
            , legend.pos = "topleft"
            , legend.title = "Residential crimes per 1000 households")
dev.off()
# 3d. Save plot

# 4. Define neighbourhoods and their weights
neighbour <- fn_neighbours(spdf = poly.sp, queen = T, zero.policy = T)
neighbour$lw$weights
saveRDS(neighbour, paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name), "_neighbours.RDS"))

# # If you get issues related to orphaned holes and geometry validity issues, run the following code to clean
# # your data. https://gis.stackexchange.com/questions/113964/fixing-orphaned-holes-in-r
# require(devtools)
# install_github("eblondel/cleangeo")
# require(cleangeo)
# # Get a report of geometry validity & issues for a sp spatial object
# report <- clgeo_CollectionReport(poly.sp)
# # Check how many problems
# report %>% dplyr::group_by(type, issue_type) %>% summarise(n = n())
# # If you want to remove the problem polygons, run the following line and re-run neighbour assignment
# poly.sp <- poly.sp[which(is.na(report$type)),]


# 5. Plot an example spatial entity and its neighbours.



# 6. Create Moran's Scatterplot
# 6a. Define variables to use
colsuse <- c("crime_vehicle_1000hh_use", "crime_res_1000hh_use")
# 6b. Run functions to create plots
crime_veh_moranscatter <- fn_moranscatter(spdf = poly.sp
                                      , lw = neighbour$lw
                                      , var = "crime_vehicle_1000hh_use"
                                      , title = "Moran's Scatterplot of Vehicle \nCrime per 1000 households"
                                      , col = crime_veh_colours$uniquecolour[length(crime_veh_colours$uniquecolour)]
                                      , zero.policy = F)


crime_res_moranscatter <- fn_moranscatter(spdf = poly.sp
                                          , lw = neighbour$lw
                                          , var = "crime_res_1000hh_use"
                                          , title = "Moran's Scatterplot of Residential \nCrime per 1000 households"
                                          , col = crime_res_colours$uniquecolour[length(crime_res_colours$uniquecolour)]
                                          , zero.policy = F)

# 6c. View plots
png(paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name) , "_2_global_moran_scatter.png"), width = 750, height = 500)
par(mfrow=c(1,2), mar = c(5,5,5,5), bg = "Grey80")

crime_veh_moranscatter$moranscatter()
crime_res_moranscatter$moranscatter()

dev.off()

# 7. Calculate Global Moran's I
# 7a. Calculate value
crime_veh_moran <- spdep::moran(x = poly.sp$crime_vehicle_1000hh_use
                                , listw = neighbour$lw
                                , n = 1000
                                , S0 = Szero(neighbour$lw)
                                , zero.policy = T)
crime_veh_moran$I
# cr
crime_res_moran <- spdep::moran(x = poly.sp$crime_res_1000hh_use
                                , listw = neighbour$lw
                                , n = 1000
                                , S0 = Szero(neighbour$lw)
                                , zero.policy = T)
crime_res_moran$I

# 7b. Combine all Moran's I's

global_moran_i <- data.frame(Variable = colsuse
                             , `Global Moran's I` = rbind(crime_veh_moran$I, crime_res_moran$I))
write.csv(global_moran_i
          , paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name) , "_3_global_moran_i.csv")
          , row.names = F)

# 7b. Check significance with Monte Carlo Method
crime_veh_mc <- spdep::moran.mc(x = poly.sp$crime_vehicle_1000hh_use
                                , listw = neighbour$lw
                                , nsim = 1000
                                , zero.policy = T)
crime_res_mc <- spdep::moran.mc(x = poly.sp$crime_res_1000hh_use
                                , listw = neighbour$lw
                                , nsim = 1000
                                , zero.policy = T)
global_moran_i_mc <- data.frame(Variable = colsuse
                                , Statistic = rbind(crime_veh_mc$statistic, crime_res_mc$statistic)
                                , `p-value` = rbind(crime_veh_mc$p.value, crime_res_mc$p.value))
write.csv(global_moran_i_mc
          , paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name) , "_4_global_moran_i_mc.csv")
          , row.names = F)


# 8. Calculate Local Moran's I
# 8a. Calculate and determine colours
crime_veh_local_moran <- fn_localmoran(spdf = poly.sp
                                       , lw = neighbour$lw
                                       , var = "crime_vehicle_1000hh_use"
                                       , zero.policy = F)
saveRDS(crime_veh_local_moran, paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name) , "_vehicle_crime_local_obj.RDS"))
crime_veh_local_moran_colours <- fn_mapcolours(var = crime_veh_local_moran$Ii
                                               , num_colours = 5
                                               , brewerpal = "YlGn"
                                               , dig.lab = 3)
crime_res_local_moran <- fn_localmoran(spdf = poly.sp
                                       , lw = neighbour$lw
                                       , var = "crime_res_1000hh_use"
                                       , zero.policy = F)
crime_res_local_moran_colours <- fn_mapcolours(var = crime_res_local_moran$Ii
                                               , num_colours = 5
                                               , brewerpal = "PuBu"
                                               , dig.lab = 3)
# 8b. Plot values
png(paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name) , "_5_local_moran_plots.png"), width = 1200, height = 900)
par(mfrow = c(2,1), mar = c(1,1,1,1), bg = "Grey80")
fn_spdfplot(poly.sp
            , colours = crime_veh_local_moran_colours
            , main = "(a), Local Moran's I of Vehicle crimes per 1000 households in Gauteng"
            , legend.pos = "topleft"
            , legend.title = "Vehicle crimes per 1000 households")
fn_spdfplot(poly.sp
            , colours = crime_res_local_moran_colours
            , main = "(a), Local Moran's I of Vehicle crimes per 1000 households in Gauteng"
            , legend.pos = "topleft"
            , legend.title = "Vehicle crimes per 1000 households")
dev.off()

# 8c. Determine cluster quadrants
crime_veh_local_quadrant <- fn_quadrant(spdf = crime_veh_local_moran
                                        , lw = neighbour$lw
                                        , var = "crime_vehicle_1000hh_use")
saveRDS(crime_veh_local_quadrant, paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name) , "_vehicle_crime_local_quadrant.RDS"))
table(crime_veh_local_quadrant$quadrant)
crime_res_local_quadrant <- fn_quadrant(spdf = crime_res_local_moran
                                        , lw = neighbour$lw
                                        , var = "crime_res_1000hh_use")
table(crime_res_local_quadrant$quadrant)

# 8d. Plot LISA cluster quadrants and proportion of spatial entities that fall in each cluster

png(paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name) , "_6_vehicle_crime_LISA.png"), width = 750, height = 500)
fn_plotLISAquad(spdf = poly.sp
                , var = "crime_vehicle_1000hh_use"
                , quadrant = crime_veh_local_quadrant
                , legend.pos = "bottomright"
                , graph.legend = 1
                , margin = 5
                , varname.print = "Vehicle Crime")
dev.off()

png(paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name) , "_7_res_crime_LISA.png"), width = 750, height = 500)
fn_plotLISAquad(spdf = poly.sp
                , var = "crime_res_1000hh_use"
                , quadrant = crime_res_local_quadrant
                , legend.pos = "bottomright"
                , graph.legend = 1
                , margin = 5
                , varname.print = "Residential Crime")
dev.off()

# 8e. Plot LISA clusters
png(paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name), "_8_LISA_Clusters.png"), width = 750, height = 500)

par(mfrow = c(1,2), mar = c(2,2,2,2), bg = "Grey80")
fn_LISAclustermap(spdf = poly.sp
                  , title = "(a), Vehicle Crime"
                  , colors = crime_veh_local_moran_colours
                  , quadrant = crime_veh_local_quadrant)

fn_LISAclustermap(spdf = poly.sp
                  , title = "(b), Residential Crime"
                  , colors = crime_res_local_moran_colours
                  , quadrant = crime_res_local_quadrant)
dev.off()

# 9. Bivariate Moran's I
# 9a. Calculate bivariate Global and Local Moran's I 
crime_veh_res_bivmoran <- fn_bivariatemoran(x = poly.sp@data$crime_vehicle_1000hh_use
                                            , y = poly.sp@data$crime_res_1000hh_use
                                            , nb = neighbour$nb)
saveRDS(crime_veh_res_bivmoran, paste0("Datasets/"
       , ds.name
       , "/esda/esda_"
       , tolower(ds.name)
       , "_crime_veh_res_bivmoran_obj.RDS"))
crime_veh_res_bivmoran_colours <- fn_mapcolours(var = crime_veh_res_bivmoran$local
                                                , num_colours = 8
                                                , brewerpal = "Reds"
                                                , dig.lab = 0)

# 9b. Plot bivariate Local Moran's I
png(paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name), "_9_vehicle_and_res_crime_bivariate_local_moran.png"), width = 750, height = 500)
par(mfrow = c(1,1), mar = c(2,2,2,2), bg = "Grey80")

fn_spdfplot(spdf = poly.sp
            , colours = crime_veh_res_bivmoran_colours
            , main = "Bivariate Local Moran's I values for \nVehicle and Residential Crimes per 1000 households"
            , bordercol = "transparent"
            , legend.pos = "topleft"
            , legend.title = "Bivariate Local Moran's I"
            )
dev.off()


#9c. Determine cluster quadrants
crime_veh_res_bivquadrant <- fn_bivariatequadrant(x = poly.sp@data$crime_vehicle_1000hh_use
                                                  , y = poly.sp@data$crime_res_1000hh_use
                                                  , sig = crime_veh_res_bivmoran$local_sig
                                                  , W = crime_veh_res_bivmoran$W)
saveRDS(crime_veh_res_bivquadrant
        , paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name) , "_vehicle_res_crime_bivariate_quadrant.RDS"))

# 9d. Plot LISA clusters and proportion of spatial entities that fall in each cluster
png(paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name), "_10_vehicle_and_res_crime_bivariate_LISA.png"), width = 750, height = 500)
fn_plotbivariateLISAquad(quadrant = crime_veh_res_bivquadrant)
dev.off()
# 9e. Plot spatial LISA clusters
png(paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name), "_10_vehicle_and_res_crime_bivariate_LISA_Clusters.png"), width = 750, height = 500)

par(mfrow = c(1,1), mar = c(2,2,2,2), bg = "Grey80")

fn_LISAclustermap(spdf = poly.sp
                  , title = "Spatial bivariate LISA clusters \n Vehicle and Residential Crime"
                  , colors = crime_veh_res_bivmoran_colours
                  , quadrant = crime_veh_res_bivquadrant)
legend("topleft"
       , title = "Cluster"
       , legend = c("insignificant (I)","low-low (L-L)","low-high (L-H)","high-low (H-L)","high-high (H-H)")
       , fill=c("grey50","blue",rgb(0,0,1,alpha=0.4),rgb(1,0,0,alpha=0.4),"red"))
dev.off()

# 10. Calculate dispersion

var(poly.sp@data$crime_vehicle_1000hh_use, na.rm = T)/mean(poly.sp@data$crime_vehicle_1000hh_use, na.rm = T)
var(poly.sp@data$crime_res_1000hh_use, na.rm = T)/mean(poly.sp@data$crime_res_1000hh_use, na.rm = T)




spdf <- poly.sp
vars <- c("crime_vehicle_1000hh_use", "crime_res_1000hh_use")
df_dispersion <- fn_dispersion(poly.sp, c("crime_vehicle_1000hh_use", "crime_res_1000hh_use"))


# 11. Summary statistics
summary(df_dispersion$dispersion)



summary(crime_res_moran$I)
summary(crime_veh_moran$I)

summary(crime_res_local_moran$Ii)
# for 05/07/2021"
# get permission and send ethics 
# 
# write story for application chapter
