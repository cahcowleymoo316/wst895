setwd("H:/Personal/Data Science/MSc/WST 895")
library(RODBC)
library(sp)
library(gstat)
library(tidyverse)
library(rgeos)
library(tmap)
library(spdep)
library(data.table)

# 0. Source functions
source("Datasets/esda_functions.R")



# 1. Load data #
ds.name <- "Lansing_Trees"
poly.sp <- gcKrig::LansingTrees

# 2. Convert data to SpatialPointsDataFrame
sp::coordinates(poly.sp) <- ~Easting+Northing
class(poly.sp)
saveRDS(poly.sp, paste0("Datasets/", ds.name, "/", ds.name, "_Data.RDS"))
section <- "esda"

# 3. Plot variables of interest:
# 3a. Create colours for the variables
maple_colours <- fn_mapcolours(var = poly.sp@data$maple, brewerpal = "YlGn", num_colours = 5, dig.lab = 0)
hickory_colours <- fn_mapcolours(var = poly.sp@data$hickory, brewerpal = "YlGn", num_colours = 5, dig.lab = 0)
blackoak_colours <- fn_mapcolours(var = poly.sp@data$blackoak, brewerpal = "YlGn", num_colours = 5, dig.lab = 0)
redoak_colours <- fn_mapcolours(var = poly.sp@data$redoak, brewerpal = "YlGn", num_colours = 5, dig.lab = 0)
whiteoak_colours <- fn_mapcolours(var = poly.sp@data$whiteoak, brewerpal = "YlGn", num_colours = 5, dig.lab = 0)
# 3b. Set plotting parameters
# Uncomment png and dev.off to save plot
# 3c. Plots
png(paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name) , "_1_variable_plots.png"), width = 750, height = 500)
par(mfrow = c(2,3), mar = c(5,5,5,5), bg = "Grey80")

fn_spdfplot(poly.sp
            , colours = maple_colours
            , main = "(a), Maple")
fn_spdfplot(poly.sp
            , colours = hickory_colours
            , main = "(b), Hickory"
            )
fn_spdfplot(poly.sp
            , colours = blackoak_colours
            , main = "(c), Black Oak"
            )
fn_spdfplot(poly.sp
            , colours = redoak_colours
            , main = "(d), Red Oak"
            )
fn_spdfplot(poly.sp
            , colours = whiteoak_colours
            , main = "(e), White Oak"
            )
dev.off()

# 4. Define neighbourhoods and their weights
neighbour <- fn_neighbours(spdf = poly.sp,k = 5, style = "W", zero.policy = T)
saveRDS(neighbour, paste0("Datasets/", ds.name, "/", section, "/",section, "_", tolower(ds.name), "_neighbours.RDS"))

neighbour$lw$weights

# 5. Plot an example spatial entity and its neighbours.



# 6. Create Moran's Scatterplot
# 6a. Define variables to use
colsuse <- c("maple", "hickory", "blackoak", "redoak", "whiteoak")
# 6b. Run functions to create plots
maple_moranscatter <- fn_moranscatter(spdf = poly.sp
                                          , lw = neighbour$lw
                                          , var = "maple"
                                          , title = "(a), Maple"
                                          , col = maple_colours$uniquecolour[length(maple_colours$uniquecolour)]
                                          , zero.policy = T)
hickory_moranscatter <- fn_moranscatter(spdf = poly.sp
                                          , lw = neighbour$lw
                                          , var = "hickory"
                                          , title = "(b), Hickory"
                                          , col = hickory_colours$uniquecolour[length(hickory_colours$uniquecolour)]
                                          , zero.policy = T)
blackoak_moranscatter <- fn_moranscatter(spdf = poly.sp
                                        , lw = neighbour$lw
                                        , var = "blackoak"
                                        , title = "(c), Black Oak"
                                        , col = blackoak_colours$uniquecolour[length(blackoak_colours$uniquecolour)]
                                        , zero.policy = T)
redoak_moranscatter <- fn_moranscatter(spdf = poly.sp
                                         , lw = neighbour$lw
                                         , var = "redoak"
                                         , title = "(d), Red Oak"
                                         , col = redoak_colours$uniquecolour[length(redoak_colours$uniquecolour)]
                                         , zero.policy = T)
whiteoak_moranscatter <- fn_moranscatter(spdf = poly.sp
                                       , lw = neighbour$lw
                                       , var = "whiteoak"
                                       , title = "(e), White Oak"
                                       , col = whiteoak_colours$uniquecolour[length(whiteoak_colours$uniquecolour)]
                                       , zero.policy = T)


maple_moranscatter$I
hickory_moranscatter$I
blackoak_moranscatter$I
redoak_moranscatter$I
whiteoak_moranscatter$I


# 6c. View plots
png(paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name), "_2_global_moran_scatter.png"), width = 750, height = 500)

par(mfrow = c(2,3), mar = c(5,5,5,5), bg = "Grey80")
maple_moranscatter$moranscatter()
hickory_moranscatter$moranscatter()
blackoak_moranscatter$moranscatter()
redoak_moranscatter$moranscatter()
whiteoak_moranscatter$moranscatter()

dev.off()

# 7. Calculate Global Moran's I
# 7a. Calculate value

# 7b. Combine all Moran's I's



# 7b. Check significance with Monte Carlo Method
maple_mc <- spdep::moran.mc(x = poly.sp$maple
                                , listw = neighbour$lw
                                , nsim = 1000
                                , zero.policy = T)
hickory_mc <- spdep::moran.mc(x = poly.sp$hickory
                            , listw = neighbour$lw
                            , nsim = 1000
                            , zero.policy = T)
blackoak_mc <- spdep::moran.mc(x = poly.sp$blackoak
                            , listw = neighbour$lw
                            , nsim = 1000
                            , zero.policy = T)
redoak_mc <- spdep::moran.mc(x = poly.sp$redoak
                            , listw = neighbour$lw
                            , nsim = 1000
                            , zero.policy = T)
whiteoak_mc <- spdep::moran.mc(x = poly.sp$whiteoak
                            , listw = neighbour$lw
                            , nsim = 1000
                            , zero.policy = T)
global_moran_i <- data.frame(Variable = colsuse
                             , `Global Moran's I` = rbind(maple_mc$statistic
                                                          , hickory_mc$statistic
                                                          , blackoak_mc$statistic
                                                          , redoak_mc$statistic
                                                          , whiteoak_mc$statistic))
write.csv(global_moran_i
          , paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name), "_3_global_moran_i.csv")
          , row.names = F)

global_moran_i_mc <- data.frame(Variable = colsuse
                                , Statistic = rbind(maple_mc$statistic
                                                    , hickory_mc$statistic
                                                    , blackoak_mc$statistic
                                                    , redoak_mc$statistic
                                                    , whiteoak_mc$statistic)
                                , `p-value` = rbind(maple_mc$p.value
                                                    , hickory_mc$p.value
                                                    , blackoak_mc$p.value
                                                    , redoak_mc$p.value
                                                    , whiteoak_mc$p.value))
write.csv(global_moran_i_mc
          , paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name), "_4_global_moran_i_mc.csv")
          , row.names = F)
# 8. Calculate Local Moran's I
# 8a. Calculate and determine colours
maple_local_moran <- fn_localmoran(spdf = poly.sp
                                       , lw = neighbour$lw
                                       , var = "maple"
                                       , zero.policy = F)
saveRDS(maple_local_moran, paste0("Datasets/", ds.name, "/esda/esda_", tolower(ds.name) , "_maple_local_obj.RDS"))
maple_local_moran_colours <- fn_mapcolours(var = maple_local_moran$Ii
                                               , num_colours = 5
                                               , brewerpal = "YlGn"
                                               , dig.lab = 3)
hickory_local_moran <- fn_localmoran(spdf = poly.sp
                                       , lw = neighbour$lw
                                       , var = "hickory"
                                       , zero.policy = F)
hickory_local_moran_colours <- fn_mapcolours(var = hickory_local_moran$Ii
                                           , num_colours = 5
                                           , brewerpal = "YlGn"
                                           , dig.lab = 3)
blackoak_local_moran <- fn_localmoran(spdf = poly.sp
                                   , lw = neighbour$lw
                                   , var = "blackoak"
                                   , zero.policy = F)
blackoak_local_moran_colours <- fn_mapcolours(var = blackoak_local_moran$Ii
                                           , num_colours = 5
                                           , brewerpal = "YlGn"
                                           , dig.lab = 3)
redoak_local_moran <- fn_localmoran(spdf = poly.sp
                                   , lw = neighbour$lw
                                   , var = "redoak"
                                   , zero.policy = F)
redoak_local_moran_colours <- fn_mapcolours(var = redoak_local_moran$Ii
                                           , num_colours = 5
                                           , brewerpal = "YlGn"
                                           , dig.lab = 3)

whiteoak_local_moran <- fn_localmoran(spdf = poly.sp
                                   , lw = neighbour$lw
                                   , var = "whiteoak"
                                   , zero.policy = F)
whiteoak_local_moran_colours <- fn_mapcolours(var = whiteoak_local_moran$Ii
                                           , num_colours = 5
                                           , brewerpal = "YlGn"
                                           , dig.lab = 3)
# 8b. Plot values
png(paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name), "_5_local_moran_plots.png"), width = 750, height = 500)
par(mfrow = c(2,3), mar = c(5,5,5,5), bg = "Grey80")
fn_spdfplot(poly.sp
            , colours = maple_local_moran_colours
            , main = "(a), Maples")
fn_spdfplot(poly.sp
            , colours = hickory_local_moran_colours
            , main = "(b), Hickory")
fn_spdfplot(poly.sp
            , colours = blackoak_local_moran_colours
            , main = "(c), Black Oak")
fn_spdfplot(poly.sp
            , colours = redoak_local_moran_colours
            , main = "(d), Red Oak")
fn_spdfplot(poly.sp
            , colours = whiteoak_local_moran_colours
            , main = "(e), White Oak")

dev.off()

# 8c. Determine cluster quadrants
maple_local_quadrant <- fn_quadrant(spdf = maple_local_moran
                                        , lw = neighbour$lw
                                        , var = "maple")
saveRDS(maple_local_quadrant, paste0("Datasets/", ds.name, "/esda/esda_", tolower(ds.name) , "_maple_local_quadrant.RDS"))
table(maple_local_quadrant$quadrant)

hickory_local_quadrant <- fn_quadrant(spdf = hickory_local_moran
                                    , lw = neighbour$lw
                                    , var = "hickory")
table(hickory_local_quadrant$quadrant)

redoak_local_quadrant <- fn_quadrant(spdf = redoak_local_moran
                                    , lw = neighbour$lw
                                    , var = "redoak")
table(redoak_local_quadrant$quadrant)

blackoak_local_quadrant <- fn_quadrant(spdf = blackoak_local_moran
                                    , lw = neighbour$lw
                                    , var = "blackoak")
table(blackoak_local_quadrant$quadrant)

whiteoak_local_quadrant <- fn_quadrant(spdf = whiteoak_local_moran
                                    , lw = neighbour$lw
                                    , var = "whiteoak")
table(whiteoak_local_quadrant$quadrant)

# 8d. Plot LISA cluster quadrants and proportion of spatial entities that fall in each cluster

png(paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name), "_6_maple_LISA.png"), width = 750, height = 500)
fn_plotLISAquad(spdf = poly.sp
                , var = "maple"
                , varname.print = "Maple"
                , quadrant = maple_local_quadrant
                , legend.pos = "bottomright"
                , graph.legend = 1
                , margin = 4.85)
dev.off()
png(paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name), "_6_hickory_LISA.png"), width = 750, height = 500)
fn_plotLISAquad(spdf = poly.sp
                , var = "hickory"
                , varname.print = "Hickory"
                , quadrant = hickory_local_quadrant
                , legend.pos = "bottomright"
                , graph.legend = 1
                , margin = 5)
dev.off()
png(paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name), "_6_redoak_LISA.png"), width = 750, height = 500)
fn_plotLISAquad(spdf = poly.sp
                , var = "redoak"
                , varname.print = "Red Oak"
                , quadrant = redoak_local_quadrant
                , legend.pos = "bottomright"
                , graph.legend = 1
                , margin = 5)
dev.off()
png(paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name), "_6_blackoak_LISA.png"), width = 750, height = 500)
fn_plotLISAquad(spdf = poly.sp
                , var = "blackoak"
                , varname.print = "Black Oak"
                , quadrant = blackoak_local_quadrant
                , legend.pos = "bottomright"
                , graph.legend = 1
                , margin = 5)
dev.off()
png(paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name), "_6_whiteoak_LISA.png"), width = 750, height = 500)
fn_plotLISAquad(spdf = poly.sp
                , var = "whiteoak"
                , varname.print = "White Oak"
                , quadrant = whiteoak_local_quadrant
                , legend.pos = "bottomright"
                , graph.legend = 1
                , margin = 5)
dev.off()

# 8e. Plot LISA clusters

png(paste0("Datasets/", ds.name, "/", section, "/", section, "_", tolower(ds.name), "_7_LISA_Clusters.png"), width = 750, height = 500)
par(mfrow = c(2,3), mar = c(5,5,5,5), bg = "Grey80")
fn_LISAclustermap(spdf = poly.sp
                  , title = "(a), Maple"
                  , colors = maple_local_moran_colours
                  , quadrant = maple_local_quadrant)

fn_LISAclustermap(spdf = poly.sp
                  , title = "(b), Hickory"
                  , colors = hickory_local_moran_colours
                  , quadrant = hickory_local_quadrant)

fn_LISAclustermap(spdf = poly.sp
                  , title = "(c), Black Oak"
                  , colors = blackoak_local_moran_colours
                  , quadrant = blackoak_local_quadrant)

fn_LISAclustermap(spdf = poly.sp
                  , title = "(d), Red Oak"
                  , colors = redoak_local_moran_colours
                  , quadrant = redoak_local_quadrant)

fn_LISAclustermap(spdf = poly.sp
                  , title = "(e), White Oak"
                  , colors = whiteoak_local_moran_colours
                  , quadrant = whiteoak_local_quadrant)
plot(0,type='n',axes=FALSE,ann=FALSE)
legend("top"
       , title = "Cluster"
       , legend = c("insignificant (I)","low-low (L-L)","low-high (L-H)","high-low (H-L)","high-high (H-H)")
       , fill=c("grey50","blue",rgb(0,0,1,alpha=0.4),rgb(1,0,0,alpha=0.4),"red"))

dev.off()



# 9. Bivariate Moran's I
# 9a. Calculate bivariate Global and Local Moran's I 
maple_hickory_bivmoran <- fn_bivariatemoran(x = poly.sp@data[,"maple"]
                                            , y = poly.sp@data[,"hickory"]
                                            , nb = neighbour$nb)

crime_veh_res_bivmoran_colours <- fn_mapcolours(var = crime_veh_res_bivmoran$local
                                                , num_colours = 8
                                                , brewerpal = "Reds"
                                                , dig.lab = 0)
# 10. Calculate dispersion
df_dispersion <- fn_dispersion(poly.sp, colsuse)

dispersion <- apply(data.frame(poly.sp) %>% dplyr::select(colsuse), 2, function(x) var(x, na.rm = T)/mean(x, na.rm = T))
df <- data.frame(colsuse, dispersion)
rownames(df) <- 1:length(colsuse)
return(df)
