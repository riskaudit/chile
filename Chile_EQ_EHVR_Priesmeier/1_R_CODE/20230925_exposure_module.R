# load libraries ####
library(readxl)
library(sp)
library(sf)
library(data.table)
library(tidyr)
library(raster)
library(rgeos)
library(ggplot2)
library(gstat)
library(xts)
library(tmap)
library(reshape)
library(scales)
library(stringr)
library(units)
library(rgdal)
library(tmaptools)
library(plyr)
library(dplyr) 

# clean R environment ####
rm(list=ls())

# set working directory ####
setwd("C:/Users/admin/Desktop/GitHub/summer2023/Chile_EQ_EHVR_Priesmeier")

# enable garbage collection ####
gc()

# select option for implementation ####
Calculate_Buildings <- "YES"                                                          
if (Calculate_Buildings == "YES") {                                                  
  nameUnit <- ""                                                                 
} else {                                                                              
  nameUnit <- "_A"                                                                    
}                                                                                     

# load exposure dataset (census) ####
Exposure_Model <- read.csv("2_INFO_INPUT/SARA_ENGINE/Chile Exposure Model Municipality Level_Peter.csv",sep=",",header= TRUE, stringsAsFactors = FALSE)  

# correct some entries
Exposure_Model[which(Exposure_Model$building_type == "CR/LWAL/HBET:3,9/RES+RES2"),"building_type"] <- "CR/LWAL/HBET:4,9/RES+RES2"   

# load grid system ####  
Grid <- st_read("3_SHAPES_INPUT/SANTIAGO_HD_GRID_small.shp") # (normaled Basis Grid)

# load working coordinate reference system ####
Working_Crs <- st_read("3_SHAPES_INPUT/COMUNA_C17.shp")                               
# load comuna admin zones ####
Comunas_InputSHP <- st_read("3_SHAPES_INPUT/COMUNA_C17_small.shp")   

# load nDSM file #### 
nDSM <- raster("3_SHAPES_INPUT/finfin_ndsm_meters.tif") 
nameRast <- "" 

# load Class H categories ####
BT_HClass_Excel <- read_excel("2_INFO_INPUT/Building_type_ClassH_HBET.xls") 

# sort Class H categories
HBET_Class_Order <- c("HBET_1_2","HBET_3", "HBET_4_5", "HBET_6_9","HBET_10_24","HBET_25_40")

# create directory #####    
dir.create(paste0(getwd(),"/4_RESULTS/",Sys.Date(),
                  "_exposure",nameUnit,nameRast))                                                       
path_abbr <- paste0(getwd(),"/4_RESULTS/",Sys.Date(),
                    "_exposure",nameUnit,nameRast)

# add subfolders
dir.create(paste0(path_abbr,"/4_1_Exposure"))
dir.create(paste0(path_abbr,"/4_2_Hazard"))
dir.create(paste0(path_abbr,"/4_3_Vulnerability"))
dir.create(paste0(path_abbr,"/4_4_Risk"))
dir.create(paste0(path_abbr,"/4_5_Population"))

# use GEM taxonomy naming convention #### 

Exposure_Model$building_type <- recode(Exposure_Model$building_type, 
                                       "CR/LWAL/HBET:1,3/RES+RES1" = "CR_1_3", 
                                       "MATO/LWAL+DNO/HBET:1,2/RES+RES6" = "MATO_1_2",
                                       "MCF+CBH+MOC/LWAL/HBET:1,2/RES+RES1(73%)andMR+CBH+RS+MOC/LWAL/HBET:1,2/RES+RES1(27%)" = "MCF_MR_1_2",
                                       "MCF+CLBRS+MOC/LWAL/HBET:1,2/RES+RES1(27%)andMCF+CLBRH+MOC/LWAL/HBET:1,2/RES+RES1(73%)" = "MCF_MCF_1_2",
                                       "MR+CLBRH+RS+MOC/LWAL/HBET:1,2/RES+RES1" = "MR_1_2",
                                       "MUR+ADO+MOM/LWAL+DNO/HBET:1,2/RES+RES1" = "MURADO_1_2",
                                       "MUR+CLBRS+MOC/LWAL/HBET:1,2/RES+RES1" = "MUR_1_2",
                                       "W/LWAL/HBET:1,3/RES+RES1" = "W_1_3",
                                       "W/LWAL+DNO/HBET:1,2/RES" = "WDNO_1_2",
                                       "MCF+CBH+MOC/LWAL/HEX:3/RES+RES2(70%)andMR+CBH+RS+MOC/LWAL/HEX:3/RES+RES2(30%)" = "MCF_MR_3",
                                       "MCF+CLBRS+MOC/LWAL/HEX:3/RES+RES2(30%)andMCF+CLBRH+MOC/LWAL/HEX:3/RES+RES2(70%)" = "MCF_MCF_3",
                                       "MR+CLBRH+RS+MOC/LWAL/HEX:3/RES+RES2" = "MR_3",
                                       "MCF+CBH+MOC/LWAL/HBET:4,5/RES+RES2(83%)andMR+CBH+RS+MOC/LWAL/HBET:4,5/RES+RES2(17%)" = "MCF_MR_4_5",
                                       "MR+CLBRH+RS+MOC/LWAL/HBET:4,5/RES+RES2" = "MR_4_5",
                                       "MCF+CLBRS+MOC/LWAL/HBET:4,5/RES+RES2(80%)andMCF+CLBRH+MOC/LWAL/HBET:4,5/RES+RES2(20%)" = "MCF_MCF_4_5",
                                       "CR/LWAL/HBET:4,9/RES+RES2" = "CR_4_9",
                                       "CR/LWAL/HBET:10,24/RES+RES2" = "CR_10_24",
                                       "CR/LWAL/HBET:25,40/RES+RES2" = "CR_25_40",)

BT_HClass_Excel$building_type <- recode(BT_HClass_Excel$building_type, "CR/LWAL/HBET:1,3/RES+RES1" = "CR_1_3", 
                                        "MATO/LWAL+DNO/HBET:1,2/RES+RES6" = "MATO_1_2",
                                        "MCF+CBH+MOC/LWAL/HBET:1,2/RES+RES1(73%)andMR+CBH+RS+MOC/LWAL/HBET:1,2/RES+RES1(27%)" = "MCF_MR_1_2",
                                        "MCF+CLBRS+MOC/LWAL/HBET:1,2/RES+RES1(27%)andMCF+CLBRH+MOC/LWAL/HBET:1,2/RES+RES1(73%)" = "MCF_MCF_1_2",
                                        "MR+CLBRH+RS+MOC/LWAL/HBET:1,2/RES+RES1" = "MR_1_2",
                                        "MUR+ADO+MOM/LWAL+DNO/HBET:1,2/RES+RES1" = "MURADO_1_2",
                                        "MUR+CLBRS+MOC/LWAL/HBET:1,2/RES+RES1" = "MUR_1_2",
                                        "W/LWAL/HBET:1,3/RES+RES1" = "W_1_3",
                                        "W/LWAL+DNO/HBET:1,2/RES" = "WDNO_1_2",
                                        "MCF+CBH+MOC/LWAL/HEX:3/RES+RES2(70%)andMR+CBH+RS+MOC/LWAL/HEX:3/RES+RES2(30%)" = "MCF_MR_3",
                                        "MCF+CLBRS+MOC/LWAL/HEX:3/RES+RES2(30%)andMCF+CLBRH+MOC/LWAL/HEX:3/RES+RES2(70%)" = "MCF_MCF_3",
                                        "MR+CLBRH+RS+MOC/LWAL/HEX:3/RES+RES2" = "MR_3",
                                        "MCF+CBH+MOC/LWAL/HBET:4,5/RES+RES2(83%)andMR+CBH+RS+MOC/LWAL/HBET:4,5/RES+RES2(17%)" = "MCF_MR_4_5",
                                        "MR+CLBRH+RS+MOC/LWAL/HBET:4,5/RES+RES2" = "MR_4_5",
                                        "MCF+CLBRS+MOC/LWAL/HBET:4,5/RES+RES2(80%)andMCF+CLBRH+MOC/LWAL/HBET:4,5/RES+RES2(20%)" = "MCF_MCF_4_5",
                                        "CR/LWAL/HBET:4,9/RES+RES2" = "CR_4_9",
                                        "CR/LWAL/HBET:10,24/RES+RES2" = "CR_10_24",
                                        "CR/LWAL/HBET:25,40/RES+RES2" = "CR_25_40",)


# pre-clean the exposure dataset (validating avg floor area values) ####
Exposure_Model <- Exposure_Model[ ,c("municipality_id","building_type",
                                     "number_buildings","avg_floor_area_building..m2.",
                                     "replace_cost_per_building_area..USD.m2.","is_urban")]
colnames(Exposure_Model)[1]<-"COMUNA"

# replace na values with 30
Exposure_Model$avg_floor_area_building..m2.[is.na(Exposure_Model$avg_floor_area_building..m2.)] <- 30

# validate avg floor area with zero number of building
for (i in 1:nrow(Exposure_Model)) {
  if (Exposure_Model[i,"number_buildings"] == 0) {
    Exposure_Model[i,"avg_floor_area_building..m2."] <- 0 
  }
}

# further consistency checks between avg floor area and number of buildings
Exposure_Model$total_floor_area <- Exposure_Model$number_buildings * Exposure_Model$avg_floor_area_building..m2.
TempZW <- Exposure_Model %>% group_by(COMUNA,building_type) %>% summarise(nbr_build_p_COM=sum(number_buildings))
Exposure_Model <- merge(Exposure_Model,TempZW, by = c("COMUNA","building_type"))
TempZW <- Exposure_Model %>% group_by(COMUNA,building_type) %>% summarise(flr_area_p_COM=sum(total_floor_area))
Exposure_Model <- merge(Exposure_Model,TempZW, by = c("COMUNA","building_type")) 
Exposure_Model$avg_floor_area_building..m2. <- Exposure_Model$flr_area_p_COM / Exposure_Model$nbr_build_p_COM
Exposure_Model[is.na(Exposure_Model)] <- 0

# combine features regardless of Is_Urban value
Exposure_Model <- Exposure_Model %>% group_by(COMUNA,building_type) %>% summarise(avg_floor_area_building..m2.=mean(avg_floor_area_building..m2.), number_buildings=sum(number_buildings),replace_cost_per_building_area..USD.m2.=mean(replace_cost_per_building_area..USD.m2.))

# export the cleaned dataset
write.csv(Exposure_Model,paste0(path_abbr,"/","Angepasstes_Exposure_Model_Santa_Maria",".txt"),sep="\t",row.names=FALSE)

# initialize exposure array variable based on area = avg floor area x no. of floors ####
ExpM_A_BT <- Exposure_Model
ExpM_A_BT$Area_BT <- ExpM_A_BT$avg_floor_area_building..m2. * ExpM_A_BT$number_buildings    
ExpM_A_BT <- ExpM_A_BT[ ,c("COMUNA","building_type","Area_BT")]

# extract building types
BT <- unique(ExpM_A_BT$building_type)
  
# spread area by building type
ExpM_A_BT_Vertical<- ExpM_A_BT %>% spread(key=building_type, value=Area_BT)

# spread no. of floors by building type
ExpM_No_BT_Vertical <- Exposure_Model[,c("COMUNA","building_type","number_buildings")] %>% spread(key=building_type, value=number_buildings)

# append the corresponding Class H category to every building type
ExpM_A_BT_HBET <- merge(ExpM_A_BT, BT_HClass_Excel, by = "building_type")

# sum floor areas by comuna and by Class H
Area_p_HBET_COM <-ExpM_A_BT_HBET %>% group_by(COMUNA,ClassH) %>% summarise(Area_BT=sum(Area_BT))
colnames(Area_p_HBET_COM)<-c("COMUNA", "ClassH", "Area_p_HBET_COM")
ExpM_A_BT_HBET <- merge(ExpM_A_BT_HBET,Area_p_HBET_COM, by = c("COMUNA","ClassH"))

# calculate the percentage of floor area of a particular building type (with respect to the summed floor areas by comuna and by Class H)
ExpM_A_BT_HBET$Area_BT_p_HBET_Perc <- ExpM_A_BT_HBET$Area_BT / ExpM_A_BT_HBET$Area_p_HBET_COM
ExpM_A_BT_HBET[is.na(ExpM_A_BT_HBET)] <- 0

# prepare grid data with percentile and proportion (share) information ####
Grid <- Grid[,c("share_of_b","share_of_1","Q90_obj","Q70_obj","Q50_obj")]

# clean trim the communa administrative zone
Comunas_InputSHP <- st_transform(Comunas_InputSHP,st_crs(Working_Crs)) 
Comunas_InputSHP <- Comunas_InputSHP[ ,c("COMUNA","geometry")]

# transform the crs of grid to be consistent with other geospatial files
Grid <- st_transform(Grid,st_crs(Comunas_InputSHP)) 
grid_input <- Grid

# add row numbering
Grid$GridNr <- seq.int(nrow(Grid))

# identify the centroid
Grid_Centroid<-st_centroid(Grid)

# add communa identification to every centroid point
Intersection_Centroid<-st_intersection(Grid_Centroid, Comunas_InputSHP[][])
st_geometry(Intersection_Centroid) <- NULL

# clean trim grid with intersected information from various layers
Intersection_Grid <- merge(Grid, Intersection_Centroid, by = c("share_of_b", "share_of_1", "Q90_obj", "Q70_obj", "Q50_obj", "GridNr"), all.x = T)
Intersection_Grid_OG <- st_drop_geometry(Intersection_Grid)
Intersection_Grid_OG <- Intersection_Grid_OG[order(Intersection_Grid_OG$GridNr),]

# correct Grid number for NA entries
for (i in Intersection_Grid_OG$GridNr) {
  if (is.na(Intersection_Grid_OG[which(Intersection_Grid_OG$GridNr == i),"COMUNA"])) {
    Intersection_Grid_OG[which(Intersection_Grid_OG$GridNr == i),"COMUNA"] <- 
      Intersection_Grid_OG[which(Intersection_Grid_OG$GridNr == i-1),"COMUNA"]
  }
}
Intersection_Grid <- merge(Grid, Intersection_Grid_OG, by = c("share_of_b", "share_of_1", "Q90_obj", "Q70_obj", "Q50_obj", "GridNr"))

# identify the communa identification within the grid extent ####
Comunas_InputSHP <- subset(Comunas_InputSHP, COMUNA %in% unique(Intersection_Grid$COMUNA))
# VARIOUS EXPOSURE DISAGGREGATION ####
# 1. Comuna ####

# prepare array with communa and redistributed floor areas (by communa, by building type)
Comuna <- merge(Comunas_InputSHP, ExpM_A_BT_Vertical, by = "COMUNA")

# add redistributed number of buildings (by communa, by building type) to the said array
if (Calculate_Buildings == "YES") {
  Comuna <- merge(Comunas_InputSHP, ExpM_No_BT_Vertical, by = "COMUNA")
}

# calculate the total number of buildings per communa across all building types
Comuna$Total_p_LvL <- rowSums(st_drop_geometry(Comuna[,BT]))


# compute the area of every communa polygon (to be expressed in km2)
Comuna$km2_p_COM <- as.numeric(st_area(Comuna)) / 1000000

# re-express the computed area using 500m-by-500m box 
Comuna$G500pCOM <- Comuna$km2_p_COM * 4
Comuna[is.na(Comuna)] <- 0

# compute how may buildings are there for every 500m-by-500m box  
Comuna$U_p_Grd500 <- Comuna$Total_p_LvL / Comuna$G500pCOM 

# export shapefile
Comuna_sp <- as_Spatial(Comuna)
writeOGR(obj=Comuna_sp,
         dsn=paste0(path_abbr,"/4_1_Exposure"),
         layer=paste0("Comuna","_ExpM"),
         driver="ESRI Shapefile",
         overwrite_layer = T)

# 2. Grid500_Comuna (G500C) ####

# prepare array with communa and redistributed floor areas (by communa, by building type)
Temp_Grid500_Comuna_G500C <- merge(Intersection_Grid, ExpM_A_BT_Vertical, by = "COMUNA")

# count the number of grids for each comuna
TempZW_G500C <- Temp_Grid500_Comuna_G500C %>% count(COMUNA)
st_geometry(TempZW_G500C) <- NULL
Temp_Grid500_Comuna_G500C <- merge(Temp_Grid500_Comuna_G500C,TempZW_G500C,by="COMUNA")

# distribute the communa-level BT-specific area to the respective number of grids present for each communa
TempZW_G500C <- Temp_Grid500_Comuna_G500C
st_geometry(TempZW_G500C) <- NULL
for (i in BT) {
  for (l in 1:nrow(TempZW_G500C)) {
    TempZW_G500C[l,i] <- TempZW_G500C[l,i]/TempZW_G500C[l,"n"]
  }
}
Temp_Grid500_Comuna_G500C[,BT] <- TempZW_G500C[,BT]
Grid500_Comuna <- Temp_Grid500_Comuna_G500C

# compute the corresponding number of buildings using ave floor area (by communa, by building type)
if (Calculate_Buildings == "YES") {
  ComunaNr_G500C_G500C <- unique(Grid500_Comuna$COMUNA)
  for (CN_G500C in ComunaNr_G500C_G500C) {
    TempZW_G500C <- subset(Exposure_Model[,c("building_type","avg_floor_area_building..m2.","COMUNA")], COMUNA == CN_G500C)
    for (BTe in BT) {
      Grid500_Comuna[which(Grid500_Comuna$COMUNA==CN_G500C),BTe] <- c(st_drop_geometry(Grid500_Comuna[which(Grid500_Comuna$COMUNA==CN_G500C),BTe])) / TempZW_G500C[which(TempZW_G500C$building_type==BTe),"avg_floor_area_building..m2."]
    }
  }
}
Grid500_Comuna[is.na(Grid500_Comuna)] <- 0

# count (by summing up) the number of buildings per communa across all building types
Grid500_Comuna$Total_p_LvL <- rowSums(st_drop_geometry(Grid500_Comuna[,BT]))

# export shapefile
Grid500_Comuna_sp <- as_Spatial(Grid500_Comuna)
writeOGR(obj=Grid500_Comuna_sp,
         dsn=paste0(path_abbr,"/4_1_Exposure"),
         layer=paste0("Grid500_Comuna","_ExpM"),
         driver="ESRI Shapefile",
         overwrite_layer = T)
 
# 3. Grid500_Absolute (G500A) ####

# define height limits for each Class H category      
Abs_Int_Heigth_G500A <- c(0, 3, 4.5, 7.5, 13.5, 36, max(Grid$Q90_obj))
Abs_HBET_G500A <- c("HBET_1_2","HBET_3", "HBET_4_5", "HBET_6_9","HBET_10_24","HBET_25_40")

# sort the grid of Q90 (increasing order) and label the Class H given respective height limit
Abs_Grid_G500A <- Intersection_Grid
Abs_Grid_G500A <- Abs_Grid_G500A[order(Abs_Grid_G500A$Q90_obj),]  
Abs_Grid_G500A$ClassH <- cut (x = Abs_Grid_G500A$Q90_obj,
                              breaks = Abs_Int_Heigth_G500A,
                              labels = Abs_HBET_G500A)

# add communa-specific and ClassH-specific area (regardless of BT) to every grid cell
Abs_Grid_G500A <- merge(Abs_Grid_G500A,
                        Area_p_HBET_COM,
                        by = c("COMUNA","ClassH"))

# compute the built-up area per grid cell using the satellite-derived "share_of_1" value
Abs_Grid_G500A$Area_WholeCell <- st_area(Abs_Grid_G500A)
Abs_Grid_G500A$Area_WholeCell <- as.numeric(Abs_Grid_G500A$Area_WholeCell)
Abs_Grid_G500A$BUA_p_Cell <- Abs_Grid_G500A$Area_WholeCell * Abs_Grid_G500A$share_of_1

# compute (by summing up across all grid cells) the total built-up area (by communa, by Class H)
BUA_p_HBET_COM_G500A <-Abs_Grid_G500A %>% group_by(COMUNA,ClassH) %>% summarise(BUA_p_Cell=sum(BUA_p_Cell))
st_geometry(BUA_p_HBET_COM_G500A) <- NULL
colnames(BUA_p_HBET_COM_G500A)<-c("COMUNA", "ClassH", "BUA_p_HBET_COM")
Abs_Grid_G500A <- merge(Abs_Grid_G500A,BUA_p_HBET_COM_G500A, by = c("COMUNA","ClassH"))

# compute the grid-cell-specific percentage (with respect to the computed total built-up area, by communa and by Class H)
Abs_Grid_G500A$BUA_p_Cell_HBET_Perc <- Abs_Grid_G500A$BUA_p_Cell / Abs_Grid_G500A$BUA_p_HBET_COM

# using the computed percentage, distribute the  communa-specific and ClassH-specific area (regardless of BT) to every grid cell
Abs_Grid_G500A$Area_HBET_p_Cell <- Abs_Grid_G500A$Area_p_HBET_COM * Abs_Grid_G500A$BUA_p_Cell_HBET_Perc

# get the Area_BT_p_HBET_Perc, which is the percentage of area of a particular BT, given communa and given Class H category
TempZW_G500A <- subset(ExpM_A_BT_HBET, select = c(COMUNA,ClassH,building_type,Area_BT_p_HBET_Perc))

# clean subset and redistribute the Area_BT_p_HBET_Perc for Class H (HBET_1_2)
BT_ExpM_HBET_1_2 <- subset(TempZW_G500A, ClassH == "HBET_1_2")
BT_ExpM_HBET_1_2 <- spread(BT_ExpM_HBET_1_2, key = "building_type", value = "Area_BT_p_HBET_Perc")
Abs_Grid_G500A <- merge(Abs_Grid_G500A, BT_ExpM_HBET_1_2, by = c("COMUNA", "ClassH"), all.x = TRUE)

# clean subset and redistribute the Area_BT_p_HBET_Perc for Class H (HBET_3)
BT_ExpM_HBET_3 <- subset(TempZW_G500A, ClassH == "HBET_3")
BT_ExpM_HBET_3 <- spread(BT_ExpM_HBET_3, key = "building_type", value = "Area_BT_p_HBET_Perc")
Abs_Grid_G500A <- merge(Abs_Grid_G500A, BT_ExpM_HBET_3, by = c("COMUNA", "ClassH"), all.x = TRUE)

# clean subset and redistribute the Area_BT_p_HBET_Perc for Class H (HBET_4_5)
BT_ExpM_HBET_4_5 <- subset(TempZW_G500A, ClassH == "HBET_4_5")
BT_ExpM_HBET_4_5 <- spread(BT_ExpM_HBET_4_5, key = "building_type", value = "Area_BT_p_HBET_Perc")
Abs_Grid_G500A <- merge(Abs_Grid_G500A, BT_ExpM_HBET_4_5, by = c("COMUNA", "ClassH"), all.x = TRUE)

# clean subset and redistribute the Area_BT_p_HBET_Perc for Class H (HBET_6_9)
BT_ExpM_HBET_6_9 <- subset(TempZW_G500A, ClassH == "HBET_6_9")
BT_ExpM_HBET_6_9 <- spread(BT_ExpM_HBET_6_9, key = "building_type", value = "Area_BT_p_HBET_Perc")
Abs_Grid_G500A <- merge(Abs_Grid_G500A, BT_ExpM_HBET_6_9, by = c("COMUNA", "ClassH"), all.x = TRUE)

# clean subset and redistribute the Area_BT_p_HBET_Perc for Class H (HBET_10_24)
BT_ExpM_HBET_10_24 <- subset(TempZW_G500A, ClassH == "HBET_10_24")
BT_ExpM_HBET_10_24 <- spread(BT_ExpM_HBET_10_24, key = "building_type", value = "Area_BT_p_HBET_Perc")
Abs_Grid_G500A <- merge(Abs_Grid_G500A, BT_ExpM_HBET_10_24, by = c("COMUNA", "ClassH"), all.x = TRUE)

# clean subset and redistribute the Area_BT_p_HBET_Perc for Class H (HBET_25_40)
BT_ExpM_HBET_25_40 <- subset(TempZW_G500A, ClassH == "HBET_25_40")
BT_ExpM_HBET_25_40 <- spread(BT_ExpM_HBET_25_40, key = "building_type", value = "Area_BT_p_HBET_Perc")
Abs_Grid_G500A <- merge(Abs_Grid_G500A, BT_ExpM_HBET_25_40, by = c("COMUNA", "ClassH"), all.x = TRUE)

# using the percentage, Area_BT_p_HBET_Perc, for every grid cell, further distribute the grid-cell-specific area (which is communa-specific and CLassH-specific already) to every BT
Abs_Grid_G500A_Calc <- Abs_Grid_G500A
Abs_Grid_G500A_Calc[,BT] <- Abs_Grid_G500A_Calc[,BT] * Abs_Grid_G500A_Calc$Area_HBET_p_Cell
Abs_Grid_G500A_Calc[is.na(Abs_Grid_G500A_Calc)] <- 0
Grid500_Absolut <- Abs_Grid_G500A_Calc

# for every grid cell, compute the total area by Class H category
Grid500_Absolut$HBET_1_2 <- rowSums(st_drop_geometry
                                    (Grid500_Absolut
                                      [,c("CR_1_3","MATO_1_2","MCF_MR_1_2","MCF_MCF_1_2",
                                          "MR_1_2","MURADO_1_2","MUR_1_2","W_1_3","WDNO_1_2")]))
Grid500_Absolut$HBET_3 <- rowSums(st_drop_geometry
                                  (Grid500_Absolut
                                    [,c("MCF_MR_3","MCF_MCF_3","MR_3")]))
Grid500_Absolut$HBET_4_5 <- rowSums(st_drop_geometry
                                    (Grid500_Absolut
                                      [,c("MCF_MR_4_5","MR_4_5","MCF_MCF_4_5")]))
Grid500_Absolut$HBET_6_9 <- rowSums(st_drop_geometry
                                    (Grid500_Absolut
                                      [,c("CR_4_9")]))
Grid500_Absolut$HBET_10_24 <- rowSums(st_drop_geometry
                                      (Grid500_Absolut
                                        [,c("CR_10_24")]))
Grid500_Absolut$HBET_25_40 <- rowSums(st_drop_geometry
                                      (Grid500_Absolut
                                        [,c("CR_25_40")]))

# for every grid cell (now with completely distributed floor areas to every BT), compute the number og buildings using the ave floor area that is specific to BT
if (Calculate_Buildings == "YES") {
 ComunaNr_G500C_G500A <- unique(Grid500_Absolut$COMUNA)
  for (CN_G500A in ComunaNr_G500C_G500A) {
    TempZW_G500A <- subset(Exposure_Model
                           [,c("building_type","avg_floor_area_building..m2.","COMUNA")],
                           COMUNA == CN_G500A)
    for (BTe in BT) {
      Grid500_Absolut[which(Grid500_Absolut$COMUNA==CN_G500A),BTe] <- 
        c(st_drop_geometry(Grid500_Absolut[which(Grid500_Absolut$COMUNA==CN_G500A),BTe])) / 
        TempZW_G500A[which(TempZW_G500A$building_type==BTe),"avg_floor_area_building..m2."]
    }
  }
}
Grid500_Absolut[is.na(Grid500_Absolut)] <- 0

# compute the total number of buildings per grid cell
Grid500_Absolut$Total_p_LvL <- rowSums(st_drop_geometry(Grid500_Absolut[,BT]))

# export shapefile
 Grid500_Absolut_sp <- as_Spatial(Grid500_Absolut)
 writeOGR(obj=Grid500_Absolut_sp,
          dsn=paste0(path_abbr,"/4_1_Exposure"),
          layer=paste0("Grid500_Absolut","_ExpM"),
          driver="ESRI Shapefile",
          overwrite_layer = T)

# 4. Grid500_Relative (G500R) ####

# using the exposure array with area by communa & by BT, convert it to by communa & by Class H category     
TempZW_G500R <- merge(ExpM_A_BT, BT_HClass_Excel, by = "building_type")
ExpM_A_HBET_COM <- TempZW_G500R %>% group_by(COMUNA,ClassH) %>% summarise(Area_BT=sum(Area_BT))
colnames(ExpM_A_HBET_COM)<-c("COMUNA","ClassH","Area_p_HBET_COM")

# add the area occupied for each communa
Area_COM_Tot <-ExpM_A_BT_HBET %>% group_by(COMUNA) %>% summarise(Area_BT=sum(Area_BT))
colnames(Area_COM_Tot)<-c("COMUNA","Area_COM_Tot")
ExpM_A_HBET_COM <- merge(ExpM_A_HBET_COM,Area_COM_Tot, by = "COMUNA")

# compute the percentage of area of each Class H as a fraction of the total area of a communa
ExpM_A_HBET_COM$Area_HBET_p_COM_Perc <- ExpM_A_HBET_COM$Area_p_HBET_COM / ExpM_A_HBET_COM$Area_COM_Tot

# prepare variables (incl. satellite-derived built-up area) for the succeeding for-loop
ComunaNr_G500R <- as.character(unique(Intersection_Grid$COMUNA))
Rel_Grid_G500R <- Intersection_Grid
Rel_Grid_G500R$Area_WholeCell <- st_area(Rel_Grid_G500R)
Rel_Grid_G500R$Area_WholeCell <- as.numeric(Rel_Grid_G500R$Area_WholeCell)
Rel_Grid_G500R$BUA_p_Cell <- Rel_Grid_G500R$Area_WholeCell * Rel_Grid_G500R$share_of_1

# initialize the array resulting from for-loop
Loop_Result_G500R <- Rel_Grid_G500R[0,]
Loop_Result_G500R$Order_ID_Q90 <- vector()
Loop_Result_G500R$Total_BUA <- vector()
Loop_Result_G500R$BUA_Perc_p_Com <- vector()
Loop_Result_G500R$BUA_Perc_p_Com_Cum <- vector()
Loop_Result_G500R$NewHClass <- vector()
st_geometry(Loop_Result_G500R) <- NULL

# start for-loop
for (l in ComunaNr_G500R) {
  
  # given communa l, get the exposure data disaggregated by Class H category
  table_ExpM_G500R <- ExpM_A_HBET_COM[ExpM_A_HBET_COM$COMUNA == l,]
  
  # sort the extracted exposure data based on Class H ordering
  table_ExpM_G500R <- table_ExpM_G500R[match(HBET_Class_Order, table_ExpM_G500R$ClassH),]
  
  # remove zero entries
  table_ExpM_G500R <- table_ExpM_G500R[which(table_ExpM_G500R$Area_HBET_p_COM_Perc != 0),]

  # add cumulative percentage following the Class H ordering
  table_ExpM_G500R$Cum_Perc_HBET_p_COM <- cumsum(table_ExpM_G500R$Area_HBET_p_COM_Perc)  
  
  # given communa l, get the grid cells and sort them based on increasing Q90 
  table_Grid_G500R <- Rel_Grid_G500R[Rel_Grid_G500R$COMUNA == l,]
  table_Grid_G500R <- table_Grid_G500R[order(table_Grid_G500R$Q90_obj),]
  table_Grid_G500R$Order_ID_Q90 <- seq.int(nrow(table_Grid_G500R))
  
  # compute the total satellite-derived BUA (based on 'share_of_1') given communa l and 
  # compute percentage cumulative following the increasing Q90 order 
  table_Grid_G500R$Total_BUA <- sum(table_Grid_G500R$BUA_p_Cell)
  table_Grid_G500R$BUA_Perc_p_Com <- table_Grid_G500R$BUA_p_Cell / table_Grid_G500R$Total_BUA
  table_Grid_G500R$BUA_Perc_p_Com_Cum <- cumsum(table_Grid_G500R$BUA_Perc_p_Com)
  
  # remove invalid geometry
  st_geometry(table_Grid_G500R) <- NULL
  
  # initialize Class H column
  table_Grid_G500R$ClassH <- NA
  
  # identify the only available Class H for our for-loop
  HBET_forloop_G500R <- table_ExpM_G500R$ClassH
  
  # start for-loop
  for (k in length(HBET_forloop_G500R):1) {

    # find the row number in our table_Grid_G500R
    # that corresponds to the cumulative percentage
    # of a particular Class H category (which is ordered)
    x_G500R <- findInterval(table_ExpM_G500R[k,"Cum_Perc_HBET_p_COM"],
                            table_Grid_G500R$BUA_Perc_p_Com_Cum)
    
    # determine which of the lower or upper limit has the closer value
    # to the indexing cumulative percentage value
    Grenzintervall <- which.min(
      c(
      table_ExpM_G500R[k,"Cum_Perc_HBET_p_COM"] - 
        table_Grid_G500R[table_Grid_G500R$Order_ID_Q90 == x_G500R,"BUA_Perc_p_Com_Cum"], 
      table_Grid_G500R[table_Grid_G500R$Order_ID_Q90 == x_G500R+1,"BUA_Perc_p_Com_Cum"] -
        table_ExpM_G500R[k,"Cum_Perc_HBET_p_COM"]
      ))
    
    # upper limit does not apply if the for-loop is in its second run
    if (HBET_forloop_G500R[k] != dplyr::last(HBET_forloop_G500R)) {
      if (length(
        table_Grid_G500R[which(table_Grid_G500R[,"ClassH"] == HBET_forloop_G500R[k+1]),"ClassH"]) 
        <= length(1:c(x_G500R-1+c(Grenzintervall)))) {
        
        x_G500R <- length(
          table_Grid_G500R[which(table_Grid_G500R[,"ClassH"] == HBET_forloop_G500R[k+1]),"ClassH"]) - 1
        table_Grid_G500R[1:x_G500R,"ClassH"] <- HBET_forloop_G500R[k]} 
      
      else {table_Grid_G500R[1:c(x_G500R-1+c(Grenzintervall)),"ClassH"] <- HBET_forloop_G500R[k]}} 
    else {table_Grid_G500R[1:c(x_G500R-1+c(Grenzintervall)),"ClassH"] <- HBET_forloop_G500R[k]
    }
  }    
  Loop_Result_G500R <- rbind(Loop_Result_G500R, table_Grid_G500R)
}

# Add the assigned Class H to the
Rel_Grid_G500R <- merge(Rel_Grid_G500R,
                        Loop_Result_G500R, 
                        by = c("share_of_b","share_of_1","Q90_obj","Q70_obj","Q50_obj",
                               "GridNr","COMUNA","Area_WholeCell","BUA_p_Cell"),
                        all.x = TRUE)
                                                 
# using the new Rel_Grid_G500R with assigned Class H categories,
# compute the total absolute BUA by communa and by Class H
BUA_p_HBET_COM <-Rel_Grid_G500R %>% group_by(COMUNA,ClassH) %>% summarise(BUA_p_Cell=sum(BUA_p_Cell))
st_geometry(BUA_p_HBET_COM) <- NULL
colnames(BUA_p_HBET_COM)<-c("COMUNA", "ClassH", "BUA_p_HBET_COM")
Rel_Grid_G500R <- merge(Rel_Grid_G500R, BUA_p_HBET_COM, by = c("COMUNA","ClassH"))

# compute the percentage of grid-cell satellite-derived BUA
# with respect to the total absolute BUA by communa and by Class H
Rel_Grid_G500R$BUA_p_Cell_HBET_Perc <- Rel_Grid_G500R$BUA_p_Cell / Rel_Grid_G500R$BUA_p_HBET_COM

# add the census-based area by comuna and by Class H category
Rel_Grid_G500R <- merge(Rel_Grid_G500R,
                        Area_p_HBET_COM,
                        by = c("COMUNA","ClassH"))

# distribute the census-based area by comuna and by Class H category
# using the satellite-based percentage of area (aka BUA)
Rel_Grid_G500R$Area_HBET_p_Cell <- Rel_Grid_G500R$Area_p_HBET_COM * Rel_Grid_G500R$BUA_p_Cell_HBET_Perc

# add the percentage of each BT for given comuna and for given Class H category
Rel_Grid_G500R <- merge(Rel_Grid_G500R, BT_ExpM_HBET_1_2, by = c("COMUNA", "ClassH"), all.x = TRUE)
Rel_Grid_G500R <- merge(Rel_Grid_G500R, BT_ExpM_HBET_3, by = c("COMUNA", "ClassH"), all.x = TRUE)
Rel_Grid_G500R <- merge(Rel_Grid_G500R, BT_ExpM_HBET_4_5, by = c("COMUNA", "ClassH"), all.x = TRUE)
Rel_Grid_G500R <- merge(Rel_Grid_G500R, BT_ExpM_HBET_6_9, by = c("COMUNA", "ClassH"), all.x = TRUE)
Rel_Grid_G500R <- merge(Rel_Grid_G500R, BT_ExpM_HBET_10_24, by = c("COMUNA", "ClassH"), all.x = TRUE)
Rel_Grid_G500R <- merge(Rel_Grid_G500R, BT_ExpM_HBET_25_40, by = c("COMUNA", "ClassH"), all.x = TRUE)

# using the percentages from the previous step,
# distribute the census-based area for each grid cell
Rel_Grid_G500R_Calc <- Rel_Grid_G500R
Rel_Grid_G500R_Calc[,BT] <- Rel_Grid_G500R_Calc[,BT] * Rel_Grid_G500R_Calc$Area_HBET_p_Cell
Rel_Grid_G500R_Calc[is.na(Rel_Grid_G500R_Calc)] <- 0
Grid500_Relative <- Rel_Grid_G500R_Calc

# compute the total area for each Class H category
# by summing up the distributed area to every building type
Grid500_Relative$HBET_1_2 <- rowSums(
  st_drop_geometry(Grid500_Relative[
  ,c("CR_1_3","MATO_1_2","MCF_MR_1_2","MCF_MCF_1_2",
     "MR_1_2","MURADO_1_2","MUR_1_2","W_1_3","WDNO_1_2")]))
Grid500_Relative$HBET_3 <- rowSums(
  st_drop_geometry(Grid500_Relative[
  ,c("MCF_MR_3","MCF_MCF_3","MR_3")]))
Grid500_Relative$HBET_4_5 <- rowSums(
  st_drop_geometry(Grid500_Relative[
  ,c("MCF_MR_4_5","MR_4_5","MCF_MCF_4_5")]))
Grid500_Relative$HBET_6_9 <- rowSums(
  st_drop_geometry(Grid500_Relative[
  ,c("CR_4_9")]))
Grid500_Relative$HBET_10_24 <- rowSums(
  st_drop_geometry(Grid500_Relative[
  ,c("CR_10_24")]))
Grid500_Relative$HBET_25_40 <- rowSums(
  st_drop_geometry(Grid500_Relative[
  ,c("CR_25_40")]))

# compute the corresponding number of buildings given average floor area
if (Calculate_Buildings == "YES") {
  ComunaNr_G500C_G500R <- unique(Grid500_Relative$COMUNA)
  
  # for each comuna
  for (CN_G500R in ComunaNr_G500C_G500R) {
    
    # subset the exposure model given selected comuna
    TempZW_G500R <- subset(Exposure_Model[
      ,c("building_type","avg_floor_area_building..m2.","COMUNA")], COMUNA == CN_G500R)
    
    # for each building type
    for (BTe in BT) {
      Grid500_Relative[which(Grid500_Relative$COMUNA==CN_G500R),BTe] <- c(
        st_drop_geometry(Grid500_Relative[which(Grid500_Relative$COMUNA==CN_G500R),BTe])) / 
        TempZW_G500R[which(TempZW_G500R$building_type==BTe),"avg_floor_area_building..m2."]
    }
  }
}
Grid500_Relative[is.na(Grid500_Relative)] <- 0

# compute the total number of buildings per grid cell
Grid500_Relative$Total_p_LvL <- rowSums(st_drop_geometry(Grid500_Relative[,BT]))

# export shapefile
Grid500_Relative_sp <- as_Spatial(Grid500_Relative)
writeOGR(obj=Grid500_Relative_sp,
         dsn=paste0(path_abbr,"/4_1_Exposure"),
         layer=paste0("Grid500_Relative","_ExpM"),
         driver="ESRI Shapefile",
         overwrite_layer = T)

# 5. GridShare_Abs (GSA) ####

# set grid projection to use nDSM projection
Grid_sp <- as_Spatial(Grid)
Grid_spTr <- spTransform(Grid_sp, crs(nDSM))

# initialize list containing clipped raster images
Raster_Clip_List <- list()

# label the raster clip using the grid ID
for (i in 1:nrow(Grid_spTr)) {
  Raster_Clip   <- crop(nDSM,Grid_spTr[i,])
  Raster_Clip_List[[i]] <- Raster_Clip
  name <- paste("Raster_Clip", Grid[i,"GridNr"], sep = "_")
  names(Raster_Clip_List)[i] <- name
}

# initialize list variable to be assigned with Class H value
Raster_HBET_List_GSA <- list()  

# reclassify the nDSM values using the following criteria
for (i in 1:length(Raster_Clip_List)) {
  Raster_HBET_GSA <- reclassify(Raster_Clip_List[[i]],
                                c(0,6,6,
                                  6,9,9,
                                  9,15,15,
                                  15,27,27,
                                  27,72,72,
                                  72,Inf,120) )
  Raster_HBET_List_GSA[[i]] <- Raster_HBET_GSA 
  name <- paste("Raster_HBET_GSA", Grid[i,"GridNr"], sep = "_")
  names(Raster_HBET_List_GSA)[i] <- name
}  

# initialize the result array
BT_Share_Perc_GSA <- setNames(data.frame(matrix(ncol = 9, nrow = 0)),
                              c("COMUNA","ClassH","HBET","Count","GridNr","Total_Pixel",
                                "Perc_Count_Tot_Pix","building_type","Area_BT_p_HBET_Perc"))

# for raster clips with no pixels in it, initialize a temporary data frame
If_RasterClip_noPixel_GSA <-  setNames(data.frame(matrix(ncol = 13, nrow = 0)), 
                                       c("COMUNA","ClassH","HBET","Count","GridNr",
                                         "Pixel_Count_ClassH_p_COM","Perc_HBETPixel_p_COM",
                                         "Area_p_HBET_COM","Area_p_HBET_Share","building_type",
                                         "Area_BT","Area_BT_p_HBET_Perc","Area_BT_GSA"))
No_Pixel_Table_GSA <- If_RasterClip_noPixel_GSA
If_RasterClip_noPixel_GSA[1,] <- NA

# start for-loop using the re-classified nDSM (takes time)
for (i in 1:length(Raster_HBET_List_GSA)) { 

  # for the selected grid ID, get the pixel values or aka as Class H category
  Pixel_Count_GSA <-as.data.frame(Raster_HBET_List_GSA[[i]])
  colnames(Pixel_Count_GSA) <- "HBET"

  # perform summary statistics (counting)
  Pixel_Count_Share_GSA <- plyr::count(Pixel_Count_GSA, vars = "HBET")
  colnames(Pixel_Count_Share_GSA) <- c("HBET","Count")
  
  # track the grid ID
  Pixel_Count_Share_GSA$GridNr <- i

  # for incomplete clips, replace NA values with zero
  if (any(is.na(Pixel_Count_Share_GSA$HBET))) {
    Pixel_Count_Share_GSA[1,"Count"] <-  Pixel_Count_Share_GSA[1,"Count"] + 
      Pixel_Count_Share_GSA[which(is.na(Pixel_Count_Share_GSA)),"Count"]
    Pixel_Count_Share_GSA <- na.omit(Pixel_Count_Share_GSA)
  }

  # remove the row of HBET being 0
  Pixel_Count_Share_GSA <- Pixel_Count_Share_GSA[-c(Pixel_Count_Share_GSA$HBET==0),]

  # if the raster clip has no pixels, identify and note the communa and grid ID (aka gridNr)
  if (empty(Pixel_Count_Share_GSA)) {
    If_RasterClip_noPixel_GSA[1,"GridNr"]<- i
    TempZW_GSA <- Intersection_Grid
    st_geometry(TempZW_GSA) <- NULL
    If_RasterClip_noPixel_GSA$COMUNA <- TempZW_GSA[TempZW_GSA$GridNr == i,"COMUNA"]
    No_Pixel_Table_GSA <- rbind(No_Pixel_Table_GSA,If_RasterClip_noPixel_GSA)
    print(paste("No Pixels_GSA in Raster_Clip/GridNr:",i))} 
  
  # if the raster clip has pixels, assign Class H category
  else {
    
    # initialize
    Pixel_Count_Share_GSA$ClassH <- 0

    if (any(Pixel_Count_Share_GSA$HBET == 6)) {
        Pixel_Count_Share_GSA[Pixel_Count_Share_GSA$HBET == 6,] <- 
          Pixel_Count_Share_GSA %>%  
          filter(Pixel_Count_Share_GSA[,"HBET"] == 6 ) %>% 
          mutate(ClassH = "HBET_1_2")}

    if (any(Pixel_Count_Share_GSA$HBET == 9)) {
        Pixel_Count_Share_GSA[Pixel_Count_Share_GSA$HBET == 9,] <- 
          Pixel_Count_Share_GSA %>%  
          filter(Pixel_Count_Share_GSA[,"HBET"] == 9 ) %>% 
          mutate(ClassH = "HBET_3")}
 
    if (any(Pixel_Count_Share_GSA$HBET == 15)) {
        Pixel_Count_Share_GSA[Pixel_Count_Share_GSA$HBET == 15,] <-
          Pixel_Count_Share_GSA %>%
          filter(Pixel_Count_Share_GSA[,"HBET"] == 15 ) %>%
          mutate(ClassH = "HBET_4_5")}

    if (any(Pixel_Count_Share_GSA$HBET == 27)) {
        Pixel_Count_Share_GSA[Pixel_Count_Share_GSA$HBET == 27,] <-
          Pixel_Count_Share_GSA %>% 
          filter(Pixel_Count_Share_GSA[,"HBET"] == 27 ) %>%
          mutate(ClassH = "HBET_6_9")} 

    if (any(Pixel_Count_Share_GSA$HBET == 72)) {
        Pixel_Count_Share_GSA[Pixel_Count_Share_GSA$HBET == 72,] <-
          Pixel_Count_Share_GSA %>%
          filter(Pixel_Count_Share_GSA[,"HBET"] == 72 ) %>%
          mutate(ClassH = "HBET_10_24")}

    if (any(Pixel_Count_Share_GSA$HBET == 120)) {
        Pixel_Count_Share_GSA[Pixel_Count_Share_GSA$HBET == 120,] <-
          Pixel_Count_Share_GSA %>%
          filter(Pixel_Count_Share_GSA[,"HBET"] == 120 ) %>%
          mutate(ClassH = "HBET_25_40")}

    # add the comuna information
    TempZW_GSA <- Intersection_Grid
    st_geometry(TempZW_GSA) <- NULL
    Pixel_Count_Share_GSA$COMUNA <- TempZW_GSA[TempZW_GSA$GridNr == i,"COMUNA"]

    # combine all together for all for-loops
    BT_Share_Perc_GSA <- rbind(BT_Share_Perc_GSA,Pixel_Count_Share_GSA)
      
  }
}

# compute the total pixel counts, by comuna and by Class H, covering all grid IDs
TempZW_GSA <- BT_Share_Perc_GSA %>% group_by(COMUNA,ClassH) %>% dplyr::summarise(Count=sum(Count))
colnames(TempZW_GSA) <- c("COMUNA","ClassH","Pixel_Count_ClassH_p_COM")

# simply add the previously calculated "total pixel counts" to our for-loop array
PD_ZW_GSA <- merge(BT_Share_Perc_GSA, TempZW_GSA, by = c("COMUNA","ClassH"))

# compute the percentage of pixel count
# (with respect to the total pixel counts by comuna and by Class H)
PD_ZW_GSA$Perc_HBETPixel_p_COM <- PD_ZW_GSA$Count / PD_ZW_GSA$Pixel_Count_ClassH_p_COM

# add the census-based area by comuna and by Class H category
PD_ZW_GSA <- merge(PD_ZW_GSA, Area_p_HBET_COM, by = c("COMUNA", "ClassH"))

# distribute the census-based area using the calculated percentage of pixel counts
PD_ZW_GSA$Area_p_HBET_Share <- PD_ZW_GSA$Area_p_HBET_COM * PD_ZW_GSA$Perc_HBETPixel_p_COM

# add the census-derived proportion of each BT so that we could further distribute to each BT 
Pixeldensity_GSA <- merge(PD_ZW_GSA,
                          ExpM_A_BT_HBET[,c("COMUNA","ClassH","building_type",
                                            "Area_BT","Area_BT_p_HBET_Perc")], 
                          by = c("COMUNA", "ClassH"))

# compute the area for each BT using the census-dervied proportion and distributed area based on nDSM
Pixeldensity_GSA$Area_BT_GSA <- Pixeldensity_GSA$Area_p_HBET_Share * Pixeldensity_GSA$Area_BT_p_HBET_Perc
Pixeldensity_GSA <- rbind(Pixeldensity_GSA,No_Pixel_Table_GSA)

# vertically spread the previously computed array (by BT)
GridShare_Abs_Complete_Result <- Pixeldensity_GSA
GridShare_Abs_HBET_Result <- Pixeldensity_GSA[
  ,c("COMUNA","ClassH","HBET","Count","GridNr",
     "Pixel_Count_ClassH_p_COM","Perc_HBETPixel_p_COM",
     "Area_p_HBET_COM","Area_p_HBET_Share","building_type",
     "Area_BT_GSA")] %>% 
  spread(key=building_type,
         value=Area_BT_GSA)
GridShare_Abs <- Pixeldensity_GSA[
  ,c("COMUNA","GridNr","building_type","Area_BT_GSA")] %>% 
  spread(key=building_type,
         value=Area_BT_GSA)
GridShare_Abs <- merge(Intersection_Grid, GridShare_Abs, by = c("COMUNA","GridNr"))
GridShare_Abs[is.na(GridShare_Abs)] <- 0

# compute the total area for each Class H category
# by summing up the distributed area to every building type
GridShare_Abs$HBET_1_2_NF <- rowSums(
  st_drop_geometry(GridShare_Abs[
    ,c("CR_1_3","MATO_1_2","MCF_MR_1_2","MCF_MCF_1_2",
       "MR_1_2","MURADO_1_2","MUR_1_2","W_1_3","WDNO_1_2")]))
GridShare_Abs$HBET_3_NF <- rowSums(
  st_drop_geometry(GridShare_Abs[
    ,c("MCF_MR_3","MCF_MCF_3","MR_3")]))
GridShare_Abs$HBET_4_5_NF <- rowSums(
  st_drop_geometry(GridShare_Abs[
    ,c("MCF_MR_4_5","MR_4_5","MCF_MCF_4_5")]))
GridShare_Abs$HBET_6_9_NF <- rowSums(
  st_drop_geometry(GridShare_Abs[
    ,c("CR_4_9")]))
GridShare_Abs$HBET_10_24_NF <- rowSums(
  st_drop_geometry(GridShare_Abs[
    ,c("CR_10_24")]))
GridShare_Abs$HBET_25_40_NF <- rowSums(
  st_drop_geometry(GridShare_Abs[
    ,c("CR_25_40")]))

# compute the corresponding number of buildings per BT
if (Calculate_Buildings == "YES") {
  ComunaNr_G500C_GSA <- unique(GridShare_Abs$COMUNA)
  
  for (CN_GSA in ComunaNr_G500C_GSA) {
    
    TempZW_GSA <- subset(Exposure_Model[,c("building_type","avg_floor_area_building..m2.","COMUNA")],
                         COMUNA == CN_GSA)
    for (BTe in BT) {
      GridShare_Abs[which(GridShare_Abs$COMUNA==CN_GSA),BTe] <-
        c(st_drop_geometry(GridShare_Abs[which(GridShare_Abs$COMUNA==CN_GSA),BTe])) / 
        TempZW_GSA[which(TempZW_GSA$building_type==BTe),"avg_floor_area_building..m2."]
    }
  }
}
GridShare_Abs[is.na(GridShare_Abs)] <- 0

# compute the total number of buildings per grid cell
GridShare_Abs$Total_p_LvL <- rowSums(st_drop_geometry(GridShare_Abs[,BT]))

# export shapefile
GridShare_Abs_sp <- as_Spatial(GridShare_Abs)
writeOGR(obj=GridShare_Abs_sp,
         dsn=paste0(path_abbr,"/4_1_Exposure"),
         layer=paste0("GridShare_Abs","_ExpM"),
         driver="ESRI Shapefile",
         overwrite_layer = T)

# 6. GridShare_Rel (GSR) ####

# initialize the re-classification of nDSM pixels into Class H categories
Hight_Thresh_GridShare_GSR <- data.frame(matrix(nrow = 0, ncol = 1 + length(HBET_Class_Order)))
colnames(Hight_Thresh_GridShare_GSR) <- c("COMUNA",HBET_Class_Order)
TempZW_Schablone_GSR <- Hight_Thresh_GridShare_GSR  

# initialize the GSR list
Raster_HBET_List_GSR <- list() 

# get a list of unique comuna labels
ComunaNr_GSR <- as.character(unique(Intersection_Grid$COMUNA))

# start for-loop for each communa i
for (i in ComunaNr_GSR) {

  # for given communa i, get all grid IDs
  GridNr_fl_GSR <- Intersection_Grid_OG[which(Intersection_Grid_OG$COMUNA == i),"GridNr"]
  
  # initialize arrays
  Pixel_p_COM_GSR <- numeric()
  remove_zero_GSR <- c(0)
  
  # for-loop to combine all pixels belonging to communa i
  for (l in GridNr_fl_GSR) {
    Pixel_p_Grid_GSR <- Raster_Clip_List[[paste("Raster_Clip_",l, sep="")]]@data@values    
    Pixel_p_Grid_GSR <- Pixel_p_Grid_GSR [! Pixel_p_Grid_GSR %in% remove_zero_GSR]
    Pixel_p_Grid_GSR[Pixel_p_Grid_GSR == 0] <- NA
    Pixel_p_Grid_GSR <- na.omit(Pixel_p_Grid_GSR)
    Pixel_p_COM_GSR <- append(Pixel_p_COM_GSR, Pixel_p_Grid_GSR)
  }
  Pixel_p_COM_GSR <- as.data.frame(Pixel_p_COM_GSR)
  setorder(Pixel_p_COM_GSR)
  
  # add cumulative count percentage
  Pixel_p_COM_GSR$Perc <-  100 / nrow(Pixel_p_COM_GSR)
  Pixel_p_COM_GSR$Perc_Cum <-cumsum(Pixel_p_COM_GSR$Perc)
  
  # for communa i, extract and sort the proportion of each Class H category
  HBET_Perc_COM_GSR <- ExpM_A_HBET_COM[which(ExpM_A_HBET_COM$COMUNA == i),]
  HBET_Perc_COM_GSR <- HBET_Perc_COM_GSR[match(HBET_Class_Order, HBET_Perc_COM_GSR$ClassH),]
  HBET_Perc_COM_GSR$Area_HBET_p_COM_Perc_Cum <- cumsum(HBET_Perc_COM_GSR$Area_HBET_p_COM_Perc)
  
  # identify the only available Class H category from the census records given communa i
  HBET_forloop_GSR <- HBET_Perc_COM_GSR[which(HBET_Perc_COM_GSR$Area_p_HBET_COM != 0),"ClassH"]
  
  # add NA entries and communa info to the initialized array 
  TempZW_Schablone_GSR[1,] <- NA
  TempZW_Schablone_GSR[1,"COMUNA"] <- i
  
  # find the interval within the pixel count where the corresponding cumsum from census records falls
  for (h in length(HBET_forloop_GSR):1) {
        
    # get the index where the percentage cumsum appears
    Interval_GSR <- findInterval(
      HBET_Perc_COM_GSR[which(HBET_Perc_COM_GSR$ClassH == HBET_forloop_GSR[h]),
                        "Area_HBET_p_COM_Perc_Cum"]*100,
      Pixel_p_COM_GSR$Perc_Cum)
    
    # extract the pixel value corresponding to the perc_cum intersection
    Pixel_Height_GSR <- Pixel_p_COM_GSR[Interval_GSR,"Pixel_p_COM_GSR"]
    
    # save the pixel value (or height)
    TempZW_Schablone_GSR[1,HBET_forloop_GSR[h]] <- Pixel_Height_GSR
    TempZW_Schablone_GSR[is.na(TempZW_Schablone_GSR)] <- 0
    
  }    
  
  # adjust if there is any overlapping limit between any two or more Class H categories
  TempZW_GSR <- 2:ncol(TempZW_Schablone_GSR)
  for (c in 2: ncol(TempZW_Schablone_GSR)) {
    if (any(TempZW_Schablone_GSR[1,TempZW_GSR[! TempZW_GSR %in% c]] %in% TempZW_Schablone_GSR[1,c])){
      TempZW_Schablone_GSR[1,c] <- TempZW_Schablone_GSR[1,c] - 0.1
      for (c in 2: ncol(TempZW_Schablone_GSR)) {
        if (any(TempZW_Schablone_GSR[1,TempZW_GSR[! TempZW_GSR %in% c]] %in% TempZW_Schablone_GSR[1,c])){
          TempZW_Schablone_GSR[1,c] <- TempZW_Schablone_GSR[1,c] - 0.1
          for (c in 2: ncol(TempZW_Schablone_GSR)) {
            if (any(TempZW_Schablone_GSR[1,TempZW_GSR[! TempZW_GSR %in% c]] %in% TempZW_Schablone_GSR[1,c])){
              TempZW_Schablone_GSR[1,c] <- TempZW_Schablone_GSR[1,c] - 0.1
              for (c in 2: ncol(TempZW_Schablone_GSR)) {
                if (any(TempZW_Schablone_GSR[1,TempZW_GSR[! TempZW_GSR %in% c]] %in% TempZW_Schablone_GSR[1,c])){
                  TempZW_Schablone_GSR[1,c] <- TempZW_Schablone_GSR[1,c] - 0.1
                  
                } 
              }  
            } 
          }  
        } 
      }  
    } 
  }
  TempZW_Schablone_GSR[TempZW_Schablone_GSR<0] <- 0 
  
  # initialize classifying array
  classify_GSR <- numeric()
  First_GSR <- TempZW_Schablone_GSR[1,HBET_forloop_GSR[1]]
  classify_GSR <- append(classify_GSR,c(0,First_GSR,as.numeric(gsub("\\D", "",HBET_forloop_GSR[1])))) 
  
  # do the same if we have 2+ Class H categories
  if (length(HBET_forloop_GSR) >1) {
    for (z in 2:length(HBET_forloop_GSR)) {
      #z <- 2
      Second_1_GSR <- TempZW_Schablone_GSR[1,HBET_forloop_GSR[z-1]]  
      Second_2_GSR <- TempZW_Schablone_GSR[1,HBET_forloop_GSR[z]]
      classify_GSR <- append(classify_GSR,c(Second_1_GSR,
                                            Second_2_GSR,
                                            as.numeric(gsub("\\D", "",HBET_forloop_GSR[z]))))
    }
  }
  
  # postprocess
  classify_GSR[length(classify_GSR)-1] <- Inf     
  
  # for each grid ID in communa i,
  # reclassify each raster clip using the previously calculated reclassification table
  # based on perc_cum of pixel values of height and
  # census-derived proportion by communa and by Class H category
  for (m in GridNr_fl_GSR) {
    Raster_HBET_GSR <- reclassify(Raster_Clip_List[[paste("Raster_Clip_",m,sep="")]],
                                  rcl = classify_GSR)  
    Raster_HBET_List_GSR[[m]] <- Raster_HBET_GSR
    name <- paste("Raster_HBET", m,"GSR", sep = "_")    
    names(Raster_HBET_List_GSR)[m] <- name
  }
  
  Hight_Thresh_GridShare_GSR <- rbind(Hight_Thresh_GridShare_GSR,TempZW_Schablone_GSR)
  
} # end of for-loop for each communa i

# initialize array for the distribution of area to every building type
BT_Share_Perc_GSR <- setNames(data.frame(matrix(ncol = 9, nrow = 0)),
                              c("COMUNA","ClassH","HBET","Count","GridNr",
                                "Total_Pixel","Perc_Count_Tot_Pix","building_type","Area_BT_p_HBET_Perc"))

# identify clips with no pixels
If_RasterClip_noPixel_GSR <-  setNames(data.frame(matrix(ncol = 13, nrow = 0)),
                                       c("COMUNA","ClassH","HBET","Count","GridNr","Pixel_Count_ClassH_p_COM",
                                         "Perc_HBETPixel_p_COM","Area_p_HBET_COM","Area_p_HBET_Share",
                                         "building_type","Area_BT","Area_BT_p_HBET_Perc","Area_BT_GSR"))
No_Pixel_Table_GSR <- If_RasterClip_noPixel_GSR
If_RasterClip_noPixel_GSR[1,] <- NA

# start for-loop for each raster clip i
for (i in 1:length(Raster_HBET_List_GSR)) { 

  # get the pixel count and rename the column with HBET (because it's already re-classified)
  Pixel_Count_GSR <- as.data.frame(Raster_HBET_List_GSR[[i]])
  colnames(Pixel_Count_GSR) <- "HBET"
  
  # get summary statistics (sum of counts) for each HBET
  Pixel_Count_Share_GSR <- plyr::count(Pixel_Count_GSR, vars = "HBET")
  colnames(Pixel_Count_Share_GSR) <- c("HBET","Count")

  # note the raster clip i (or the grid ID, equivalently)
  Pixel_Count_Share_GSR$GridNr <- i
  
  # if there are any NA pixel values, add it to the "HBET = 0"
  if (any(is.na(Pixel_Count_Share_GSR$HBET))) {
    Pixel_Count_Share_GSR[1,"Count"] <-  Pixel_Count_Share_GSR[1,"Count"] + 
      Pixel_Count_Share_GSR[which(is.na(Pixel_Count_Share_GSR)),"Count"]
    Pixel_Count_Share_GSR <- na.omit(Pixel_Count_Share_GSR)
  }
  
  # remove the row with the "HBET = 0"
  Pixel_Count_Share_GSR <- Pixel_Count_Share_GSR[-c(Pixel_Count_Share_GSR$HBET==0),]
  
  # 
  if (empty(Pixel_Count_Share_GSR)) {
    If_RasterClip_noPixel_GSR[1,"GridNr"]<- i
    TempZW_GSR <- Intersection_Grid
    st_geometry(TempZW_GSR) <- NULL
    If_RasterClip_noPixel_GSR$COMUNA <- TempZW_GSR[TempZW_GSR$GridNr == i,"COMUNA"]
    No_Pixel_Table_GSR <- rbind(No_Pixel_Table_GSR,If_RasterClip_noPixel_GSR)
    print(paste("No Pixels_GSR in Raster_Clip/GridNr:",i))} 
  
  # assign the string version of Class H category instead of numeric input
  else {
    Pixel_Count_Share_GSR$ClassH <- 0
    
    if (any(Pixel_Count_Share_GSR$HBET == 12)) {
      Pixel_Count_Share_GSR[Pixel_Count_Share_GSR$HBET == 12,]  <- 
        Pixel_Count_Share_GSR %>%  
        filter(Pixel_Count_Share_GSR[,"HBET"] == 12 ) %>% 
        mutate(ClassH = "HBET_1_2")  
    }
    
    if (any(Pixel_Count_Share_GSR$HBET == 3)) {
      Pixel_Count_Share_GSR[Pixel_Count_Share_GSR$HBET == 3,] <- 
        Pixel_Count_Share_GSR %>% 
        filter(Pixel_Count_Share_GSR[,"HBET"] == 3 ) %>%
        mutate(ClassH = "HBET_3")  
    }
    
    if (any(Pixel_Count_Share_GSR$HBET == 45)) {
      Pixel_Count_Share_GSR[Pixel_Count_Share_GSR$HBET == 45,] <- 
        Pixel_Count_Share_GSR %>%  
        filter(Pixel_Count_Share_GSR[,"HBET"] == 45 ) %>% 
        mutate(ClassH = "HBET_4_5")  
    }
    
    if (any(Pixel_Count_Share_GSR$HBET == 69)) {
      Pixel_Count_Share_GSR[Pixel_Count_Share_GSR$HBET == 69,] <- 
        Pixel_Count_Share_GSR %>%  
        filter(Pixel_Count_Share_GSR[,"HBET"] == 69 ) %>% 
        mutate(ClassH = "HBET_6_9")
    } 
    
    if (any(Pixel_Count_Share_GSR$HBET == 1024)) {
      Pixel_Count_Share_GSR[Pixel_Count_Share_GSR$HBET == 1024,] <- 
        Pixel_Count_Share_GSR %>%  
        filter(Pixel_Count_Share_GSR[,"HBET"] == 1024 ) %>% 
        mutate(ClassH = "HBET_10_24") 
    }
    
    if (any(Pixel_Count_Share_GSR$HBET == 2540)) {
      Pixel_Count_Share_GSR[Pixel_Count_Share_GSR$HBET == 2540,] <- 
        Pixel_Count_Share_GSR %>%  
        filter(Pixel_Count_Share_GSR[,"HBET"] == 2540 ) %>% 
        mutate(ClassH = "HBET_25_40") 
    }
    
    # add the communa information
    TempZW_GSR <- Intersection_Grid
    st_geometry(TempZW_GSR) <- NULL
    Pixel_Count_Share_GSR$COMUNA <- TempZW_GSR[TempZW_GSR$GridNr == i,"COMUNA"]
    
    # combine all results of for-loop
    BT_Share_Perc_GSR <- rbind(BT_Share_Perc_GSR,Pixel_Count_Share_GSR)
    
  }
  
} # end of for-loop for each raster clip i

# compute the total number of pixel counts by comuna and by Class H category
TempZW_GSR <- BT_Share_Perc_GSR %>% group_by(COMUNA,ClassH) %>% dplyr::summarise(Count=sum(Count))
colnames(TempZW_GSR) <- c("COMUNA","ClassH","Pixel_Count_ClassH_p_COM")

# add the previously computed to the for-loop array
PD_ZW_GSR <- merge(BT_Share_Perc_GSR, TempZW_GSR, by = c("COMUNA","ClassH"))

# compute the percentage based on pixel counts
PD_ZW_GSR$Perc_HBETPixel_p_COM <- PD_ZW_GSR$Count / PD_ZW_GSR$Pixel_Count_ClassH_p_COM

# add the area occupied (by communa and by Class H category) from census records
PD_ZW_GSR <- merge(PD_ZW_GSR, Area_p_HBET_COM, by = c("COMUNA", "ClassH"))

# multiply the pixel-based or satellite-derived percentage to the area to distribute it grid-wise
PD_ZW_GSR$Area_p_HBET_Share <- PD_ZW_GSR$Area_p_HBET_COM * PD_ZW_GSR$Perc_HBETPixel_p_COM

# add the information on proportion of each building type from census records
Pixeldensity_GSR <- merge(PD_ZW_GSR,
                          ExpM_A_BT_HBET[,c("COMUNA","ClassH","building_type",
                                            "Area_BT","Area_BT_p_HBET_Perc")],
                          by = c("COMUNA", "ClassH"))

# compute the area of each building type
Pixeldensity_GSR$Area_BT_GSR <- Pixeldensity_GSR$Area_p_HBET_Share * Pixeldensity_GSR$Area_BT_p_HBET_Perc

# combine the resulting array with no_pixel_data to be complete
Pixeldensity_GSR <- rbind(Pixeldensity_GSR,No_Pixel_Table_GSR) 

# vertically spread the area by building type
GridShare_Rel_Complete_Result <- Pixeldensity_GSR
GridShare_Rel_HBET_Result <- Pixeldensity_GSR[,c("COMUNA","ClassH","HBET","Count","GridNr",
                                                 "Pixel_Count_ClassH_p_COM","Perc_HBETPixel_p_COM",
                                                 "Area_p_HBET_COM","Area_p_HBET_Share",
                                                 "building_type","Area_BT_GSR")] %>% 
  spread(key=building_type,
         value=Area_BT_GSR)
GridShare_Rel <- Pixeldensity_GSR[,c("COMUNA","GridNr","building_type","Area_BT_GSR")] %>%
  spread(key=building_type,
         value=Area_BT_GSR)
GridShare_Rel <- merge(Intersection_Grid, GridShare_Rel, by = c("COMUNA","GridNr"))
GridShare_Rel[is.na(GridShare_Rel)] <- 0

# compute the total area for each building type
GridShare_Rel$HBET_1_2 <- rowSums(
  st_drop_geometry(GridShare_Rel[
    ,c("CR_1_3","MATO_1_2","MCF_MR_1_2","MCF_MCF_1_2",
       "MR_1_2","MURADO_1_2","MUR_1_2","W_1_3","WDNO_1_2")]))
GridShare_Rel$HBET_3 <- rowSums(
  st_drop_geometry(GridShare_Rel[
    ,c("MCF_MR_3","MCF_MCF_3","MR_3")]))
GridShare_Rel$HBET_4_5 <- rowSums(
  st_drop_geometry(GridShare_Rel[
    ,c("MCF_MR_4_5","MR_4_5","MCF_MCF_4_5")]))
GridShare_Rel$HBET_6_9 <- rowSums(
  st_drop_geometry(GridShare_Rel[
    ,c("CR_4_9")]))
GridShare_Rel$HBET_10_24 <- rowSums(
  st_drop_geometry(GridShare_Rel[
  ,c("CR_10_24")]))
GridShare_Rel$HBET_25_40 <- rowSums(
  st_drop_geometry(GridShare_Rel[
    ,c("CR_25_40")]))

# compute the number of buildings per building type
if (Calculate_Buildings == "YES") {
  ComunaNr_G500C_GSR <- unique(GridShare_Rel$COMUNA)
  for (CN_GSR in ComunaNr_G500C_GSR) {
    TempZW_GSR <- subset(Exposure_Model[,c("building_type","avg_floor_area_building..m2.","COMUNA")],
                         COMUNA == CN_GSR)
    for (BTe in BT) {
      GridShare_Rel[which(GridShare_Rel$COMUNA==CN_GSR),BTe] <- 
        c(st_drop_geometry(GridShare_Rel[which(GridShare_Rel$COMUNA==CN_GSR),BTe])) / 
        TempZW_GSR[which(TempZW_GSR$building_type==BTe),"avg_floor_area_building..m2."]
    }
  }
}
GridShare_Rel[is.na(GridShare_Rel)] <- 0

# compute the total number of buildings grid-wise
GridShare_Rel$Total_p_LvL <- rowSums(st_drop_geometry(GridShare_Rel[,BT]))

# export shapefile
GridShare_Rel_sp <- as_Spatial(GridShare_Rel)
writeOGR(obj=GridShare_Rel_sp,
         dsn=paste0(path_abbr,"/4_1_Exposure"),
         layer=paste0("GridShare_Rel","_ExpM"),
         driver="ESRI Shapefile",
         overwrite_layer = T)

# 7. WSF ####

# load WSF3D building height
WSF3D_bldgheight <- raster("C:/Users/admin/Desktop/GitHub/summer2023/WSF3D/WSF3D_V02_BuildingHeight_clipped_chile_nodata0.tif") 

# initialize the re-classification of WSF pixels into Class H categories
Hight_Thresh_GridShare_WSF <- data.frame(matrix(nrow = 0, ncol = 1 + length(HBET_Class_Order)))
colnames(Hight_Thresh_GridShare_WSF) <- c("COMUNA",HBET_Class_Order)
TempZW_Schablone_WSF <- Hight_Thresh_GridShare_WSF

# initialize the GSR list
Raster_HBET_List_WSF <- list() 

# get a list of unique comuna labels
ComunaNr_WSF <- as.character(unique(Intersection_Grid$COMUNA))

# set grid projection to use nDSM projection
Grid_sp <- as_Spatial(Grid)
Grid_spTr <- spTransform(Grid_sp, crs(WSF3D_bldgheight))

# initialize list containing clipped raster images
Raster_Clip_List_WSF <- list()

# label the raster clip using the grid ID
for (i in 1:nrow(Grid_spTr)) {
  Raster_Clip   <- crop(WSF3D_bldgheight,Grid_spTr[i,])
  Raster_Clip_List_WSF[[i]] <- Raster_Clip
  name <- paste("Raster_Clip", Grid[i,"GridNr"], sep = "_")
  names(Raster_Clip_List_WSF)[i] <- name
}

# start for-loop for each communa i
for (i in ComunaNr_WSF) {
  
  # for given communa i, get all grid IDs
  GridNr_fl_WSF <- Intersection_Grid_OG[which(Intersection_Grid_OG$COMUNA == i),"GridNr"]
  
  # initialize arrays
  Pixel_p_COM_WSF <- numeric()
  remove_zero_WSF <- c(0)
  
  # for-loop to combine all pixels belonging to communa i
  for (l in GridNr_fl_WSF) {
    Pixel_p_Grid_WSF <- Raster_Clip_List_WSF[[paste("Raster_Clip_",l, sep="")]]@data@values    
    Pixel_p_Grid_WSF <- Pixel_p_Grid_WSF [! Pixel_p_Grid_WSF %in% remove_zero_WSF]
    Pixel_p_Grid_WSF[Pixel_p_Grid_WSF == 0] <- NA
    Pixel_p_Grid_WSF <- na.omit(Pixel_p_Grid_WSF)
    Pixel_p_COM_WSF <- append(Pixel_p_COM_WSF, Pixel_p_Grid_WSF)
  }
  Pixel_p_COM_WSF <- as.data.frame(Pixel_p_COM_WSF)
  setorder(Pixel_p_COM_WSF)
  
  # add cumulative count percentage
  Pixel_p_COM_WSF$Perc <-  100 / nrow(Pixel_p_COM_WSF)
  Pixel_p_COM_WSF$Perc_Cum <- cumsum(Pixel_p_COM_WSF$Perc)
  
  # for communa i, extract and sort the proportion of each Class H category
  HBET_Perc_COM_WSF <- ExpM_A_HBET_COM[which(ExpM_A_HBET_COM$COMUNA == i),]
  HBET_Perc_COM_WSF <- HBET_Perc_COM_WSF[match(HBET_Class_Order, HBET_Perc_COM_WSF$ClassH),]
  HBET_Perc_COM_WSF$Area_HBET_p_COM_Perc_Cum <- cumsum(HBET_Perc_COM_WSF$Area_HBET_p_COM_Perc)
  
  # identify the only available Class H category from the census records given communa i
  HBET_forloop_WSF <- HBET_Perc_COM_WSF[which(HBET_Perc_COM_WSF$Area_p_HBET_COM != 0),"ClassH"]
  
  # add NA entries and communa info to the initialized array 
  TempZW_Schablone_WSF[1,] <- NA
  TempZW_Schablone_WSF[1,"COMUNA"] <- i
  
  # find the interval within the pixel count where the corresponding cumsum from census records falls
  for (h in length(HBET_forloop_WSF):1) {
    
    # get the index where the percentage cumsum appears
    Interval_WSF <- findInterval(
      HBET_Perc_COM_WSF[which(HBET_Perc_COM_WSF$ClassH == HBET_forloop_WSF[h]),
                        "Area_HBET_p_COM_Perc_Cum"]*100,
      Pixel_p_COM_WSF$Perc_Cum)
    
    # extract the pixel value corresponding to the perc_cum intersection
    Pixel_Height_WSF <- Pixel_p_COM_WSF[Interval_WSF,"Pixel_p_COM_WSF"]
    
    # save the pixel value (or height)
    TempZW_Schablone_WSF[1,HBET_forloop_WSF[h]] <- Pixel_Height_WSF
    TempZW_Schablone_WSF[is.na(TempZW_Schablone_WSF)] <- 0
    
  }    
  
  # adjust if there is any overlapping limit between any two or more Class H categories
  TempZW_WSF <- 2:ncol(TempZW_Schablone_WSF)
  for (c in 2: ncol(TempZW_Schablone_WSF)) {
    if (any(TempZW_Schablone_WSF[1,TempZW_WSF[! TempZW_WSF %in% c]] %in% TempZW_Schablone_WSF[1,c])){
      TempZW_Schablone_WSF[1,c] <- TempZW_Schablone_WSF[1,c] - 0.1
      for (c in 2: ncol(TempZW_Schablone_WSF)) {
        if (any(TempZW_Schablone_WSF[1,TempZW_WSF[! TempZW_WSF %in% c]] %in% TempZW_Schablone_WSF[1,c])){
          TempZW_Schablone_WSF[1,c] <- TempZW_Schablone_WSF[1,c] - 0.1
          for (c in 2: ncol(TempZW_Schablone_WSF)) {
            if (any(TempZW_Schablone_WSF[1,TempZW_WSF[! TempZW_WSF %in% c]] %in% TempZW_Schablone_WSF[1,c])){
              TempZW_Schablone_WSF[1,c] <- TempZW_Schablone_WSF[1,c] - 0.1
              for (c in 2: ncol(TempZW_Schablone_WSF)) {
                if (any(TempZW_Schablone_WSF[1,TempZW_WSF[! TempZW_WSF %in% c]] %in% TempZW_Schablone_WSF[1,c])){
                  TempZW_Schablone_WSF[1,c] <- TempZW_Schablone_WSF[1,c] - 0.1
                  
                } 
              }  
            } 
          }  
        } 
      }  
    } 
  }
  TempZW_Schablone_WSF[TempZW_Schablone_WSF<0] <- 0 
  
  # initialize classifying array
  classify_WSF <- numeric()
  First_WSF <- TempZW_Schablone_WSF[1,HBET_forloop_WSF[1]]
  classify_WSF <- append(classify_WSF,c(0,First_WSF,as.numeric(gsub("\\D", "",HBET_forloop_WSF[1])))) 
  
  # do the same if we have 2+ Class H categories
  if (length(HBET_forloop_WSF) >1) {
    for (z in 2:length(HBET_forloop_WSF)) {
      #z <- 2
      Second_1_WSF <- TempZW_Schablone_WSF[1,HBET_forloop_WSF[z-1]]  
      Second_2_WSF <- TempZW_Schablone_WSF[1,HBET_forloop_WSF[z]]
      classify_WSF <- append(classify_WSF,c(Second_1_WSF,
                                            Second_2_WSF,
                                            as.numeric(gsub("\\D", "",HBET_forloop_WSF[z]))))
    }
  }
  
  # postprocess
  classify_WSF[length(classify_WSF)-1] <- Inf     
  
  # for each grid ID in communa i,
  # reclassify each raster clip using the previously calculated reclassification table
  # based on perc_cum of pixel values of height and
  # census-derived proportion by communa and by Class H category
  for (m in GridNr_fl_WSF) {
    Raster_HBET_WSF <- reclassify(Raster_Clip_List_WSF[[paste("Raster_Clip_",m,sep="")]],
                                  rcl = classify_WSF)  
    Raster_HBET_List_WSF[[m]] <- Raster_HBET_WSF
    name <- paste("Raster_HBET", m,"GSR", sep = "_")    
    names(Raster_HBET_List_WSF)[m] <- name
  }
  
  Hight_Thresh_GridShare_WSF <- rbind(Hight_Thresh_GridShare_WSF,TempZW_Schablone_WSF)
  
} # end of for-loop for each communa i

# initialize array for the distribution of area to every building type
BT_Share_Perc_WSF <- setNames(data.frame(matrix(ncol = 9, nrow = 0)),
                              c("COMUNA","ClassH","HBET","Count","GridNr",
                                "Total_Pixel","Perc_Count_Tot_Pix","building_type","Area_BT_p_HBET_Perc"))

# identify clips with no pixels
If_RasterClip_noPixel_WSF <-  setNames(data.frame(matrix(ncol = 13, nrow = 0)),
                                       c("COMUNA","ClassH","HBET","Count","GridNr","Pixel_Count_ClassH_p_COM",
                                         "Perc_HBETPixel_p_COM","Area_p_HBET_COM","Area_p_HBET_Share",
                                         "building_type","Area_BT","Area_BT_p_HBET_Perc","Area_BT_WSF"))
No_Pixel_Table_WSF <- If_RasterClip_noPixel_WSF
If_RasterClip_noPixel_WSF[1,] <- NA

# start for-loop for each raster clip i
for (i in 1:length(Raster_HBET_List_WSF)) { 
  
  # get the pixel count and rename the column with HBET (because it's already re-classified)
  Pixel_Count_WSF <- as.data.frame(Raster_HBET_List_WSF[[i]])
  colnames(Pixel_Count_WSF) <- "HBET"
  
  # get summary statistics (sum of counts) for each HBET
  Pixel_Count_Share_WSF <- plyr::count(Pixel_Count_WSF, vars = "HBET")
  colnames(Pixel_Count_Share_WSF) <- c("HBET","Count")
  
  # note the raster clip i (or the grid ID, equivalently)
  Pixel_Count_Share_WSF$GridNr <- i
  
  # if there are any NA pixel values, add it to the "HBET = 0"
  if (any(is.na(Pixel_Count_Share_WSF$HBET))) {
    Pixel_Count_Share_WSF[1,"Count"] <-  Pixel_Count_Share_WSF[1,"Count"] + 
      Pixel_Count_Share_WSF[which(is.na(Pixel_Count_Share_WSF)),"Count"]
    Pixel_Count_Share_WSF <- na.omit(Pixel_Count_Share_WSF)
  }
  
  # 
  if (empty(Pixel_Count_Share_WSF)) {
    If_RasterClip_noPixel_WSF[1,"GridNr"]<- i
    TempZW_WSF <- Intersection_Grid
    st_geometry(TempZW_WSF) <- NULL
    If_RasterClip_noPixel_WSF$COMUNA <- TempZW_WSF[TempZW_WSF$GridNr == i,"COMUNA"]
    No_Pixel_Table_WSF <- rbind(No_Pixel_Table_WSF,If_RasterClip_noPixel_WSF)
    print(paste("No Pixels_WSF in Raster_Clip/GridNr:",i))} 
  
  # assign the string version of Class H category instead of numeric input
  else {
    Pixel_Count_Share_WSF$ClassH <- 0
    
    if (any(Pixel_Count_Share_WSF$HBET == 12)) {
      Pixel_Count_Share_WSF[Pixel_Count_Share_WSF$HBET == 12,]  <- 
        Pixel_Count_Share_WSF %>%  
        filter(Pixel_Count_Share_WSF[,"HBET"] == 12 ) %>% 
        mutate(ClassH = "HBET_1_2")  
    }
    
    if (any(Pixel_Count_Share_WSF$HBET == 3)) {
      Pixel_Count_Share_WSF[Pixel_Count_Share_WSF$HBET == 3,] <- 
        Pixel_Count_Share_WSF %>% 
        filter(Pixel_Count_Share_WSF[,"HBET"] == 3 ) %>%
        mutate(ClassH = "HBET_3")  
    }
    
    if (any(Pixel_Count_Share_WSF$HBET == 45)) {
      Pixel_Count_Share_WSF[Pixel_Count_Share_WSF$HBET == 45,] <- 
        Pixel_Count_Share_WSF %>%  
        filter(Pixel_Count_Share_WSF[,"HBET"] == 45 ) %>% 
        mutate(ClassH = "HBET_4_5")  
    }
    
    if (any(Pixel_Count_Share_WSF$HBET == 69)) {
      Pixel_Count_Share_WSF[Pixel_Count_Share_WSF$HBET == 69,] <- 
        Pixel_Count_Share_WSF %>%  
        filter(Pixel_Count_Share_WSF[,"HBET"] == 69 ) %>% 
        mutate(ClassH = "HBET_6_9")
    } 
    
    if (any(Pixel_Count_Share_WSF$HBET == 1024)) {
      Pixel_Count_Share_WSF[Pixel_Count_Share_WSF$HBET == 1024,] <- 
        Pixel_Count_Share_WSF %>%  
        filter(Pixel_Count_Share_WSF[,"HBET"] == 1024 ) %>% 
        mutate(ClassH = "HBET_10_24") 
    }
    
    if (any(Pixel_Count_Share_WSF$HBET == 2540)) {
      Pixel_Count_Share_WSF[Pixel_Count_Share_WSF$HBET == 2540,] <- 
        Pixel_Count_Share_WSF %>%  
        filter(Pixel_Count_Share_WSF[,"HBET"] == 2540 ) %>% 
        mutate(ClassH = "HBET_25_40") 
    }
    
    # add the communa information
    TempZW_WSF <- Intersection_Grid
    st_geometry(TempZW_WSF) <- NULL
    Pixel_Count_Share_WSF$COMUNA <- TempZW_WSF[TempZW_WSF$GridNr == i,"COMUNA"]
    
    # combine all results of for-loop
    BT_Share_Perc_WSF <- rbind(BT_Share_Perc_WSF,Pixel_Count_Share_WSF)
    
  }
  
} # end of for-loop for each raster clip i

# compute the total number of pixel counts by comuna and by Class H category
TempZW_WSF <- BT_Share_Perc_WSF %>% group_by(COMUNA,ClassH) %>% dplyr::summarise(Count=sum(Count))
colnames(TempZW_WSF) <- c("COMUNA","ClassH","Pixel_Count_ClassH_p_COM")

# add the previously computed to the for-loop array
PD_ZW_WSF <- merge(BT_Share_Perc_WSF, TempZW_WSF, by = c("COMUNA","ClassH"))

# compute the percentage based on pixel counts
PD_ZW_WSF$Perc_HBETPixel_p_COM <- PD_ZW_WSF$Count / PD_ZW_WSF$Pixel_Count_ClassH_p_COM

# add the area occupied (by communa and by Class H category) from census records
PD_ZW_WSF <- merge(PD_ZW_WSF, Area_p_HBET_COM, by = c("COMUNA", "ClassH"))

# multiply the pixel-based or satellite-derived percentage to the area to distribute it grid-wise
PD_ZW_WSF$Area_p_HBET_Share <- PD_ZW_WSF$Area_p_HBET_COM * PD_ZW_WSF$Perc_HBETPixel_p_COM

# add the information on proportion of each building type from census records
Pixeldensity_WSF <- merge(PD_ZW_WSF,
                          ExpM_A_BT_HBET[,c("COMUNA","ClassH","building_type",
                                            "Area_BT","Area_BT_p_HBET_Perc")],
                          by = c("COMUNA", "ClassH"))

# compute the area of each building type
Pixeldensity_WSF$Area_BT_WSF <- Pixeldensity_WSF$Area_p_HBET_Share * Pixeldensity_WSF$Area_BT_p_HBET_Perc

# combine the resulting array with no_pixel_data to be complete
Pixeldensity_WSF <- rbind(Pixeldensity_WSF,No_Pixel_Table_WSF) 

# vertically spread the area by building type
GridShare_Rel_WSF_Complete_Result <- Pixeldensity_WSF
GridShare_Rel_WSF_HBET_Result <- Pixeldensity_WSF[,c("COMUNA","ClassH","HBET","Count","GridNr",
                                                 "Pixel_Count_ClassH_p_COM","Perc_HBETPixel_p_COM",
                                                 "Area_p_HBET_COM","Area_p_HBET_Share",
                                                 "building_type","Area_BT_WSF")] %>% 
  spread(key=building_type,
         value=Area_BT_WSF)
GridShare_Rel_WSF <- Pixeldensity_WSF[,c("COMUNA","GridNr","building_type","Area_BT_WSF")] %>%
  spread(key=building_type,
         value=Area_BT_WSF)
GridShare_Rel_WSF <- merge(Intersection_Grid, GridShare_Rel_WSF, by = c("COMUNA","GridNr"))
GridShare_Rel_WSF[is.na(GridShare_Rel_WSF)] <- 0

# compute the total area for each building type
GridShare_Rel_WSF$HBET_1_2 <- rowSums(
  st_drop_geometry(GridShare_Rel_WSF[
    ,c("CR_1_3","MATO_1_2","MCF_MR_1_2","MCF_MCF_1_2",
       "MR_1_2","MURADO_1_2","MUR_1_2","W_1_3","WDNO_1_2")]))
GridShare_Rel_WSF$HBET_3 <- rowSums(
  st_drop_geometry(GridShare_Rel_WSF[
    ,c("MCF_MR_3","MCF_MCF_3","MR_3")]))
GridShare_Rel_WSF$HBET_4_5 <- rowSums(
  st_drop_geometry(GridShare_Rel_WSF[
    ,c("MCF_MR_4_5","MR_4_5","MCF_MCF_4_5")]))
GridShare_Rel_WSF$HBET_6_9 <- rowSums(
  st_drop_geometry(GridShare_Rel_WSF[
    ,c("CR_4_9")]))
GridShare_Rel_WSF$HBET_10_24 <- rowSums(
  st_drop_geometry(GridShare_Rel_WSF[
    ,c("CR_10_24")]))
GridShare_Rel_WSF$HBET_25_40 <- rowSums(
  st_drop_geometry(GridShare_Rel_WSF[
    ,c("CR_25_40")]))

# compute the number of buildings per building type
if (Calculate_Buildings == "YES") {
  ComunaNr_G500C_WSF <- unique(GridShare_Rel_WSF$COMUNA)
  for (CN_WSF in ComunaNr_G500C_WSF) {
    TempZW_WSF <- subset(Exposure_Model[,c("building_type","avg_floor_area_building..m2.","COMUNA")],
                         COMUNA == CN_WSF)
    for (BTe in BT) {
      GridShare_Rel_WSF[which(GridShare_Rel_WSF$COMUNA==CN_WSF),BTe] <- 
        c(st_drop_geometry(GridShare_Rel_WSF[which(GridShare_Rel_WSF$COMUNA==CN_WSF),BTe])) / 
        TempZW_WSF[which(TempZW_WSF$building_type==BTe),"avg_floor_area_building..m2."]
    }
  }
}
GridShare_Rel_WSF[is.na(GridShare_Rel_WSF)] <- 0

# compute the total number of buildings grid-wise
GridShare_Rel_WSF$Total_p_LvL <- rowSums(st_drop_geometry(GridShare_Rel_WSF[,BT]))

# export shapefile
GridShare_Rel_WSF_sp <- as_Spatial(GridShare_Rel_WSF)
writeOGR(obj=GridShare_Rel_WSF_sp,
         dsn=paste0(path_abbr,"/4_1_Exposure"),
         layer=paste0("GridShare_Rel_WSF","_ExpM"),
         driver="ESRI Shapefile",
         overwrite_layer = T)
# 8. HighRes ####

# load WSF3D building height
HighRes_bldgheight <- raster("C:/Users/admin/Desktop/GitHub/summer2023/Chile_EQ_EHVR_Priesmeier/3_SHAPES_INPUT/nDSM_SPOT_guf.tif") 

# initialize the re-classification of HighRes pixels into Class H categories
Hight_Thresh_GridShare_HighRes <- data.frame(matrix(nrow = 0, ncol = 1 + length(HBET_Class_Order)))
colnames(Hight_Thresh_GridShare_HighRes) <- c("COMUNA",HBET_Class_Order)
TempZW_Schablone_HighRes <- Hight_Thresh_GridShare_HighRes

# initialize the GSR list
Raster_HBET_List_HighRes <- list() 

# get a list of unique comuna labels
ComunaNr_HighRes <- as.character(unique(Intersection_Grid$COMUNA))

# set grid projection to use nDSM projection
Grid_sp <- as_Spatial(Grid)
Grid_spTr <- spTransform(Grid_sp, crs(HighRes_bldgheight))

# initialize list containing clipped raster images
Raster_Clip_List_HighRes <- list()

# label the raster clip using the grid ID
for (i in 1:nrow(Grid_spTr)) {
  Raster_Clip   <- crop(HighRes_bldgheight,Grid_spTr[i,])
  Raster_Clip_List_HighRes[[i]] <- Raster_Clip
  name <- paste("Raster_Clip", Grid[i,"GridNr"], sep = "_")
  names(Raster_Clip_List_HighRes)[i] <- name
}

# start for-loop for each communa i
for (i in ComunaNr_HighRes) {
  
  # for given communa i, get all grid IDs
  GridNr_fl_HighRes <- Intersection_Grid_OG[which(Intersection_Grid_OG$COMUNA == i),"GridNr"]
  
  # initialize arrays
  Pixel_p_COM_HighRes <- numeric()
  remove_zero_HighRes <- c(0)
  
  # for-loop to combine all pixels belonging to communa i
  for (l in GridNr_fl_HighRes) {
    Pixel_p_Grid_HighRes <- Raster_Clip_List_HighRes[[paste("Raster_Clip_",l, sep="")]]@data@values    
    Pixel_p_Grid_HighRes <- Pixel_p_Grid_HighRes [! Pixel_p_Grid_HighRes %in% remove_zero_HighRes]
    Pixel_p_Grid_HighRes[Pixel_p_Grid_HighRes == 0] <- NA
    Pixel_p_Grid_HighRes <- na.omit(Pixel_p_Grid_HighRes)
    Pixel_p_COM_HighRes <- append(Pixel_p_COM_HighRes, Pixel_p_Grid_HighRes)
  }
  Pixel_p_COM_HighRes <- as.data.frame(Pixel_p_COM_HighRes)
  setorder(Pixel_p_COM_HighRes)
  
  # add cumulative count percentage
  Pixel_p_COM_HighRes$Perc <-  100 / nrow(Pixel_p_COM_HighRes)
  Pixel_p_COM_HighRes$Perc_Cum <- cumsum(Pixel_p_COM_HighRes$Perc)
  
  # for communa i, extract and sort the proportion of each Class H category
  HBET_Perc_COM_HighRes <- ExpM_A_HBET_COM[which(ExpM_A_HBET_COM$COMUNA == i),]
  HBET_Perc_COM_HighRes <- HBET_Perc_COM_HighRes[match(HBET_Class_Order, HBET_Perc_COM_HighRes$ClassH),]
  HBET_Perc_COM_HighRes$Area_HBET_p_COM_Perc_Cum <- cumsum(HBET_Perc_COM_HighRes$Area_HBET_p_COM_Perc)
  
  # identify the only available Class H category from the census records given communa i
  HBET_forloop_HighRes <- HBET_Perc_COM_HighRes[which(HBET_Perc_COM_HighRes$Area_p_HBET_COM != 0),"ClassH"]
  
  # add NA entries and communa info to the initialized array 
  TempZW_Schablone_HighRes[1,] <- NA
  TempZW_Schablone_HighRes[1,"COMUNA"] <- i
  
  # find the interval within the pixel count where the corresponding cumsum from census records falls
  for (h in length(HBET_forloop_HighRes):1) {
    
    # get the index where the percentage cumsum appears
    Interval_HighRes <- findInterval(
      HBET_Perc_COM_HighRes[which(HBET_Perc_COM_HighRes$ClassH == HBET_forloop_HighRes[h]),
                        "Area_HBET_p_COM_Perc_Cum"]*100,
      Pixel_p_COM_HighRes$Perc_Cum)
    
    # extract the pixel value corresponding to the perc_cum intersection
    Pixel_Height_HighRes <- Pixel_p_COM_HighRes[Interval_HighRes,"Pixel_p_COM_HighRes"]
    
    # save the pixel value (or height)
    TempZW_Schablone_HighRes[1,HBET_forloop_HighRes[h]] <- Pixel_Height_HighRes
    TempZW_Schablone_HighRes[is.na(TempZW_Schablone_HighRes)] <- 0
    
  }    
  
  # adjust if there is any overlapping limit between any two or more Class H categories
  TempZW_HighRes <- 2:ncol(TempZW_Schablone_HighRes)
  for (c in 2: ncol(TempZW_Schablone_HighRes)) {
    if (any(TempZW_Schablone_HighRes[1,TempZW_HighRes[! TempZW_HighRes %in% c]] %in% TempZW_Schablone_HighRes[1,c])){
      TempZW_Schablone_HighRes[1,c] <- TempZW_Schablone_HighRes[1,c] - 0.1
      for (c in 2: ncol(TempZW_Schablone_HighRes)) {
        if (any(TempZW_Schablone_HighRes[1,TempZW_HighRes[! TempZW_HighRes %in% c]] %in% TempZW_Schablone_HighRes[1,c])){
          TempZW_Schablone_HighRes[1,c] <- TempZW_Schablone_HighRes[1,c] - 0.1
          for (c in 2: ncol(TempZW_Schablone_HighRes)) {
            if (any(TempZW_Schablone_HighRes[1,TempZW_HighRes[! TempZW_HighRes %in% c]] %in% TempZW_Schablone_HighRes[1,c])){
              TempZW_Schablone_HighRes[1,c] <- TempZW_Schablone_HighRes[1,c] - 0.1
              for (c in 2: ncol(TempZW_Schablone_HighRes)) {
                if (any(TempZW_Schablone_HighRes[1,TempZW_HighRes[! TempZW_HighRes %in% c]] %in% TempZW_Schablone_HighRes[1,c])){
                  TempZW_Schablone_HighRes[1,c] <- TempZW_Schablone_HighRes[1,c] - 0.1
                  
                } 
              }  
            } 
          }  
        } 
      }  
    } 
  }
  TempZW_Schablone_HighRes[TempZW_Schablone_HighRes<0] <- 0 
  
  # initialize classifying array
  classify_HighRes <- numeric()
  First_HighRes <- TempZW_Schablone_HighRes[1,HBET_forloop_HighRes[1]]
  classify_HighRes <- append(classify_HighRes,c(0,First_HighRes,as.numeric(gsub("\\D", "",HBET_forloop_HighRes[1])))) 
  
  # do the same if we have 2+ Class H categories
  if (length(HBET_forloop_HighRes) >1) {
    for (z in 2:length(HBET_forloop_HighRes)) {
      #z <- 2
      Second_1_HighRes <- TempZW_Schablone_HighRes[1,HBET_forloop_HighRes[z-1]]  
      Second_2_HighRes <- TempZW_Schablone_HighRes[1,HBET_forloop_HighRes[z]]
      classify_HighRes <- append(classify_HighRes,c(Second_1_HighRes,
                                            Second_2_HighRes,
                                            as.numeric(gsub("\\D", "",HBET_forloop_HighRes[z]))))
    }
  }
  
  # postprocess
  classify_HighRes[length(classify_HighRes)-1] <- Inf     
  
  # for each grid ID in communa i,
  # reclassify each raster clip using the previously calculated reclassification table
  # based on perc_cum of pixel values of height and
  # census-derived proportion by communa and by Class H category
  for (m in GridNr_fl_HighRes) {
    Raster_HBET_HighRes <- reclassify(Raster_Clip_List_HighRes[[paste("Raster_Clip_",m,sep="")]],
                                  rcl = classify_HighRes)  
    Raster_HBET_List_HighRes[[m]] <- Raster_HBET_HighRes
    name <- paste("Raster_HBET", m,"GSR", sep = "_")    
    names(Raster_HBET_List_HighRes)[m] <- name
  }
  
  Hight_Thresh_GridShare_HighRes <- rbind(Hight_Thresh_GridShare_HighRes,TempZW_Schablone_HighRes)
  
} # end of for-loop for each communa i

# initialize array for the distribution of area to every building type
BT_Share_Perc_HighRes <- setNames(data.frame(matrix(ncol = 9, nrow = 0)),
                              c("COMUNA","ClassH","HBET","Count","GridNr",
                                "Total_Pixel","Perc_Count_Tot_Pix","building_type","Area_BT_p_HBET_Perc"))

# identify clips with no pixels
If_RasterClip_noPixel_HighRes <-  setNames(data.frame(matrix(ncol = 13, nrow = 0)),
                                       c("COMUNA","ClassH","HBET","Count","GridNr","Pixel_Count_ClassH_p_COM",
                                         "Perc_HBETPixel_p_COM","Area_p_HBET_COM","Area_p_HBET_Share",
                                         "building_type","Area_BT","Area_BT_p_HBET_Perc","Area_BT_HighRes"))
No_Pixel_Table_HighRes <- If_RasterClip_noPixel_HighRes
If_RasterClip_noPixel_HighRes[1,] <- NA

# start for-loop for each raster clip i
for (i in 1:length(Raster_HBET_List_HighRes)) { 
  
  # get the pixel count and rename the column with HBET (because it's already re-classified)
  Pixel_Count_HighRes <- as.data.frame(Raster_HBET_List_HighRes[[i]])
  colnames(Pixel_Count_HighRes) <- "HBET"
  
  # get summary statistics (sum of counts) for each HBET
  Pixel_Count_Share_HighRes <- plyr::count(Pixel_Count_HighRes, vars = "HBET")
  colnames(Pixel_Count_Share_HighRes) <- c("HBET","Count")
  
  # note the raster clip i (or the grid ID, equivalently)
  Pixel_Count_Share_HighRes$GridNr <- i
  
  # if there are any NA pixel values, add it to the "HBET = 0"
  if (any(is.na(Pixel_Count_Share_HighRes$HBET))) {
    Pixel_Count_Share_HighRes[1,"Count"] <-  Pixel_Count_Share_HighRes[1,"Count"] + 
      Pixel_Count_Share_HighRes[which(is.na(Pixel_Count_Share_HighRes)),"Count"]
    Pixel_Count_Share_HighRes <- na.omit(Pixel_Count_Share_HighRes)
  }
  
  # 
  if (empty(Pixel_Count_Share_HighRes)) {
    If_RasterClip_noPixel_HighRes[1,"GridNr"]<- i
    TempZW_HighRes <- Intersection_Grid
    st_geometry(TempZW_HighRes) <- NULL
    If_RasterClip_noPixel_HighRes$COMUNA <- TempZW_HighRes[TempZW_HighRes$GridNr == i,"COMUNA"]
    No_Pixel_Table_HighRes <- rbind(No_Pixel_Table_HighRes,If_RasterClip_noPixel_HighRes)
    print(paste("No Pixels_HighRes in Raster_Clip/GridNr:",i))} 
  
  # assign the string version of Class H category instead of numeric input
  else {
    Pixel_Count_Share_HighRes$ClassH <- 0
    
    if (any(Pixel_Count_Share_HighRes$HBET == 12)) {
      Pixel_Count_Share_HighRes[Pixel_Count_Share_HighRes$HBET == 12,]  <- 
        Pixel_Count_Share_HighRes %>%  
        filter(Pixel_Count_Share_HighRes[,"HBET"] == 12 ) %>% 
        mutate(ClassH = "HBET_1_2")  
    }
    
    if (any(Pixel_Count_Share_HighRes$HBET == 3)) {
      Pixel_Count_Share_HighRes[Pixel_Count_Share_HighRes$HBET == 3,] <- 
        Pixel_Count_Share_HighRes %>% 
        filter(Pixel_Count_Share_HighRes[,"HBET"] == 3 ) %>%
        mutate(ClassH = "HBET_3")  
    }
    
    if (any(Pixel_Count_Share_HighRes$HBET == 45)) {
      Pixel_Count_Share_HighRes[Pixel_Count_Share_HighRes$HBET == 45,] <- 
        Pixel_Count_Share_HighRes %>%  
        filter(Pixel_Count_Share_HighRes[,"HBET"] == 45 ) %>% 
        mutate(ClassH = "HBET_4_5")  
    }
    
    if (any(Pixel_Count_Share_HighRes$HBET == 69)) {
      Pixel_Count_Share_HighRes[Pixel_Count_Share_HighRes$HBET == 69,] <- 
        Pixel_Count_Share_HighRes %>%  
        filter(Pixel_Count_Share_HighRes[,"HBET"] == 69 ) %>% 
        mutate(ClassH = "HBET_6_9")
    } 
    
    if (any(Pixel_Count_Share_HighRes$HBET == 1024)) {
      Pixel_Count_Share_HighRes[Pixel_Count_Share_HighRes$HBET == 1024,] <- 
        Pixel_Count_Share_HighRes %>%  
        filter(Pixel_Count_Share_HighRes[,"HBET"] == 1024 ) %>% 
        mutate(ClassH = "HBET_10_24") 
    }
    
    if (any(Pixel_Count_Share_HighRes$HBET == 2540)) {
      Pixel_Count_Share_HighRes[Pixel_Count_Share_HighRes$HBET == 2540,] <- 
        Pixel_Count_Share_HighRes %>%  
        filter(Pixel_Count_Share_HighRes[,"HBET"] == 2540 ) %>% 
        mutate(ClassH = "HBET_25_40") 
    }
    
    # add the communa information
    TempZW_HighRes <- Intersection_Grid
    st_geometry(TempZW_HighRes) <- NULL
    Pixel_Count_Share_HighRes$COMUNA <- TempZW_HighRes[TempZW_HighRes$GridNr == i,"COMUNA"]
    
    # combine all results of for-loop
    BT_Share_Perc_HighRes <- rbind(BT_Share_Perc_HighRes,Pixel_Count_Share_HighRes)
    
  }
  
} # end of for-loop for each raster clip i

# compute the total number of pixel counts by comuna and by Class H category
TempZW_HighRes <- BT_Share_Perc_HighRes %>% group_by(COMUNA,ClassH) %>% dplyr::summarise(Count=sum(Count))
colnames(TempZW_HighRes) <- c("COMUNA","ClassH","Pixel_Count_ClassH_p_COM")

# add the previously computed to the for-loop array
PD_ZW_HighRes <- merge(BT_Share_Perc_HighRes, TempZW_HighRes, by = c("COMUNA","ClassH"))

# compute the percentage based on pixel counts
PD_ZW_HighRes$Perc_HBETPixel_p_COM <- PD_ZW_HighRes$Count / PD_ZW_HighRes$Pixel_Count_ClassH_p_COM

# add the area occupied (by communa and by Class H category) from census records
PD_ZW_HighRes <- merge(PD_ZW_HighRes, Area_p_HBET_COM, by = c("COMUNA", "ClassH"))

# multiply the pixel-based or satellite-derived percentage to the area to distribute it grid-wise
PD_ZW_HighRes$Area_p_HBET_Share <- PD_ZW_HighRes$Area_p_HBET_COM * PD_ZW_HighRes$Perc_HBETPixel_p_COM

# add the information on proportion of each building type from census records
Pixeldensity_HighRes <- merge(PD_ZW_HighRes,
                          ExpM_A_BT_HBET[,c("COMUNA","ClassH","building_type",
                                            "Area_BT","Area_BT_p_HBET_Perc")],
                          by = c("COMUNA", "ClassH"))

# compute the area of each building type
Pixeldensity_HighRes$Area_BT_HighRes <- Pixeldensity_HighRes$Area_p_HBET_Share * Pixeldensity_HighRes$Area_BT_p_HBET_Perc

# combine the resulting array with no_pixel_data to be complete
Pixeldensity_HighRes <- rbind(Pixeldensity_HighRes,No_Pixel_Table_HighRes) 

# vertically spread the area by building type
GridShare_Rel_HighRes_Complete_Result <- Pixeldensity_HighRes
GridShare_Rel_HighRes_HBET_Result <- Pixeldensity_HighRes[,c("COMUNA","ClassH","HBET","Count","GridNr",
                                                     "Pixel_Count_ClassH_p_COM","Perc_HBETPixel_p_COM",
                                                     "Area_p_HBET_COM","Area_p_HBET_Share",
                                                     "building_type","Area_BT_HighRes")] %>% 
  spread(key=building_type,
         value=Area_BT_HighRes)
GridShare_Rel_HighRes <- Pixeldensity_HighRes[,c("COMUNA","GridNr","building_type","Area_BT_HighRes")] %>%
  spread(key=building_type,
         value=Area_BT_HighRes)
GridShare_Rel_HighRes <- merge(Intersection_Grid, GridShare_Rel_HighRes, by = c("COMUNA","GridNr"))
GridShare_Rel_HighRes[is.na(GridShare_Rel_HighRes)] <- 0

# compute the total area for each building type
GridShare_Rel_HighRes$HBET_1_2 <- rowSums(
  st_drop_geometry(GridShare_Rel_HighRes[
    ,c("CR_1_3","MATO_1_2","MCF_MR_1_2","MCF_MCF_1_2",
       "MR_1_2","MURADO_1_2","MUR_1_2","W_1_3","WDNO_1_2")]))
GridShare_Rel_HighRes$HBET_3 <- rowSums(
  st_drop_geometry(GridShare_Rel_HighRes[
    ,c("MCF_MR_3","MCF_MCF_3","MR_3")]))
GridShare_Rel_HighRes$HBET_4_5 <- rowSums(
  st_drop_geometry(GridShare_Rel_HighRes[
    ,c("MCF_MR_4_5","MR_4_5","MCF_MCF_4_5")]))
GridShare_Rel_HighRes$HBET_6_9 <- rowSums(
  st_drop_geometry(GridShare_Rel_HighRes[
    ,c("CR_4_9")]))
GridShare_Rel_HighRes$HBET_10_24 <- rowSums(
  st_drop_geometry(GridShare_Rel_HighRes[
    ,c("CR_10_24")]))
GridShare_Rel_HighRes$HBET_25_40 <- rowSums(
  st_drop_geometry(GridShare_Rel_HighRes[
    ,c("CR_25_40")]))

# compute the number of buildings per building type
if (Calculate_Buildings == "YES") {
  ComunaNr_G500C_HighRes <- unique(GridShare_Rel_HighRes$COMUNA)
  for (CN_HighRes in ComunaNr_G500C_HighRes) {
    TempZW_HighRes <- subset(Exposure_Model[,c("building_type","avg_floor_area_building..m2.","COMUNA")],
                         COMUNA == CN_HighRes)
    for (BTe in BT) {
      GridShare_Rel_HighRes[which(GridShare_Rel_HighRes$COMUNA==CN_HighRes),BTe] <- 
        c(st_drop_geometry(GridShare_Rel_HighRes[which(GridShare_Rel_HighRes$COMUNA==CN_HighRes),BTe])) / 
        TempZW_HighRes[which(TempZW_HighRes$building_type==BTe),"avg_floor_area_building..m2."]
    }
  }
}
GridShare_Rel_HighRes[is.na(GridShare_Rel_HighRes)] <- 0

# compute the total number of buildings grid-wise
GridShare_Rel_HighRes$Total_p_LvL <- rowSums(st_drop_geometry(GridShare_Rel_HighRes[,BT]))

# export shapefile
GridShare_Rel_HighRes_sp <- as_Spatial(GridShare_Rel_HighRes)
writeOGR(obj=GridShare_Rel_HighRes_sp,
         dsn=paste0(path_abbr,"/4_1_Exposure"),
         layer=paste0("GridShare_Rel_WSF","_ExpM"),
         driver="ESRI Shapefile",
         overwrite_layer = T)