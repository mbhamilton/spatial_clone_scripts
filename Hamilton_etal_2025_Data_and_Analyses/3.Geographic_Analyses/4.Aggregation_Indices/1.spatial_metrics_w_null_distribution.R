#########################################
#                                       #
#   Title: Spatial Metrics and Null     #
#   Simulation, w/ Plots                #
#   Authors: Charli D. Minsavage-Davis, #
#   Matthew B. Hamilton, and Sean Cole  #
#   Date last edited: 26 June 2025      #
#   Manuscript: S. patens landscape     # 
#   genetics - Long Island & NYC        #
#                                       #
#########################################

##requisite libraries if necessary
##if you do not have them installed, simply install.packages() or use the tool
##Some of these are likely unnecessary!!!
library (readxl)
library (Rcpp)
library (sp)
library (raster)
library (rstudioapi)
library (ggpubr)
library (nlme)
library (MASS)
library (rcompanion)
library (VGAM)
library (vegan)
library (ggplot2)
library (dplyr)
library (tidyr)
library (forcats)
library (hrbrthemes)
library (RColorBrewer)
library (car)
library (Rmisc)
library (Hmisc)
library (boot)
library (car)
library (nlme)
library (EnvStats)
library (emmeans)
library (utils)
library (multcomp)
library (viridis)
library (forecast)
library (grid)
library (landscapemetrics)
library (tidyverse)
library (tidybayes) 
library (brms)
library (cowplot)

#Getting the path of your current open file
current_path <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
print(getwd())

###############################################
#FIND AGGREGATION INDICIES FOR ALL POPULATIONS#
###############################################

##Read in data from all populations
PB <- read.csv("2.PB_Patch_Matrices.csv", header = FALSE)
CS <- read.csv("3.CS_Patch_Matrices.csv", header = FALSE)
JA <- read.csv("4.JA_Patch_Matrices.csv", header = FALSE)
WR <- read.csv("5.WR_Patch_Matrices.csv", header = FALSE)
SA <- read.csv("6.SA_Patch_Matrices.csv", header = FALSE)
WT <- read.csv("7.WT_Patch_Matrices.csv", header = FALSE)
AC <- read.csv("8.AC_Patch_Matrices.csv", header = FALSE)
JB <- read.csv("9.JB_Patch_Matrices.csv", header = FALSE)
GI <- read.csv("10.GI_Patch_Matrices.csv", header = FALSE)
CB <- read.csv("11.CB_Patch_Matrices.csv", header = FALSE)

#Segment data for each separate patch
PB_A <- PB[2:5,2:8]
PB_B <- PB[8:11,2:4]
PB_C <- PB[14:17,2:11]
PB_D <- PB[20:23,2:6]
PB_E <- PB[26:29,2:9]
 
CS_A <- CS[2:5,2:9]
CS_B <- CS[8:11,2:9]
CS_C <- CS[14:17,2:5]
CS_D <- CS[20:23,2:6]
CS_E <- CS[26:29,2:5]

JA_A <- JA[2:5,2:7]
JA_B <- JA[8:11,2:5]
JA_C <- JA[14:17,2:6]
JA_D <- JA[20:23,2:7]
JA_E <- JA[26:29,2:5]

WR_A <- WR[2:5,2:5]
WR_B <- WR[8:11,2:5]
WR_C <- WR[14:17,2:5]
WR_D <- WR[20:23,2:6]
WR_E <- WR[26:29,2:3]

SA_A <- SA[2:5,2:6]
SA_B <- SA[8:11,2:6]
SA_C <- SA[14:17,2:6]
SA_D <- SA[20:23,2:6]
SA_E <- SA[26:29,2:7]

WT_A <- WT[2:5,2:7]
WT_B <- WT[8:11,2:8]
WT_C <- WT[14:17,2:5]
WT_D <- WT[20:23,2:8]
WT_E <- WT[26:29,2:8]

AC_A <- AC[2:5,2:8]
AC_B <- AC[8:11,2:8]
AC_C <- AC[14:17,2:8]
AC_D <- AC[20:23,2:5]
AC_E <- AC[26:29,2:6]

JB_A <- JB[2:5,2:6]
JB_B <- JB[8:11,2:7]
JB_C <- JB[14:17,2:6]
JB_D <- JB[20:23,2:5]
JB_E <- JB[26:29,2:7]

GI_A <- GI[2:5,2:7]
GI_B <- GI[8:11,2:9]
GI_C <- GI[14:17,2:8]
GI_D <- GI[20:23,2:7]
GI_E <- GI[26:29,2:8]

CB_A <- CB[2:5,2:7]
CB_B <- CB[8:11,2:7]
CB_C <- CB[14:17,2:5]
CB_D <- CB[20:23,2:6]
CB_E <- CB[26:29,2:7]

##Create Rasters from Above
#PB
PB_A_list <- list()
PB_A_list$x = seq(1, by = 1, len = nrow(PB_A))
PB_A_list$y = seq(1, by = 1, len = ncol(PB_A))
PB_A_list$z = as.matrix(PB_A)

PB_A_rast <- raster(PB_A_list)
plot(PB_A_rast)

PB_B_list <- list()
PB_B_list$x = seq(1, by = 1, len = nrow(PB_B))
PB_B_list$y = seq(1, by = 1, len = ncol(PB_B))
PB_B_list$z = as.matrix(PB_B)

PB_B_rast <- raster(PB_B_list)
plot(PB_B_rast)

PB_C_list <- list()
PB_C_list$x = seq(1, by = 1, len = nrow(PB_C))
PB_C_list$y = seq(1, by = 1, len = ncol(PB_C))
PB_C_list$z = as.matrix(PB_C)

PB_C_rast <- raster(PB_C_list)
plot(PB_C_rast)

PB_D_list <- list()
PB_D_list$x = seq(1, by = 1, len = nrow(PB_D))
PB_D_list$y = seq(1, by = 1, len = ncol(PB_D))
PB_D_list$z = as.matrix(PB_D)

PB_D_rast <- raster(PB_D_list)
plot(PB_D_rast)

PB_E_list <- list()
PB_E_list$x = seq(1, by = 1, len = nrow(PB_E))
PB_E_list$y = seq(1, by = 1, len = ncol(PB_E))
PB_E_list$z = as.matrix(PB_E)

PB_E_rast <- raster(PB_E_list)
plot(PB_E_rast)

#CS
CS_A_list <- list()
CS_A_list$x = seq(1, by = 1, len = nrow(CS_A))
CS_A_list$y = seq(1, by = 1, len = ncol(CS_A))
CS_A_list$z = as.matrix(CS_A)

CS_A_rast <- raster(CS_A_list)
plot(CS_A_rast)

CS_B_list <- list()
CS_B_list$x = seq(1, by = 1, len = nrow(CS_B))
CS_B_list$y = seq(1, by = 1, len = ncol(CS_B))
CS_B_list$z = as.matrix(CS_B)

CS_B_rast <- raster(CS_B_list)
plot(CS_B_rast)

CS_C_list <- list()
CS_C_list$x = seq(1, by = 1, len = nrow(CS_C))
CS_C_list$y = seq(1, by = 1, len = ncol(CS_C))
CS_C_list$z = as.matrix(CS_C)

CS_C_rast <- raster(CS_C_list)
plot(CS_C_rast)

CS_D_list <- list()
CS_D_list$x = seq(1, by = 1, len = nrow(CS_D))
CS_D_list$y = seq(1, by = 1, len = ncol(CS_D))
CS_D_list$z = as.matrix(CS_D)

CS_D_rast <- raster(CS_D_list)
plot(CS_D_rast)

CS_E_list <- list()
CS_E_list$x = seq(1, by = 1, len = nrow(CS_E))
CS_E_list$y = seq(1, by = 1, len = ncol(CS_E))
CS_E_list$z = as.matrix(CS_E)

CS_E_rast <- raster(CS_E_list)
plot(CS_E_rast)

#JA
JA_A_list <- list()
JA_A_list$x = seq(1, by = 1, len = nrow(JA_A))
JA_A_list$y = seq(1, by = 1, len = ncol(JA_A))
JA_A_list$z = as.matrix(JA_A)

JA_A_rast <- raster(JA_A_list)
plot(JA_A_rast)

JA_B_list <- list()
JA_B_list$x = seq(1, by = 1, len = nrow(JA_B))
JA_B_list$y = seq(1, by = 1, len = ncol(JA_B))
JA_B_list$z = as.matrix(JA_B)

JA_B_rast <- raster(JA_B_list)
plot(JA_B_rast)

JA_C_list <- list()
JA_C_list$x = seq(1, by = 1, len = nrow(JA_C))
JA_C_list$y = seq(1, by = 1, len = ncol(JA_C))
JA_C_list$z = as.matrix(JA_C)

JA_C_rast <- raster(JA_C_list)
plot(JA_C_rast)

JA_D_list <- list()
JA_D_list$x = seq(1, by = 1, len = nrow(JA_D))
JA_D_list$y = seq(1, by = 1, len = ncol(JA_D))
JA_D_list$z = as.matrix(JA_D)

JA_D_rast <- raster(JA_D_list)
plot(JA_D_rast)

JA_E_list <- list()
JA_E_list$x = seq(1, by = 1, len = nrow(JA_E))
JA_E_list$y = seq(1, by = 1, len = ncol(JA_E))
JA_E_list$z = as.matrix(JA_E)

JA_E_rast <- raster(JA_E_list)
plot(JA_E_rast)

#WR
WR_A_list <- list()
WR_A_list$x = seq(1, by = 1, len = nrow(WR_A))
WR_A_list$y = seq(1, by = 1, len = ncol(WR_A))
WR_A_list$z = as.matrix(WR_A)

WR_A_rast <- raster(WR_A_list)
plot(WR_A_rast)

WR_B_list <- list()
WR_B_list$x = seq(1, by = 1, len = nrow(WR_B))
WR_B_list$y = seq(1, by = 1, len = ncol(WR_B))
WR_B_list$z = as.matrix(WR_B)

WR_B_rast <- raster(WR_B_list)
plot(WR_B_rast)

WR_C_list <- list()
WR_C_list$x = seq(1, by = 1, len = nrow(WR_C))
WR_C_list$y = seq(1, by = 1, len = ncol(WR_C))
WR_C_list$z = as.matrix(WR_C)

WR_C_rast <- raster(WR_C_list)
plot(WR_C_rast)

WR_D_list <- list()
WR_D_list$x = seq(1, by = 1, len = nrow(WR_D))
WR_D_list$y = seq(1, by = 1, len = ncol(WR_D))
WR_D_list$z = as.matrix(WR_D)

WR_D_rast <- raster(WR_D_list)
plot(WR_D_rast)

WR_E_list <- list()
WR_E_list$x = seq(1, by = 1, len = nrow(WR_E))
WR_E_list$y = seq(1, by = 1, len = ncol(WR_E))
WR_E_list$z = as.matrix(WR_E)

WR_E_rast <- raster(WR_E_list)
plot(WR_E_rast)

#SA
SA_A_list <- list()
SA_A_list$x = seq(1, by = 1, len = nrow(SA_A))
SA_A_list$y = seq(1, by = 1, len = ncol(SA_A))
SA_A_list$z = as.matrix(SA_A)

SA_A_rast <- raster(SA_A_list)
plot(SA_A_rast)

SA_B_list <- list()
SA_B_list$x = seq(1, by = 1, len = nrow(SA_B))
SA_B_list$y = seq(1, by = 1, len = ncol(SA_B))
SA_B_list$z = as.matrix(SA_B)

SA_B_rast <- raster(SA_B_list)
plot(SA_B_rast)

SA_C_list <- list()
SA_C_list$x = seq(1, by = 1, len = nrow(SA_C))
SA_C_list$y = seq(1, by = 1, len = ncol(SA_C))
SA_C_list$z = as.matrix(SA_C)

SA_C_rast <- raster(SA_C_list)
plot(SA_C_rast)

SA_D_list <- list()
SA_D_list$x = seq(1, by = 1, len = nrow(SA_D))
SA_D_list$y = seq(1, by = 1, len = ncol(SA_D))
SA_D_list$z = as.matrix(SA_D)

SA_D_rast <- raster(SA_D_list)
plot(SA_D_rast)

SA_E_list <- list()
SA_E_list$x = seq(1, by = 1, len = nrow(SA_E))
SA_E_list$y = seq(1, by = 1, len = ncol(SA_E))
SA_E_list$z = as.matrix(SA_E)

SA_E_rast <- raster(SA_E_list)
plot(SA_E_rast)

#WT
WT_A_list <- list()
WT_A_list$x = seq(1, by = 1, len = nrow(WT_A))
WT_A_list$y = seq(1, by = 1, len = ncol(WT_A))
WT_A_list$z = as.matrix(WT_A)

WT_A_rast <- raster(WT_A_list)
plot(WT_A_rast)

WT_B_list <- list()
WT_B_list$x = seq(1, by = 1, len = nrow(WT_B))
WT_B_list$y = seq(1, by = 1, len = ncol(WT_B))
WT_B_list$z = as.matrix(WT_B)

WT_B_rast <- raster(WT_B_list)
plot(WT_B_rast)

WT_C_list <- list()
WT_C_list$x = seq(1, by = 1, len = nrow(WT_C))
WT_C_list$y = seq(1, by = 1, len = ncol(WT_C))
WT_C_list$z = as.matrix(WT_C)

WT_C_rast <- raster(WT_C_list)
plot(WT_C_rast)

WT_D_list <- list()
WT_D_list$x = seq(1, by = 1, len = nrow(WT_D))
WT_D_list$y = seq(1, by = 1, len = ncol(WT_D))
WT_D_list$z = as.matrix(WT_D)

WT_D_rast <- raster(WT_D_list)
plot(WT_D_rast)

WT_E_list <- list()
WT_E_list$x = seq(1, by = 1, len = nrow(WT_E))
WT_E_list$y = seq(1, by = 1, len = ncol(WT_E))
WT_E_list$z = as.matrix(WT_E)

WT_E_rast <- raster(WT_E_list)
plot(WT_E_rast)

#AC
AC_A_list <- list()
AC_A_list$x = seq(1, by = 1, len = nrow(AC_A))
AC_A_list$y = seq(1, by = 1, len = ncol(AC_A))
AC_A_list$z = as.matrix(AC_A)

AC_A_rast <- raster(AC_A_list)
plot(AC_A_rast)

AC_B_list <- list()
AC_B_list$x = seq(1, by = 1, len = nrow(AC_B))
AC_B_list$y = seq(1, by = 1, len = ncol(AC_B))
AC_B_list$z = as.matrix(AC_B)

AC_B_rast <- raster(AC_B_list)
plot(AC_B_rast)

AC_C_list <- list()
AC_C_list$x = seq(1, by = 1, len = nrow(AC_C))
AC_C_list$y = seq(1, by = 1, len = ncol(AC_C))
AC_C_list$z = as.matrix(AC_C)

AC_C_rast <- raster(AC_C_list)
plot(AC_C_rast)

AC_D_list <- list()
AC_D_list$x = seq(1, by = 1, len = nrow(AC_D))
AC_D_list$y = seq(1, by = 1, len = ncol(AC_D))
AC_D_list$z = as.matrix(AC_D)

AC_D_rast <- raster(AC_D_list)
plot(AC_D_rast)

AC_E_list <- list()
AC_E_list$x = seq(1, by = 1, len = nrow(AC_E))
AC_E_list$y = seq(1, by = 1, len = ncol(AC_E))
AC_E_list$z = as.matrix(AC_E)

AC_E_rast <- raster(AC_E_list)
plot(AC_E_rast)

#JB
JB_A_list <- list()
JB_A_list$x = seq(1, by = 1, len = nrow(JB_A))
JB_A_list$y = seq(1, by = 1, len = ncol(JB_A))
JB_A_list$z = as.matrix(JB_A)

JB_A_rast <- raster(JB_A_list)
plot(JB_A_rast)

JB_B_list <- list()
JB_B_list$x = seq(1, by = 1, len = nrow(JB_B))
JB_B_list$y = seq(1, by = 1, len = ncol(JB_B))
JB_B_list$z = as.matrix(JB_B)

JB_B_rast <- raster(JB_B_list)
plot(JB_B_rast)

JB_C_list <- list()
JB_C_list$x = seq(1, by = 1, len = nrow(JB_C))
JB_C_list$y = seq(1, by = 1, len = ncol(JB_C))
JB_C_list$z = as.matrix(JB_C)

JB_C_rast <- raster(JB_C_list)
plot(JB_C_rast)

JB_D_list <- list()
JB_D_list$x = seq(1, by = 1, len = nrow(JB_D))
JB_D_list$y = seq(1, by = 1, len = ncol(JB_D))
JB_D_list$z = as.matrix(JB_D)

JB_D_rast <- raster(JB_D_list)
plot(JB_D_rast)

JB_E_list <- list()
JB_E_list$x = seq(1, by = 1, len = nrow(JB_E))
JB_E_list$y = seq(1, by = 1, len = ncol(JB_E))
JB_E_list$z = as.matrix(JB_E)

JB_E_rast <- raster(JB_E_list)
plot(JB_E_rast)

#GI
GI_A_list <- list()
GI_A_list$x = seq(1, by = 1, len = nrow(GI_A))
GI_A_list$y = seq(1, by = 1, len = ncol(GI_A))
GI_A_list$z = as.matrix(GI_A)

GI_A_rast <- raster(GI_A_list)
plot(GI_A_rast)

GI_B_list <- list()
GI_B_list$x = seq(1, by = 1, len = nrow(GI_B))
GI_B_list$y = seq(1, by = 1, len = ncol(GI_B))
GI_B_list$z = as.matrix(GI_B)

GI_B_rast <- raster(GI_B_list)
plot(GI_B_rast)

GI_C_list <- list()
GI_C_list$x = seq(1, by = 1, len = nrow(GI_C))
GI_C_list$y = seq(1, by = 1, len = ncol(GI_C))
GI_C_list$z = as.matrix(GI_C)

GI_C_rast <- raster(GI_C_list)
plot(GI_C_rast)

GI_D_list <- list()
GI_D_list$x = seq(1, by = 1, len = nrow(GI_D))
GI_D_list$y = seq(1, by = 1, len = ncol(GI_D))
GI_D_list$z = as.matrix(GI_D)

GI_D_rast <- raster(GI_D_list)
plot(GI_D_rast)

GI_E_list <- list()
GI_E_list$x = seq(1, by = 1, len = nrow(GI_E))
GI_E_list$y = seq(1, by = 1, len = ncol(GI_E))
GI_E_list$z = as.matrix(GI_E)

GI_E_rast <- raster(GI_E_list)
plot(GI_E_rast)

#CB
CB_A_list <- list()
CB_A_list$x = seq(1, by = 1, len = nrow(CB_A))
CB_A_list$y = seq(1, by = 1, len = ncol(CB_A))
CB_A_list$z = as.matrix(CB_A)

CB_A_rast <- raster(CB_A_list)
plot(CB_A_rast)

CB_B_list <- list()
CB_B_list$x = seq(1, by = 1, len = nrow(CB_B))
CB_B_list$y = seq(1, by = 1, len = ncol(CB_B))
CB_B_list$z = as.matrix(CB_B)

CB_B_rast <- raster(CB_B_list)
plot(CB_B_rast)

CB_C_list <- list()
CB_C_list$x = seq(1, by = 1, len = nrow(CB_C))
CB_C_list$y = seq(1, by = 1, len = ncol(CB_C))
CB_C_list$z = as.matrix(CB_C)

CB_C_rast <- raster(CB_C_list)
plot(CB_C_rast)

CB_D_list <- list()
CB_D_list$x = seq(1, by = 1, len = nrow(CB_D))
CB_D_list$y = seq(1, by = 1, len = ncol(CB_D))
CB_D_list$z = as.matrix(CB_D)

CB_D_rast <- raster(CB_D_list)
plot(CB_D_rast)

CB_E_list <- list()
CB_E_list$x = seq(1, by = 1, len = nrow(CB_E))
CB_E_list$y = seq(1, by = 1, len = ncol(CB_E))
CB_E_list$z = as.matrix(CB_E)

CB_E_rast <- raster(CB_E_list)
plot(CB_E_rast)

#Calculate cohesion index for each patch
coh_PB_A <- lsm_c_cohesion(PB_A_rast)
coh_PB_B <- lsm_c_cohesion(PB_B_rast)
coh_PB_C <- lsm_c_cohesion(PB_C_rast)
coh_PB_D <- lsm_c_cohesion(PB_D_rast)
coh_PB_E <- lsm_c_cohesion(PB_E_rast)

coh_CS_A <- lsm_c_cohesion(CS_A_rast)
coh_CS_B <- lsm_c_cohesion(CS_B_rast)
coh_CS_C <- lsm_c_cohesion(CS_C_rast)
coh_CS_D <- lsm_c_cohesion(CS_D_rast)
coh_CS_E <- lsm_c_cohesion(CS_E_rast)

coh_JA_A <- lsm_c_cohesion(JA_A_rast)
coh_JA_B <- lsm_c_cohesion(JA_B_rast)
coh_JA_C <- lsm_c_cohesion(JA_C_rast)
coh_JA_D <- lsm_c_cohesion(JA_D_rast)
coh_JA_E <- lsm_c_cohesion(JA_E_rast)

coh_WR_A <- lsm_c_cohesion(WR_A_rast)
coh_WR_B <- lsm_c_cohesion(WR_B_rast)
coh_WR_C <- lsm_c_cohesion(WR_C_rast)
coh_WR_D <- lsm_c_cohesion(WR_D_rast)
coh_WR_E <- lsm_c_cohesion(WR_E_rast)

coh_SA_A <- lsm_c_cohesion(SA_A_rast)
coh_SA_B <- lsm_c_cohesion(SA_B_rast)
coh_SA_C <- lsm_c_cohesion(SA_C_rast)
coh_SA_D <- lsm_c_cohesion(SA_D_rast)
coh_SA_E <- lsm_c_cohesion(SA_E_rast)

coh_WT_A <- lsm_c_cohesion(WT_A_rast)
coh_WT_B <- lsm_c_cohesion(WT_B_rast)
coh_WT_C <- lsm_c_cohesion(WT_C_rast)
coh_WT_D <- lsm_c_cohesion(WT_D_rast)
coh_WT_E <- lsm_c_cohesion(WT_E_rast)

coh_AC_A <- lsm_c_cohesion(AC_A_rast)
coh_AC_B <- lsm_c_cohesion(AC_B_rast)
coh_AC_C <- lsm_c_cohesion(AC_C_rast)
coh_AC_D <- lsm_c_cohesion(AC_D_rast)
coh_AC_E <- lsm_c_cohesion(AC_E_rast)

coh_JB_A <- lsm_c_cohesion(JB_A_rast)
coh_JB_B <- lsm_c_cohesion(JB_B_rast)
coh_JB_C <- lsm_c_cohesion(JB_C_rast)
coh_JB_D <- lsm_c_cohesion(JB_D_rast)
coh_JB_E <- lsm_c_cohesion(JB_E_rast)

coh_GI_A <- lsm_c_cohesion(GI_A_rast)
coh_GI_B <- lsm_c_cohesion(GI_B_rast)
coh_GI_C <- lsm_c_cohesion(GI_C_rast)
coh_GI_D <- lsm_c_cohesion(GI_D_rast)
coh_GI_E <- lsm_c_cohesion(GI_E_rast)

coh_CB_A <- lsm_c_cohesion(CB_A_rast)
coh_CB_B <- lsm_c_cohesion(CB_B_rast)
coh_CB_C <- lsm_c_cohesion(CB_C_rast)
coh_CB_D <- lsm_c_cohesion(CB_D_rast)
coh_CB_E <- lsm_c_cohesion(CB_E_rast)

#Calculate class area for each patch
ca_PB_A <- lsm_c_ca(PB_A_rast)
ca_PB_B <- lsm_c_ca(PB_B_rast)
ca_PB_C <- lsm_c_ca(PB_C_rast)
ca_PB_D <- lsm_c_ca(PB_D_rast)
ca_PB_E <- lsm_c_ca(PB_E_rast)

ca_CS_A <- lsm_c_ca(CS_A_rast)
ca_CS_B <- lsm_c_ca(CS_B_rast)
ca_CS_C <- lsm_c_ca(CS_C_rast)
ca_CS_D <- lsm_c_ca(CS_D_rast)
ca_CS_E <- lsm_c_ca(CS_E_rast)

ca_JA_A <- lsm_c_ca(JA_A_rast)
ca_JA_B <- lsm_c_ca(JA_B_rast)
ca_JA_C <- lsm_c_ca(JA_C_rast)
ca_JA_D <- lsm_c_ca(JA_D_rast)
ca_JA_E <- lsm_c_ca(JA_E_rast)

ca_WR_A <- lsm_c_ca(WR_A_rast)
ca_WR_B <- lsm_c_ca(WR_B_rast)
ca_WR_C <- lsm_c_ca(WR_C_rast)
ca_WR_D <- lsm_c_ca(WR_D_rast)
ca_WR_E <- lsm_c_ca(WR_E_rast)

ca_SA_A <- lsm_c_ca(SA_A_rast)
ca_SA_B <- lsm_c_ca(SA_B_rast)
ca_SA_C <- lsm_c_ca(SA_C_rast)
ca_SA_D <- lsm_c_ca(SA_D_rast)
ca_SA_E <- lsm_c_ca(SA_E_rast)

ca_WT_A <- lsm_c_ca(WT_A_rast)
ca_WT_B <- lsm_c_ca(WT_B_rast)
ca_WT_C <- lsm_c_ca(WT_C_rast)
ca_WT_D <- lsm_c_ca(WT_D_rast)
ca_WT_E <- lsm_c_ca(WT_E_rast)

ca_AC_A <- lsm_c_ca(AC_A_rast)
ca_AC_B <- lsm_c_ca(AC_B_rast)
ca_AC_C <- lsm_c_ca(AC_C_rast)
ca_AC_D <- lsm_c_ca(AC_D_rast)
ca_AC_E <- lsm_c_ca(AC_E_rast)

ca_JB_A <- lsm_c_ca(JB_A_rast)
ca_JB_B <- lsm_c_ca(JB_B_rast)
ca_JB_C <- lsm_c_ca(JB_C_rast)
ca_JB_D <- lsm_c_ca(JB_D_rast)
ca_JB_E <- lsm_c_ca(JB_E_rast)

ca_GI_A <- lsm_c_ca(GI_A_rast)
ca_GI_B <- lsm_c_ca(GI_B_rast)
ca_GI_C <- lsm_c_ca(GI_C_rast)
ca_GI_D <- lsm_c_ca(GI_D_rast)
ca_GI_E <- lsm_c_ca(GI_E_rast)

ca_CB_A <- lsm_c_ca(CB_A_rast)
ca_CB_B <- lsm_c_ca(CB_B_rast)
ca_CB_C <- lsm_c_ca(CB_C_rast)
ca_CB_D <- lsm_c_ca(CB_D_rast)
ca_CB_E <- lsm_c_ca(CB_E_rast)

#Calculate interespersion-juxtaposition index
iji_PB_A <- lsm_c_iji(PB_A_rast)
iji_PB_B <- lsm_c_iji(PB_B_rast)
iji_PB_C <- lsm_c_iji(PB_C_rast)
iji_PB_D <- lsm_c_iji(PB_D_rast)
iji_PB_E <- lsm_c_iji(PB_E_rast)

iji_CS_A <- lsm_c_iji(CS_A_rast)
iji_CS_B <- lsm_c_iji(CS_B_rast)
iji_CS_C <- lsm_c_iji(CS_C_rast)
iji_CS_D <- lsm_c_iji(CS_D_rast)
iji_CS_E <- lsm_c_iji(CS_E_rast)

iji_JA_A <- lsm_c_iji(JA_A_rast)
iji_JA_B <- lsm_c_iji(JA_B_rast)
iji_JA_C <- lsm_c_iji(JA_C_rast)
iji_JA_D <- lsm_c_iji(JA_D_rast)
iji_JA_E <- lsm_c_iji(JA_E_rast)

iji_WR_A <- lsm_c_iji(WR_A_rast)
iji_WR_B <- lsm_c_iji(WR_B_rast)
iji_WR_C <- lsm_c_iji(WR_C_rast)
iji_WR_D <- lsm_c_iji(WR_D_rast)
iji_WR_E <- lsm_c_iji(WR_E_rast)

iji_SA_A <- lsm_c_iji(SA_A_rast)
iji_SA_B <- lsm_c_iji(SA_B_rast)
iji_SA_C <- lsm_c_iji(SA_C_rast)
iji_SA_D <- lsm_c_iji(SA_D_rast)
iji_SA_E <- lsm_c_iji(SA_E_rast)

iji_WT_A <- lsm_c_iji(WT_A_rast)
iji_WT_B <- lsm_c_iji(WT_B_rast)
iji_WT_C <- lsm_c_iji(WT_C_rast)
iji_WT_D <- lsm_c_iji(WT_D_rast)
iji_WT_E <- lsm_c_iji(WT_E_rast)

iji_AC_A <- lsm_c_iji(AC_A_rast)
iji_AC_B <- lsm_c_iji(AC_B_rast)
iji_AC_C <- lsm_c_iji(AC_C_rast)
iji_AC_D <- lsm_c_iji(AC_D_rast)
iji_AC_E <- lsm_c_iji(AC_E_rast)

iji_JB_A <- lsm_c_iji(JB_A_rast)
iji_JB_B <- lsm_c_iji(JB_B_rast)
iji_JB_C <- lsm_c_iji(JB_C_rast)
iji_JB_D <- lsm_c_iji(JB_D_rast)
iji_JB_E <- lsm_c_iji(JB_E_rast)

iji_GI_A <- lsm_c_iji(GI_A_rast)
iji_GI_B <- lsm_c_iji(GI_B_rast)
iji_GI_C <- lsm_c_iji(GI_C_rast)
iji_GI_D <- lsm_c_iji(GI_D_rast)
iji_GI_E <- lsm_c_iji(GI_E_rast)

iji_CB_A <- lsm_c_iji(CB_A_rast)
iji_CB_B <- lsm_c_iji(CB_B_rast)
iji_CB_C <- lsm_c_iji(CB_C_rast)
iji_CB_D <- lsm_c_iji(CB_D_rast)
iji_CB_E <- lsm_c_iji(CB_E_rast)

#export all of the above as csvs
##This is for cohesion!
# coh_PB_A.f <- coh_PB_A[, c(1:3,6)]
# colnames(coh_PB_A.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_PB_A.f$LCP <- "PB"
# coh_PB_A.f$Patch <- "A"
# 
# coh_PB_B.f <- coh_PB_B[, c(1:3,6)]
# colnames(coh_PB_B.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_PB_B.f$LCP <- "PB"
# coh_PB_B.f$Patch <- "B"
# 
# coh_PB_C.f <- coh_PB_C[, c(1:3,6)]
# colnames(coh_PB_C.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_PB_C.f$LCP <- "PB"
# coh_PB_C.f$Patch <- "C"
# 
# coh_PB_D.f <- coh_PB_D[, c(1:3,6)]
# colnames(coh_PB_D.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_PB_D.f$LCP <- "PB"
# coh_PB_D.f$Patch <- "D"
# 
# coh_PB_E.f <- coh_PB_E[, c(1:3,6)]
# colnames(coh_PB_E.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_PB_E.f$LCP <- "PB"
# coh_PB_E.f$Patch <- "E"
# 
# coh_PB.f <- rbind(coh_PB_A.f, coh_PB_B.f, coh_PB_C.f, coh_PB_D.f, coh_PB_E.f)
# 
# coh_CS_A.f <- coh_CS_A[, c(1:3,6)]
# colnames(coh_CS_A.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_CS_A.f$LCP <- "CS"
# coh_CS_A.f$Patch <- "A"
# 
# coh_CS_B.f <- coh_CS_B[, c(1:3,6)]
# colnames(coh_CS_B.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_CS_B.f$LCP <- "CS"
# coh_CS_B.f$Patch <- "B"
# 
# coh_CS_C.f <- coh_CS_C[, c(1:3,6)]
# colnames(coh_CS_C.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_CS_C.f$LCP <- "CS"
# coh_CS_C.f$Patch <- "C"
# 
# coh_CS_D.f <- coh_CS_D[, c(1:3,6)]
# colnames(coh_CS_D.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_CS_D.f$LCP <- "CS"
# coh_CS_D.f$Patch <- "D"
# 
# coh_CS_E.f <- coh_CS_E[, c(1:3,6)]
# colnames(coh_CS_E.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_CS_E.f$LCP <- "CS"
# coh_CS_E.f$Patch <- "E"
# 
# coh_CS.f <- rbind(coh_CS_A.f, coh_CS_B.f, coh_CS_C.f, coh_CS_D.f, coh_CS_E.f)
# 
# coh_JA_A.f <- coh_JA_A[, c(1:3,6)]
# colnames(coh_JA_A.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_JA_A.f$LCP <- "JA"
# coh_JA_A.f$Patch <- "A"
# 
# coh_JA_B.f <- coh_JA_B[, c(1:3,6)]
# colnames(coh_JA_B.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_JA_B.f$LCP <- "JA"
# coh_JA_B.f$Patch <- "B"
# 
# coh_JA_C.f <- coh_JA_C[, c(1:3,6)]
# colnames(coh_JA_C.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_JA_C.f$LCP <- "JA"
# coh_JA_C.f$Patch <- "C"
# 
# coh_JA_D.f <- coh_JA_D[, c(1:3,6)]
# colnames(coh_JA_D.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_JA_D.f$LCP <- "JA"
# coh_JA_D.f$Patch <- "D"
# 
# coh_JA_E.f <- coh_JA_E[, c(1:3,6)]
# colnames(coh_JA_E.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_JA_E.f$LCP <- "JA"
# coh_JA_E.f$Patch <- "E"
# 
# coh_JA.f <- rbind(coh_JA_A.f, coh_JA_B.f, coh_JA_C.f, coh_JA_D.f, coh_JA_E.f)
# 
# coh_WR_A.f <- coh_WR_A[, c(1:3,6)]
# colnames(coh_WR_A.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_WR_A.f$LCP <- "WR"
# coh_WR_A.f$Patch <- "A"
# 
# coh_WR_B.f <- coh_WR_B[, c(1:3,6)]
# colnames(coh_WR_B.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_WR_B.f$LCP <- "WR"
# coh_WR_B.f$Patch <- "B"
# 
# coh_WR_C.f <- coh_WR_C[, c(1:3,6)]
# colnames(coh_WR_C.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_WR_C.f$LCP <- "WR"
# coh_WR_C.f$Patch <- "C"
# 
# coh_WR_D.f <- coh_WR_D[, c(1:3,6)]
# colnames(coh_WR_D.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_WR_D.f$LCP <- "WR"
# coh_WR_D.f$Patch <- "D"
# 
# coh_WR_E.f <- coh_WR_E[, c(1:3,6)]
# colnames(coh_WR_E.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_WR_E.f$LCP <- "WR"
# coh_WR_E.f$Patch <- "E"
# 
# coh_WR.f <- rbind(coh_WR_A.f, coh_WR_B.f, coh_WR_C.f, coh_WR_D.f, coh_WR_E.f)
# 
# coh_SA_A.f <- coh_SA_A[, c(1:3,6)]
# colnames(coh_SA_A.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_SA_A.f$LCP <- "SA"
# coh_SA_A.f$Patch <- "A"
# 
# coh_SA_B.f <- coh_SA_B[, c(1:3,6)]
# colnames(coh_SA_B.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_SA_B.f$LCP <- "SA"
# coh_SA_B.f$Patch <- "B"
# 
# coh_SA_C.f <- coh_SA_C[, c(1:3,6)]
# colnames(coh_SA_C.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_SA_C.f$LCP <- "SA"
# coh_SA_C.f$Patch <- "C"
# 
# coh_SA_D.f <- coh_SA_D[, c(1:3,6)]
# colnames(coh_SA_D.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_SA_D.f$LCP <- "SA"
# coh_SA_D.f$Patch <- "D"
# 
# coh_SA_E.f <- coh_SA_E[, c(1:3,6)]
# colnames(coh_SA_E.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_SA_E.f$LCP <- "SA"
# coh_SA_E.f$Patch <- "E"
# 
# coh_SA.f <- rbind(coh_SA_A.f, coh_SA_B.f, coh_SA_C.f, coh_SA_D.f, coh_SA_E.f)
# 
# coh_WT_A.f <- coh_WT_A[, c(1:3,6)]
# colnames(coh_WT_A.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_WT_A.f$LCP <- "WT"
# coh_WT_A.f$Patch <- "A"
# 
# coh_WT_B.f <- coh_WT_B[, c(1:3,6)]
# colnames(coh_WT_B.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_WT_B.f$LCP <- "WT"
# coh_WT_B.f$Patch <- "B"
# 
# coh_WT_C.f <- coh_WT_C[, c(1:3,6)]
# colnames(coh_WT_C.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_WT_C.f$LCP <- "WT"
# coh_WT_C.f$Patch <- "C"
# 
# coh_WT_D.f <- coh_WT_D[, c(1:3,6)]
# colnames(coh_WT_D.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_WT_D.f$LCP <- "WT"
# coh_WT_D.f$Patch <- "D"
# 
# coh_WT_E.f <- coh_WT_E[, c(1:3,6)]
# colnames(coh_WT_E.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_WT_E.f$LCP <- "WT"
# coh_WT_E.f$Patch <- "E"
# 
# coh_WT.f <- rbind(coh_WT_A.f, coh_WT_B.f, coh_WT_C.f, coh_WT_D.f, coh_WT_E.f)
# 
# coh_AC_A.f <- coh_AC_A[, c(1:3,6)]
# colnames(coh_AC_A.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_AC_A.f$LCP <- "AC"
# coh_AC_A.f$Patch <- "A"
# 
# coh_AC_B.f <- coh_AC_B[, c(1:3,6)]
# colnames(coh_AC_B.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_AC_B.f$LCP <- "AC"
# coh_AC_B.f$Patch <- "B"
# 
# coh_AC_C.f <- coh_AC_C[, c(1:3,6)]
# colnames(coh_AC_C.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_AC_C.f$LCP <- "AC"
# coh_AC_C.f$Patch <- "C"
# 
# coh_AC_D.f <- coh_AC_D[, c(1:3,6)]
# colnames(coh_AC_D.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_AC_D.f$LCP <- "AC"
# coh_AC_D.f$Patch <- "D"
# 
# coh_AC_E.f <- coh_AC_E[, c(1:3,6)]
# colnames(coh_AC_E.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_AC_E.f$LCP <- "AC"
# coh_AC_E.f$Patch <- "E"
# 
# coh_AC.f <- rbind(coh_AC_A.f, coh_AC_B.f, coh_AC_C.f, coh_AC_D.f, coh_AC_E.f)
# 
# coh_JB_A.f <- coh_JB_A[, c(1:3,6)]
# colnames(coh_JB_A.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_JB_A.f$LCP <- "JB"
# coh_JB_A.f$Patch <- "A"
# 
# coh_JB_B.f <- coh_JB_B[, c(1:3,6)]
# colnames(coh_JB_B.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_JB_B.f$LCP <- "JB"
# coh_JB_B.f$Patch <- "B"
# 
# coh_JB_C.f <- coh_JB_C[, c(1:3,6)]
# colnames(coh_JB_C.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_JB_C.f$LCP <- "JB"
# coh_JB_C.f$Patch <- "C"
# 
# coh_JB_D.f <- coh_JB_D[, c(1:3,6)]
# colnames(coh_JB_D.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_JB_D.f$LCP <- "JB"
# coh_JB_D.f$Patch <- "D"
# 
# coh_JB_E.f <- coh_JB_E[, c(1:3,6)]
# colnames(coh_JB_E.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_JB_E.f$LCP <- "JB"
# coh_JB_E.f$Patch <- "E"
# 
# coh_JB.f <- rbind(coh_JB_A.f, coh_JB_B.f, coh_JB_C.f, coh_JB_D.f, coh_JB_E.f)
# 
# coh_GI_A.f <- coh_GI_A[, c(1:3,6)]
# colnames(coh_GI_A.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_GI_A.f$LCP <- "GI"
# coh_GI_A.f$Patch <- "A"
# 
# coh_GI_B.f <- coh_GI_B[, c(1:3,6)]
# colnames(coh_GI_B.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_GI_B.f$LCP <- "GI"
# coh_GI_B.f$Patch <- "B"
# 
# coh_GI_C.f <- coh_GI_C[, c(1:3,6)]
# colnames(coh_GI_C.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_GI_C.f$LCP <- "GI"
# coh_GI_C.f$Patch <- "C"
# 
# coh_GI_D.f <- coh_GI_D[, c(1:3,6)]
# colnames(coh_GI_D.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_GI_D.f$LCP <- "GI"
# coh_GI_D.f$Patch <- "D"
# 
# coh_GI_E.f <- coh_GI_E[, c(1:3,6)]
# colnames(coh_GI_E.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_GI_E.f$LCP <- "GI"
# coh_GI_E.f$Patch <- "E"
# 
# coh_GI.f <- rbind(coh_GI_A.f, coh_GI_B.f, coh_GI_C.f, coh_GI_D.f, coh_GI_E.f)
# 
# coh_CB_A.f <- coh_CB_A[, c(1:3,6)]
# colnames(coh_CB_A.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_CB_A.f$LCP <- "CB"
# coh_CB_A.f$Patch <- "A"
# 
# coh_CB_B.f <- coh_CB_B[, c(1:3,6)]
# colnames(coh_CB_B.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_CB_B.f$LCP <- "CB"
# coh_CB_B.f$Patch <- "B"
# 
# coh_CB_C.f <- coh_CB_C[, c(1:3,6)]
# colnames(coh_CB_C.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_CB_C.f$LCP <- "CB"
# coh_CB_C.f$Patch <- "C"
# 
# coh_CB_D.f <- coh_CB_D[, c(1:3,6)]
# colnames(coh_CB_D.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_CB_D.f$LCP <- "CB"
# coh_CB_D.f$Patch <- "D"
# 
# coh_CB_E.f <- coh_CB_E[, c(1:3,6)]
# colnames(coh_CB_E.f) <- c("LCP", "Patch", "Genet", "Cohesion")
# coh_CB_E.f$LCP <- "CB"
# coh_CB_E.f$Patch <- "E"
# 
# coh_CB.f <- rbind(coh_CB_A.f, coh_CB_B.f, coh_CB_C.f, coh_CB_D.f, coh_CB_E.f)
# 
# coh_all <- rbind(coh_PB.f, coh_CS.f, coh_JA.f, coh_WR.f, coh_SA.f, coh_WT.f, coh_AC.f, coh_JB.f, coh_GI.f, coh_CB.f)
# 
# write.csv(coh_all, "coh_all.csv")

##This is for core area!
# ca_PB_A.f <- ca_PB_A[, c(1:3,6)]
# colnames(ca_PB_A.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_PB_A.f$LCP <- "PB"
# ca_PB_A.f$Patch <- "A"
# 
# ca_PB_B.f <- ca_PB_B[, c(1:3,6)]
# colnames(ca_PB_B.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_PB_B.f$LCP <- "PB"
# ca_PB_B.f$Patch <- "B"
# 
# ca_PB_C.f <- ca_PB_C[, c(1:3,6)]
# colnames(ca_PB_C.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_PB_C.f$LCP <- "PB"
# ca_PB_C.f$Patch <- "C"
# 
# ca_PB_D.f <- ca_PB_D[, c(1:3,6)]
# colnames(ca_PB_D.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_PB_D.f$LCP <- "PB"
# ca_PB_D.f$Patch <- "D"
# 
# ca_PB_E.f <- ca_PB_E[, c(1:3,6)]
# colnames(ca_PB_E.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_PB_E.f$LCP <- "PB"
# ca_PB_E.f$Patch <- "E"
# 
# ca_PB.f <- rbind(ca_PB_A.f, ca_PB_B.f, ca_PB_C.f, ca_PB_D.f, ca_PB_E.f)
# 
# ca_CS_A.f <- ca_CS_A[, c(1:3,6)]
# colnames(ca_CS_A.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_CS_A.f$LCP <- "CS"
# ca_CS_A.f$Patch <- "A"
# 
# ca_CS_B.f <- ca_CS_B[, c(1:3,6)]
# colnames(ca_CS_B.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_CS_B.f$LCP <- "CS"
# ca_CS_B.f$Patch <- "B"
# 
# ca_CS_C.f <- ca_CS_C[, c(1:3,6)]
# colnames(ca_CS_C.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_CS_C.f$LCP <- "CS"
# ca_CS_C.f$Patch <- "C"
# 
# ca_CS_D.f <- ca_CS_D[, c(1:3,6)]
# colnames(ca_CS_D.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_CS_D.f$LCP <- "CS"
# ca_CS_D.f$Patch <- "D"
# 
# ca_CS_E.f <- ca_CS_E[, c(1:3,6)]
# colnames(ca_CS_E.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_CS_E.f$LCP <- "CS"
# ca_CS_E.f$Patch <- "E"
# 
# ca_CS.f <- rbind(ca_CS_A.f, ca_CS_B.f, ca_CS_C.f, ca_CS_D.f, ca_CS_E.f)
# 
# ca_JA_A.f <- ca_JA_A[, c(1:3,6)]
# colnames(ca_JA_A.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_JA_A.f$LCP <- "JA"
# ca_JA_A.f$Patch <- "A"
# 
# ca_JA_B.f <- ca_JA_B[, c(1:3,6)]
# colnames(ca_JA_B.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_JA_B.f$LCP <- "JA"
# ca_JA_B.f$Patch <- "B"
# 
# ca_JA_C.f <- ca_JA_C[, c(1:3,6)]
# colnames(ca_JA_C.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_JA_C.f$LCP <- "JA"
# ca_JA_C.f$Patch <- "C"
# 
# ca_JA_D.f <- ca_JA_D[, c(1:3,6)]
# colnames(ca_JA_D.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_JA_D.f$LCP <- "JA"
# ca_JA_D.f$Patch <- "D"
# 
# ca_JA_E.f <- ca_JA_E[, c(1:3,6)]
# colnames(ca_JA_E.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_JA_E.f$LCP <- "JA"
# ca_JA_E.f$Patch <- "E"
# 
# ca_JA.f <- rbind(ca_JA_A.f, ca_JA_B.f, ca_JA_C.f, ca_JA_D.f, ca_JA_E.f)
# 
# ca_WR_A.f <- ca_WR_A[, c(1:3,6)]
# colnames(ca_WR_A.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_WR_A.f$LCP <- "WR"
# ca_WR_A.f$Patch <- "A"
# 
# ca_WR_B.f <- ca_WR_B[, c(1:3,6)]
# colnames(ca_WR_B.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_WR_B.f$LCP <- "WR"
# ca_WR_B.f$Patch <- "B"
# 
# ca_WR_C.f <- ca_WR_C[, c(1:3,6)]
# colnames(ca_WR_C.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_WR_C.f$LCP <- "WR"
# ca_WR_C.f$Patch <- "C"
# 
# ca_WR_D.f <- ca_WR_D[, c(1:3,6)]
# colnames(ca_WR_D.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_WR_D.f$LCP <- "WR"
# ca_WR_D.f$Patch <- "D"
# 
# ca_WR_E.f <- ca_WR_E[, c(1:3,6)]
# colnames(ca_WR_E.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_WR_E.f$LCP <- "WR"
# ca_WR_E.f$Patch <- "E"
# 
# ca_WR.f <- rbind(ca_WR_A.f, ca_WR_B.f, ca_WR_C.f, ca_WR_D.f, ca_WR_E.f)
# 
# ca_SA_A.f <- ca_SA_A[, c(1:3,6)]
# colnames(ca_SA_A.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_SA_A.f$LCP <- "SA"
# ca_SA_A.f$Patch <- "A"
# 
# ca_SA_B.f <- ca_SA_B[, c(1:3,6)]
# colnames(ca_SA_B.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_SA_B.f$LCP <- "SA"
# ca_SA_B.f$Patch <- "B"
# 
# ca_SA_C.f <- ca_SA_C[, c(1:3,6)]
# colnames(ca_SA_C.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_SA_C.f$LCP <- "SA"
# ca_SA_C.f$Patch <- "C"
# 
# ca_SA_D.f <- ca_SA_D[, c(1:3,6)]
# colnames(ca_SA_D.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_SA_D.f$LCP <- "SA"
# ca_SA_D.f$Patch <- "D"
# 
# ca_SA_E.f <- ca_SA_E[, c(1:3,6)]
# colnames(ca_SA_E.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_SA_E.f$LCP <- "SA"
# ca_SA_E.f$Patch <- "E"
# 
# ca_SA.f <- rbind(ca_SA_A.f, ca_SA_B.f, ca_SA_C.f, ca_SA_D.f, ca_SA_E.f)
# 
# ca_WT_A.f <- ca_WT_A[, c(1:3,6)]
# colnames(ca_WT_A.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_WT_A.f$LCP <- "WT"
# ca_WT_A.f$Patch <- "A"
# 
# ca_WT_B.f <- ca_WT_B[, c(1:3,6)]
# colnames(ca_WT_B.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_WT_B.f$LCP <- "WT"
# ca_WT_B.f$Patch <- "B"
# 
# ca_WT_C.f <- ca_WT_C[, c(1:3,6)]
# colnames(ca_WT_C.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_WT_C.f$LCP <- "WT"
# ca_WT_C.f$Patch <- "C"
# 
# ca_WT_D.f <- ca_WT_D[, c(1:3,6)]
# colnames(ca_WT_D.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_WT_D.f$LCP <- "WT"
# ca_WT_D.f$Patch <- "D"
# 
# ca_WT_E.f <- ca_WT_E[, c(1:3,6)]
# colnames(ca_WT_E.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_WT_E.f$LCP <- "WT"
# ca_WT_E.f$Patch <- "E"
# 
# ca_WT.f <- rbind(ca_WT_A.f, ca_WT_B.f, ca_WT_C.f, ca_WT_D.f, ca_WT_E.f)
# 
# ca_AC_A.f <- ca_AC_A[, c(1:3,6)]
# colnames(ca_AC_A.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_AC_A.f$LCP <- "AC"
# ca_AC_A.f$Patch <- "A"
# 
# ca_AC_B.f <- ca_AC_B[, c(1:3,6)]
# colnames(ca_AC_B.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_AC_B.f$LCP <- "AC"
# ca_AC_B.f$Patch <- "B"
# 
# ca_AC_C.f <- ca_AC_C[, c(1:3,6)]
# colnames(ca_AC_C.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_AC_C.f$LCP <- "AC"
# ca_AC_C.f$Patch <- "C"
# 
# ca_AC_D.f <- ca_AC_D[, c(1:3,6)]
# colnames(ca_AC_D.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_AC_D.f$LCP <- "AC"
# ca_AC_D.f$Patch <- "D"
# 
# ca_AC_E.f <- ca_AC_E[, c(1:3,6)]
# colnames(ca_AC_E.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_AC_E.f$LCP <- "AC"
# ca_AC_E.f$Patch <- "E"
# 
# ca_AC.f <- rbind(ca_AC_A.f, ca_AC_B.f, ca_AC_C.f, ca_AC_D.f, ca_AC_E.f)
# 
# ca_JB_A.f <- ca_JB_A[, c(1:3,6)]
# colnames(ca_JB_A.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_JB_A.f$LCP <- "JB"
# ca_JB_A.f$Patch <- "A"
# 
# ca_JB_B.f <- ca_JB_B[, c(1:3,6)]
# colnames(ca_JB_B.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_JB_B.f$LCP <- "JB"
# ca_JB_B.f$Patch <- "B"
# 
# ca_JB_C.f <- ca_JB_C[, c(1:3,6)]
# colnames(ca_JB_C.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_JB_C.f$LCP <- "JB"
# ca_JB_C.f$Patch <- "C"
# 
# ca_JB_D.f <- ca_JB_D[, c(1:3,6)]
# colnames(ca_JB_D.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_JB_D.f$LCP <- "JB"
# ca_JB_D.f$Patch <- "D"
# 
# ca_JB_E.f <- ca_JB_E[, c(1:3,6)]
# colnames(ca_JB_E.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_JB_E.f$LCP <- "JB"
# ca_JB_E.f$Patch <- "E"
# 
# ca_JB.f <- rbind(ca_JB_A.f, ca_JB_B.f, ca_JB_C.f, ca_JB_D.f, ca_JB_E.f)
# 
# ca_GI_A.f <- ca_GI_A[, c(1:3,6)]
# colnames(ca_GI_A.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_GI_A.f$LCP <- "GI"
# ca_GI_A.f$Patch <- "A"
# 
# ca_GI_B.f <- ca_GI_B[, c(1:3,6)]
# colnames(ca_GI_B.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_GI_B.f$LCP <- "GI"
# ca_GI_B.f$Patch <- "B"
# 
# ca_GI_C.f <- ca_GI_C[, c(1:3,6)]
# colnames(ca_GI_C.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_GI_C.f$LCP <- "GI"
# ca_GI_C.f$Patch <- "C"
# 
# ca_GI_D.f <- ca_GI_D[, c(1:3,6)]
# colnames(ca_GI_D.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_GI_D.f$LCP <- "GI"
# ca_GI_D.f$Patch <- "D"
# 
# ca_GI_E.f <- ca_GI_E[, c(1:3,6)]
# colnames(ca_GI_E.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_GI_E.f$LCP <- "GI"
# ca_GI_E.f$Patch <- "E"
# 
# ca_GI.f <- rbind(ca_GI_A.f, ca_GI_B.f, ca_GI_C.f, ca_GI_D.f, ca_GI_E.f)
# 
# ca_CB_A.f <- ca_CB_A[, c(1:3,6)]
# colnames(ca_CB_A.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_CB_A.f$LCP <- "CB"
# ca_CB_A.f$Patch <- "A"
# 
# ca_CB_B.f <- ca_CB_B[, c(1:3,6)]
# colnames(ca_CB_B.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_CB_B.f$LCP <- "CB"
# ca_CB_B.f$Patch <- "B"
# 
# ca_CB_C.f <- ca_CB_C[, c(1:3,6)]
# colnames(ca_CB_C.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_CB_C.f$LCP <- "CB"
# ca_CB_C.f$Patch <- "C"
# 
# ca_CB_D.f <- ca_CB_D[, c(1:3,6)]
# colnames(ca_CB_D.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_CB_D.f$LCP <- "CB"
# ca_CB_D.f$Patch <- "D"
# 
# ca_CB_E.f <- ca_CB_E[, c(1:3,6)]
# colnames(ca_CB_E.f) <- c("LCP", "Patch", "Genet", "Ramets")
# ca_CB_E.f$LCP <- "CB"
# ca_CB_E.f$Patch <- "E"
# 
# ca_CB.f <- rbind(ca_CB_A.f, ca_CB_B.f, ca_CB_C.f, ca_CB_D.f, ca_CB_E.f)
# 
# ca_all <- rbind(ca_PB.f, ca_CS.f, ca_JA.f, ca_WR.f, ca_SA.f, ca_WT.f, ca_AC.f, ca_JB.f, ca_GI.f, ca_CB.f)
# 
# write.csv(ca_all, "ca_all.csv")

#This is for the interspersion-juxtaposition index!
# iji_PB_A.f <- iji_PB_A[, c(1:3,6)]
# colnames(iji_PB_A.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_PB_A.f$LCP <- "PB"
# iji_PB_A.f$Patch <- "A"
# 
# iji_PB_B.f <- iji_PB_B[, c(1:3,6)]
# colnames(iji_PB_B.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_PB_B.f$LCP <- "PB"
# iji_PB_B.f$Patch <- "B"
# 
# iji_PB_C.f <- iji_PB_C[, c(1:3,6)]
# colnames(iji_PB_C.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_PB_C.f$LCP <- "PB"
# iji_PB_C.f$Patch <- "C"
# 
# iji_PB_D.f <- iji_PB_D[, c(1:3,6)]
# colnames(iji_PB_D.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_PB_D.f$LCP <- "PB"
# iji_PB_D.f$Patch <- "D"
# 
# iji_PB_E.f <- iji_PB_E[, c(1:3,6)]
# colnames(iji_PB_E.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_PB_E.f$LCP <- "PB"
# iji_PB_E.f$Patch <- "E"
# 
# iji_PB.f <- rbind(iji_PB_A.f, iji_PB_B.f, iji_PB_C.f, iji_PB_D.f, iji_PB_E.f)
# 
# iji_CS_A.f <- iji_CS_A[, c(1:3,6)]
# colnames(iji_CS_A.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_CS_A.f$LCP <- "CS"
# iji_CS_A.f$Patch <- "A"
# 
# iji_CS_B.f <- iji_CS_B[, c(1:3,6)]
# colnames(iji_CS_B.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_CS_B.f$LCP <- "CS"
# iji_CS_B.f$Patch <- "B"
# 
# iji_CS_C.f <- iji_CS_C[, c(1:3,6)]
# colnames(iji_CS_C.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_CS_C.f$LCP <- "CS"
# iji_CS_C.f$Patch <- "C"
# 
# iji_CS_D.f <- iji_CS_D[, c(1:3,6)]
# colnames(iji_CS_D.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_CS_D.f$LCP <- "CS"
# iji_CS_D.f$Patch <- "D"
# 
# iji_CS_E.f <- iji_CS_E[, c(1:3,6)]
# colnames(iji_CS_E.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_CS_E.f$LCP <- "CS"
# iji_CS_E.f$Patch <- "E"
# 
# iji_CS.f <- rbind(iji_CS_A.f, iji_CS_B.f, iji_CS_C.f, iji_CS_D.f, iji_CS_E.f)
# 
# iji_JA_A.f <- iji_JA_A[, c(1:3,6)]
# colnames(iji_JA_A.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_JA_A.f$LCP <- "JA"
# iji_JA_A.f$Patch <- "A"
# 
# iji_JA_B.f <- iji_JA_B[, c(1:3,6)]
# colnames(iji_JA_B.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_JA_B.f$LCP <- "JA"
# iji_JA_B.f$Patch <- "B"
# 
# iji_JA_C.f <- iji_JA_C[, c(1:3,6)]
# colnames(iji_JA_C.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_JA_C.f$LCP <- "JA"
# iji_JA_C.f$Patch <- "C"
# 
# iji_JA_D.f <- iji_JA_D[, c(1:3,6)]
# colnames(iji_JA_D.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_JA_D.f$LCP <- "JA"
# iji_JA_D.f$Patch <- "D"
# 
# iji_JA_E.f <- iji_JA_E[, c(1:3,6)]
# colnames(iji_JA_E.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_JA_E.f$LCP <- "JA"
# iji_JA_E.f$Patch <- "E"
# 
# iji_JA.f <- rbind(iji_JA_A.f, iji_JA_B.f, iji_JA_C.f, iji_JA_D.f, iji_JA_E.f)
# 
# iji_WR_A.f <- iji_WR_A[, c(1:3,6)]
# colnames(iji_WR_A.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_WR_A.f$LCP <- "WR"
# iji_WR_A.f$Patch <- "A"
# 
# iji_WR_B.f <- iji_WR_B[, c(1:3,6)]
# colnames(iji_WR_B.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_WR_B.f$LCP <- "WR"
# iji_WR_B.f$Patch <- "B"
# 
# iji_WR_C.f <- iji_WR_C[, c(1:3,6)]
# colnames(iji_WR_C.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_WR_C.f$LCP <- "WR"
# iji_WR_C.f$Patch <- "C"
# 
# iji_WR_D.f <- iji_WR_D[, c(1:3,6)]
# colnames(iji_WR_D.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_WR_D.f$LCP <- "WR"
# iji_WR_D.f$Patch <- "D"
# 
# iji_WR_E.f <- iji_WR_E[, c(1:3,6)]
# colnames(iji_WR_E.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_WR_E.f$LCP <- "WR"
# iji_WR_E.f$Patch <- "E"
# 
# iji_WR.f <- rbind(iji_WR_A.f, iji_WR_B.f, iji_WR_C.f, iji_WR_D.f, iji_WR_E.f)
# 
# iji_SA_A.f <- iji_SA_A[, c(1:3,6)]
# colnames(iji_SA_A.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_SA_A.f$LCP <- "SA"
# iji_SA_A.f$Patch <- "A"
# 
# iji_SA_B.f <- iji_SA_B[, c(1:3,6)]
# colnames(iji_SA_B.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_SA_B.f$LCP <- "SA"
# iji_SA_B.f$Patch <- "B"
# 
# iji_SA_C.f <- iji_SA_C[, c(1:3,6)]
# colnames(iji_SA_C.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_SA_C.f$LCP <- "SA"
# iji_SA_C.f$Patch <- "C"
# 
# iji_SA_D.f <- iji_SA_D[, c(1:3,6)]
# colnames(iji_SA_D.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_SA_D.f$LCP <- "SA"
# iji_SA_D.f$Patch <- "D"
# 
# iji_SA_E.f <- iji_SA_E[, c(1:3,6)]
# colnames(iji_SA_E.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_SA_E.f$LCP <- "SA"
# iji_SA_E.f$Patch <- "E"
# 
# iji_SA.f <- rbind(iji_SA_A.f, iji_SA_B.f, iji_SA_C.f, iji_SA_D.f, iji_SA_E.f)
# 
# iji_WT_A.f <- iji_WT_A[, c(1:3,6)]
# colnames(iji_WT_A.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_WT_A.f$LCP <- "WT"
# iji_WT_A.f$Patch <- "A"
# 
# iji_WT_B.f <- iji_WT_B[, c(1:3,6)]
# colnames(iji_WT_B.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_WT_B.f$LCP <- "WT"
# iji_WT_B.f$Patch <- "B"
# 
# iji_WT_C.f <- iji_WT_C[, c(1:3,6)]
# colnames(iji_WT_C.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_WT_C.f$LCP <- "WT"
# iji_WT_C.f$Patch <- "C"
# 
# iji_WT_D.f <- iji_WT_D[, c(1:3,6)]
# colnames(iji_WT_D.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_WT_D.f$LCP <- "WT"
# iji_WT_D.f$Patch <- "D"
# 
# iji_WT_E.f <- iji_WT_E[, c(1:3,6)]
# colnames(iji_WT_E.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_WT_E.f$LCP <- "WT"
# iji_WT_E.f$Patch <- "E"
# 
# iji_WT.f <- rbind(iji_WT_A.f, iji_WT_B.f, iji_WT_C.f, iji_WT_D.f, iji_WT_E.f)
# 
# iji_AC_A.f <- iji_AC_A[, c(1:3,6)]
# colnames(iji_AC_A.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_AC_A.f$LCP <- "AC"
# iji_AC_A.f$Patch <- "A"
# 
# iji_AC_B.f <- iji_AC_B[, c(1:3,6)]
# colnames(iji_AC_B.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_AC_B.f$LCP <- "AC"
# iji_AC_B.f$Patch <- "B"
# 
# iji_AC_C.f <- iji_AC_C[, c(1:3,6)]
# colnames(iji_AC_C.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_AC_C.f$LCP <- "AC"
# iji_AC_C.f$Patch <- "C"
# 
# iji_AC_D.f <- iji_AC_D[, c(1:3,6)]
# colnames(iji_AC_D.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_AC_D.f$LCP <- "AC"
# iji_AC_D.f$Patch <- "D"
# 
# iji_AC_E.f <- iji_AC_E[, c(1:3,6)]
# colnames(iji_AC_E.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_AC_E.f$LCP <- "AC"
# iji_AC_E.f$Patch <- "E"
# 
# iji_AC.f <- rbind(iji_AC_A.f, iji_AC_B.f, iji_AC_C.f, iji_AC_D.f, iji_AC_E.f)
# 
# iji_JB_A.f <- iji_JB_A[, c(1:3,6)]
# colnames(iji_JB_A.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_JB_A.f$LCP <- "JB"
# iji_JB_A.f$Patch <- "A"
# 
# iji_JB_B.f <- iji_JB_B[, c(1:3,6)]
# colnames(iji_JB_B.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_JB_B.f$LCP <- "JB"
# iji_JB_B.f$Patch <- "B"
# 
# iji_JB_C.f <- iji_JB_C[, c(1:3,6)]
# colnames(iji_JB_C.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_JB_C.f$LCP <- "JB"
# iji_JB_C.f$Patch <- "C"
# 
# iji_JB_D.f <- iji_JB_D[, c(1:3,6)]
# colnames(iji_JB_D.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_JB_D.f$LCP <- "JB"
# iji_JB_D.f$Patch <- "D"
# 
# iji_JB_E.f <- iji_JB_E[, c(1:3,6)]
# colnames(iji_JB_E.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_JB_E.f$LCP <- "JB"
# iji_JB_E.f$Patch <- "E"
# 
# iji_JB.f <- rbind(iji_JB_A.f, iji_JB_B.f, iji_JB_C.f, iji_JB_D.f, iji_JB_E.f)
# 
# iji_GI_A.f <- iji_GI_A[, c(1:3,6)]
# colnames(iji_GI_A.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_GI_A.f$LCP <- "GI"
# iji_GI_A.f$Patch <- "A"
# 
# iji_GI_B.f <- iji_GI_B[, c(1:3,6)]
# colnames(iji_GI_B.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_GI_B.f$LCP <- "GI"
# iji_GI_B.f$Patch <- "B"
# 
# iji_GI_C.f <- iji_GI_C[, c(1:3,6)]
# colnames(iji_GI_C.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_GI_C.f$LCP <- "GI"
# iji_GI_C.f$Patch <- "C"
# 
# iji_GI_D.f <- iji_GI_D[, c(1:3,6)]
# colnames(iji_GI_D.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_GI_D.f$LCP <- "GI"
# iji_GI_D.f$Patch <- "D"
# 
# iji_GI_E.f <- iji_GI_E[, c(1:3,6)]
# colnames(iji_GI_E.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_GI_E.f$LCP <- "GI"
# iji_GI_E.f$Patch <- "E"
# 
# iji_GI.f <- rbind(iji_GI_A.f, iji_GI_B.f, iji_GI_C.f, iji_GI_D.f, iji_GI_E.f)
# 
# iji_CB_A.f <- iji_CB_A[, c(1:3,6)]
# colnames(iji_CB_A.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_CB_A.f$LCP <- "CB"
# iji_CB_A.f$Patch <- "A"
# 
# iji_CB_B.f <- iji_CB_B[, c(1:3,6)]
# colnames(iji_CB_B.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_CB_B.f$LCP <- "CB"
# iji_CB_B.f$Patch <- "B"
# 
# iji_CB_C.f <- iji_CB_C[, c(1:3,6)]
# colnames(iji_CB_C.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_CB_C.f$LCP <- "CB"
# iji_CB_C.f$Patch <- "C"
# 
# iji_CB_D.f <- iji_CB_D[, c(1:3,6)]
# colnames(iji_CB_D.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_CB_D.f$LCP <- "CB"
# iji_CB_D.f$Patch <- "D"
# 
# iji_CB_E.f <- iji_CB_E[, c(1:3,6)]
# colnames(iji_CB_E.f) <- c("LCP", "Patch", "Genet", "IJI")
# iji_CB_E.f$LCP <- "CB"
# iji_CB_E.f$Patch <- "E"
# 
# iji_CB.f <- rbind(iji_CB_A.f, iji_CB_B.f, iji_CB_C.f, iji_CB_D.f, iji_CB_E.f)
# 
# iji_all <- rbind(iji_PB.f, iji_CS.f, iji_JA.f, iji_WR.f, iji_SA.f, iji_WT.f, iji_AC.f, iji_JB.f, iji_GI.f, iji_CB.f)
# 
# write.csv(iji_all, "iji_all.csv")

################################
#SIMULATING RANDOM DISTRIBUTION#
################################

#set random number generator seed, can use an integer to test code
set.seed(Sys.time())

# tail probability for confidence intervals
alpha = 0.05

# number of randomization iterations for null distributions
iterations <- 1000

# read in data
original.df <- read.csv("12.Long_Island_MLL_pop_patch_x_y_Feb_2025.csv")

split.df <- original.df %>% group_by(Population) %>% group_split()

split.df.1 <- split.df[[1]] %>% group_by(Patch) %>% group_split()
split.df.2 <- split.df[[2]] %>% group_by(Patch) %>% group_split()
split.df.3 <- split.df[[3]] %>% group_by(Patch) %>% group_split()
split.df.4 <- split.df[[4]] %>% group_by(Patch) %>% group_split()
split.df.5 <- split.df[[5]] %>% group_by(Patch) %>% group_split()
split.df.6 <- split.df[[6]] %>% group_by(Patch) %>% group_split()
split.df.7 <- split.df[[7]] %>% group_by(Patch) %>% group_split()
split.df.8 <- split.df[[8]] %>% group_by(Patch) %>% group_split()
split.df.9 <- split.df[[9]] %>% group_by(Patch) %>% group_split()
split.df.10 <- split.df[[10]] %>% group_by(Patch) %>% group_split()

all.split <- list(split.df.1, split.df.2, split.df.3, split.df.4, split.df.5,
  split.df.6, split.df.7, split.df.8, split.df.9, split.df.10)

n.genotypes <- c()

# get the number of genotypes in patch
for (i in 1:length(unique(original.df$Patch))){
  n.genotypes.1 <- length(unique(split.df.1[[i]]$clone))
  n.genotypes.2 <- length(unique(split.df.2[[i]]$clone))
  n.genotypes.3 <- length(unique(split.df.3[[i]]$clone))
  n.genotypes.4 <- length(unique(split.df.4[[i]]$clone))
  n.genotypes.5 <- length(unique(split.df.5[[i]]$clone))
  n.genotypes.6 <- length(unique(split.df.6[[i]]$clone))
  n.genotypes.7 <- length(unique(split.df.7[[i]]$clone))
  n.genotypes.8 <- length(unique(split.df.8[[i]]$clone))
  n.genotypes.9 <- length(unique(split.df.9[[i]]$clone))
  n.genotypes.10 <- length(unique(split.df.10[[i]]$clone))
  
  n.genotypes <- rbind(n.genotypes, n.genotypes.1, n.genotypes.2, n.genotypes.3, 
  n.genotypes.4, n.genotypes.5, n.genotypes.6, n.genotypes.7, n.genotypes.8,
  n.genotypes.9, n.genotypes.10)
  }

# reformat data as a raster
###
# Randomize Genotypes to generate null distributions
# copy original data

k <- 0

end.index <- 0

cohesion.vector <- vector(length = iterations*sum(n.genotypes))
#aggregation.vector <- vector(length = iterations*sum(n.genotypes))
iji.vector <- vector(length = iterations*sum(n.genotypes))

for (i in 1:length(unique(original.df$Population))){
  
  for (j in 1:length(unique(original.df$Patch))){
  
    k <- k+1
      
    for (v in 1:iterations){
    
      curr.pop <- all.split[[i]]
      curr.patch <- curr.pop[[j]]
      
      curr.genetonly <- curr.patch[, 4:6]
      
      randomized.df <- curr.patch[, 4:6]
      
      random.rows <- sample(nrow(curr.genetonly))
      
      randomized.df[3] <- curr.genetonly[random.rows, 3]
      
      rand.data.rast <- rasterFromXYZ(randomized.df)
      
      end.index <- end.index + n.genotypes[k]
      
      start.index <- end.index - n.genotypes[k] + 1
      sprintf("start.index = %i", start.index)
      sprintf("end.index = %i", end.index) 
  
      cohesion.null.estimate <- lsm_c_cohesion(rand.data.rast)
      #cohesion.landscape.null.estimates[[start.index:end.index]] <- cohesion.null.estimate$value
      cohesion.vector[start.index:end.index] <- cohesion.null.estimate$value
    
      #aggregation.null.estimate <- lsm_c_ai(data.rast)
      #aggregation.landscape.null.estimates <- aggregation.null.estimate$value
      #aggregation.vector[start.index:end.index] <- aggregation.null.estimate$value
    
      iji.null.estimate <- lsm_c_iji(rand.data.rast)
      #iji.landscape.null.estimates[[i]] <- iji.null.estimate$value
      iji.vector[start.index:end.index] <- iji.null.estimate$value
      
      print(end.index)
      
      }

    }
  }

# make the vector a data frame
cohesion.df <- data.frame(cohesion.vector)
cohesion.df$metric <- 'rand cohesion'
# assigning new names to the columns of the data frame 
colnames(cohesion.df) <- c('value','metric') 

# aggregation.df <- data.frame(aggregation.vector)
# aggregation.df$metric <- 'rand aggregation'
# # assigning new names to the columns of the data frame 
# colnames(aggregation.df) <- c('value','metric') 

iji.df <- data.frame(iji.vector)
iji.df$metric <- 'rand iji'
# assigning new names to the columns of the data frame 
colnames(iji.df) <- c('value','metric') 

# combine into into a single data frame (stacking rows)
null.metrics <- rbind(cohesion.df, iji.df)

write.csv(null.metrics, "nullmetrics.csv")

###################################
#PLOT AGGREGATION AND NULL METRICS#
###################################

Agg.Data <- read.csv("13.Aggregation_Full_LI_SP2.csv", header = TRUE)
#Agg.Data2 <- read.csv("Aggregation2_Full_LI_SP.csv", header = TRUE)

Rand.Data <- read.csv("nullmetrics.csv", header = TRUE)

agg.obs <- ggplot(Agg.Data, aes(x = Value, color = Index, fill = Index)) +
  geom_freqpoly(bins=25, size = 0.5, aes(y = after_stat(density))) +
  scale_color_manual(values=c("#d45973", "#3c155d")) +  
  facet_zoom(ylim = c(0, 0.0735), zoom.size = 4.05) + 
  scale_x_continuous(expand = c(0, 0), limits = c(-5.4, 105.5)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + theme_light() +
  xlab("") + ylab("Density") +
  theme(axis.title.x = element_text(hjust = 0.34)) +
  theme(plot.title = element_text(hjust = 0.075)) + 
  theme(legend.position="none", axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), plot.margin=unit(c(0.1,0.1,0.5,0.1), "cm"))

Agg.Data2 <- Agg.Data[Agg.Data$Value > 0,]

Rand.iji <- Rand.Data[Rand.Data$metric == "rand iji",]
Rand.coh <- Rand.Data[Rand.Data$metric == "rand cohesion",]

Rand.Data2 <- Rand.Data[Rand.Data$value > 0,]

Rand.iji2 <- Rand.iji[Rand.iji$value > 0,]
Rand.coh2 <- Rand.coh[Rand.coh$value > 0,]

Rand.Data$zeroes <- "With Zeroes"
Rand.Data2$zeroes <- "Without Zeroes"

Rand.Data.f <- rbind(Rand.Data, Rand.Data2)

Rand.Data.f$zeroes <- factor(Rand.Data.f$zeroes, levels = c("With Zeroes", "Without Zeroes"))

agg.rand <- ggplot(Rand.Data, aes(x = value, fill = metric, color = metric)) +
  geom_histogram(alpha = 0.5, aes(y = after_stat(density)), position = 'identity', bins = 10)+
  scale_fill_manual(values=c("#d45973", "#3c155d"), name = "Index\n(Randomized)", labels = c("Cohesion", "IJI"))+ 
  scale_color_manual(values=c("#d45973", "#3c155d"), name = "Index\n(Randomized)", labels = c("Cohesion", "IJI"))+
  scale_x_continuous(expand = c(0, 0), limits = c(NA, NA)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + theme_light() +
  xlab("Index Value") + ylab("Density") +
  #facet_wrap(~ zeroes, nrow = 2) +
  theme(strip.background=element_rect(color=NA, fill=NA), 
    strip.text = element_text(color = "white", size = 14)) +
  theme(legend.title = element_text(size = 16.85), plot.margin=unit(c(0.5,0.1,0.1,0.1), "cm"))


grid.arrange(agg.obs, agg.rand, ncol = 1)

agg.obs2 <- ggplot(Agg.Data2, aes(x = Value, color = Index, fill = Index)) +
  geom_freqpoly(bins=25, size = 0.5, aes(y = after_stat(density))) +
  scale_color_manual(values=c("#d45973", "#3c155d")) +  
  #facet_zoom(ylim = c(0, 0.035), zoom.size = 4.05) + 
  scale_x_continuous(expand = c(0, 0), limits = c(0, 105.5)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.075)) + theme_light() +
  xlab("") + ylab("Density") +
  theme(axis.title.x = element_text(hjust = 0.34)) +
  theme(plot.title = element_text(hjust = 0.075)) + 
  theme(legend.position="bottom", axis.title.x = element_blank(), axis.ticks.x = element_blank(), plot.margin=unit(c(0.1,0.1,0.5,0.1), "cm"))

agg.rand2 <- ggplot(Rand.Data2, aes(x = value, fill = metric, color = metric)) +
  geom_histogram(alpha = 0.5, aes(y = after_stat(density)), position = 'identity', bins = 10)+
  scale_fill_manual(values=c("#d45973", "#3c155d"), name = "Index\n(Randomized)", labels = c("Cohesion", "IJI"))+ 
  scale_color_manual(values=c("#d45973", "#3c155d"), name = "Index\n(Randomized)", labels = c("Cohesion", "IJI"))+
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.076)) + theme_light() +
  xlab("Index Value") + ylab("Density") +
  #facet_wrap(~ zeroes, nrow = 2) +
  theme(strip.background=element_rect(color=NA, fill=NA), 
    strip.text = element_text(color = "white", size = 14)) +
  theme(legend.position="none", plot.margin=unit(c(0.5,0.1,0.1,0.1), "cm"))

grid.arrange(agg.obs2, agg.rand2, ncol = 1)

agg.rand2

agg.obs2

