# script to plot clone identity on x-y axes, requires input file with x, y and clone identity columns

#########################################
#                                       #
#   Title: Plot clone identity on       #
#   x-y axes                            # 
#   Author: Matthew B. Hamilton and     #
#   Charli D. Minsavage-Davis           #
#   Date last edited: 26 June 2025      #
#   Manuscript: S. patens landscape     # 
#   genetics - Long Island & NYC        #
#                                       #
#########################################

##Read in all necessary packages
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
library(rstudioapi) # Load rstudioapi package

setwd(dirname(getActiveDocumentContext()$path)) # Set working directory to source file location

getwd() # Check updated working directory

#Read in data
full_data <- read.csv("2.Long_Island_MLL_pop_patch_x_y_Feb_2025.csv")

##The below are only examples, you can choose any patch from any population to plot

# get values for one patch within one population to use in plot
plot_data1 <- subset(full_data, Population == "1", Patch == "A")

plot_data2 <- subset(full_data, Population == "2", Patch == "E")

plot_data3 <- subset(full_data, Population == "3", Patch == "C")

plot_data4 <- subset(full_data, Population == "4", Patch == "B")

#Plot selected patches
#MUST CHANGE X AND Y LIMITS TO MATCH THE ABOVE DATA!

p1 <- ggplot(plot_data1, aes(x=y, y=x, color=clone)) +
  ggtitle("E2: Tuckahoe, NJ MLLs")+
  scale_y_discrete(name = "x-transect (m)", limits = c("1","2","3","4", "5"), expand = c(0, 0.5))+
  scale_x_discrete(name = "y-coordinate (m)", limits = c("1","2","3","4", "5", "6", "7", "8", "9", "10", "11","12","13","14", "15", "16", "17", "18", "19", "20", "21","22","23","24", "25", "26", "27", "28", "29", "30"), expand = c(0, 0.8))+
  geom_point(size=0, show.legend = FALSE)+
  geom_label(label = plot_data1$clone, aes(fill=factor(clone)), color = "white", fontface="bold", label.padding = unit(0.5, "lines"), label.r = unit(0.5, "lines"), size = 4, show.legend = FALSE)+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title.y=element_text(size=16), axis.title.x=element_text(size=16), title = element_text(size = 18))

p1

p2 <- ggplot(plot_data2, aes(x=y, y=x, color=clone)) +
  ggtitle("D2: Tuckerton, NJ MLLs")+
  scale_y_discrete(name = "x-transect (m)", limits = c("1","2","3","4", "5"), expand = c(0, 0.5))+
  scale_x_discrete(name = "", limits = c("1","2","3","4", "5", "6", "7", "8", "9", "10", "11","12","13","14", "15", "16", "17", "18", "19", "20", "21","22","23","24", "25", "26", "27", "28", "29", "30", "31", "32"), expand = c(0, 0.8))+
  geom_point(size=0, show.legend = FALSE)+
  geom_label(label = plot_data2$clone, aes(fill=factor(clone)), color = "white", fontface="bold", label.padding = unit(0.5, "lines"), label.r = unit(0.5, "lines"), size = 4, show.legend = FALSE)+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title.y=element_text(size=16), axis.title.x=element_text(size=0), title = element_text(size = 18))

#p2

p3 <- ggplot(plot_data3, aes(x=y+1, y=x, color=clone)) +
  ggtitle("C2: Cedar Beach, NY MLLs")+
  scale_x_discrete(name = "", limits = c("0", "1","2","3","4", "5", "6", "7", "8", "9", "10", "11","12","13","14", "15", "16", "17", "18", "19", "20"), expand = c(0, 0.8))+
  scale_y_discrete(name = "x-transect (m)", limits = c("1","2","3","4"), expand = c(0, 0.5))+
  geom_point(size=0, show.legend = FALSE)+
  geom_label(label = plot_data3$clone, aes(fill=factor(clone)), color = "white", fontface="bold", label.padding = unit(0.5, "lines"), label.r = unit(0.5, "lines"), size = 4, show.legend = FALSE)+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title.y=element_text(size=16), axis.title.x=element_text(size=0), title = element_text(size = 18))

#p3

p4 <- ggplot(plot_data4, aes(x=y, y=x, color=clone)) +
  ggtitle("B2: Cape Cod, MA MLLs")+
  scale_y_discrete(name = "x-transect (m)", limits = c("1","2","3","4", "5"), expand = c(0, 0.5))+
  scale_x_discrete(name = "", limits = c("1","2","3","4", "5", "6", "7", "8", "9", "10", "11","12","13","14", "15", "16", "17", "18", "19", "20", "21","22","23","24", "25", "26", "27", "28", "29", "30"), expand = c(0, 0.8))+
  geom_point(size=0, show.legend = FALSE)+
  geom_label(label = plot_data4$clone, aes(fill=factor(clone)), color = "white", fontface="bold", label.padding = unit(0.5, "lines"), label.r = unit(0.5, "lines"), size = 4, show.legend = FALSE)+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title.y=element_text(size=16), axis.title.x=element_text(size=0), title = element_text(size = 18))

#p4

