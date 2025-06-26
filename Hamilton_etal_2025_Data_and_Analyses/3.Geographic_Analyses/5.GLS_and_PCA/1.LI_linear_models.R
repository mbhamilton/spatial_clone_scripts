#########################################
#                                       #
#   Title: Linear Models to Analyze     #
#   Genetic Differentiation             #
#   Author: Charli D. Minsavage-Davis   #
#   Date last edited: 26 June 2025      #
#   Manuscript: S. patens landscape     # 
#   genetics - Long Island & NYC        #
#                                       #
#########################################

#Library for getting path
library(rstudioapi)
library(nlme)
library(rcompanion)
library(vegan)
library(FactoMineR)
library(factoextra)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(scales)

#Getting the path of your current open file
current_path <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
print(getwd())

#All independent and dependent data
dat <- read.csv("2.data_for_linear_models.csv")

datstand <- as.data.frame(dat[, c(1:24)])

datstand <- datstand[!duplicated(datstand$i.noherb),]

for (i in 3:17){
  datstand[, i] <- abs(datstand[, i])/max(abs(datstand[, i]))
}

#Run PCA to understand the weight of variables
corr_matrix <- cor(datstand[, c(3, 5:9, 13:14)])
#ggcorrplot(corr_matrix)

data.pca <- princomp(corr_matrix)
summary(data.pca)

fviz_pca_var(data.pca, col.var = "black", xlim = c(-0.65, 0.65))
fviz_cos2(data.pca, choice = "var", axes = 1:2)

##Run suite of models for each independent and several ecologically relevant 
##combinations against linearized G'st 

g.r1 <- gls(GpST/(1-GpST)~d+(pat2019-pat1974)+(marsh2019-marsh1974)+i.estab, data = datstand)
summary(g.r1) 

g.r2 <- gls(GpST/(1-GpST)~d+(pat2019-pat1974)+(marsh2019-marsh1974)+i.noherb, data = datstand)
summary(g.r2) 

g.r3 <- gls(GpST/(1-GpST)~d+(pat2019-pat1974)+(marsh2019-marsh1974)+i.1974, data = datstand)
summary(g.r3) 

g.r4 <- gls(GpST/(1-GpST)~d+i.estab, data = datstand)
summary(g.r4)

g.r5 <- gls(GpST/(1-GpST)~d+i.noherb, data = datstand)
summary(g.r5)

g.r6 <- gls(GpST/(1-GpST)~d+i.1974, data = datstand)
summary(g.r6)

g.r7 <- gls(GpST/(1-GpST)~d, data = datstand)
summary(g.r7)

g.r8 <- gls(GpST/(1-GpST)~coast+i.estab, data = datstand)
summary(g.r8)

g.r9 <- gls(GpST/(1-GpST)~coast+i.noherb, data = datstand)
summary(g.r9)

g.r10 <- gls(GpST/(1-GpST)~coast+i.1974, data = datstand)
summary(g.r10)

##Read in reformatted data for ease of plotting
datdist <- read.csv("3.datdist.csv", header = TRUE)
datconnect <- read.csv("4.datconnect.csv", header = TRUE)
datland <- read.csv("5.datland.csv", header = TRUE)

#find means based on distance
means.dist <- aggregate(linear~type, data = datdist, mean)

#Plot distance against differentiation
p <- ggplot(datdist, aes(x=value, y=G.st/(1-G.st), color=type)) + 
  #xlim(0, 0.5)+
  #ylim(0, 0.25)+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlim(0, 1)+
  geom_point()+
  scale_color_manual(values = c("#E85362ff", "#400f73ff"), labels = c("Coastal", "Euclidean")) +
  guides(col=guide_legend(title="")) +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  #geom_text(aes(x = -0.05, y = linear+0.115), label = "*", color = "#400f73ff", size = 6, data = means.dist, show.legend=FALSE)+
  #geom_text(aes(x = -0.05, y = linear+0.02), label = "ns", color = "#E85362ff", size = 4, data = means.dist, show.legend=FALSE)+
  labs(x="Standardized Distance", y = expression(italic("G'"["ST"])~"/(1 — "~italic("G'"["ST"])~")"))+
  theme_classic() +
  theme(legend.margin = margin(c(0, 0, -10, 0)), legend.position = "top")+
  theme(axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), legend.text = element_text(size=12))

p

#uni <- unique(datconnect$type)#uni G.st<- unique(datconnect$type)
#linedat1 <- datconnect %>% filter(type == uni[1])
#linedat2 <- datconnect %>% filter(type == uni[2])
#linedat3 <- datconnect %>% filter(type == uni[3])

#Find means for connectivity
means.connect <- aggregate(linear~type, data = datconnect, mean)

#Plot connectivity against differentiation
g <- ggplot(datconnect, aes(x=value, y=G.st/(1-G.st), color=type)) + 
  #ylim(0, 0.25)+
  geom_point()+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlim(0, 1)+
  scale_color_manual(values = c("#E85362ff", "#952c80ff", "#400f73ff"), labels = c("1974", "2019A", "2019B")) +
  guides(col=guide_legend(title="")) +
  geom_smooth(method = "lm", se=FALSE, fullrange=TRUE)+
  #geom_text(aes(x = 1, y = linear+0.115), label = "*", color = "#400f73ff", size = 6, data = means.connect, show.legend=FALSE)+
 # geom_text(aes(x = 1, y = linear+0.115), label = "*", color = "#952c80ff", size = 6, data = means.connect, show.legend=FALSE)+
 # geom_text(aes(x = 1, y = linear+0.02), label = "ns", color = "#E85362ff", size = 4, data = means.connect, show.legend=FALSE)+ 
  #stat_function(fun=function(x) {1/x}, args = list(mean = mean(datconnect$value), sd = sd(datconnect$value)))+
  labs(x="Standardized Current", y = expression(italic("G'"["ST"])~"/(1 — "~italic("G'"["ST"])~")"))+
  theme_classic()+
  theme(legend.margin = margin(c(0, 0, -10, 0)), legend.position = "top")+
  theme(axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), legend.text = element_text(size=12))

g
#g2 + geom_line(data = linedat1, size = 1, stat = "function", fun = function(x) {1/((x-0.35)*26)}) + 
 # geom_line(data = linedat2, size = 1, stat = "function", fun = function(x) {1/((x-0.45)*26)}) +
 # geom_line(data = linedat3, size = 1, stat = "function", fun = function(x) {1/((x-0.15)*26)})

#Plot land comparisons against differentiation
b <- ggplot(datland, aes(x=value, y=G.st/(1-G.st), color=type)) + 
  #xlim(0, 1)+
  #ylim(0, 0.25)+
  geom_point()+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + 
  xlim(0, 1)+
  scale_color_manual(values = c("#E85362ff", "#400f73ff"), labels = c("All Marsh", expression(italic("S. patens")))) +
  guides(col=guide_legend(title="")) +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  labs(x="Standardized Area Change", y = expression(italic("G'"["ST"])~"/(1 — "~italic("G'"["ST"])~")"))+
  theme_classic()+
  theme(legend.margin = margin(c(0, 0, -10, 0)), legend.position = "top")+
  theme(axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), legend.text = element_text(size=12))
  
b


p
g
b
