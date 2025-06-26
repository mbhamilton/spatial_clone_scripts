#########################################
#                                       #
#   Title: MLL Phylogenies              # 
#   Author: Charli D. Minsavage-Davis   #
#   Date last edited: 26 June 2025      #
#   Manuscript: S. patens landscape     # 
#   genetics - Long Island & NYC        #
#                                       #
#########################################

#library("BiocManager")

#BiocManager::install("remotes")
#BiocManager::install("YuLab-SMU/treedataverse")
#BiocManager::install("treeio") 

library(ape)
library(ggtree)
library(remotes)
library(treeio)
library(ggplot2)
library(RColorBrewer)
library(aplot)
library(rstudioapi)
library(viridis)
library(tidytree)
library(ggpubr)
library(grid)

#Getting the path of your current open file
current_path <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
print(getwd())

#Read in trees generated with Bruvos distance in Matlab
newick_tree1 <- ape::read.tree(file = "2.BruvodistNJtreesamplenametiplabel.txt")

newick_tree2 <- ape::read.tree(file = "3.BruvodistNJtreeclonenumbertiplabel.txt")

pops <- read.csv("4.LongIslandfull10popsMLLassignments2.csv")

#Organize data
newdat <- data.frame(label = newick_tree2$tip.label, numpops = pops$numpops)

#Join data into new tree
tree2 <- full_join(newick_tree2, newdat, by = "label")

#Set random seed
set.seed(123)

#Subsdet tree
tree <- ape::rtree(6)
sub_tree <- tree_subset(tree, node = "t1", levels_back = 2)

#Plot full tree
q <-
  ggtree(tree2, options(ignore.negative.edge=TRUE)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 0.75)) + scale_y_continuous(expand = c(0.01, 0), limits = c(0, 676))+
  #xlab("Evolutionary Distance") +
  #ylab("Unique MLL")+
  geom_tippoint(aes(color = numpops), size = 3) + 
  scale_color_viridis_c(option = "magma", begin = 0.2, end = 0.8, direction = 1, name = "Population\nCount", limits = c(0, 15), guide = "none") +
  theme_update()+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.5), panel.background = element_rect(colour = "black", fill=NA), axis.text.x=element_text(size = 18), axis.text.y=element_text(size = 18),
    legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=18), #change legend title font size
        legend.text = element_text(size =14))+ 
    geom_hilight(node=772, fill="steelblue", alpha=.4)+
    geom_hilight(node=933, fill="darkmagenta", alpha=.4)+
  coord_cartesian(clip = "off")

  #geom_tiplab(size=-10, align=TRUE, linesize=0.3)

q

#Create new sub tree to plot slices from phylogeny
subtree <- tree_subset(tree2, node = 772, levels_back = 0)

#Plot subtree
p <-
  ggtree(subtree, options(ignore.negative.edge=TRUE)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 0.375)) + scale_y_continuous(expand = c(0.01, 0), limits = c(0, 21))+
  #xlab("Evolutionary Distance") +
  #ylab("Unique MLL")+
  geom_tippoint(aes(color = numpops), size = 3) + 
  scale_color_viridis_c(option = "magma", begin = 0.2, end = 0.8, direction = 1, name = "Population\nCount\n", limits = c(0, 15)) +
  theme_update()+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.5), panel.background = element_rect(colour = "black", fill=NA), axis.text.y=element_text(size = 18),axis.text.x=element_text(size = 18),
    legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=18), #change legend title font size
        legend.text = element_text(size =14),
        plot.margin = unit(c(20,5.5,5.5,5.5), "pt"))+
  geom_tiplab(size=5, hjust = -0.1)+ 
  annotate("rect", xmin=0.35, xmax=0.3725, ymin=0, ymax=21, alpha=0.4, fill="steelblue")+
  coord_cartesian(clip = "off")

p

#Create new subtree
subtree2 <- tree_subset(tree2, node = 933, levels_back = 0)

#Plot this subtree
p2 <-
  ggtree(subtree2, options(ignore.negative.edge=TRUE)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 0.375)) + scale_y_continuous(expand = c(0.01, 0), limits = c(0, 23))+
  #xlab("Evolutionary Distance") +
  #ylab("Unique MLL")+
  geom_tippoint(aes(color = numpops), size = 3) + 
  scale_color_viridis_c(option = "magma", begin = 0.2, end = 0.8, direction = 1, name = "Population\nCount\n", limits = c(0, 15)) +
  theme_update()+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.5), panel.background = element_rect(colour = "black", fill=NA), axis.text.y=element_text(size = 18),axis.text.x=element_text(size = 18),
    legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=18), #change legend title font size
        legend.text = element_text(size =14),
        plot.margin = unit(c(5.5,5.5,20,5.5), "pt"))+
  geom_tiplab(size=5, hjust = -0.1)+ 
  annotate("rect", xmin=0.35, xmax=0.3725, ymin=0, ymax=23, alpha=0.4, fill="darkmagenta")+
  coord_cartesian(clip = "off")

p2

#Arrange the above plots into facets
p.p2 <- ggarrange(p2, p, common.legend = TRUE, legend = "right", nrow = 2)
p.p2

#Again
p.q <- ggarrange(q, p.p2)

#Annotate the full figure
annotate_figure(p.q, left = textGrob("Unique MLLs", rot = 90, vjust = 1, gp = gpar(cex = 1.65)),
                    bottom = textGrob("Evolutionary Distance", hjust = 0.425, gp = gpar(cex = 1.65)))


#Plot full tree
ttt1 <-
  ggtree(tree2, options(ignore.negative.edge=TRUE)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 0.71)) + scale_y_continuous(expand = c(0.01, 0), limits = c(1, 169))+
  xlab("Evolutionary Distance") +
  ylab("Unique MLL")+
  #geom_tippoint(aes(color = numpops), size = 3) + 
  scale_color_viridis_c(option = "magma", begin = 0.2, end = 0.8, direction = 1, name = "Population\nCount", limits = c(0, 15), guide = "none") +
  theme_update()+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.5), panel.background = element_rect(colour = "black", fill=NA), axis.text.x=element_text(size = 18), axis.text.y=element_text(size = 18), axis.title.x=element_text(size = 20), axis.title.y=element_text(size = 20),
    legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=18), #change legend title font size
        legend.text = element_text(size =14))+ 
    #geom_hilight(node=772, fill="steelblue", alpha=.4)+
    #geom_hilight(node=933, fill="darkmagenta", alpha=.4)+
  coord_cartesian(clip = "off")+

  geom_tiplab(size=4, align=TRUE, linesize=0.3)

ttt1

#Different format
ttt2 <-
  ggtree(tree2, options(ignore.negative.edge=TRUE)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 0.71)) + scale_y_continuous(expand = c(0.001, 0), limits = c(170, 338))+
  xlab("Evolutionary Distance") +
  ylab("Unique MLL")+
  #geom_tippoint(aes(color = numpops), size = 3) + 
  scale_color_viridis_c(option = "magma", begin = 0.2, end = 0.8, direction = 1, name = "Population\nCount", limits = c(0, 15), guide = "none") +
  theme_update()+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.5), panel.background = element_rect(colour = "black", fill=NA), axis.text.x=element_text(size = 18), axis.text.y=element_text(size = 18), axis.title.x=element_text(size = 20), axis.title.y=element_text(size = 20),
    legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=18), #change legend title font size
        legend.text = element_text(size =14))+ 
    #geom_hilight(node=772, fill="steelblue", alpha=.4)+
    #geom_hilight(node=933, fill="darkmagenta", alpha=.4)+
  coord_cartesian(clip = "off")+

  geom_tiplab(size=4, align=TRUE, linesize=0.3)

ttt2

#Different format
ttt3 <-
  ggtree(tree2, options(ignore.negative.edge=TRUE)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 0.71)) + scale_y_continuous(expand = c(0.001, 0), limits = c(339, 507))+
  xlab("Evolutionary Distance") +
  ylab("Unique MLL")+
  #geom_tippoint(aes(color = numpops), size = 3) + 
  scale_color_viridis_c(option = "magma", begin = 0.2, end = 0.8, direction = 1, name = "Population\nCount", limits = c(0, 15), guide = "none") +
  theme_update()+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.5), panel.background = element_rect(colour = "black", fill=NA), axis.text.x=element_text(size = 18), axis.text.y=element_text(size = 18), axis.title.x=element_text(size = 20), axis.title.y=element_text(size = 20),
    legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=18), #change legend title font size
        legend.text = element_text(size =14))+ 
    #geom_hilight(node=772, fill="steelblue", alpha=.4)+
    #geom_hilight(node=933, fill="darkmagenta", alpha=.4)+
  coord_cartesian(clip = "off")+

  geom_tiplab(size=4, align=TRUE, linesize=0.3)

ttt3

#Different format
ttt4 <-
  ggtree(tree2, options(ignore.negative.edge=TRUE)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 0.71)) + scale_y_continuous(expand = c(0.001, 0), limits = c(508, 676))+
  xlab("Evolutionary Distance") +
  ylab("Unique MLL")+
  #geom_tippoint(aes(color = numpops), size = 3) + 
  scale_color_viridis_c(option = "magma", begin = 0.2, end = 0.8, direction = 1, name = "Population\nCount", limits = c(0, 15), guide = "none") +
  theme_update()+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.5), panel.background = element_rect(colour = "black", fill=NA), axis.text.x=element_text(size = 18), axis.text.y=element_text(size = 18), axis.title.x=element_text(size = 20), axis.title.y=element_text(size = 20),
    legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=18), #change legend title font size
        legend.text = element_text(size =14))+ 
    #geom_hilight(node=772, fill="steelblue", alpha=.4)+
    #geom_hilight(node=933, fill="darkmagenta", alpha=.4)+
  coord_cartesian(clip = "off")+

  geom_tiplab(size=4, align=TRUE, linesize=0.3)

ttt4





