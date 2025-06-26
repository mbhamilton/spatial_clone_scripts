#########################################
#                                       #
#   Title: Expected genotypes           # 
#   Author: Matthew B. Hamilton         #
#   Date last edited: 26 June 2025      #
#   Manuscript: S. patens landscape     # 
#   genetics - Long Island & NYC        #
#                                       #
#########################################

# Load rstudioapi package to be able to get directory of source file
library("rstudioapi")

setwd(dirname(getActiveDocumentContext()$path)) # Set working directory to source file location

getwd() # Check updated working directory

library(devtools)

# # RClone is a R package version of GenClone program (Arnaud-Haond & Belkhir 2007)
devtools::install_github("dbailleul/RClone")

library(RClone)

# read in five locus data set 
#five_locus_data <- read.csv(file = "/Users/hamilton_lab/Dropbox/Organized SP Data (Moved from lab computer)/Long Island/Rclone scripts for Long Island data/H_L_microsat_data_5_dip_loci.csv", header=TRUE, fileEncoding="UTF-8-BOM")

orig.genotype.data <- read.table("2.Long_Island_Ten_Pops_07_22_2024.txt", sep="\t", header=TRUE)

# look at the header and top few lines of the data file
head(orig.genotype.data)

# determine number of loci in genotype data
n.loci <- ncol(orig.genotype.data)/2

#make a copy of the genotype data to modify
genotype.data <- orig.genotype.data
  
# need to impute alleles for any missing values
# missing is scored as zero

#library(dplyr)
#library(tidy)
#library(purrr)

# one_locus_alleles <- orig.genotype.data[,1:2]
# 
# # find rows with zero values
# #list <- one_locus_alleles[which(one_locus_alleles == 0), ]
# 
# list <- which(one_locus_alleles == 0, arr.ind=TRUE)
# 
# #replace.size <- length(list)
# 
# replace.size <- as.numeric(length(one_locus_alleles[(one_locus_alleles == 0)]))
# 
# # from alleles not equal zero sample two with replacement
# #rand_allele_index <- sample(which(one_locus_alleles != 0), nrow(list), replace = TRUE)
# 
# # from alleles not equal zero sample two with replacement
# new_alleles <- sample(one_locus_alleles[!(one_locus_alleles == 0)], replace.size, replace = T)
# 
# #new_one_locus_alleles <- replace(one_locus_alleles, list, new_alleles)
# 
# genotype.data[,1:2] <- replace(orig.genotype.data[,1:2], list, new_alleles)

# Impute missing data
# for each locus, find zeros (or the missing value) and then replace them with alleles
# sampled at random from among the alleles in both columns for that locus. This "imputes" 
# missing data using alleles based on the observed allele frequencies
missing.value = 0

for (col in seq(from=1, to=(2*n.loci), by=2)) {
  
  #print(col)

  one_locus_alleles <- orig.genotype.data[,col:(col+1)]
  
  list <- which(one_locus_alleles == missing.value, arr.ind=TRUE)
  
  #print(list)
  
  replace.size <- as.numeric(length(one_locus_alleles[(one_locus_alleles == 0)]))
  
  #print(replace.size)
  
  if (replace.size != 0){
    # from alleles not equal zero sample two with replacement
    new_alleles <- sample(one_locus_alleles[!(one_locus_alleles == 0)], replace.size, replace = T)
    
    genotype.data[,col:(col+1)] <- replace(one_locus_alleles, list, new_alleles)
  }
}

# check to be sure that all zeros were detected and replaced
#list <- which(genotype.data == 0, arr.ind=TRUE)

# save the imputed genotype data to tab delimited text file if you want
write.table(genotype.data, "Long_Island_Ten_Pops_07_22_2024 imputed missing.txt", append = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)

# see https://www.rdocumentation.org/packages/RClone/versions/1.0.3
vignette("RClone_quickmanual")

# make a data table to summarize the alleles at each locus
list_all_tab(genotype.data)

# make a table showing the multilocus genotypes - could take a while to compute
MLG_tab(genotype.data)

# make allele frequency table, round-robin estimation method
freq_RR(genotype.data, RR = TRUE) #on ramets

# 
MLLlist <- MLL_generator(genotype.data, manh = TRUE)


genet_result <- psex(genotype.data, genet = TRUE) # psex on genets, gives the expected genotype frequency of each genet under psex assumptions
plot(genet_result, main='Expected frequency of each genet multilocus genotype\nunder strict sexual reproduction', xlab='Genet', ylab='Expected frequency')

# compute the probability that repeated genotypes originate from distinct sexual events (i.e. being different genets and not ramets of the same MLG) assuming H-W equilibrium
table <- psex(genotype.data, RR = TRUE, nbrepeat = 100, bar = TRUE) # with p-values and a progress bar

# display the results table

table[[1]]

# compute the probability that repeated genotypes originate from distinct sexual events (i.e. being different genets and not ramets of the same MLG) assuming adjusting for Fis
# *** this is VERY SLOW to execute for the H_L_microsat_data_5_dip_loci.csv data set, requiring +/- one day to run ***
table_Fisadj <- psex_Fis(genotype.data, RR = TRUE, nbrepeat = 100, bar = TRUE) # with p-values and a progress bar

# display the results table
table_Fisadj[[1]]

# make a scatter plot of the MLG probabilities for each genet
plot(table_Fisadj[[1]]$genet, table_Fisadj[[1]]$psexFis, xlim = c(1, 1683), ylim = c(0, 0.4), main='Distribution of probabilities of identical MLG for each genet \n (round-robin allele frequencies, fixation-index adjusted genotype freqs)', xlab='genet', ylab='prob identical MLG with sexual reproduction')


# make a histogram of the MLG probabilities for each genet
# remove any rows that are blank (the ramets)
psex <- table_Fisadj[[1]] [(! table_Fisadj[[1]]$psexFis==""), ]

class(psex$psexFis)
psex_num <- as.numeric(psex$psexFis)

hist(psex_num, breaks=16, xlim=c(0,0.4), ylim=c(0,1400), main='Distribution of probabilities of identical MLG for each genet \n (round-robin allele frequencies, fixation-index adjusted genotype freqs)', xlab='Probability of multilocus genotype (MLG)', ylab='Count of genets')

length(psex_num)

