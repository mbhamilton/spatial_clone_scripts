#########################################
#                                       #
#   Title: Reverse Cumulative           #
#   Distributions for MLL               #
#   Size Distributions                  # 
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


# read in MLL assignment data for all populations and generate reverse cumulative distributions
# of ramets per genet for each population
library(readxl)

library(dplyr) # for count function, bind_rows function

library(pracma) # for histc function

library(spatstat.utils) # for revcumsum function

# Excel file has columns: Population	Individual	Patch	x	y	clone
full_data <- read.csv("2.Long_Island_MLL_pop_patch_x_y_Feb_2025.csv")

# determine how many populations in the data file
n.pops = max(full_data$Population)

# make an empty dataframe to hold the reverse cumulative distributions for all the populations
full_results_df <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(full_results_df) <- c('ramet.count', 'obs.genets', 'rev.cum.count', 'rev.cum.proportion', 'pop')

ramet_values <- data.frame(matrix(ncol = 1, nrow = n.pops))
colnames(ramet_values) <- c('maxramets')

#one_pop <- subset(full_data, Population == "1")

#function to count clones and return reverse cumulative distribution proportions for one population
# formatted as four columns - ramet count, genets obs with ramet count, rev cum proportion, and population
pop_rev_dist <- function(one_pop) {
  
  clone_counts <- one_pop %>% count(clone, sort = TRUE)
  
  maxcount = as.numeric(max(clone_counts$n))
  
  dist_bins <- seq(1, maxcount, by = 1)
  
  # dist.ctn gives the counts, each row represents the same numbered bin
  counts = histc(clone_counts$n, dist_bins)
  
  reverse_cumulative_counts = revcumsum(counts$cnt)
  
  reverse_cumulative_proportions = reverse_cumulative_counts/nrow(one_pop)
  
  result <- data.frame(matrix(ncol = 5, nrow = maxcount))
  colnames(result) <- c('ramet.count', 'obs.genets', 'rev.cum.count', 'rev.cum.proportion', 'pop')

  result$rev.cum.proportion <- reverse_cumulative_proportions #fill in reverse cumulative proportions column
  
  result$ramet.count <- c(1:maxcount)
  
  result$rev.cum.count <- reverse_cumulative_counts
  
  result$obs.genets <- counts$cnt
  
  result$pop = rep.int(one_pop[1,1], 1) #fill in population column, all same value
  
  return (result)
}

for (i in (1:n.pops)) {
  
  one_pop <- subset(full_data, Population == i)
  
  result <- pop_rev_dist(one_pop)
  
  full_results_df <- rbind(full_results_df, result)
  
  #record the maximum x count for each population
  ramet_values$maxramets[i] = nrow(result)
  
  ### need to compute pooled tails for each population too
  
}

# determine the maximum ramet value across all pops
max.ramet = max(ramet_values)

# determine the maximum ramet value with pooled tails across all pops
#max.ramet.pooled.tails = max(ramet_values_pooled)

# coerce to dataframe; see https://stackoverflow.com/questions/74702850/exporting-dataframe-by-write-csv
full_results_df  <- as.data.frame(full_results_df)
full_results_df[] <- lapply(full_results_df, unlist)

write.csv(full_results_df,file = "Long_Island_MLLs.csv", row.names = FALSE)

