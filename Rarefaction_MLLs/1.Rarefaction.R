#########################################
#                                       #
#   Title: MLL Rarefaction              # 
#   Author: Charli D. Minsavage-Davis   #
#   Date last edited: 23 June 2025      #
#   Manuscript: S. patens landscape     # 
#   genetics - Long Island & NYC        #
#                                       #
#########################################

#Read in requisite package
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rstudioapi)

#Getting the path of your current open file and setting working directory
current_path <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
print(getwd())

#Set seed!
set.seed(123123123)

#Reformat data into tibble with genet observations per patch, population,
#or among all populations tallied


#This is for population data
data_pop <- read_tsv("3b.Rarefaction_pop.txt") %>%
  select(Group, starts_with("clone")) %>%
  pivot_longer(-Group) %>%
  group_by(Group) %>%
  mutate(total = sum(value)) %>%
  group_by(name) %>%
  mutate(total = sum(value)) %>%
  filter(total != 0) %>%
  ungroup() %>%
  select(-total)

#This is for among populations; all data
data_all <- read_tsv("2b.Rarefaction_all.txt") %>%
  select(Group, starts_with("clone")) %>%
  pivot_longer(-Group) %>%
  group_by(Group) %>%
  mutate(total = sum(value)) %>%
  group_by(name) %>%
  mutate(total = sum(value)) %>%
  filter(total != 0) %>%
  ungroup() %>%
  select(-total)

#User-defined function to produce collector's sampling curves
collect <- function(data, group){
  
  data %>%
    filter(Group == group) %>% #Group is an aggregator; i.e., population or patch
    uncount(value) %>%
    sample_n(n()) %>% #Here we sample from all observations to randomize observation data; i.e., collecting
    mutate(observation  = row_number()) %>% #Return row no.
    arrange(name, observation) %>% 
    group_by(name) %>%
    mutate(distinct = row_number() == 1) %>% #For each group, only take 1 distinct MLL
    ungroup() %>% 
    arrange(observation) %>% #Reorder to original
    mutate (s = cumsum(distinct)) %>% #sum all unique MLLs only
    select(observation, s) #take s alongside the current observation!

}

#map_dfr will run the collector's curve x times and aggregate
collect_curves <- map_dfr(1:100, ~collect(data_all, "all"), .id = "iteration")  

#This takes a mean of the collector's curves; i.e., rarefaction
rarefaction_curve <- collect_curves %>%
  group_by(observation) %>%
  summarize(r = mean(s))

#Plot the collector's curves and overlay the rarefaction curve
Rare_all <- collect_curves %>%
  ggplot(aes(x = observation, y = s, group = iteration)) + 
  xlab("Sample Size (N)") +
  ylab("No. Observed MLLs") +
  labs(title = "MLL Rarefaction Among Populations") +
  geom_line(color = "gray") + 
  geom_line(data = rarefaction_curve, 
    aes(x = observation, y = r), inherit.aes = FALSE,
    color = "black", linewidth = 1) +
  theme_classic()
  
#Run the collection curve function 100 times for each population
collect_curves1 <- map_dfr(1:100, ~collect(data_pop, "pop1"), .id = "iteration")  
collect_curves2 <- map_dfr(1:100, ~collect(data_pop, "pop2"), .id = "iteration")
collect_curves3 <- map_dfr(1:100, ~collect(data_pop, "pop3"), .id = "iteration")
collect_curves4 <- map_dfr(1:100, ~collect(data_pop, "pop4"), .id = "iteration")
collect_curves5 <- map_dfr(1:100, ~collect(data_pop, "pop5"), .id = "iteration")
collect_curves6 <- map_dfr(1:100, ~collect(data_pop, "pop6"), .id = "iteration")
collect_curves7 <- map_dfr(1:100, ~collect(data_pop, "pop7"), .id = "iteration")
collect_curves8 <- map_dfr(1:100, ~collect(data_pop, "pop8"), .id = "iteration")
collect_curves9 <- map_dfr(1:100, ~collect(data_pop, "pop9"), .id = "iteration")
collect_curves10 <- map_dfr(1:100, ~collect(data_pop, "pop10"), .id = "iteration")

#Label within each population data frame
collect_curves1[, 4] <- "1"
names(collect_curves1)[names(collect_curves1) == "...4"] <- "Population"

collect_curves2[, 4] <- "2"
names(collect_curves2)[names(collect_curves2) == "...4"] <- "Population"

collect_curves3[, 4] <- "3"
names(collect_curves3)[names(collect_curves3) == "...4"] <- "Population"

collect_curves4[, 4] <- "4"
names(collect_curves4)[names(collect_curves4) == "...4"] <- "Population"

collect_curves5[, 4] <- "5"
names(collect_curves5)[names(collect_curves5) == "...4"] <- "Population"

collect_curves6[, 4] <- "6"
names(collect_curves6)[names(collect_curves6) == "...4"] <- "Population"

collect_curves7[, 4] <- "7"
names(collect_curves7)[names(collect_curves7) == "...4"] <- "Population"

collect_curves8[, 4] <- "8"
names(collect_curves8)[names(collect_curves8) == "...4"] <- "Population"

collect_curves9[, 4] <- "9"
names(collect_curves9)[names(collect_curves9) == "...4"] <- "Population"

collect_curves10[, 4] <- "10"
names(collect_curves10)[names(collect_curves10) == "...4"] <- "Population"

#Generate rarefaction curve for each population
rarefaction_curve1 <- collect_curves1 %>%
  group_by(observation) %>%
  summarize(r = mean(s))

rarefaction_curve2 <- collect_curves2 %>%
  group_by(observation) %>%
  summarize(r = mean(s))

rarefaction_curve3 <- collect_curves3 %>%
  group_by(observation) %>%
  summarize(r = mean(s))

rarefaction_curve4 <- collect_curves4 %>%
  group_by(observation) %>%
  summarize(r = mean(s))

rarefaction_curve5 <- collect_curves5 %>%
  group_by(observation) %>%
  summarize(r = mean(s))

rarefaction_curve6 <- collect_curves6 %>%
  group_by(observation) %>%
  summarize(r = mean(s))

rarefaction_curve7 <- collect_curves7 %>%
  group_by(observation) %>%
  summarize(r = mean(s))

rarefaction_curve8 <- collect_curves8 %>%
  group_by(observation) %>%
  summarize(r = mean(s))

rarefaction_curve9 <- collect_curves9 %>%
  group_by(observation) %>%
  summarize(r = mean(s))

rarefaction_curve10 <- collect_curves10 %>%
  group_by(observation) %>%
  summarize(r = mean(s))

#Bind curves for each population together with their unique IDs
collect_curves_pop <- rbind(collect_curves1, collect_curves2, collect_curves3, 
  collect_curves4, collect_curves5, collect_curves6, collect_curves7, 
  collect_curves8, collect_curves9, collect_curves10)

#Reorder because R is dumb and puts 10 after 1...
collect_curves_pop$Population <- factor(collect_curves_pop$Population, levels = c("1", "2", "3",
  "4", "5", "6", "7", "8", "9", "10"))

#Plot all population rarefaction curves in same figure
Pop_Rare <- collect_curves_pop %>%
  ggplot(aes(x = observation, y = s, group = interaction(Population, iteration))) + 
  xlab("Sample Size (n)") +
  ylab("No. Observed MLLs") +
  labs(title = "MLL Rarefaction Within Populations") + 
  geom_line(data = rarefaction_curve1, 
    aes(x = observation, y = r, color = "Pelham Bay (PB)"), linetype = 1, inherit.aes = FALSE, linewidth = 1) +
  geom_line(data = rarefaction_curve2, 
    aes(x = observation, y = r, color = "Caumsett SP (CS)"), linetype = 2, inherit.aes = FALSE, linewidth = 1) +
  geom_line(data = rarefaction_curve3, 
    aes(x = observation, y = r, color = "Ambro Preserve (JA)"), linetype = 3, inherit.aes = FALSE, linewidth = 1) +
  geom_line(data = rarefaction_curve4, 
    aes(x = observation, y = r, color = "Wading River (WR)"), linetype = 4, inherit.aes = FALSE, linewidth = 1) +
  geom_line(data = rarefaction_curve5, 
    aes(x = observation, y = r, color = "Sammys Beach (SB)"), linetype = 5, inherit.aes = FALSE, linewidth = 1) +
  geom_line(data = rarefaction_curve6, 
    aes(x = observation, y = r, color = "William T Davis (WT)"), linetype = 6, inherit.aes = FALSE, linewidth = 1) +
  geom_line(data = rarefaction_curve7, 
    aes(x = observation, y = r, color = "Alley Creek (AC)"), linetype = 7, inherit.aes = FALSE, linewidth = 1) +
  geom_line(data = rarefaction_curve8, 
    aes(x = observation, y = r, color = "Jones Beach (JB)"), linetype = 8, inherit.aes = FALSE, linewidth = 1) +
  geom_line(data = rarefaction_curve9, 
    aes(x = observation, y = r, color = "Cedar Beach (GI)"), linetype = 9, inherit.aes = FALSE, linewidth = 1) +
  geom_line(data = rarefaction_curve10, 
    aes(x = observation, y = r, color = "Cupsogue Beach (CB)"), linetype = 10, inherit.aes = FALSE, linewidth = 1) +
  scale_color_manual(name = "", values = 1:10) +
  theme_classic()

#Arrange the total and population level rarefaction curves side-by-side
ggarrange(Rare_all, Pop_Rare, widths = c(1, 1.5))

