library(tidyverse)

filepaths <- c("output/0.1kgen_additive_landscape_20250722_15'54.csv",
               "output/0.2kgen_additive_landscape_20250722_14'45.csv",
               "output/0.8kgen_additive_landscape_20250722_16'13.csv",
               "output/1kgen_additive_landscape_20250722_15'55.csv",
               "output/2kgen_additive_landscape_20250722_14'47.csv",
               "output/3kgen_additive_landscape_20250722_16'03.csv",
               "output/20kgen_additive_landscape_20250722_12'43.csv")
timepoint <- c(100, 200, 800, 1000, 2000, 3000, 20000)

data <- data.frame(Loci = c(), 
                   Replicate = c(), 
                   loci_diff = c(), 
                   Final_fitness = c(),
                   Max_fitness = c(),
                   timepoint = c())
for (i in 1:length(filepaths)){
  temp_data <- read_csv(filepaths[[i]]) %>% 
    add_column(`timepoint` = timepoint[[i]])
  data <- rbind(data, temp_data)
}

ggplot(data, aes(y = loci_diff, x = log(timepoint), color = as.factor(Loci))) +
  geom_point()
