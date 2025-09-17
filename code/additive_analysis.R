library(tidyverse)
library(ggridges)
setwd("~/mastersarbeit/thesis/code")

# Visualizing the adaptation time of uniform landscapes of various sizes  

filepaths <- c("output/0.1kgen_additive_landscape_20250723_1625.csv",
               "output/0.2kgen_additive_landscape_20250723_1626.csv",
               "output/0.3kgen_additive_landscape_20250723_1627.csv",
               "output/0.5kgen_additive_landscape_20250723_1628.csv",
               "output/0.8kgen_additive_landscape_20250723_1629.csv",
               "output/1kgen_additive_landscape_20250723_1631.csv",
               "output/2kgen_additive_landscape_20250723_1633.csv")
timepoint <- c(100, 200, 300, 500, 800, 1000, 2000)

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

ggplot(data, aes(x = loci_diff, y = as.factor(timepoint), fill = as.factor(Loci), color = as.factor(Loci))) +
  facet_wrap(data$Loci) +
  geom_density_ridges(scale = 3) +
  scale_fill_viridis_d(option = "C", begin = 0.2, end = 0.85, alpha = 0.65) +
  scale_color_viridis_d(option = "C", begin = 0.2, end= 0.85) +
  labs(x = "Hamming Distance from Global Peak", y = "Generations", 
       title = "Time to Equilibrium on Uniform Additive Landscapes (n = 33/Condition)",
       fill = "Loci",
       color = "Loci") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_minimal()
ggsave("outputs/figures/uniform_additive_t_eq.pdf", bg = "transparent")

# Now determine if the variable additive landscapes are randomly walking at a 50-25-25% ratio. 

df02 <- read.csv("outputs/data/04_expanding_additive_replicate6_20250901.csv") %>%
  rename(`1` = X1, `10` = X10, `20` = X20) %>% 
  pivot_longer(cols = c(`1`, `10`, `20`), names_to = "loci")

ggplot(df02, aes(y = value, group = loci, fill = loci, color = loci)) +
  geom_density(bw = 0.5, alpha = 0.5) +
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.8) +
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.8) +
  ylim(c(0,20)) +
  labs(y = "Genome Size", 
       x = "Density", 
#       title = "Histogram of genome sizes (1e5 Generations)",
       fill = "Initial genome size",
       color = "Initial genome size"
       ) + 
  geom_hline(yintercept = mean(df02$value), color= "black", lty = 2) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_minimal()
ggsave("outputs/figures/additive_densities_rep6.pdf", bg = "transparent")
