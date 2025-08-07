library(tidyverse)
library(ggridges)
setwd("~/mastersarbeit/thesis/code")

# Begin with the static landscape
filepath_df01 <- "output/epistatic_landscape_20250805_1654.csv"
df01 <- read.csv(filepath_df01)

ggplot(df01, aes(x = Final_fitness, y = as.factor(Generation), fill = as.factor(Loci), color = as.factor(Loci))) +
  facet_wrap(df01$Loci) +
  geom_density_ridges(scale = 3, from = 0) +
  scale_fill_viridis_d(option = "C", begin = 0.2, end = 0.85, alpha = 0.65) +
  scale_color_viridis_d(option = "C", begin = 0.2, end= 0.85) +
  labs(x = "Fitness", y = "Generations", 
       title = "Final fitness achieved on Static Epistatic Landscapes (n = 33/Condition)",
       fill = "Loci",
       color = "Loci") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_minimal()

# Moving on to expanding landscape
filepath_df02 <- "output/expanding_epistatic_landscape_20250806_1850.csv"
df02<- read.csv(filepath_df02) %>% 
  mutate("Final_genome_size" = round(Final_genome_size, 0), "Loci" = as.factor(Loci))

# Does initial size correlate with final fitness?
ggplot(df02) +
  geom_boxplot(aes(x = Loci, y = Final_fitness, group = Loci, fill = Loci, color = Loci)) +
  scale_fill_viridis_d(option = "B", begin = 0.25, end = 0.9, alpha = 0.7) +
  scale_color_viridis_d(option = "B", begin = 0.25, end = 0.9, ) + 
  labs(x = "Initial Genome Size", y = "Fitness",
       fill = "Initial Size",
       color = "Initial Size") +
  theme_minimal()
model_size1 <- lm(Final_fitness ~ Loci, df02)
summary(model_size1)
anova(model_size1)

# What about final size?
ggplot(df02, aes(x = Final_genome_size, y = Final_fitness, color = Final_genome_size)) +
  geom_point() +
  stat_smooth(method = "lm", color = "#000000BB", se = FALSE) +
  scale_color_viridis_c(option = "B", begin = 0.25, end = 0.85, alpha = 0.7) +
  labs(x = "Final Genome Size", y = "Fitness",
       fill = "Final Size") +
  theme_minimal()
model_size2 <- lm(Final_fitness ~ Final_genome_size, df02)
summary(model_size2)

# And what about between Initial and final genome size?
ggplot(df02, aes(x = Final_genome_size, y = Loci, fill = Loci, color = Loci)) +
  geom_density_ridges(stat = "binline", binwidth = 1, from = 0, to = 20, scale = 3) + 
  scale_fill_viridis_d(option = "A", begin = 0.25, end = 0.95, alpha = 0.7) +
  scale_color_viridis_d(option = "A", begin = 0.25, end = 0.95) +
  labs(y = "Initial Genome Size", x = "Final Genome Size",
       color = "Initial Size",
       fill = "Initial Size") + 
  theme_minimal()

# It doesn't seem right to me that the largest (20) all conform to the same value. 
filepath_df03 <- "output/expanding_epistatic_landscape_20250807_1234.csv"
df03 <- read.csv(filepath_df03) %>% 
  mutate("Final_genome_size" = round(Final_genome_size, 0), "Loci" = as.factor(Loci))

ggplot(df03, aes(x = Final_genome_size, y = Loci, fill = Loci, color = Loci)) +
  geom_density_ridges(stat = "binline", binwidth = 1, from = 0, to = 20, scale = 7) + 
  scale_fill_viridis_d(option = "A", begin = 0.25, end = 0.95, alpha = 0.7) +
  scale_color_viridis_d(option = "A", begin = 0.25, end = 0.95) +
  labs(y = "Initial Genome Size", x = "Final Genome Size",
       color = "Initial Size",
       fill = "Initial Size") + 
  theme_minimal()
# It doesn't seem to be scale-dependent, so I think it's an artifact of the algorithm somehow?



