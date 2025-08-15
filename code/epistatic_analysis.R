library(tidyverse)
library(ggridges)
setwd("~/mastersarbeit/thesis/code")

# Function for the solution from Suman that yields the expectation for the distribution
# of fitness peaks on a landscape
expected_max_fitness <- function(L, add, epi){
  return(exp((2 * L * (((L / 2) * add^2) + epi^2) * log(2))^0.5))
}

################################################################################
# Begin with the static landscape
################################################################################

filepath_df01 <- "outputs/data/05_static_epistatic_landscape_20250808_1624.csv"
df01 <- read.csv(filepath_df01) %>% 
  mutate("percent_max" = Final_fitness * 100 / expected_max_fitness(Loci, 0, 1))

ggplot(df01, aes(x = Final_fitness, y = as.factor(Generation), fill = as.factor(Loci), color = as.factor(Loci))) +
  facet_wrap(df01$Loci) +
  geom_density_ridges(scale = 3, from = 0) +
  scale_fill_viridis_d(option = "C", begin = 0.2, end = 0.85, alpha = 0.65) +
  scale_color_viridis_d(option = "C", begin = 0.2, end= 0.85) +
  labs(x = "Fitness", y = "Generations", 
       title = "Fitness Achieved on Static Uncorrelated Landscapes over time (99 Replicates)",
       fill = "Loci",
       color = "Loci") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_minimal()
ggsave("outputs/figures/05-1_static_timeseries.png")

# Also plotting the final fitnesses vs size for comparison
# Filter only for endpoints as to compare to later simulations
ggplot(filter(df01, Generation == 50000), aes(x = as.factor(Loci), y = Final_fitness, group = Loci,
                                              fill = as.factor(Loci), color = as.factor(Loci))) +
  geom_boxplot() +
  lims(y = c(0,200)) +
  scale_fill_viridis_d(option = "B", begin = 0.25, end = 0.9, alpha = 0.7) +
  scale_color_viridis_d(option = "B", begin = 0.25, end = 0.9) + 
  labs(x = "Genome Size", y = "Fitness",
       title = "Final Fitnesses of Static Uncorrelated Landscapes (99 Replicates)",
       fill = "Genome Size (Fixed)",
       color = "Genome Size (Fixed)") +
  theme_minimal()
ggsave("outputs/figures/05-2_static_finalfitness.png")

################################################################################
# Moving on to expanding landscape (max 20, c(1,5,10,15,20))
################################################################################

filepath_df02 <- "outputs/data/06_expanding_epistatic_landscape_20250807_2245.csv"
df02<- read.csv(filepath_df02) %>% 
  mutate("Final_genome_size" = round(Final_genome_size, 0), "Loci" = as.factor(Loci))

# Does initial size correlate with final fitness?
ggplot(df02) +
  geom_boxplot(aes(x = Loci, y = Final_fitness, group = Loci, fill = Loci, color = Loci)) +
  lims(y = c(0,200)) +
  scale_fill_viridis_d(option = "B", begin = 0.25, end = 0.9, alpha = 0.7) +
  scale_color_viridis_d(option = "B", begin = 0.25, end = 0.9) + 
  labs(x = "Initial Genome Size", y = "Fitness",
       title = "Initial Genome Size vs. Final Fitness (200 Replicates)",
       fill = "Initial Genome Size",
       color = "Initial Genome Size") +
  theme_minimal()
ggsave("outputs/figures/06-1_expanding_initialsize_fitness.png")
model_size1 <- lm(Final_fitness ~ Loci, df02)
summary(model_size1)
anova(model_size1)

# What about final size?
# Unclear. Boxplot shows the correlation best
ggplot(df02, aes(x = as.factor(Final_genome_size), y = Final_fitness, 
                 color = Final_genome_size, fill = Final_genome_size,
                 group = Final_genome_size)) +
  geom_boxplot() +
  #stat_summary(fun = mean, geom = "crossbar", width = 0.2) + 
  lims(y = c(0,200)) +
  #stat_smooth(method = "lm", color = "#000000BB", se = FALSE) +
  scale_fill_viridis_c(option = "B", begin = 0.25, end = 0.9, alpha = 0.7) +
  scale_color_viridis_c(option = "B", begin = 0.25, end = 0.9, alpha = 0.7) +
  labs(x = "Final Genome Size", y = "Fitness", 
       fill = "Final Genome Size", color = "Final Genome Size",
       title = "Final Genome Size vs. Final Fitness (n = 1000)") +
  theme_minimal()
ggsave("outputs/figures/06-2_expanding_finalsizefitness.png")
model_size2 <- lm(Final_fitness ~ Final_genome_size, filter(df02, !(Final_genome_size %in% c(1,19))))
summary(model_size2) # even if you remove 1 and 19 as outliers it still shows significance
anova(model_size2)

# And what about between Initial and final genome size?
ggplot(df02, aes(x = Final_genome_size, y = Loci, fill = Loci, color = Loci)) +
  geom_density_ridges(stat = "binline", binwidth = 1, from = 0, to = 20, scale = 3) + 
  stat_summary(fun = mean) +
  scale_fill_viridis_d(option = "A", begin = 0.25, end = 0.95, alpha = 0.7) +
  scale_color_viridis_d(option = "A", begin = 0.25, end = 0.95) +
  labs(y = "Initial Genome Size", x = "Final Genome Size",
       title = "Initial vs. Final Genome Size (200 Replicates)",
       color = "Initial Genome Size",
       fill = "Initial Genome Size") + 
  theme_minimal()
ggsave("outputs/figures/06-3_expanding_initialfinalsize.png")

################################################################################
# Now trying with a slightly different set of sizes (max 100, c(1,25,50,75,100))
################################################################################

filepath_df03 <- "outputs/data/06_expanding_epistatic_landscape_20250811_0947.csv"
df03 <- read.csv(filepath_df03) %>% 
  mutate("Final_genome_size" = round(Final_genome_size, 0), "Loci" = as.factor(Loci))

# Does initial size correlate with final fitness?
ggplot(df03) +
  geom_boxplot(aes(x = Loci, y = Final_fitness, group = Loci, fill = Loci, color = Loci)) +
  lims(y = c(0,200)) +
  scale_fill_viridis_d(option = "B", begin = 0.25, end = 0.9, alpha = 0.7) +
  scale_color_viridis_d(option = "B", begin = 0.25, end = 0.9) + 
  labs(x = "Initial Genome Size", y = "Fitness",
       title = "Initial Genome Size vs. Final Fitness (200 Replicates)",
       fill = "Initial Genome Size",
       color = "Initial Genome Size") +
  theme_minimal()
ggsave("outputs/figures/06-4_expanding_initialsizefitness_2.png")
model_size1 <- lm(Final_fitness ~ Loci, df03)
summary(model_size1)
anova(model_size1)

# What about final size?
# Unclear. Boxplot shows the correlation best
ggplot(df03, aes(x = as.factor(Final_genome_size), y = Final_fitness, 
                 color = Final_genome_size, fill = Final_genome_size,
                 group = Final_genome_size)) +
  geom_boxplot() +
  #stat_summary(fun = mean, geom = "crossbar", width = 0.2) + 
  lims(y = c(0,200)) +
  #stat_smooth(method = "lm", color = "#000000BB", se = FALSE) +
  scale_fill_viridis_c(option = "B", begin = 0.25, end = 0.9, alpha = 0.7) +
  scale_color_viridis_c(option = "B", begin = 0.25, end = 0.9, alpha = 0.7) +
  labs(x = "Final Genome Size", y = "Fitness", 
       fill = "Final Genome Size", color = "Final Genome Size",
       title = "Final Genome Size vs. Final Fitness (n = 1000)") +
  theme_minimal()
ggsave("outputs/figures/06-5_expanding_finalsizefitness_2.png")
model_size2 <- lm(Final_fitness ~ Final_genome_size, df03)
summary(model_size2) # even if you remove 1 and 19 as outliers it still shows significance
anova(model_size2)

# And what about between initial and final genome size?
ggplot(df03, aes(x = Final_genome_size, y = Loci, fill = Loci, color = Loci)) +
  geom_density_ridges(stat = "binline", binwidth = 1, from = 0, to = 20, scale = 3) + 
  stat_summary(fun = mean) +
  scale_fill_viridis_d(option = "A", begin = 0.25, end = 0.95, alpha = 0.7) +
  scale_color_viridis_d(option = "A", begin = 0.25, end = 0.95) +
  labs(y = "Initial Genome Size", x = "Final Genome Size",
       title = "Initial vs. Final Genome Size (200 Replicates)",
       color = "Initial Genome Size",
       fill = "Initial Genome Size") + 
  theme_minimal()
ggsave("outputs/figures/06-6_initialfinalsize_2.png")

# Running the 75 and 100 loci simulations for 200k generations because I have a
# suspicion that 50k isn't enough, and this may explain why the upper ranges have
# lower fitnesses than expected from the HoC
# filepath_df04 <- "outputs/data/06_expanding_epistatic_landscape_20250811_2253.csv"
# df04 <- read.csv(filepath_df04) %>% 
#   mutate("Final_genome_size" = round(Final_genome_size, 0), "Loci" = as.factor(Loci))
# 
# # Replace values in df03 with those from df04
# df04 <- rbind(filter(df03, !(Loci %in% c(75, 100))), df04)
# # Does initial size correlate with final fitness?
# ggplot(df04) +
#   geom_boxplot(aes(x = Loci, y = Final_fitness, group = Loci, fill = Loci, color = Loci)) +
#   lims(y = c(0,200)) +
#   scale_fill_viridis_d(option = "B", begin = 0.25, end = 0.9, alpha = 0.7) +
#   scale_color_viridis_d(option = "B", begin = 0.25, end = 0.9) + 
#   labs(x = "Initial Genome Size", y = "Fitness",
#        fill = "Initial Size",
#        color = "Initial Size",
#        title = "Initial Genome Size vs. Final Fitness (200k generations)") +
#   theme_minimal()
# Yes, but this is no more robust that the original analysis. Since in both cases
# we only ran for 50k generations, the static landscapes could also benefit from more fitness gains at longer timescales.
# we can justify this in saying that - the distribution of fitness increases vs time is likely to be similar between landscapes,
# so running for 50k generations is sufficient to see the pattern of interest, since most of the adaptation happens early

################################################################################
# Returning to the skewing question:
################################################################################

max20_sizes <- c(1, 5, 10, 15, 20)
max20_means <- sapply(max20_sizes, function(x) mean(filter(df02, Loci == x)$Final_genome_size))
max20_skews <- (max20_means - max20_sizes)# / 20

max100_sizes <- c(1, 25, 50, 75, 100)
max100_means <- sapply(max100_sizes, function(x) mean(filter(df03, Loci == x)$Final_genome_size))
max100_skews <- (max100_means - max100_sizes)# / 100

plot(x = max20_sizes, y = max20_skews)
summary(lm(max20_skews ~ max20_means))
plot(x = max100_sizes, y = max100_skews)
summary(lm(max100_skews ~ max100_means))
# the slope on 100 loci is suspiciously close to 1/5 of that on 20 loci

################################################################################
# Are different fitnesses achieved on static vs. expanding landscapes?
################################################################################

df_static <- filter(df01, Generation == 50000)
df03
# need to rerun 10 for the expanding landscape?

