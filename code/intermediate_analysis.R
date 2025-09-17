library(tidyverse)
library(ggridges)
setwd("~/mastersarbeit/thesis/code")

# Function for the solution from Suman that yields the expectation for the distribution
# of fitness peaks on a landscape
expected_max_fitness <- function(L, add, epi){
  return(exp((2 * L * (((L / 2) * add^2) + epi^2) * log(2))^0.5))
}

################################################################################
# Now the expanding landscape
################################################################################

filepath_df02 <- "outputs/data/09_intermediate_expanding_landscape_20250827_2339.csv"
df02 <- read.csv(filepath_df02) %>% 
  mutate("percent_max" = Final_fitness * 100 / expected_max_fitness(Final_genome_size, 0.1, 0.1),
         "Final_genome_size" = round(Final_genome_size, 0))

# On a max 20 loci landscape, we started replicates at 5, 6, 7, 8, 9, & 10 active loci.
# What's the distribution of final sizes?

ggplot(df02, aes(x = Final_genome_size)) +
  geom_histogram(binwidth = 1)

# 7, 8, 9 are all good candidates. 10 will let us analyze the nonchanging genomes as well.

ggplot(filter(df02, Final_genome_size == 7), aes(x = as.factor(Loci), y = Final_fitness, group = as.factor(Loci))) +
  geom_boxplot()

################################################################################
# Looking at dynamics (through the series of selective sweeps) of expanding landscapes
################################################################################

filepath_df03.1 <- "outputs/data/10_dynamics_intermediate_expanding_landscape_20250828_1702_1kreps.csv"
filepath_df03.2 <- "outputs/data/10_dynamics_intermediate_expanding_landscape_20250831_2329_2kreps.csv"
filepath_df03.3 <- "outputs/data/10_dynamics_intermediate_expanding_landscape_20250903_1422_3kreps.csv"
filepath_df03.4 <- "outputs/data/10_dynamics_intermediate_expanding_landscape_20250904_0818_5kreps.csv"

df03 <- read.csv(filepath_df03.1) %>% 
  rbind(read.csv(filepath_df03.2)) %>% 
  rbind(read.csv(filepath_df03.3)) %>% 
  rbind(read.csv(filepath_df03.4)) %>% 
  mutate("FinalGenomeSize" = round(FinalGenomeSize, 0),
         "CurrentGenomeSize" = round(CurrentGenomeSize, 0)) %>% 
  filter(FinalGenomeSize >= 5) %>% # Get rid of everything that decreased genome size
  group_by(Replicate) %>% 
  mutate("Endpoint" = c(rep("No", n()-1), "Yes"))

print(length(unique(df03$Replicate))) # 4460 / 5000 replicates increased in size (89.2%)

# we now create a filter for simulations that display nonmonotonic genome size change
# over the course of the simulation
nonmonotonic_reps <- df03 %>%
  group_by(Replicate) %>%
  arrange(Step, .by_group = TRUE) %>%
  mutate(diff = sign(CurrentGenomeSize - lag(CurrentGenomeSize)),
         fitness_diff = sign(Fitness - lag(Fitness))) %>% 
  filter(diff < 0 & fitness_diff > 0) # tunneling is when fitness increases from decreasing genome size in these data

print(length(unique(nonmonotonic_reps$Replicate))) # 1402 / 4460 exhibit tunneling behavior (31.4%)

# Exclude these from the dataset (it just makes things harder. I include them later to test)
df03_no_tunneling <- df03 %>%
  filter(!(Replicate %in% nonmonotonic_reps$Replicate)) %>% 
  group_by(Replicate, CurrentGenomeSize) %>% 
  filter(row_number() == n()) # filter out local landscape adaptation (only the last instance of a selective sweep with that genome size)

# Ok now are all the distributions of fitnesses at the same size equal, regardless of endpoint?
ggplot(data = filter(df03_no_tunneling, CurrentGenomeSize != 4), aes(y = Fitness, x = Endpoint, fill = CurrentGenomeSize, color = CurrentGenomeSize)) +
  geom_boxplot() +
  facet_wrap(vars(CurrentGenomeSize)) +
  scale_fill_viridis_c(option = "A", begin = 0.2, end = 0.9, alpha = 0.65, guide = "none") +
  scale_color_viridis_c(option = "A", begin = 0.2, end= 0.9) +
  labs(x = "At endpoint at current genome size?",
       fill = "Genome size",
       color = "Genome size") +
  theme_minimal()
ggsave("outputs/figures/endpoint_boxplots.pdf")
# NO THEY ARE NOT!

# What about timing?
ggplot(data = filter(df03_no_tunneling, CurrentGenomeSize != 4), aes(y = Step, x = Endpoint, fill = CurrentGenomeSize, color = CurrentGenomeSize)) +
  geom_boxplot() +
  facet_wrap(vars(CurrentGenomeSize)) +
  scale_fill_viridis_c(option = "A", begin = 0.2, end = 0.9, alpha = 0.65, guide = "none") +
  scale_color_viridis_c(option = "A", begin = 0.2, end= 0.9) +
  labs(x = "At endpoint at current genome size?",
       fill = "Genome size",
       color = "Genome size") +
  theme_minimal()
ggsave("outputs/figures/timings.pdf")

# What about when including the tunneling ones (which do not have a monotonically increasing genome size)?
# Actually, the same?!
ggplot(df03, aes(y = Fitness, x = Endpoint, fill = CurrentGenomeSize, color = CurrentGenomeSize)) +
  geom_boxplot() +
  facet_wrap(vars(CurrentGenomeSize)) +
  scale_fill_viridis_c(option = "A", begin = 0.2, end = 0.9, alpha = 0.65, guide = "none") +
  scale_color_viridis_c(option = "A", begin = 0.2, end= 0.9) +
  labs(x = "At endpoint at current genome size?",
       fill = "Genome size",
       color = "Genome size") +
  theme_minimal()
ggsave("outputs/figures/endpoint_boxplots_unfiltered.pdf")


# alternative plot showing genomic trajectories across different sizes, with fitness values
ggplot(filter(df03_no_tunneling, CurrentGenomeSize != 4), aes(y = Fitness, x = CurrentGenomeSize,
        group = Replicate, color = as.factor(FinalGenomeSize))) +
  geom_line() +
  scale_color_viridis_d(option = "F", begin = 0.05, end = 0.85) +
  labs(x = "Current Genome Size",
       #title = "Fitness Trajectories of Expanding Landscapes",
       color = "Final genome size") +
  theme(legend.position = "none") +
  theme_minimal()
ggsave("outputs/figures/genome_trajectories.pdf")

ggplot(df03, aes(y = Fitness, x = CurrentGenomeSize,
        group = Replicate, color = as.factor(FinalGenomeSize))) +
  geom_line() +
  scale_color_viridis_d(option = "F", begin = 0.05, end = 0.85) +
  labs(x = "Current Genome Size",
       #title = "Fitness Trajectories of Expanding Landscapes",
       color = "Final genome size") +
  theme(legend.position = "none") +
  theme_minimal()
ggsave("outputs/figures/genome_trajectories_unfiltered.pdf")

################################################################################
# Comparing with the same timeseries but on HoC
################################################################################

filepath_df04 <- "outputs/data/10-1_dynamics_hoc_expanding_landscape_20250905_2218.csv"

df04 <- read.csv(filepath_df04) %>% 
  mutate("FinalGenomeSize" = round(FinalGenomeSize, 0),
         "CurrentGenomeSize" = round(CurrentGenomeSize, 0)) %>% 
  filter(FinalGenomeSize >= 5) %>% # Get rid of everything that decreased genome size
  group_by(Replicate) %>% 
  mutate("Endpoint" = c(rep("No", n()-1), "Yes"))

print(length(unique(df04$Replicate))) # 4343 / 5000 replicates increased in size (86.7%)

# we now create a filter for simulations that display nonmonotonic genome size change
# over the course of the simulation
nonmonotonic_reps <- df04 %>%
  group_by(Replicate) %>%
  arrange(Step, .by_group = TRUE) %>%
  mutate(diff = sign(CurrentGenomeSize - lag(CurrentGenomeSize)),
         fitness_diff = sign(Fitness - lag(Fitness))) %>% 
  filter(diff < 0 & fitness_diff > 0) # tunneling is when fitness increases from decreasing genome size in these data

print(length(unique(nonmonotonic_reps$Replicate))) # 1021 / 4343 exhibit tunneling behavior (23.5%)

# Exclude these from the dataset (it just makes things harder. I include them later to test)
df04_no_tunneling <- df04 %>%
  filter(!(Replicate %in% nonmonotonic_reps$Replicate)) %>% 
  group_by(Replicate, CurrentGenomeSize) %>% 
  filter(row_number() == n()) # filter out local landscape adaptation (only the last instance of a selective sweep with that genome size)

ggplot(data = filter(df04_no_tunneling, CurrentGenomeSize != 4), aes(y = Fitness, x = Endpoint, fill = CurrentGenomeSize, color = CurrentGenomeSize)) +
  geom_boxplot() +
  facet_wrap(vars(CurrentGenomeSize)) +
  scale_fill_viridis_c(option = "A", begin = 0.2, end = 0.9, alpha = 0.65, guide = "none") +
  scale_color_viridis_c(option = "A", begin = 0.2, end= 0.9) +
  labs(x = "At endpoint at current genome size?",
       fill = "Genome size",
       color = "Genome size") +
  theme_minimal()
ggsave("outputs/figures/hoc_boxplots.pdf")
# effect magnitude seems greater, actually