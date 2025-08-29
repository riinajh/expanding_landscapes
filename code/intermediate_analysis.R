library(tidyverse)
library(ggridges)
library(ggpmisc)
setwd("~/mastersarbeit/thesis/code")

# Function for the solution from Suman that yields the expectation for the distribution
# of fitness peaks on a landscape
expected_max_fitness <- function(L, add, epi){
  return(exp((2 * L * (((L / 2) * add^2) + epi^2) * log(2))^0.5))
}

################################################################################
# Begin with the static landscape
################################################################################

filepath_df01 <- "outputs/data/"
df01 <- read.csv(filepath_df01) %>% 
  mutate("percent_max" = Final_fitness * 100 / expected_max_fitness(Loci, 0, 1))

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

ggplot(filter(df02, Final_genome_size == 6), aes(x = as.factor(Loci), y = Final_fitness, group = as.factor(Loci))) +
  geom_boxplot()

################################################################################
# Looking at dynamics (through the series of selective sweeps) of expanding landscapes
################################################################################

filepath_df03 <- "outputs/data/10_dynamics_intermediate_expanding_landscape_20250828_1702.csv"
df03 <- read.csv(filepath_df03) %>% 
  mutate("FinalGenomeSize" = round(FinalGenomeSize, 0),
         "CurrentGenomeSize" = round(CurrentGenomeSize, 0)) %>% 
  filter(FinalGenomeSize >= 5) %>% # Get rid of everything that decreased genome size
  group_by(Replicate) %>% 
  mutate("Endpoint" = c(rep("No", n()-1), "Yes"))

print(length(unique(df03$Replicate))) # 888 / 1000 replicates increased in size (88.8%)

nonmonotonic_reps <- df03 %>%
  group_by(Replicate) %>%
  arrange(Step, .by_group = TRUE) %>%
  mutate(
    diff = sign(CurrentGenomeSize - lag(CurrentGenomeSize)),
    fitness_diff = sign(Fitness - lag(Fitness))
  ) %>% 
  filter(diff < 0 & fitness_diff > 0) # tunneling is when fitness increases from decreasing genome size in these data

print(length(unique(nonmonotonic_reps$Replicate))) # 279 / 888 exhibit tunneling behavior (31.4%)

# Exclude these from the dataset
df03_no_tunneling <- df03 %>%
  filter(!(Replicate %in% nonmonotonic_reps$Replicate)) %>% 
  group_by(Replicate, CurrentGenomeSize) %>% 
  filter(row_number() == n()) # filter out local landscape adaptation (only the last instance of a selective sweep with that genome size)

# Ok now are all the distributions of fitnesses at the same size equal, regardless of endpoint?

ggplot(df03_no_tunneling, aes(y = Fitness, x = Endpoint)) +
  geom_boxplot() +
  facet_wrap(df03_no_tunneling$CurrentGenomeSize)

# NO THEY ARE NOT!
# (discounting 9, 10, 11 here due to low sample size)

# What about when including the tunneling ones?

ggplot(df03, aes(y = Fitness, x = Endpoint)) +
  geom_boxplot() + 
  facet_wrap(df03$CurrentGenomeSize)

# alternative plot showing fitness paths over genome sizes

params <- tibble(Parameter = c("Epistatic σ", "Additive σ", "μ", "M", "Population", "Generations"), 
                 c(0.1, 0.1,1e-3, 1e-4, 1000, 5e4),
                 tb = list(tb))
ggplot(df03_no_tunneling, aes(y = Fitness, x = CurrentGenomeSize, group = Replicate, color = as.factor(FinalGenomeSize))) +
  geom_line() +
  scale_color_viridis_d(option = "F", begin = 0.05, end = 0.85) +
  labs(x = "Current Genome Size",
       title = "Fitness Increasing Paths of Expanding Landscapes") +
  theme(legend.position = "none") +
  geom_table(data = params, mapping = aes(x = x, y = y, label = "tb"), hjust = 0, vjust = 1) +
  theme_minimal()

