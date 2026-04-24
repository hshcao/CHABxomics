rm(list = ls())
library(tidyverse)

# Read the data
df <- read.delim("kaiju_eqR4m_totalNum_genus.tsv", check.names = FALSE)

# Extract domain from taxon_name
df <- df %>%
  mutate(Domain = sub(";.*", "", taxon_name))

# Reshape data: wide -> long
df_long <- df %>%
  pivot_longer(
    cols = -c(taxon_name, Domain),
    names_to = "Sample",
    values_to = "Abundance"
  )

# Extract timepoint and replicate from Sample name
df_long <- df_long %>%
  mutate(
    Time = sub("mG19_(.*)-[0-9]$", "\\1", Sample),
    Replicate = sub(".*-([0-9])$", "\\1", Sample)
  )

# Summarize abundance per domain, time, replicate
df_sum <- df_long %>%
  group_by(Domain, Time, Replicate) %>%
  summarise(TotalAbundance = sum(Abundance), .groups = "drop")

# Order timepoints properly
time_levels <- c("06-15", "07-23", "09-10", "09-29", "11-3", "12-4")
df_sum$Time <- factor(df_sum$Time, levels = time_levels)

# Change replicate names from 1,2,3 → A,B,C
df_sum$Replicate <- recode(df_sum$Replicate, "1"="A", "2"="B", "3"="C")

# Adjustable plotting parameters
line_width <- 0.6   # adjust line thickness
point_size <- 0.9     # adjust point size
facet_label_size <- 6
axis_label_size <- 6
axis_text_size <- 6
axis_title_size <- 5
legend_text_size <- 5
legend_title_size <- 5

# Plot
p <- ggplot(df_sum, aes(x = Time, y = TotalAbundance, 
                        color = Replicate, group = Replicate)) +
  geom_line(size = line_width) +
  geom_point(size = point_size) +
  facet_wrap(~Domain, nrow = 1, scales = "free_y") +   # all domains in one row
  theme_bw() +
  theme(
    strip.text = element_text(size = facet_label_size, color = "black"),   
    axis.text.x = element_text(size = axis_text_size, color = "black", 
                               angle = 90, hjust = 1),
    axis.text.y = element_text(size = axis_text_size, color = "black"),
    axis.title = element_text(size = axis_label_size, color = "black"),
    legend.text = element_text(size = legend_text_size, color = "black"),
    legend.title = element_text(size = legend_title_size, color = "black")
  ) +
  labs(x = "Date", y = "Abundance", color = "Replicate")

# Save
pdf("Abundance_4domains_byDateSite.pdf", width = 10, height = 2)
print(p)
dev.off()


##########################
rm(list = ls())
library(tidyverse)

# Read the data
df <- read.delim("kaiju_eqR4m_totalNum_species.tsv", check.names = FALSE)

# Extract domain and species (species = substring before last ";")
df <- df %>%
  mutate(
    Domain = sub(";.*", "", taxon_name),
    Species = sub(".*;(.*);$", "\\1", taxon_name)
  )

# Reshape to long format
df_long <- df %>%
  pivot_longer(
    cols = -c(taxon_name, Domain, Species),
    names_to = "Sample",
    values_to = "Abundance"
  )

# Extract timepoint and replicate from Sample
df_long <- df_long %>%
  mutate(
    Time = sub("mG19_(.*)-[0-9]$", "\\1", Sample),
    Replicate = sub(".*-([0-9])$", "\\1", Sample)
  )

# Keep only abundance > 0
df_present <- df_long %>% filter(Abundance > 0)

# Count unique species per Domain, Time, Replicate
df_count <- df_present %>%
  group_by(Domain, Time, Replicate) %>%
  summarise(UniqueSpecies = n_distinct(Species), .groups = "drop")

# Order timepoints
time_levels <- c("06-15", "07-23", "09-10", "09-29", "11-3", "12-4")
df_count$Time <- factor(df_count$Time, levels = time_levels)

# Rename replicates (1,2,3 → A,B,C)
df_count$Replicate <- recode(df_count$Replicate, "1"="A", "2"="B", "3"="C")

# ---- Adjustable plotting parameters ----
line_width <- 0.6   # adjust line thickness
point_size <- 0.9     # adjust point size
facet_label_size <- 6
axis_text_size <- 6
axis_title_size <- 5
legend_text_size <- 5
legend_title_size <- 5

# Plot
p <- ggplot(df_count, aes(x = Time, y = UniqueSpecies, 
                          color = Replicate, group = Replicate)) +
  geom_line(size = line_width) +
  geom_point(size = point_size) +
  facet_wrap(~Domain, nrow = 1, scales = "free_y") +   # all domains in one row
  theme_bw() +
  theme(
    strip.text = element_text(size = facet_label_size, color = "black"), # facet labels
    axis.text.x = element_text(size = axis_text_size, color = "black", angle = 90, hjust = 1),
    axis.text.y = element_text(size = axis_text_size, color = "black"),
    axis.title  = element_text(size = axis_title_size, color = "black"),
    legend.text = element_text(size = legend_text_size, color = "black"),
    legend.title= element_text(size = legend_title_size, color = "black")
  ) +
  labs(x = "Date", y = "Species count", color = "Replicate")

# Save
pdf("Number_species_4domains_byDateSite.pdf", 10, 2)
print(p)
dev.off()

##############################################
## draw the abundance of 4 Microcystis species

rm(list = ls())
library(tidyverse)

# Read data
df <- read.delim("Microcystis_4species_abundance_4plot.tsv", check.names = FALSE)

# Reshape into long format
df_long <- df %>%
  pivot_longer(
    -taxon_name,
    names_to = "sample",
    values_to = "abundance"
  ) %>%
  mutate(
    # Extract timepoint (before last "-")
    Date = str_replace(sample, "mG19_", "") %>%
      str_remove("-[123]$"),
    # Extract replicate number
    Replicate = str_extract(sample, "[123]$"),
    # Label replicates as sites A, B, C
    Site = recode(Replicate, "1" = "A", "2" = "B", "3" = "C"),
    # Order dates in chronological order
    Date = factor(Date, levels = c("06-15","07-23","09-10","09-29","11-3","12-4"))
  )

# Plot function with adjustable sizes
plot_taxa <- function(line_width = 0.5, text_size = 6) {
  ggplot(df_long, aes(x = Date, y = abundance, group = Site, color = Site)) +
    geom_line(linewidth = line_width) +
    geom_point(size = line_width + 0.5) +
    facet_wrap(~ taxon_name, scales = "free_y", ncol = 4) +
    theme_bw(base_size = text_size) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      strip.text = element_text(size = text_size)
    ) +
    labs(x = "Date", y = "Abundance", color = "Replicate")
}

# Example: run with thicker lines and larger text
pdf("Microcystis_4species_abundance_lineplot_9x15.pdf", 9, 1.5)
plot_taxa(line_width = 0.8, text_size = 6)
dev.off()
##############################################
## draw the abundance of epibiotic bacteria

rm(list = ls())
library(tidyverse)
library(gridExtra)

# -------------------------------
# Load species abundance table
# -------------------------------
df <- read.delim("kaiju_eqR4m_totalNum_species.tsv", sep = "\t", check.names = FALSE)
colnames(df)[1] <- "taxon"

# -------------------------------
# Define epibiont genera
# -------------------------------
genera <- c("Limnobacter", "Variovorax", "Bosea", "Pseudoxanthomonas",
            "Rhizobium", "Agrobacterium", "Xanthomonadaceae", "Brevundimonas",
            "Comamonadaceae", "Phyllobacteriaceae", "Pseudomonas",
            "Sphingopyxis", "Flavobacterium", "Hyphomonas", "Flavihumibacter")

# -------------------------------
# Collapse species -> genus level sums
# -------------------------------
df_genus <- map_dfr(genera, function(gen) {
  df %>%
    filter(str_detect(taxon, gen)) %>%
    select(-taxon) %>%
    summarise(across(everything(), sum)) %>%
    mutate(genus = gen)
})

# Reshape to long format
df_long <- df_genus %>%
  pivot_longer(-genus, names_to = "sample", values_to = "abundance") %>%
  mutate(
    Date = str_replace(sample, "mG19_", "") %>% str_remove("-[123]$"),
    Replicate = str_extract(sample, "[123]$"),
    Site = recode(Replicate, "1" = "A", "2" = "B", "3" = "C"),
    Date = factor(Date, levels = c("06-15","07-23","09-10","09-29","11-3","12-4"))
  )

# -------------------------------
# Create individual plots
# -------------------------------
plot_list <- map(genera, function(gen) {
  df_plot <- df_long %>% filter(genus == gen)
  
  ggplot(df_plot, aes(x = Date, y = abundance, group = Site, color = Site)) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1) +
    theme_bw(base_size = 6) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.title.x = element_text(size = 6),
      axis.title.y = element_text(size = 6)
    ) +
    labs(x = "Date", y = gen, color = "Site")
})

# -------------------------------
# Arrange plots in 4x4 grid
# -------------------------------
pdf("epibiotic_bacteria_abundance_lineplot_11x6.pdf",11,6)
do.call(grid.arrange, c(plot_list, ncol = 4, nrow = 4))
dev.off()
