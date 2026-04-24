################################################################################
# ----------------------------------------------------------------
# 1. TEMPERATURE PLOT (5 x 3 inches) - Bar plot of monthly means ± SD
##### plot 1-year Microcystis aeruginosa abundance in Taihu

# Required packages (install once if needed)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# ----------------------------------------------------------------
# 1. TEMPERATURE PLOT (5 x 3 inches)
# ----------------------------------------------------------------
temp_data <- read_excel("/Users/hc284/Documents/_0taihuXomicsMicrobiome/fieldCellCounts/algae-Temp_taihu-1yr.xlsx", 
                        sheet = "temperature")

temp_summary <- temp_data %>%
  filter(Month %in% 6:10) %>%
  group_by(Month) %>%
  summarise(mean_temp = mean(Temperature, na.rm = TRUE),
            sd_temp = sd(Temperature, na.rm = TRUE)) %>%
  mutate(Month_label = factor(Month, levels = 6:10,
                              labels = c("06-15", "07-15", "08-15", "09-15", "10-15")))

theme_nature <- function(base_size = 10) {
  theme(text = element_text(family = "Arial", size = base_size),
        axis.title = element_text(size = base_size + 1),
        axis.text = element_text(size = base_size, color = "black"),
        axis.line = element_line(linewidth = 0.5, color = "black"),
        axis.ticks = element_line(linewidth = 0.5, color = "black"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        legend.position = "none")
}

p_temp <- ggplot(temp_summary, aes(x = Month_label, y = mean_temp)) +
  geom_col(fill = "#56B4E9", color = "black", linewidth = 0, width = 0.2) +
  geom_errorbar(aes(ymin = mean_temp - sd_temp, ymax = mean_temp + sd_temp),
                width = 0.1, linewidth = 0.25, color = "black") +
  labs(x = NULL, y = "Temperature (°C)") +
  scale_y_continuous(breaks = pretty_breaks(n = 5)) +
  theme_nature(base_size = 11)
p_temp
ggsave("1-yr_taihu_temperature_bar5x2.pdf", p_temp, width = 5, height = 2, units = "in", dpi = 600, device = cairo_pdf)

# ----------------------------------------------------------------
# 2. MICROCYSTIS AERUGINOSA PLOT (5 x 5 inches)
# ----------------------------------------------------------------
algae_data <- read_excel("/Users/hc284/Documents/_0taihuXomicsMicrobiome/fieldCellCounts/algae-Temp_taihu-1yr.xlsx", 
                         sheet = "algal_cell_count")

if (!inherits(algae_data$Date, "Date")) {
  algae_data <- algae_data %>% mutate(Date = as.Date(Date, origin = "1899-12-30"))
}

micro_data <- algae_data %>%
  filter(Species == "铜绿微囊藻Microcystis aeruginosa Kütz.",
         Date %in% as.Date(c("2019-06-15", "2019-07-15", "2019-08-15", 
                             "2019-09-15", "2019-10-15"))) %>%
  select(Date, S1:S12) %>%
  pivot_longer(cols = S1:S12, names_to = "Site", values_to = "Cell_count") %>%
  # === THE FIX: force numerical order S1–S12 ===
  mutate(Site = factor(Site, levels = paste0("S", 1:12))) %>%
  mutate(Cell_density = log(Cell_count))

theme_nature_lines <- function(base_size = 10) {
  theme(text = element_text(family = "Arial", size = base_size),
        axis.title = element_text(size = base_size + 1),
        axis.text = element_text(size = base_size, color = "black"),
        axis.line = element_line(linewidth = 0.5, color = "black"),
        axis.ticks = element_line(linewidth = 0.5, color = "black"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = base_size),
        legend.text = element_text(size = base_size - 1),
        legend.key = element_blank())
}

p_algae <- ggplot(micro_data, aes(x = Date, y = Cell_density, color = Site, group = Site)) +
  geom_line(linewidth = 0.85) +
  geom_point(size = 1.2) +
  scale_x_date(breaks = as.Date(c("2019-06-15", "2019-07-15", "2019-08-15", 
                                  "2019-09-15", "2019-10-15")),
               labels = c("06-15", "07-15", "08-15", "09-15", "10-15")) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, 
       y = expression(paste("Log(cell density) (cells L"^"-1"*")")),
       color = "Site") +
  theme_nature_lines(base_size = 11)
p_algae
ggsave("1-yr_taihu_Microcystis_bar5x4.pdf", p_algae, width = 5, height = 4, 
       units = "in", dpi = 600, device = cairo_pdf)

####### bar and lines together one on top of the other 
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# ----------------------------------------------------------------
# 1. LOAD BOTH SHEETS
# ----------------------------------------------------------------
file_path <- "/Users/hc284/Documents/_0taihuXomicsMicrobiome/fieldCellCounts/algae-Temp_taihu-1yr.xlsx"

temp_data <- read_excel(file_path, sheet = "temperature")

algae_data <- read_excel(file_path, sheet = "algal_cell_count")
if (!inherits(algae_data$Date, "Date")) {
  algae_data <- algae_data %>% mutate(Date = as.Date(Date, origin = "1899-12-30"))
}

# ----------------------------------------------------------------
# 2. MICROCYSTIS DATA (12 sites, 5 dates, log scale)
# ----------------------------------------------------------------
target_dates <- as.Date(c("2019-06-15", "2019-07-15", "2019-08-15",
                          "2019-09-15", "2019-10-15"))

micro_data <- algae_data %>%
  filter(Species == "铜绿微囊藻Microcystis aeruginosa Kütz.",
         Date %in% target_dates) %>%
  select(Date, S1:S12) %>%
  pivot_longer(cols = S1:S12, names_to = "Site", values_to = "Cell_count") %>%
  mutate(Site = factor(Site, levels = paste0("S", 1:12)),
         Cell_density = log(Cell_count))

# ----------------------------------------------------------------
# 3. TEMPERATURE SUMMARY + DUAL-AXIS SCALING
# ----------------------------------------------------------------
temp_summary <- temp_data %>%
  filter(Month %in% 6:10) %>%
  group_by(Month) %>%
  summarise(mean_temp = mean(Temperature, na.rm = TRUE),
            sd_temp     = sd(Temperature, na.rm = TRUE)) %>%
  mutate(Date = target_dates)

y_max <- max(micro_data$Cell_density, na.rm = TRUE) * 1.02
p_min <- 12
s_min <- 20
t_max <- max(temp_summary$mean_temp, na.rm = TRUE)
s_max <- t_max

factor    <- (s_max - s_min) / (y_max - p_min)
intercept <- s_min - p_min * factor

temp_summary <- temp_summary %>%
  mutate(scaled_mean = p_min + (mean_temp - s_min) * (y_max - p_min) / (s_max - s_min),
         scaled_sd   = sd_temp * (y_max - p_min) / (s_max - s_min))

sec_transform <- function(x) x * factor + intercept

# ----------------------------------------------------------------
# 4. NATURE-STYLE THEME
# ----------------------------------------------------------------
theme_nature <- function(base_size = 12) {
  theme(text = element_text(family = "Arial", size = base_size),
        axis.title = element_text(size = base_size + 1, face = "bold"),
        axis.text = element_text(size = base_size, color = "black"),
        axis.line = element_line(linewidth = 0.6, color = "black"),
        axis.ticks = element_line(linewidth = 0.6, color = "black"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = base_size),
        legend.text = element_text(size = base_size - 1),
        legend.key = element_blank(),
        plot.margin = margin(5, 15, 8, 5),   # extra bottom margin for x-title
        axis.title.y.right = element_text(color = "red"),
        axis.text.y.right  = element_text(color = "red"),
        axis.ticks.y.right = element_line(color = "red"))
}

# ----------------------------------------------------------------
# 5. COMBINED PLOT (x-axis title added)
# ----------------------------------------------------------------
p_combined <- ggplot() +
  
  # Primary: 12 Microcystis lines + points
  geom_line(data = micro_data,
            aes(x = Date, y = Cell_density, color = Site, group = Site),
            linewidth = 1.0) +
  geom_point(data = micro_data,
             aes(x = Date, y = Cell_density, color = Site),
             size = 1.5) +
  
  # Secondary: Temperature bars (no borders, transparent blue)
  geom_col(data = temp_summary,
           aes(x = Date, y = scaled_mean),
           fill = "red", alpha = 0.2, color = NA, linewidth = 0, width = 10) +
  
  # SD error bars (same color as bars)
  geom_errorbar(data = temp_summary,
                aes(x = Date, ymin = scaled_mean - scaled_sd,
                    ymax = scaled_mean + scaled_sd),
                width = 2.5, linewidth = 0.2, alpha = 0.5, color = "red") +
  
  # Axes
  scale_x_date(breaks = target_dates,
               labels = c("06-15", "07-15", "08-15", "09-15", "10-15")) +
  scale_y_continuous(
    name = expression(paste("log (cell density) (cells L"^-1*")")),
    sec.axis = sec_axis(transform = sec_transform, name = "Temperature (°C)")
  ) +
  
  # View clipping (primary starts at 12, bars stay visible)
  coord_cartesian(ylim = c(12, NA)) +
  
  labs(x = "Sampling date (2019)", color = "Site") +
  theme_nature(base_size = 12) +
  theme(legend.key.size = unit(0.5, "lines"))

print(p_combined)
# ----------------------------------------------------------------
# 6. SAVE AS 5 × 5 INCH PDF
# ----------------------------------------------------------------
ggsave("1-yr_taihu_Microcystis_with_temperature_final_5x5.pdf",
       p_combined, width = 5, height = 5, units = "in",
       dpi = 600, device = cairo_pdf)

################################################################################
##### 1/2 year temperature vs cyanobacterial abundance
# ============================================================

library(readr)
library(dplyr)
library(ggplot2)

# ---------------------------
# Read updated data
# ---------------------------
df <- read_csv("Taihu_Shi_7yr/PCmonthly_7yr.csv", show_col_types = FALSE)

# ---------------------------
# Set x-axis order exactly as it appears in file
# ---------------------------
df <- df %>%
  mutate(
    `Temperature interval` = factor(
      `Temperature interval`,
      levels = unique(`Temperature interval`)
    )
  )

# ---------------------------
# Colors
# ---------------------------
bar_col  <- "cyan4"
line_col <- "cyan2"

# ---------------------------
# Nature-like theme
# Slightly larger for legibility after shrinking
# ---------------------------
theme_nature_big <- function(base_size = 11, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.text.x = element_text(
        angle = 45, hjust = 1, vjust = 1,
        size = base_size, color = "black"
      ),
      axis.text.y = element_text(
        size = base_size, color = "black"
      ),
      axis.title.x = element_text(
        size = base_size + 1, color = "black"
      ),
      axis.title.y = element_text(
        size = base_size + 1, color = "black"
      ),
      axis.line = element_line(linewidth = 0.7, color = "black"),
      axis.ticks = element_line(linewidth = 0.7, color = "black"),
      axis.ticks.length = unit(2, "mm"),
      plot.margin = margin(6, 8, 6, 8),
      legend.position = "none"
    )
}

# ---------------------------
# Build plot
# ---------------------------
p <- ggplot(df, aes(x = `Temperature interval`)) +
  geom_col(
    aes(y = Mean),
    fill = bar_col,
    color = "cyan4",
    linewidth = 0.3,
    width = 0.7
  ) +
  geom_errorbar(
    aes(ymin = pmax(0, Mean - SD), ymax = pmin(425, Mean + SD)),
    width = 0.15,
    linewidth = 0.7,
    color = "cyan4"
  ) +
  geom_line(
    aes(y = Maximum, group = 1),
    color = line_col,
    linewidth = 1,
    lineend = "round"
  ) +
  geom_point(
    aes(y = Maximum),
    color = line_col,
    size = 2.5
  ) +
  scale_y_continuous(
    limits = c(0, 425),
    breaks = c(0, 100, 200, 300, 400),
    labels = c("0", "100", "200", "300", "400"),
    expand = c(0, 0)
  ) +
  labs(
    x = "Temperature interval",
    y = "PC"
  ) +
  theme_nature_big(base_size = 11)
p
# ---------------------------
# Save PDF
# ---------------------------
ggsave(
  filename = "PCmonthly_7yr_mean_sd_max_5x3.pdf",
  plot = p,
  width = 5,
  height = 3,
  units = "in",
  device = cairo_pdf
)
cat("Saved: PCmonthly_7yr_mean_sd_max.pdf\n")
################################################################################
##### scatter plot of cyanobacteria (PC, cell counts, or mG) vs temperature
###-------------- 7-year PC vs T---------------------------------------
library(tidyverse)
library(mgcv)
rm(list = ls())

# load data
df <- read.csv("Taihu_Shi_7yr/PC_T_7year.csv")

# sort by temperature
df <- df %>%
  arrange(T)

# classify only increasing / decreasing
# first point is set to NA because there is no previous point for comparison
df <- df %>%
  mutate(
    dPC = c(NA, diff(PC)),
    phase = case_when(
      dPC > 0 ~ "Increasing",
      dPC < 0 ~ "Decreasing",
      TRUE ~ NA_character_
    )
  )

# ======================
# CORRELATION ANALYSIS (T vs PC)
# ======================
cat("\n\n=== Correlation between Temperature (T) and PC ===\n")

# Overall (full dataset)
pearson <- cor.test(df$T, df$PC, method = "pearson")
spearman <- cor.test(df$T, df$PC, method = "spearman")

cat("Method: Pearson\n")
cat("  Coefficient (r)  =", round(pearson$estimate, 4), "\n")
cat("  p-value          =", format.pval(pearson$p.value, digits = 4, eps = 0.0001), "\n\n")

cat("Method: Spearman\n")
cat("  Coefficient (rho) =", round(spearman$estimate, 4), "\n")
cat("  p-value           =", format.pval(spearman$p.value, digits = 4, eps = 0.0001), "\n\n")

# Per phase (for completeness — matches the color coding in the plot)
cat("--- Increasing phase only ---\n")
df_inc <- df %>% filter(phase == "Increasing")
if (nrow(df_inc) > 2) {
  pearson_inc <- cor.test(df_inc$T, df_inc$PC, method = "pearson")
  spearman_inc <- cor.test(df_inc$T, df_inc$PC, method = "spearman")
  cat("  Pearson r  =", round(pearson_inc$estimate, 4), 
      " p =", format.pval(pearson_inc$p.value, digits = 4, eps = 0.0001), "\n")
  cat("  Spearman rho =", round(spearman_inc$estimate, 4), 
      " p =", format.pval(spearman_inc$p.value, digits = 4, eps = 0.0001), "\n\n")
}

cat("--- Decreasing phase only ---\n")
df_dec <- df %>% filter(phase == "Decreasing")
if (nrow(df_dec) > 2) {
  pearson_dec <- cor.test(df_dec$T, df_dec$PC, method = "pearson")
  spearman_dec <- cor.test(df_dec$T, df_dec$PC, method = "spearman")
  cat("  Pearson r  =", round(pearson_dec$estimate, 4), 
      " p =", format.pval(pearson_dec$p.value, digits = 4, eps = 0.0001), "\n")
  cat("  Spearman rho =", round(spearman_dec$estimate, 4), 
      " p =", format.pval(spearman_dec$p.value, digits = 4, eps = 0.0001), "\n\n")
}

cat("Note: Pearson assumes linearity; Spearman is rank-based (monotonic). ",
    "The GAM smooth (non-linear) is shown in the plot.\n\n")

# GAM model
gam_model <- gam(PC ~ s(T, k = 5), data = df)

# prediction grid + 95% CI
T_seq <- seq(min(df$T), max(df$T), length.out = 200)
pred <- data.frame(T = T_seq)
pred_gam <- predict(gam_model, newdata = pred, se.fit = TRUE)
pred$PC_pred <- pred_gam$fit
pred$PC_lwr <- pred_gam$fit - 1.96 * pred_gam$se.fit
pred$PC_upr <- pred_gam$fit + 1.96 * pred_gam$se.fit

# -------------------------
# PLOT A ONLY
# -------------------------
pA <- ggplot(df, aes(x = T, y = PC)) +
  geom_point(aes(color = phase), size = 2.5, na.rm = TRUE) +
  geom_ribbon(
    data = pred,
    aes(x = T, ymin = PC_lwr, ymax = PC_upr),
    inherit.aes = FALSE,
    fill = "grey80",
    alpha = 0.6
  ) +
  geom_line(
    data = pred,
    aes(x = T, y = PC_pred),
    color = "black",
    linewidth = 1.2
  ) +
  geom_vline(
    xintercept = 30.35,
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.5
  ) +
  scale_color_manual(
    values = c(
      "Increasing" = "#D55E00",
      "Decreasing" = "#0072B2"
    ),
    na.translate = FALSE
  ) +
  labs(
    title = "A",
    x = "Temperature (°C)",
    y = "PC"
  ) +
  theme_classic(base_size = 14)
pA

ggsave(
  "PC_T_7year_regression_5x4.pdf",
  plot = pA,
  width = 5,
  height = 4
)

##@@@@@ removes the temperature above 30C
###-------------- 7-year PC vs T (filtered T ≤ 30) ---------------------------------------
library(tidyverse)
library(mgcv)
rm(list = ls())

# load data
df <- read.csv("Taihu_Shi_7yr/PC_T_7year.csv")

# <<< NEW: REMOVE ALL POINTS WITH T > 30 >>>
df <- df %>% filter(T <= 30)

# sort by temperature (only the remaining points)
df <- df %>%
  arrange(T)

# classify only increasing / decreasing
# first point is set to NA because there is no previous point for comparison
df <- df %>%
  mutate(
    dPC = c(NA, diff(PC)),
    phase = case_when(
      dPC > 0 ~ "Increasing",
      dPC < 0 ~ "Decreasing",
      TRUE ~ NA_character_
    )
  )

# ======================
# CORRELATION ANALYSIS (T vs PC) — on filtered data only
# ======================
cat("\n\n=== Correlation between Temperature (T) and PC (T ≤ 30 only) ===\n")

# Overall (filtered dataset)
pearson <- cor.test(df$T, df$PC, method = "pearson")
spearman <- cor.test(df$T, df$PC, method = "spearman")

cat("Method: Pearson\n")
cat("  Coefficient (r)  =", round(pearson$estimate, 4), "\n")
cat("  p-value          =", format.pval(pearson$p.value, digits = 4, eps = 0.0001), "\n\n")

cat("Method: Spearman\n")
cat("  Coefficient (rho) =", round(spearman$estimate, 4), "\n")
cat("  p-value           =", format.pval(spearman$p.value, digits = 4, eps = 0.0001), "\n\n")

# Per phase (for completeness — matches the color coding in the plot)
cat("--- Increasing phase only ---\n")
df_inc <- df %>% filter(phase == "Increasing")
if (nrow(df_inc) > 2) {
  pearson_inc <- cor.test(df_inc$T, df_inc$PC, method = "pearson")
  spearman_inc <- cor.test(df_inc$T, df_inc$PC, method = "spearman")
  cat("  Pearson r  =", round(pearson_inc$estimate, 4), 
      " p =", format.pval(pearson_inc$p.value, digits = 4, eps = 0.0001), "\n")
  cat("  Spearman rho =", round(spearman_inc$estimate, 4), 
      " p =", format.pval(spearman_inc$p.value, digits = 4, eps = 0.0001), "\n\n")
}

cat("--- Decreasing phase only ---\n")
df_dec <- df %>% filter(phase == "Decreasing")
if (nrow(df_dec) > 2) {
  pearson_dec <- cor.test(df_dec$T, df_dec$PC, method = "pearson")
  spearman_dec <- cor.test(df_dec$T, df_dec$PC, method = "spearman")
  cat("  Pearson r  =", round(pearson_dec$estimate, 4), 
      " p =", format.pval(pearson_dec$p.value, digits = 4, eps = 0.0001), "\n")
  cat("  Spearman rho =", round(spearman_dec$estimate, 4), 
      " p =", format.pval(spearman_dec$p.value, digits = 4, eps = 0.0001), "\n\n")
}

cat("Note: Pearson assumes linearity; Spearman is rank-based (monotonic). ",
    "The GAM smooth (non-linear) is shown in the plot.\n\n")

# GAM model (fitted on filtered data only)
gam_model <- gam(PC ~ s(T, k = 5), data = df)

# prediction grid + 95% CI (on filtered range)
T_seq <- seq(min(df$T), max(df$T), length.out = 200)
pred <- data.frame(T = T_seq)
pred_gam <- predict(gam_model, newdata = pred, se.fit = TRUE)
pred$PC_pred <- pred_gam$fit
pred$PC_lwr <- pred_gam$fit - 1.96 * pred_gam$se.fit
pred$PC_upr <- pred_gam$fit + 1.96 * pred_gam$se.fit

# -------------------------
# PLOT A ONLY (filtered data + new GAM)
# -------------------------
pA <- ggplot(df, aes(x = T, y = PC)) +
  geom_point(aes(color = phase), size = 2.5, na.rm = TRUE) +
  geom_ribbon(
    data = pred,
    aes(x = T, ymin = PC_lwr, ymax = PC_upr),
    inherit.aes = FALSE,
    fill = "grey80",
    alpha = 0.6
  ) +
  geom_line(
    data = pred,
    aes(x = T, y = PC_pred),
    color = "black",
    linewidth = 1.2
  ) +
  geom_vline(
    xintercept = 30.35,
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.5
  ) +
  scale_color_manual(
    values = c(
      "Increasing" = "#D55E00",
      "Decreasing" = "#0072B2"
    ),
    na.translate = FALSE
  ) +
  labs(
    title = "A",
    x = "Temperature (°C)",
    y = "PC"
  ) +
  theme_classic(base_size = 14)
pA

ggsave(
  "PC_T_7year_regression_noMaxT_5x4.pdf",
  plot = pA,
  width = 5,
  height = 4
)
#@@@@

##### ============================================================
#####-------1 year MA cell counts vs Temperature --------------------
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
rm(list = ls())
# ---------------------------
# Read data
# ---------------------------
ma_raw   <- read_excel("algae-Temp_taihu-1yr.xlsx", sheet = "MAcellcount")
temp_raw <- read_excel("algae-Temp_taihu-1yr.xlsx", sheet = "temperature")

# ---------------------------
# Reshape MAcellcount
# ---------------------------
ma_long <- ma_raw %>%
  pivot_longer(
    cols = starts_with("S"),
    names_to = "Site",
    values_to = "MAcellcount"
  ) %>%
  mutate(
    Month = as.integer(Month),
    Site = as.character(Site),
    MAcellcount = as.numeric(MAcellcount)
  ) %>%
  select(Date, Month, Site, MAcellcount)

# ---------------------------
# Temperature + Phase
# ---------------------------
temp_df <- temp_raw %>%
  mutate(
    Month = as.integer(Month),
    Site = as.character(Site),
    Temperature = as.numeric(Temperature),
    Phase = as.character(Phase)
  ) %>%
  select(Year, Month, Site, Temperature, Phase)

# ---------------------------
# Join + transform
# ---------------------------
df_all <- ma_long %>%
  inner_join(temp_df, by = c("Month", "Site")) %>%
  filter(!is.na(MAcellcount), !is.na(Temperature), !is.na(Phase)) %>%
  mutate(
    logMA = log2(MAcellcount + 1),
    Phase = factor(Phase, levels = c("Increasing", "Decreasing"))
  )

# ---------------------------
# Nature-like theme
# ---------------------------
theme_nature_big <- function(base_size = 11, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.text = element_text(size = base_size, color = "black"),
      axis.title = element_text(size = base_size + 1, color = "black"),
      axis.line = element_line(linewidth = 0.7, color = "black"),
      axis.ticks = element_line(linewidth = 0.7, color = "black"),
      axis.ticks.length = unit(2, "mm"),
      legend.position = c(0.18, 0.82),
      legend.title = element_blank(),
      legend.text = element_text(size = base_size),
      legend.background = element_blank(),
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0),
      plot.margin = margin(6, 8, 6, 8)
    )
}

# ---------------------------
# Plot function
# ---------------------------
make_linear_plot <- function(dat, out_file, plot_title = NULL) {
  
  model <- lm(logMA ~ Temperature, data = dat)
  sm <- summary(model)
  
  intercept <- coef(model)[1]
  slope     <- coef(model)[2]
  p_value   <- sm$coefficients[2, 4]
  r2        <- sm$r.squared
  
  # 95% confidence interval for slope
  slope_ci <- confint(model, level = 0.95)["Temperature", ]
  
  cat("====================================\n")
  cat(plot_title, "\n")
  cat("Intercept:", intercept, "\n")
  cat("Slope:", slope, "\n")
  cat("Slope 95% CI:", slope_ci[1], "to", slope_ci[2], "\n")
  cat("p-value:", p_value, "\n")
  cat("R²:", r2, "\n")
  
  eq_label <- paste0(
    "y = ", round(slope, 3), "x + ", round(intercept, 2),
    "\nR² = ", round(r2, 3),
    "\np = ", signif(p_value, 3),
    "\n95% CI = [", round(slope_ci[1], 3), ", ", round(slope_ci[2], 3), "]"
  )
  
  p <- ggplot(dat, aes(x = Temperature, y = logMA)) +
    geom_point(
      aes(color = Phase),
      size = 2.8,
      alpha = 0.95
    ) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      level = 0.95,
      color = "black",
      fill = "grey80",
      alpha = 0.6,
      linewidth = 1.2
    ) +
    geom_vline(
      xintercept = 30,
      linetype = "dashed",
      color = "grey40",
      linewidth = 0.6
    ) +
    scale_color_manual(values = c(
      "Increasing" = "#D55E00",
      "Decreasing" = "#0072B2"
    )) +
    annotate(
      "text",
      x = min(dat$Temperature) + 0.05 * diff(range(dat$Temperature)),
      y = max(dat$logMA),
      label = eq_label,
      hjust = 0,
      vjust = 1,
      size = 4
    ) +
    labs(
      title = plot_title,
      x = "Temperature (°C)",
      y = expression(log[2]*"(MA cell count + 1)")
    ) +
    theme_nature_big(base_size = 11)
  
  print(p)
  
  # ---------------------------
  # Save (5 x 5 inches)
  # ---------------------------
  ggsave(
    filename = out_file,
    plot = p,
    width = 5,
    height = 5,
    units = "in",
    device = cairo_pdf
  )
  
  cat("Saved:", out_file, "\n\n")
}

# ---------------------------
# Plot 1: all data
# ---------------------------
make_linear_plot(
  dat = df_all,
  out_file = "MA_log2_vs_T_linear_all.pdf",
  plot_title = "All months"
)

# ---------------------------
# Plot 2: without Month 1 and 12
# ---------------------------
df_no_1_12 <- df_all %>%
  filter(!Month %in% c(1, 12))

make_linear_plot(
  dat = df_no_1_12,
  out_file = "MA_log2_vs_T_linear_no_M1_M12.pdf",
  plot_title = "Without Months 1 and 12"
)

## @@@@@ >>
#####-------1 year MA cell counts vs Temperature (T ≤ 30 °C only) --------------------
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
rm(list = ls())

# ---------------------------
# Read data
# ---------------------------
ma_raw <- read_excel("algae-Temp_taihu-1yr.xlsx", sheet = "MAcellcount")
temp_raw <- read_excel("algae-Temp_taihu-1yr.xlsx", sheet = "temperature")

# ---------------------------
# Reshape MAcellcount
# ---------------------------
ma_long <- ma_raw %>%
  pivot_longer(
    cols = starts_with("S"),
    names_to = "Site",
    values_to = "MAcellcount"
  ) %>%
  mutate(
    Month = as.integer(Month),
    Site = as.character(Site),
    MAcellcount = as.numeric(MAcellcount)
  ) %>%
  select(Date, Month, Site, MAcellcount)

# ---------------------------
# Temperature + Phase
# ---------------------------
temp_df <- temp_raw %>%
  mutate(
    Month = as.integer(Month),
    Site = as.character(Site),
    Temperature = as.numeric(Temperature),
    Phase = as.character(Phase)
  ) %>%
  select(Year, Month, Site, Temperature, Phase)

# ---------------------------
# Join + transform
# ---------------------------
df_all <- ma_long %>%
  inner_join(temp_df, by = c("Month", "Site")) %>%
  filter(!is.na(MAcellcount), !is.na(Temperature), !is.na(Phase)) %>%
  mutate(
    logMA = log2(MAcellcount + 1),
    Phase = factor(Phase, levels = c("Increasing", "Decreasing"))
  )

# <<< NEW: REMOVE ALL POINTS WITH Temperature > 30 °C >>>
# (applied to BOTH plots and BOTH regressions)
df_all <- df_all %>% filter(Temperature <= 30)

# ---------------------------
# Nature-like theme
# ---------------------------
theme_nature_big <- function(base_size = 11, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.text = element_text(size = base_size, color = "black"),
      axis.title = element_text(size = base_size + 1, color = "black"),
      axis.line = element_line(linewidth = 0.7, color = "black"),
      axis.ticks = element_line(linewidth = 0.7, color = "black"),
      axis.ticks.length = unit(2, "mm"),
      legend.position = c(0.18, 0.82),
      legend.title = element_blank(),
      legend.text = element_text(size = base_size),
      legend.background = element_blank(),
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0),
      plot.margin = margin(6, 8, 6, 8)
    )
}

# ---------------------------
# Plot function (regression & plot on filtered data only)
# ---------------------------
make_linear_plot <- function(dat, out_file, plot_title = NULL) {
  
  model <- lm(logMA ~ Temperature, data = dat)
  sm <- summary(model)
  
  intercept <- coef(model)[1]
  slope <- coef(model)[2]
  p_value <- sm$coefficients[2, 4]
  r2 <- sm$r.squared
  
  # 95% confidence interval for slope
  slope_ci <- confint(model, level = 0.95)["Temperature", ]
  
  cat("====================================\n")
  cat(plot_title, " (T ≤ 30 °C only)\n")
  cat("Intercept:", round(intercept, 3), "\n")
  cat("Slope:", round(slope, 3), "\n")
  cat("Slope 95% CI:", round(slope_ci[1], 3), "to", round(slope_ci[2], 3), "\n")
  cat("p-value:", format.pval(p_value, digits = 4, eps = 0.0001), "\n")
  cat("R²:", round(r2, 3), "\n\n")
  
  eq_label <- paste0(
    "y = ", round(slope, 3), "x + ", round(intercept, 2),
    "\nR² = ", round(r2, 3),
    "\np = ", signif(p_value, 3),
    "\n95% CI = [", round(slope_ci[1], 3), ", ", round(slope_ci[2], 3), "]"
  )
  
  p <- ggplot(dat, aes(x = Temperature, y = logMA)) +
    geom_point(
      aes(color = Phase),
      size = 2.8,
      alpha = 0.95
    ) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      level = 0.95,
      color = "black",
      fill = "grey80",
      alpha = 0.6,
      linewidth = 1.2
    ) +
    geom_vline(
      xintercept = 30,
      linetype = "dashed",
      color = "grey40",
      linewidth = 0.6
    ) +
    scale_color_manual(values = c(
      "Increasing" = "#D55E00",
      "Decreasing" = "#0072B2"
    )) +
    annotate(
      "text",
      x = min(dat$Temperature) + 0.05 * diff(range(dat$Temperature)),
      y = max(dat$logMA),
      label = eq_label,
      hjust = 0,
      vjust = 1,
      size = 4
    ) +
    labs(
      title = plot_title,
      x = "Temperature (°C)",
      y = expression(log[2]*"(MA cell count + 1)")
    ) +
    theme_nature_big(base_size = 11)
  
  print(p)
  
  # ---------------------------
  # Save (5 x 5 inches)
  # ---------------------------
  ggsave(
    filename = out_file,
    plot = p,
    width = 5,
    height = 5,
    units = "in",
    device = cairo_pdf
  )
  
  cat("Saved:", out_file, "\n\n")
}

# ---------------------------
# Plot 1: all months (T ≤ 30 °C only)
# ---------------------------
make_linear_plot(
  dat = df_all,
  out_file = "MA_log2_vs_T_linear_all_Tle30.pdf",
  plot_title = "All months (T ≤ 30 °C)"
)

# ---------------------------
# Plot 2: without Month 1 and 12 (T ≤ 30 °C only)
# ---------------------------
df_no_1_12 <- df_all %>%
  filter(!Month %in% c(1, 12))
make_linear_plot(
  dat = df_no_1_12,
  out_file = "MA_log2_vs_T_linear_no_M1_M12_Tle30.pdf",
  plot_title = "Without Months 1 and 12 (T ≤ 30 °C)"
)
## @@@@@ <<


##### ============================================================
###--------- 1/2 year MA relative abundance vs Temperature --------------------
library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
rm(list = ls())
# ---------------------------
# Read data
# ---------------------------
df <- read_csv("water_chemistry.csv", show_col_types = FALSE)

# ---------------------------
# Clean column names for easier coding
# ---------------------------
df <- df %>%
  rename(
    Chla = `Chl a`
  )

# ---------------------------
# Nature-like theme
# Larger fonts/lines for later shrinking
# ---------------------------
theme_nature_big <- function(base_size = 11, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.text = element_text(size = base_size, color = "black"),
      axis.title = element_text(size = base_size + 1, color = "black"),
      axis.line = element_line(linewidth = 0.7, color = "black"),
      axis.ticks = element_line(linewidth = 0.7, color = "black"),
      axis.ticks.length = unit(2, "mm"),
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0),
      plot.margin = margin(6, 8, 6, 8),
      legend.position = "none"
    )
}

# ---------------------------
# Function to make one scatter plot
# ---------------------------
make_scatter_plot <- function(data, yvar, ylab, panel_label,
                              point_color = "#0072B2") {
  
  x <- data$T
  y <- data[[yvar]]
  
  # correlation
  ct <- cor.test(x, y, method = "pearson")
  r_val <- unname(ct$estimate)
  p_val <- ct$p.value
  
  # annotation text
  cor_label <- paste0(
    "r = ", round(r_val, 3),
    "\np = ", signif(p_val, 3)
  )
  
  # coordinates for annotation
  x_pos <- min(x, na.rm = TRUE) + 0.05 * diff(range(x, na.rm = TRUE))
  y_pos <- max(y, na.rm = TRUE) - 0.05 * diff(range(y, na.rm = TRUE))
  
  p <- ggplot(data, aes(x = T, y = .data[[yvar]])) +
    geom_point(
      color = point_color,
      size = 2.8,
      alpha = 0.95
    ) +
    geom_smooth(
      method = "lm",
      se = TRUE,
      color = "black",
      fill = "grey80",
      linewidth = 1.2
    ) +
    
    # -------- ADD THIS --------
  geom_vline(
    xintercept = 30,
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.6
  ) +
    # --------------------------
  
  annotate(
    "text",
    x = x_pos,
    y = y_pos,
    label = cor_label,
    hjust = 0,
    vjust = 1,
    size = 4
  ) +
    labs(
      title = panel_label,
      x = "Temperature (°C)",
      y = ylab
    ) +
    theme_nature_big(base_size = 11)
  
  return(p)
}
# ---------------------------
# Make plots
# ---------------------------
p1 <- make_scatter_plot(
  data = df,
  yvar = "Chla",
  ylab = "Chl a",
  panel_label = "A"
)
print(p1)
p2 <- make_scatter_plot(
  data = df,
  yvar = "PC",
  ylab = "PC",
  panel_label = "B"
)
print(p2)
p3 <- make_scatter_plot(
  data = df,
  yvar = "MAmg",
  ylab = "MAmg",
  panel_label = "C"
)
print(p3)
# ---------------------------
# Combined 3-panel figure
# ---------------------------
combined <- p1 + p2 + p3 + plot_layout(ncol = 3)

print(combined)

ggsave(
  filename = "mG_T_vs_Chla_PC_MAmg_combined.pdf",
  plot = combined,
  width = 12,
  height = 4.2,
  units = "in",
  device = cairo_pdf
)

# ---------------------------
# Save individual panels
# ---------------------------
ggsave(
  filename = "mG19_T_vs_Chla.pdf",
  plot = p1,
  width = 5,
  height = 5,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = "mG19_T_vs_PC.pdf",
  plot = p2,
  width = 5,
  height = 5,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = "mG19_T_vs_MAmG_relAbundance.pdf",
  plot = p3,
  width = 5,
  height = 5,
  units = "in",
  device = cairo_pdf
)

## @@ >> 
# 1/2 year MA relative abundance vs Temperature (T ≤ 30 °C only) -------
library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
rm(list = ls())

# ---------------------------
# Read data
# ---------------------------
df <- read_csv("water_chemistry.csv", show_col_types = FALSE)

# <<< NEW: REMOVE ALL POINTS WITH T > 30 °C >>>
df <- df %>% filter(T <= 30)

# ---------------------------
# Clean column names for easier coding
# ---------------------------
df <- df %>%
  rename(
    Chla = `Chl a`
  )

# ---------------------------
# Nature-like theme
# Larger fonts/lines for later shrinking
# ---------------------------
theme_nature_big <- function(base_size = 11, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.text = element_text(size = base_size, color = "black"),
      axis.title = element_text(size = base_size + 1, color = "black"),
      axis.line = element_line(linewidth = 0.7, color = "black"),
      axis.ticks = element_line(linewidth = 0.7, color = "black"),
      axis.ticks.length = unit(2, "mm"),
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0),
      plot.margin = margin(6, 8, 6, 8),
      legend.position = "none"
    )
}

# ---------------------------
# Function to make one scatter plot (now includes full lm regression output)
# ---------------------------
make_scatter_plot <- function(data, yvar, ylab, panel_label,
                              point_color = "#0072B2") {
  
  x <- data$T
  y <- data[[yvar]]
  
  # Linear regression (same as geom_smooth)
  model <- lm(as.formula(paste(yvar, "~ T")), data = data)
  sm <- summary(model)
  
  intercept <- coef(model)[1]
  slope     <- coef(model)[2]
  p_value   <- sm$coefficients[2, 4]
  r2        <- sm$r.squared
  
  # Pearson correlation (for consistency with previous version)
  ct <- cor.test(x, y, method = "pearson")
  r_val <- unname(ct$estimate)
  
  # 95% CI for slope
  slope_ci <- confint(model, level = 0.95)["T", ]
  
  # Console output (full regression results)
  cat("====================================\n")
  cat(panel_label, " — ", yvar, " vs T (T ≤ 30 °C only)\n")
  cat("Slope (b)      :", round(slope, 4), "\n")
  cat("Intercept (a)  :", round(intercept, 4), "\n")
  cat("Slope 95% CI   :", round(slope_ci[1], 4), "to", round(slope_ci[2], 4), "\n")
  cat("R²             :", round(r2, 4), "\n")
  cat("Pearson r      :", round(r_val, 4), "\n")
  cat("p-value        :", format.pval(p_value, digits = 4, eps = 0.0001), "\n\n")
  
  # Annotation text for plot (now includes slope + R² + p)
  eq_label <- paste0(
    "y = ", round(slope, 3), "x + ", round(intercept, 2),
    "\nR² = ", round(r2, 3),
    "\nr = ", round(r_val, 3),
    "\np = ", signif(p_value, 3)
  )
  
  # Coordinates for annotation (dynamic)
  x_pos <- min(x, na.rm = TRUE) + 0.05 * diff(range(x, na.rm = TRUE))
  y_pos <- max(y, na.rm = TRUE) - 0.05 * diff(range(y, na.rm = TRUE))
  
  p <- ggplot(data, aes(x = T, y = .data[[yvar]])) +
    geom_point(
      color = point_color,
      size = 2.8,
      alpha = 0.95
    ) +
    geom_smooth(
      method = "lm",
      se = TRUE,
      color = "black",
      fill = "grey80",
      linewidth = 1.2
    ) +
    geom_vline(
      xintercept = 30,
      linetype = "dashed",
      color = "grey40",
      linewidth = 0.6
    ) +
    annotate(
      "text",
      x = x_pos,
      y = y_pos,
      label = eq_label,
      hjust = 0,
      vjust = 1,
      size = 4
    ) +
    labs(
      title = panel_label,
      x = "Temperature (°C)",
      y = ylab
    ) +
    theme_nature_big(base_size = 11)
  
  return(p)
}

# ---------------------------
# Make plots (all on filtered data T ≤ 30 °C)
# ---------------------------
p1 <- make_scatter_plot(
  data = df,
  yvar = "Chla",
  ylab = "Chl a",
  panel_label = "A"
)
print(p1)

p2 <- make_scatter_plot(
  data = df,
  yvar = "PC",
  ylab = "PC",
  panel_label = "B"
)
print(p2)

p3 <- make_scatter_plot(
  data = df,
  yvar = "MAmg",
  ylab = "MAmg",
  panel_label = "C"
)
print(p3)

# ---------------------------
# Combined 3-panel figure
# ---------------------------
combined <- p1 + p2 + p3 + plot_layout(ncol = 3)
print(combined)
ggsave(
  filename = "mG_T_vs_Chla_PC_MAmg_combined_Tle30.pdf",
  plot = combined,
  width = 12,
  height = 4.2,
  units = "in",
  device = cairo_pdf
)

# ---------------------------
# Save individual panels (filtered version)
# ---------------------------
ggsave(
  filename = "mG19_T_vs_Chla_Tle30.pdf",
  plot = p1,
  width = 5,
  height = 5,
  units = "in",
  device = cairo_pdf
)
ggsave(
  filename = "mG19_T_vs_PC_Tle30.pdf",
  plot = p2,
  width = 5,
  height = 5,
  units = "in",
  device = cairo_pdf
)
ggsave(
  filename = "mG19_T_vs_MAmG_relAbundance_Tle30.pdf",
  plot = p3,
  width = 5,
  height = 5,
  units = "in",
  device = cairo_pdf
)

## @@ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
###--------- 1/2 year MA relative abundance vs Temperature --------------------
library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
rm(list = ls())

# ---------------------------
# Read data
# ---------------------------
df <- read_csv("water_chemistry.csv", show_col_types = FALSE)

# ---------------------------
# Hard-coded assignment (exact SampleIDs you specified):
# These 8 timepoints = "Increasing"; ALL other rows = "Decreasing"
# ---------------------------
increasing_samples <- c(
  "mG19_06.15.1",
  "mG19_06.15.2",
  "mG19_06.15.3",
  "mG19_07.23.1",
  "mG19_07.23.2",
  "mG19_07.23.3",
  "mG19_09.10.1",
  "mG19_09.10.2"
)

df <- df %>%
  mutate(
    Phase = ifelse(SampleID %in% increasing_samples, "Increasing", "Decreasing")
  )

# ---------------------------
# Clean column names for easier coding
# ---------------------------
df <- df %>%
  rename(
    Chla = `Chl a`
  )

# ---------------------------
# Nature-like theme (legend now enabled at bottom for the new phase colors)
# Larger fonts/lines for later shrinking
# ---------------------------
theme_nature_big <- function(base_size = 11, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.text = element_text(size = base_size, color = "black"),
      axis.title = element_text(size = base_size + 1, color = "black"),
      axis.line = element_line(linewidth = 0.7, color = "black"),
      axis.ticks = element_line(linewidth = 0.7, color = "black"),
      axis.ticks.length = unit(2, "mm"),
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0),
      plot.margin = margin(6, 8, 6, 8),
      legend.position = "bottom",          # changed from "none" to show phase legend
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 11)
    )
}

# ---------------------------
# Function to make one scatter plot (now colors points by Phase, single overall regression line)
# ---------------------------
make_scatter_plot <- function(data, yvar, ylab, panel_label) {
  
  x <- data$T
  y <- data[[yvar]]
  
  # correlation (overall, across both phases)
  ct <- cor.test(x, y, method = "pearson")
  r_val <- unname(ct$estimate)
  p_val <- ct$p.value
  
  # annotation text
  cor_label <- paste0(
    "r = ", round(r_val, 3),
    "\np = ", signif(p_val, 3)
  )
  
  # coordinates for annotation
  x_pos <- min(x, na.rm = TRUE) + 0.05 * diff(range(x, na.rm = TRUE))
  y_pos <- max(y, na.rm = TRUE) - 0.05 * diff(range(y, na.rm = TRUE))
  
  p <- ggplot(data, aes(x = T, y = .data[[yvar]])) +
    geom_point(
      aes(color = Phase),          # color by Increasing/Decreasing
      size = 2.8,
      alpha = 0.95
    ) +
    geom_smooth(
      method = "lm",
      se = TRUE,
      color = "black",
      fill = "grey80",
      linewidth = 1.2
    ) +
    geom_vline(
      xintercept = 30,
      linetype = "dashed",
      color = "grey40",
      linewidth = 0.6
    ) +
    annotate(
      "text",
      x = x_pos,
      y = y_pos,
      label = cor_label,
      hjust = 0,
      vjust = 1,
      size = 4
    ) +
    labs(
      title = panel_label,
      x = "Temperature (°C)",
      y = ylab
    ) +
    scale_color_manual(
      values = c(
        "Increasing" = "#D55E00",
        "Decreasing" = "#0072B2"
      ),
      na.translate = FALSE,
      name = "Phase"
    ) +
    theme_nature_big(base_size = 11)
  
  return(p)
}

# ---------------------------
# Make plots
# ---------------------------
p1 <- make_scatter_plot(
  data = df,
  yvar = "Chla",
  ylab = "Chl a",
  panel_label = "A"
)
print(p1)
p2 <- make_scatter_plot(
  data = df,
  yvar = "PC",
  ylab = "PC",
  panel_label = "B"
)
print(p2)
p3 <- make_scatter_plot(
  data = df,
  yvar = "MAmg",
  ylab = "MAmg",
  panel_label = "C"
)
print(p3)

# ---------------------------
# Combined 3-panel figure
# ---------------------------
combined <- p1 + p2 + p3 + plot_layout(ncol = 3)
print(combined)
ggsave(
  filename = "mG_T_vs_Chla_PC_MAmg_combined.pdf",
  plot = combined,
  width = 12,
  height = 4.2,
  units = "in",
  device = cairo_pdf
)

# ---------------------------
# Save individual panels
# ---------------------------
ggsave(
  filename = "mG19_T_vs_Chla.pdf",
  plot = p1,
  width = 5,
  height = 5,
  units = "in",
  device = cairo_pdf
)
ggsave(
  filename = "mG19_T_vs_PC.pdf",
  plot = p2,
  width = 5,
  height = 5,
  units = "in",
  device = cairo_pdf
)
ggsave(
  filename = "mG19_T_vs_MAmG_relAbundance.pdf",
  plot = p3,
  width = 5,
  height = 5,
  units = "in",
  device = cairo_pdf
)

## @@ >>
# 1/2 year MA relative abundance vs Temperature (T ≤ 30 °C only) -------
library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
rm(list = ls())

# ---------------------------
# Read data
# ---------------------------
df <- read_csv("water_chemistry.csv", show_col_types = FALSE)

# ---------------------------
# Hard-coded assignment (exact same list as above):
# These 8 timepoints = "Increasing"; ALL other rows = "Decreasing"
# (assignment done on the FULL dataset first)
# ---------------------------
increasing_samples <- c(
  "mG19_06.15.1",
  "mG19_06.15.2",
  "mG19_06.15.3",
  "mG19_07.23.1",
  "mG19_07.23.2",
  "mG19_07.23.3",
  "mG19_09.10.1",
  "mG19_09.10.2"
)

df <- df %>%
  mutate(
    Phase = ifelse(SampleID %in% increasing_samples, "Increasing", "Decreasing")
  )

# <<< REMOVE ALL POINTS WITH T > 30 °C >>>
df <- df %>% filter(T <= 30)

# ---------------------------
# Clean column names for easier coding
# ---------------------------
df <- df %>%
  rename(
    Chla = `Chl a`
  )

# ---------------------------
# Nature-like theme (legend now enabled at bottom for the new phase colors)
# Larger fonts/lines for later shrinking
# ---------------------------
theme_nature_big <- function(base_size = 11, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.text = element_text(size = base_size, color = "black"),
      axis.title = element_text(size = base_size + 1, color = "black"),
      axis.line = element_line(linewidth = 0.7, color = "black"),
      axis.ticks = element_line(linewidth = 0.7, color = "black"),
      axis.ticks.length = unit(2, "mm"),
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0),
      plot.margin = margin(6, 8, 6, 8),
      legend.position = "bottom",          # changed from "none" to show phase legend
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 11)
    )
}

# ---------------------------
# Function to make one scatter plot (now colors points by Phase, single overall regression line;
# also prints full regression statistics to console as before)
# ---------------------------
make_scatter_plot <- function(data, yvar, ylab, panel_label) {
  
  x <- data$T
  y <- data[[yvar]]
  
  # Linear regression (same as geom_smooth)
  model <- lm(as.formula(paste(yvar, "~ T")), data = data)
  sm <- summary(model)
  
  intercept <- coef(model)[1]
  slope <- coef(model)[2]
  p_value <- sm$coefficients[2, 4]
  r2 <- sm$r.squared
  
  # Pearson correlation (for consistency with previous version)
  ct <- cor.test(x, y, method = "pearson")
  r_val <- unname(ct$estimate)
  
  # 95% CI for slope
  slope_ci <- confint(model, level = 0.95)["T", ]
  
  # Console output (full regression results — overall across both phases)
  cat("====================================\n")
  cat(panel_label, " — ", yvar, " vs T (T ≤ 30 °C only)\n")
  cat("Slope (b) :", round(slope, 4), "\n")
  cat("Intercept (a) :", round(intercept, 4), "\n")
  cat("Slope 95% CI :", round(slope_ci[1], 4), "to", round(slope_ci[2], 4), "\n")
  cat("R² :", round(r2, 4), "\n")
  cat("Pearson r :", round(r_val, 4), "\n")
  cat("p-value :", format.pval(p_value, digits = 4, eps = 0.0001), "\n\n")
  
  # Annotation text for plot (now includes slope + R² + p)
  eq_label <- paste0(
    "y = ", round(slope, 3), "x + ", round(intercept, 2),
    "\nR² = ", round(r2, 3),
    "\nr = ", round(r_val, 3),
    "\np = ", signif(p_value, 3)
  )
  
  # Coordinates for annotation (dynamic)
  x_pos <- min(x, na.rm = TRUE) + 0.05 * diff(range(x, na.rm = TRUE))
  y_pos <- max(y, na.rm = TRUE) - 0.05 * diff(range(y, na.rm = TRUE))
  
  p <- ggplot(data, aes(x = T, y = .data[[yvar]])) +
    geom_point(
      aes(color = Phase),          # color by Increasing/Decreasing
      size = 2.8,
      alpha = 0.95
    ) +
    geom_smooth(
      method = "lm",
      se = TRUE,
      color = "black",
      fill = "grey80",
      linewidth = 1.2
    ) +
    geom_vline(
      xintercept = 30,
      linetype = "dashed",
      color = "grey40",
      linewidth = 0.6
    ) +
    annotate(
      "text",
      x = x_pos,
      y = y_pos,
      label = eq_label,
      hjust = 0,
      vjust = 1,
      size = 4
    ) +
    labs(
      title = panel_label,
      x = "Temperature (°C)",
      y = ylab
    ) +
    scale_color_manual(
      values = c(
        "Increasing" = "#D55E00",
        "Decreasing" = "#0072B2"
      ),
      na.translate = FALSE,
      name = "Phase"
    ) +
    theme_nature_big(base_size = 11)
  
  return(p)
}

# ---------------------------
# Make plots (all on filtered data T ≤ 30 °C)
# ---------------------------
p1 <- make_scatter_plot(
  data = df,
  yvar = "Chla",
  ylab = "Chl a",
  panel_label = "A"
)
print(p1)
p2 <- make_scatter_plot(
  data = df,
  yvar = "PC",
  ylab = "PC",
  panel_label = "B"
)
print(p2)
p3 <- make_scatter_plot(
  data = df,
  yvar = "MAmg",
  ylab = "MAmg",
  panel_label = "C"
)
print(p3)

# ---------------------------
# Combined 3-panel figure
# ---------------------------
combined <- p1 + p2 + p3 + plot_layout(ncol = 3)
print(combined)
ggsave(
  filename = "mG_T_vs_Chla_PC_MAmg_combined_Tle30.pdf",
  plot = combined,
  width = 12,
  height = 4.2,
  units = "in",
  device = cairo_pdf
)

# ---------------------------
# Save individual panels (filtered version)
# ---------------------------
ggsave(
  filename = "mG19_T_vs_Chla_Tle30.pdf",
  plot = p1,
  width = 5,
  height = 5,
  units = "in",
  device = cairo_pdf
)
ggsave(
  filename = "mG19_T_vs_PC_Tle30.pdf",
  plot = p2,
  width = 5,
  height = 5,
  units = "in",
  device = cairo_pdf
)
ggsave(
  filename = "mG19_T_vs_MAmG_relAbundance_Tle30.pdf",
  plot = p3,
  width = 5,
  height = 5,
  units = "in",
  device = cairo_pdf
)

