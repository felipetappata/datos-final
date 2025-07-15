# ========================================================================= 
# FIGURE 1: BIAS COMPARISON ACROSS SAMPLE SIZES - VISUALIZATION
# =========================================================================
# 
# Usage:
# Rscript fig1.R
# 
# Author: Felipe I. Tappata
# Date: July 2025
# =========================================================================

library(tidyverse)
library(here)
library(tikzDevice)

# Set up here package to locate ourselves in the project
here::i_am("code/fig1.R")

# Function to save plot with tikz for LaTeX rendering
save_tikz_plot <- function(plot, filename, width = 6, height = 4) {
  # Set options for math mode formatting
  options(tikzDefaultEngine = "pdftex")
  options(tikzSanitizeCharacters = c("%", "$", "}", "{", "^", "_", "#", "&", "~"))
  options(tikzReplacementCharacters = c("\\%", "\\$", "\\}", "\\{", "\\^{}", "\\_{}", "\\#", "\\&", "\\~{}"))

  tikz(filename,
    width = width, height = height,
    sanitize = FALSE,
    documentDeclaration = "\\documentclass[8pt]{article}",
    packages = c("\\usepackage{tikz}", "\\usepackage{amsmath}", "\\usepackage{amssymb}")
  )
  print(plot)
  dev.off()
}

# Create output directory if it doesn't exist
if (!dir.exists(here("code", "output", "fig1"))) {
  dir.create(here("code", "output", "fig1"), recursive = TRUE)
}

# Load all CSV files from output/partial directory
csv_files <- list.files(here("code", "output", "partial"),
  pattern = "*.csv",
  full.names = TRUE
)

# Filter out master files that might exist
csv_files <- csv_files[!grepl("master", csv_files)]

# Read and combine all CSV files
simulation_data <- map_dfr(csv_files, read_csv, show_col_types = FALSE)

# Check the structure
glimpse(simulation_data)

# Prepare data for Panel 1 (alpha = 0.25, Model A)
# Filter for Model A, rho = 0.25, and sample sizes needed for the plot
panel1_data <- simulation_data %>%
  filter(
    model == "A",
    rho == 0.25,
    N %in% c(200, 400, 600, 800, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000)
  ) %>%
  # Create long format for plotting
  pivot_longer(
    cols = c(AB_bias, SYS_bias),
    names_to = "estimator_type",
    values_to = "bias"
  ) %>%
  # Create the series labels as in the figure
  mutate(series = case_when(
    estimator_type == "AB_bias" & selection_type == "full_sample" ~ "AB (all)",
    estimator_type == "AB_bias" & selection_type == "endogenous" ~ "AB (select)",
    estimator_type == "SYS_bias" & selection_type == "full_sample" ~ "System (all)",
    estimator_type == "SYS_bias" & selection_type == "endogenous" ~ "System (select)"
  )) %>%
  filter(!is.na(series)) # Remove any rows that don't match our criteria

# Check the prepared data
print(panel1_data)

# Create Panel 1 plot (mean bias, alpha = 0.25)
panel1 <- panel1_data %>%
  ggplot(aes(x = N, y = bias, color = series, shape = series)) +
  geom_line(size = 0.8) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c(
    "AB (all)" = "#2166ac",
    "AB (select)" = "#5aae61",
    "System (all)" = "#d73027",
    "System (select)" = "#fc8d59"
  )) +
  scale_shape_manual(values = c(
    "AB (all)" = 16, # circle
    "AB (select)" = 17, # triangle
    "System (all)" = 18, # diamond
    "System (select)" = 15
  )) + # square
  labs(
    x = "$N$",
    y = "",
    color = "",
    shape = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    aspect.ratio = 0.6
  ) +
  scale_x_continuous(
    breaks = c(0, 1000, 2000, 3000, 4000, 5000),
    labels = c("$0$", "$1000$", "$2000$", "$3000$", "$4000$", "$5000$")
  ) +
  scale_y_continuous(labels = function(x) paste0("$", x, "$")) +
  coord_cartesian(xlim = c(0, 5000))

# Display the plot
print(panel1)

# Prepare data for Panel 2 (alpha = 0.50, Model A)
# Filter for Model A, rho = 0.50, and sample sizes needed for the plot
panel2_data <- simulation_data %>%
  filter(
    model == "A",
    rho == 0.50,
    N %in% c(200, 400, 600, 800, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000)
  ) %>%
  # Create long format for plotting
  pivot_longer(
    cols = c(AB_bias, SYS_bias),
    names_to = "estimator_type",
    values_to = "bias"
  ) %>%
  # Create the series labels as in the figure
  mutate(series = case_when(
    estimator_type == "AB_bias" & selection_type == "full_sample" ~ "AB (all)",
    estimator_type == "AB_bias" & selection_type == "endogenous" ~ "AB (select)",
    estimator_type == "SYS_bias" & selection_type == "full_sample" ~ "System (all)",
    estimator_type == "SYS_bias" & selection_type == "endogenous" ~ "System (select)"
  )) %>%
  filter(!is.na(series)) # Remove any rows that don't match our criteria

# Check the prepared data for panel 2
print("Panel 2 data:")
print(panel2_data)

# Create Panel 2 plot (mean bias, alpha = 0.50)
panel2 <- panel2_data %>%
  ggplot(aes(x = N, y = bias, color = series, shape = series)) +
  geom_line(size = 0.8) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c(
    "AB (all)" = "#2166ac",
    "AB (select)" = "#5aae61",
    "System (all)" = "#d73027",
    "System (select)" = "#fc8d59"
  )) +
  scale_shape_manual(values = c(
    "AB (all)" = 16, # circle
    "AB (select)" = 17, # triangle
    "System (all)" = 18, # diamond
    "System (select)" = 15
  )) + # square
  labs(
    x = "$N$",
    y = "",
    color = "",
    shape = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    aspect.ratio = 0.6
  ) +
  scale_x_continuous(
    breaks = c(0, 1000, 2000, 3000, 4000, 5000),
    labels = c("$0$", "$1000$", "$2000$", "$3000$", "$4000$", "$5000$")
  ) +
  scale_y_continuous(labels = function(x) paste0("$", x, "$")) +
  coord_cartesian(xlim = c(0, 5000))

# Display Panel 2
print(panel2)

# Prepare data for Panel 3 (alpha = 0.75, Model A)
# Filter for Model A, rho = 0.75, and sample sizes needed for the plot
panel3_data <- simulation_data %>%
  filter(
    model == "A",
    rho == 0.75,
    N %in% c(200, 400, 600, 800, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000)
  ) %>%
  # Create long format for plotting
  pivot_longer(
    cols = c(AB_bias, SYS_bias),
    names_to = "estimator_type",
    values_to = "bias"
  ) %>%
  # Create the series labels as in the figure
  mutate(series = case_when(
    estimator_type == "AB_bias" & selection_type == "full_sample" ~ "AB (all)",
    estimator_type == "AB_bias" & selection_type == "endogenous" ~ "AB (select)",
    estimator_type == "SYS_bias" & selection_type == "full_sample" ~ "System (all)",
    estimator_type == "SYS_bias" & selection_type == "endogenous" ~ "System (select)"
  )) %>%
  filter(!is.na(series)) # Remove any rows that don't match our criteria

# Check the prepared data for panel 3
print("Panel 3 data:")
print(panel3_data)

# Create Panel 3 plot (mean bias, alpha = 0.75)
panel3 <- panel3_data %>%
  ggplot(aes(x = N, y = bias, color = series, shape = series)) +
  geom_line(size = 0.8) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c(
    "AB (all)" = "#2166ac",
    "AB (select)" = "#5aae61",
    "System (all)" = "#d73027",
    "System (select)" = "#fc8d59"
  )) +
  scale_shape_manual(values = c(
    "AB (all)" = 16, # circle
    "AB (select)" = 17, # triangle
    "System (all)" = 18, # diamond
    "System (select)" = 15
  )) + # square
  labs(
    x = "$N$",
    y = "",
    color = "",
    shape = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    aspect.ratio = 0.6
  ) +
  scale_x_continuous(
    breaks = c(0, 1000, 2000, 3000, 4000, 5000),
    labels = c("$0$", "$1000$", "$2000$", "$3000$", "$4000$", "$5000$")
  ) +
  scale_y_continuous(labels = function(x) paste0("$", x, "$")) +
  coord_cartesian(xlim = c(0, 5000))

# Display Panel 3
print(panel3)

# Save plots as TikZ files
save_tikz_plot(panel1, here("code", "output", "fig1", "panel1_bias_rho025.tex"), width = 6, height = 4)
save_tikz_plot(panel2, here("code", "output", "fig1", "panel2_bias_rho050.tex"), width = 6, height = 4)
save_tikz_plot(panel3, here("code", "output", "fig1", "panel3_bias_rho075.tex"), width = 6, height = 4)
