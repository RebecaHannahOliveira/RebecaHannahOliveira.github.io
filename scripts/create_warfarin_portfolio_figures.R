# ============================================================
# Warfarin PopPK portfolio figures
# ============================================================

required_packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "patchwork",
  "scales"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Install the following packages first: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

final_model_dir <- paste0(
  "C:/Users/bequi/OneDrive/Documentos/PopPK_Portfolio/",
  "01_warfarin_training/monolix/",
  "warfarin_14_1cpt_fixed_allometry_ageCl"
)

output_dir <- "assets/warfarin"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Read Monolix results
# ------------------------------------------------------------

read_monolix_csv <- function(path) {
  read.csv(
    path,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

predictions <- read_monolix_csv(
  file.path(final_model_dir, "predictions.txt")
)

random_effects <- read_monolix_csv(
  file.path(
    final_model_dir,
    "IndividualParameters",
    "estimatedRandomEffects.txt"
  )
)

individual_parameters <- read_monolix_csv(
  file.path(
    final_model_dir,
    "IndividualParameters",
    "estimatedIndividualParameters.txt"
  )
)

bootstrap_summary <- read_monolix_csv(
  file.path(
    final_model_dir,
    "Bootstrap",
    "populationSummary.txt"
  )
)

# ------------------------------------------------------------
# Portfolio theme
# ------------------------------------------------------------

rose <- "#A64D70"
rose_dark <- "#71354D"
sage <- "#56806A"
sage_light <- "#E5F0E9"
cream <- "#FCFBFD"
ink <- "#332E32"
muted <- "#6E666C"
light_grid <- "#E5DCE2"

theme_portfolio <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(
        family = "sans",
        color = ink
      ),
      plot.background = element_rect(
        fill = cream,
        color = NA
      ),
      panel.background = element_rect(
        fill = cream,
        color = NA
      ),
      plot.title = element_text(
        color = rose_dark,
        face = "bold",
        size = rel(1.35),
        margin = margin(b = 7)
      ),
      plot.subtitle = element_text(
        color = muted,
        size = rel(0.95),
        margin = margin(b = 12)
      ),
      axis.title = element_text(
        color = ink,
        face = "bold"
      ),
      axis.text = element_text(
        color = muted
      ),
      panel.grid.major = element_line(
        color = light_grid,
        linewidth = 0.35
      ),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(
        fill = sage_light,
        color = NA
      ),
      strip.text = element_text(
        color = ink,
        face = "bold"
      ),
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.caption = element_text(
        color = muted,
        hjust = 0,
        size = rel(0.78),
        margin = margin(t = 10)
      ),
      plot.margin = margin(15, 18, 15, 15)
    )
}

# ============================================================
# Figure 1: concentration-time summary
# ============================================================

time_summary <- predictions %>%
  group_by(time) %>%
  summarise(
    observed_median = median(dv, na.rm = TRUE),
    observed_q25 = quantile(dv, 0.25, na.rm = TRUE),
    observed_q75 = quantile(dv, 0.75, na.rm = TRUE),
    population_median = median(popPred, na.rm = TRUE),
    individual_median = median(indivPred_SAEM, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

p_time <- ggplot(time_summary, aes(x = time)) +
  geom_ribbon(
    aes(
      ymin = observed_q25,
      ymax = observed_q75
    ),
    fill = rose,
    alpha = 0.16
  ) +
  geom_point(
    aes(y = observed_median),
    color = rose_dark,
    size = 2.5
  ) +
  geom_line(
    aes(
      y = observed_median,
      color = "Observed median"
    ),
    linewidth = 0.9
  ) +
  geom_line(
    aes(
      y = population_median,
      color = "Population prediction"
    ),
    linewidth = 1.05
  ) +
  geom_line(
    aes(
      y = individual_median,
      color = "Individual prediction"
    ),
    linewidth = 1.05
  ) +
  scale_color_manual(
    values = c(
      "Observed median" = rose_dark,
      "Population prediction" = sage,
      "Individual prediction" = "#C98C6B"
    )
  ) +
  scale_x_continuous(
    breaks = c(0, 12, 24, 48, 72, 96, 120)
  ) +
  labs(
    title = "Observed and predicted warfarin concentrations",
    subtitle = paste(
      "Median profiles across 32 participants;",
      "shaded region shows the observed interquartile range"
    ),
    x = "Time after dose (h)",
    y = "Warfarin concentration",
    caption = paste(
      "Educational analysis using the public Monolix warfarin dataset.",
      "This summary is not a simulation-based VPC."
    )
  ) +
  theme_portfolio()

ggsave(
  filename = file.path(
    output_dir,
    "warfarin_concentration_time.png"
  ),
  plot = p_time,
  width = 8.5,
  height = 5.3,
  dpi = 320,
  bg = cream
)

# ============================================================
# Figure 2: observed versus predicted
# ============================================================

max_prediction <- max(
  predictions$dv,
  predictions$popPred,
  predictions$indivPred_SAEM,
  na.rm = TRUE
)

p_pop <- ggplot(
  predictions,
  aes(x = popPred, y = dv)
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    color = muted,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  geom_point(
    color = rose,
    alpha = 0.72,
    size = 2
  ) +
  coord_equal(
    xlim = c(0, max_prediction),
    ylim = c(0, max_prediction)
  ) +
  labs(
    title = "Population predictions",
    x = "Population prediction",
    y = "Observed concentration"
  ) +
  theme_portfolio()

p_ind <- ggplot(
  predictions,
  aes(x = indivPred_SAEM, y = dv)
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    color = muted,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  geom_point(
    color = sage,
    alpha = 0.72,
    size = 2
  ) +
  coord_equal(
    xlim = c(0, max_prediction),
    ylim = c(0, max_prediction)
  ) +
  labs(
    title = "Individual predictions",
    x = "Individual prediction",
    y = "Observed concentration"
  ) +
  theme_portfolio()

p_gof <- (
  p_pop + p_ind +
    plot_annotation(
      title = "Goodness-of-fit assessment",
      subtitle = paste(
        "Observed concentrations compared with",
        "population and individual predictions"
      ),
      caption = "Dashed lines represent perfect agreement."
    )
) &
  theme(
    plot.title = element_text(
      color = rose_dark,
      face = "bold"
    )
  )

ggsave(
  filename = file.path(
    output_dir,
    "warfarin_observed_vs_predicted.png"
  ),
  plot = p_gof,
  width = 10,
  height = 5,
  dpi = 320,
  bg = cream
)

# ============================================================
# Figure 3: residual diagnostics
# ============================================================

p_res_time <- ggplot(
  predictions,
  aes(x = time, y = indWRes_SAEM)
) +
  geom_hline(
    yintercept = 0,
    color = muted,
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = c(-2, 2),
    color = rose,
    linetype = "dotted",
    alpha = 0.8
  ) +
  geom_point(
    color = sage,
    alpha = 0.72,
    size = 1.9
  ) +
  geom_smooth(
    method = "loess",
    formula = y ~ x,
    se = FALSE,
    color = rose_dark,
    linewidth = 0.9
  ) +
  labs(
    title = "IWRES versus time",
    x = "Time after dose (h)",
    y = "Individual weighted residual"
  ) +
  theme_portfolio()

p_res_pred <- ggplot(
  predictions,
  aes(x = indivPred_SAEM, y = indWRes_SAEM)
) +
  geom_hline(
    yintercept = 0,
    color = muted,
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = c(-2, 2),
    color = rose,
    linetype = "dotted",
    alpha = 0.8
  ) +
  geom_point(
    color = sage,
    alpha = 0.72,
    size = 1.9
  ) +
  geom_smooth(
    method = "loess",
    formula = y ~ x,
    se = FALSE,
    color = rose_dark,
    linewidth = 0.9
  ) +
  labs(
    title = "IWRES versus individual prediction",
    x = "Individual prediction",
    y = "Individual weighted residual"
  ) +
  theme_portfolio()

p_residuals <- (
  p_res_time + p_res_pred +
    plot_annotation(
      title = "Residual diagnostics",
      subtitle = paste(
        "Residuals remained generally centered near zero",
        "with mild curvature and several influential observations"
      )
    )
) &
  theme(
    plot.title = element_text(
      color = rose_dark,
      face = "bold"
    )
  )

ggsave(
  filename = file.path(
    output_dir,
    "warfarin_residual_diagnostics.png"
  ),
  plot = p_residuals,
  width = 10,
  height = 5,
  dpi = 320,
  bg = cream
)

# ============================================================
# Figure 4: bootstrap uncertainty relative to reference
# ============================================================

parameter_labels <- c(
  "Tlag_pop" = "Tlag",
  "ka_pop" = "ka",
  "V_pop" = "V",
  "Cl_pop" = "Cl",
  "beta_Cl_age10_c40" = "Age effect on Cl",
  "omega_ka" = "ω ka",
  "omega_V" = "ω V",
  "omega_Cl" = "ω Cl",
  "a" = "Additive error",
  "b" = "Proportional error"
)

parameter_order <- c(
  "Tlag",
  "ka",
  "V",
  "Cl",
  "Age effect on Cl",
  "ω ka",
  "ω V",
  "ω Cl",
  "Additive error",
  "Proportional error"
)

bootstrap_plot_data <- bootstrap_summary %>%
  mutate(
    label = unname(
      parameter_labels[as.character(Parameters)]
    ),
    label = ifelse(
      is.na(label),
      as.character(Parameters),
      label
    ),
    relative_median = Median / reference,
    relative_lower = P2.5 / reference,
    relative_upper = P97.5 / reference,
    parameter_group = case_when(
      grepl("^omega", Parameters) ~
        "Interindividual variability",
      Parameters %in% c("a", "b") ~
        "Residual error",
      TRUE ~
        "Fixed effects"
    ),
    label = factor(
      label,
      levels = rev(parameter_order)
    )
  )

p_bootstrap <- ggplot(
  bootstrap_plot_data,
  aes(
    y = label,
    x = relative_median,
    xmin = relative_lower,
    xmax = relative_upper
  )
) +
  geom_vline(
    xintercept = 1,
    color = rose_dark,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  geom_errorbarh(
    height = 0.18,
    color = sage,
    linewidth = 0.95
  ) +
  geom_point(
    color = sage,
    fill = cream,
    shape = 21,
    size = 3
  ) +
  facet_grid(
    parameter_group ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_x_continuous(
    breaks = seq(0, 4, by = 0.5),
    labels = label_number(accuracy = 0.1),
    expand = expansion(mult = c(0.03, 0.06))
  ) +
  labs(
    title = "Bootstrap parameter uncertainty",
    subtitle = paste(
      "Bootstrap medians and 95% intervals",
      "relative to the original estimates"
    ),
    x = "Ratio to reference estimate",
    y = NULL,
    caption = paste(
      "The dashed line indicates the original model estimate.",
      "Values below or above 1 indicate lower or higher bootstrap estimates."
    )
  ) +
  theme_portfolio() +
  theme(
    panel.grid.major.y = element_blank(),
    strip.text.y = element_text(
      size = 10
    )
  )

ggsave(
  filename = file.path(
    output_dir,
    "warfarin_bootstrap_uncertainty.png"
  ),
  plot = p_bootstrap,
  width = 8.5,
  height = 7.2,
  dpi = 320,
  bg = cream
)
# ------------------------------------------------------------
# Completion message
# ------------------------------------------------------------

created_files <- list.files(
  output_dir,
  pattern = "\\.png$",
  full.names = TRUE
)

cat(
  "\nCreated portfolio figures:\n",
  paste(created_files, collapse = "\n"),
  "\n"
)