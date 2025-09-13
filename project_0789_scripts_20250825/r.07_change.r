rm(list = ls()); gc()

ORIGINAL_DIR <- ""
output <- file.path(ORIGINAL_DIR, "07_change")

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

setwd(ORIGINAL_DIR)

set.seed(42)

library(ggplot2)
library(rworldmap)
library(meta)
library(brms)
library(dplyr)
library(tidyr)
library(countrycode)
library(readxl)

ssb_countries <- read.csv('00_rawdata/Country-level estimates/v15_cnty.csv')

ssb_countries_90 <- ssb_countries %>%
  filter(
    year == 1990,
    age!= 999,
    female != 999,
    urban != 999,
    edu != 999,
    !is.na(median),
    !is.na(lowerci_95),
    !is.na(upperci_95)
  )

ssb_countries_90 <- ssb_countries_90 %>%
  mutate(
    age = factor(age, levels = c(0.5, 1.5, 3.5, seq(7.5, 97.5, by = 5))),
    gender = factor(female, levels = c(0, 1), labels = c("Male", "Female")),      
    education = factor(edu, levels = c(1, 2, 3), labels = c("0-6", "6-12", "≥12")),
    residence = factor(urban, levels = c(0, 1), labels = c("Rural", "Urban")),
    superregion = factor(superregion2, levels = c("Asia", "FSU", "HIC", "LAC", "MENA", "SAARC", "SSA")),
    country = factor(iso3, levels = unique(ssb_countries_90$iso3))
  )

if (file.exists(file.path(output, "ssb_all_country_bm_1990.rds"))) {
  model_cnty_90 <- readRDS(file.path(output, "ssb_all_country_bm_1990.rds"))
} else {
  model_cnty_90 <- brm(
    formula = log(median) ~ age + gender + education + residence + (1 | superregion / country),
    data = ssb_countries_90,
    family = gaussian(),
    chains = 8,
    cores = 8,
    iter = 2000,
    warmup = 1000,
    control = list(
      adapt_delta = 0.95,      
      max_treedepth = 15       
    ),
    seed = 42
  )
  saveRDS(model_cnty_90, file.path(output, "ssb_all_country_bm_1990.rds"))
}

if (file.exists(file.path(ORIGINAL_DIR, "02_ssbintk", "ssb_all_country_bm.rds"))) {
  model_cnty <- readRDS(file.path(ORIGINAL_DIR, "02_ssbintk", "ssb_all_country_bm.rds"))
}

predict_data_global <- expand.grid(
  age = factor(seq(22.5, 97.5, by = 5)),
  gender = factor(c("Male", "Female")),
  education = factor(c("0-6", "6-12", "≥12")),
  residence = factor(c("Rural", "Urban")),
  superregion = NA,
  country = NA
)

posterior_preds_global_90 <- posterior_epred(
  model_cnty_90,
  newdata = predict_data_global,
  draws = 1000,
  allow_new_levels = TRUE,
  re_formula = NA
)
posterior_preds_global_21 <- posterior_epred(
  model_cnty,
  newdata = predict_data_global,
  draws = 1000,
  allow_new_levels = TRUE,
  re_formula = NA
)

posterior_long_global_90 <- as.data.frame(posterior_preds_global_90) %>%
  mutate(draw_id = 1:nrow(.)) %>%
  pivot_longer(-draw_id, names_to = "row_id", values_to = "log_intake") %>%
  mutate(
    row_id = as.integer(gsub("V", "", row_id)),
    intake = exp(log_intake)
  ) %>%
  left_join(predict_data_global %>% mutate(row_id = row_number()), by = "row_id")
posterior_long_global_21 <- as.data.frame(posterior_preds_global_21) %>%
  mutate(draw_id = 1:nrow(.)) %>%
  pivot_longer(-draw_id, names_to = "row_id", values_to = "log_intake") %>%
  mutate(
    row_id = as.integer(gsub("V", "", row_id)),
    intake = exp(log_intake)
  ) %>%
  left_join(predict_data_global %>% mutate(row_id = row_number()), by = "row_id")

bins <- c(0, 25, 50, 100, 150, 200, Inf)

posterior_long_global_90 <- posterior_long_global_90 %>%
  mutate(
    intake_level = cut(intake, breaks = bins, right = FALSE, include.lowest = TRUE,
                       labels = paste0("level_", 1:(length(bins) - 1)))
  )
posterior_long_global_21 <- posterior_long_global_21 %>%
  mutate(
    intake_level = cut(intake, breaks = bins, right = FALSE, include.lowest = TRUE,
                       labels = paste0("level_", 1:(length(bins) - 1)))
  )

pi_draws_global_90 <- posterior_long_global_90 %>%
  group_by(draw_id) %>%
  mutate(weight = 1 / n()) %>%
  group_by(draw_id, intake_level) %>%
  summarise(pi = sum(weight), .groups = "drop") %>%
  pivot_wider(names_from = intake_level, values_from = pi, values_fill = 0) %>% 
  dplyr::select(draw_id, level_1, level_2, level_3, level_4, level_5, level_6)
pi_draws_global_21 <- posterior_long_global_21 %>%
  group_by(draw_id) %>%
  mutate(weight = 1 / n()) %>%
  group_by(draw_id, intake_level) %>%
  summarise(pi = sum(weight), .groups = "drop") %>%
  pivot_wider(names_from = intake_level, values_from = pi, values_fill = 0) %>%
  dplyr::select(draw_id, level_1, level_2, level_3, level_4, level_5, level_6)

calculate_paf <- function(pi_vec, rr_vec) {
  numerator <- sum(pi_vec * (rr_vec - 1))
  denominator <- sum(pi_vec * rr_vec)

  if(denominator <= 0) return(0)

  paf <- numerator / denominator
  return(max(0, min(1, paf)))
}

studies <- data.frame(
  study = c("HPFS", "NHS"),
  rr = c(1.84, 2.39),
  lower_ci = c(1.08, 1.34),
  upper_ci = c(3.15, 4.26)
)
studies$log_rr <- log(studies$rr)
studies$se_log_rr <- (log(studies$upper_ci) - log(studies$lower_ci)) / (2 * 1.96)

meta_result <- metagen(
  TE = log_rr,
  seTE = se_log_rr,
  data = studies,
  studlab = study,
  common = TRUE,
  random = FALSE
)

RR_fixed <- exp(meta_result$TE.common)
RR_matrix <- matrix(1, nrow = 1000, ncol = 6)
RR_matrix[, 6] <- RR_fixed

paf_draws_90 <- numeric(1000)
for(i in 1:1000) {
  pi_vec <- as.numeric(pi_draws_global_90[i, -1])
  rr_vec <- RR_matrix[i, ]
  paf_draws_90[i] <- calculate_paf(pi_vec, rr_vec)
}
paf_draws_21 <- numeric(1000)
for(i in 1:1000) {
  pi_vec <- as.numeric(pi_draws_global_21[i, -1])
  rr_vec <- RR_matrix[i, ]
  paf_draws_21[i] <- calculate_paf(pi_vec, rr_vec)
}

paf_global_90 <- data.frame(
  mean_paf = mean(paf_draws_90) * 100,
  lower_paf = quantile(paf_draws_90, 0.025) * 100,
  upper_paf = quantile(paf_draws_90, 0.975) * 100
)
paf_global_21 <- data.frame(
  mean_paf = mean(paf_draws_21) * 100,
  lower_paf = quantile(paf_draws_21, 0.025) * 100,
  upper_paf = quantile(paf_draws_21, 0.975) * 100
)

print(paf_global_90)
print(paf_global_21)

gbd_gout_1990 <- read.csv("00_rawdata/IHME-GBD_2021_DATA-d92a7383-1.csv")

asir_1990 <- gbd_gout_1990 %>%
  filter(measure == "Incidence") %>%
  pull(val)

aspr_1990 <- gbd_gout_1990 %>%
  filter(measure == "Prevalence") %>%
  pull(val)

asdr_1990 <- gbd_gout_1990 %>%
  filter(measure == "DALYs (Disability-Adjusted Life Years)") %>%
  pull(val)

global_attributable_1990 <- data.frame(
  year = 1990,


  paf_mean = paf_global_90$mean_paf,
  paf_lower = paf_global_90$lower_paf,
  paf_upper = paf_global_90$upper_paf,

  asir_attributable = (paf_global_90$mean_paf/100) * asir_1990,
  asir_attributable_lower = (paf_global_90$lower_paf/100) * asir_1990,
  asir_attributable_upper = (paf_global_90$upper_paf/100) * asir_1990,

  aspr_attributable = (paf_global_90$mean_paf/100) * aspr_1990,
  aspr_attributable_lower = (paf_global_90$lower_paf/100) * aspr_1990,
  aspr_attributable_upper = (paf_global_90$upper_paf/100) * aspr_1990,

  asdr_attributable = (paf_global_90$mean_paf/100) * asdr_1990,
  asdr_attributable_lower = (paf_global_90$lower_paf/100) * asdr_1990,
  asdr_attributable_upper = (paf_global_90$upper_paf/100) * asdr_1990
)

gbd_gout_2021 <- read.csv("00_rawdata/IHME-GBD_2021_DATA-ec7daba0-1.csv")

asir_2021 <- gbd_gout_2021 %>%
  filter(
    location_name == "Global",
    sex_name == "Both", 
    age_name == "Age-standardized",
    metric_name == "Rate",
    measure_name == "Incidence",
    year == "2021"
  ) %>%
  pull(val)

aspr_2021 <- gbd_gout_2021 %>%
  filter(
    location_name == "Global",
    sex_name == "Both", 
    age_name == "Age-standardized",
    metric_name == "Rate",
    measure_name == "Prevalence",
    year == "2021"
  ) %>%
  pull(val)

asdr_2021 <- gbd_gout_2021 %>%
  filter(
    location_name == "Global",
    sex_name == "Both", 
    age_name == "Age-standardized",
    metric_name == "Rate",
    measure_name == "DALYs (Disability-Adjusted Life Years)",
    year == "2021"
  ) %>%
  pull(val)

global_attributable_2021 <- data.frame(
  year = 2021,

  paf_mean = paf_global_21$mean_paf,
  paf_lower = paf_global_21$lower_paf,
  paf_upper = paf_global_21$upper_paf,

  asir_attributable = (paf_global_21$mean_paf/100) * asir_2021,
  asir_attributable_lower = (paf_global_21$lower_paf/100) * asir_2021,
  asir_attributable_upper = (paf_global_21$upper_paf/100) * asir_2021,

  aspr_attributable = (paf_global_21$mean_paf/100) * aspr_2021,
  aspr_attributable_lower = (paf_global_21$lower_paf/100) * aspr_2021,
  aspr_attributable_upper = (paf_global_21$upper_paf/100) * aspr_2021,

  asdr_attributable = (paf_global_21$mean_paf/100) * asdr_2021,
  asdr_attributable_lower = (paf_global_21$lower_paf/100) * asdr_2021,
  asdr_attributable_upper = (paf_global_21$upper_paf/100) * asdr_2021
)

comparison_data <- data.frame(
  year = c(1990, 2021),
  asir_attributable = c(global_attributable_1990$asir_attributable,
                        global_attributable_2021$asir_attributable),
  aspr_attributable = c(global_attributable_1990$aspr_attributable,
                        global_attributable_2021$aspr_attributable),
  asdr_attributable = c(global_attributable_1990$asdr_attributable,
                        global_attributable_2021$asdr_attributable)
)

change_rates <- data.frame(
  indicator = c("ASIR", "ASPR", "ASDR"),
  value_1990 = c(global_attributable_1990$asir_attributable,
                 global_attributable_1990$aspr_attributable,
                 global_attributable_1990$asdr_attributable),
  value_2021 = c(global_attributable_2021$asir_attributable,
                 global_attributable_2021$aspr_attributable,
                 global_attributable_2021$asdr_attributable)
) %>%
  mutate(
    absolute_change = value_2021 - value_1990,
    percent_change = ((value_2021 - value_1990) / value_1990) * 100
  )

change_plot_data <- change_rates %>%
  mutate(
    color = ifelse(percent_change > 0, "Increase", "Decrease"),
    indicator = factor(indicator, levels = c("ASIR", "ASPR", "ASDR"))
  )

p_change <- ggplot(change_plot_data, aes(x = indicator, y = percent_change, fill = color)) +
  geom_bar(stat = "identity", width = 0.6, alpha = 0.8) +
  scale_fill_manual(values = c("Increase" = "#E31A1C", "Decrease" = "#1F78B4"),
                    name = "Change Direction") +
  labs(
    title = "Percentage Change in SSBs-Attributable Gout Burden (1990-2021)",
    x = "Indicators",
    y = "Percentage Change (%)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  geom_text(aes(label = sprintf("%+.1f%%", percent_change)),
            vjust = ifelse(change_plot_data$percent_change > 0, -0.3, 1.3),
            size = 4, fontface = "bold") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15)))

print(p_change)

print(change_rates)

final_results <- data.frame(
  year = c(1990, 2021),
  paf_mean = c(paf_global_90$mean_paf, paf_global_21$mean_paf),
  paf_lower = c(paf_global_90$lower_paf, paf_global_21$lower_paf),
  paf_upper = c(paf_global_90$upper_paf, paf_global_21$upper_paf),
  asir_attributable = c(global_attributable_1990$asir_attributable, global_attributable_2021$asir_attributable),
  aspr_attributable = c(global_attributable_1990$aspr_attributable, global_attributable_2021$aspr_attributable),
  asdr_attributable = c(global_attributable_1990$asdr_attributable, global_attributable_2021$asdr_attributable)
)

final_results$asir_change_percent <- c(NA, change_rates$percent_change[1])
final_results$aspr_change_percent <- c(NA, change_rates$percent_change[2])
final_results$asdr_change_percent <- c(NA, change_rates$percent_change[3])

write.csv(final_results, file.path(output, "01.SSB_gout_burden_change.csv"), row.names = FALSE)

ggsave(filename = file.path(output, "02.SSB_gout_burden_change.png"), 
       plot = p_change, 
       width = 10, 
       height = 8, 
       dpi = 300, 
       bg = "white")
ggsave(filename = file.path(output, "02.SSB_gout_burden_change.pdf"), 
       plot = p_change, 
       width = 10, 
       height = 8, 
       bg = "white")