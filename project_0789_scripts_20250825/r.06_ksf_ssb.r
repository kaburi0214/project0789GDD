rm(list = ls()); gc()

ORIGINAL_DIR <- ""
output <- file.path(ORIGINAL_DIR, "06_ksf_ssb")


if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}


setwd(ORIGINAL_DIR)


set.seed(42)


library(brms)
library(dplyr)
library(tidyr)
library(meta)
library(ggplot2)


inv_data <- read.csv("00_rawdata/IHME-GBD_2021_DATA-7010f9b0-1.csv")
if (file.exists(file.path(ORIGINAL_DIR, "02_ssbintk", "ssb_all_country_bm.rds"))) {
  model <- readRDS(file.path(ORIGINAL_DIR, "02_ssbintk", "ssb_all_country_bm.rds"))
}
ssb_cnty <- read.csv("00_rawdata/Country-level estimates/v15_cnty.csv")


predict_data_global <- expand.grid(
  age = factor(seq(22.5, 97.5, by = 5)),
  gender = factor(c("Male", "Female")),
  education = factor(c("0-6", "6-12", "≥12")),
  residence = factor(c("Rural", "Urban")),
  superregion = NA,    
  country = NA
)


posterior_preds_global <- posterior_epred(
  model,
  newdata = predict_data_global,
  draws = 1000,
  allow_new_levels = TRUE
)


posterior_long_global <- as.data.frame(posterior_preds_global) %>%
  mutate(draw_id = 1:nrow(.)) %>%
  pivot_longer(-draw_id, names_to = "row_id", values_to = "log_intake") %>%
  mutate(
    row_id = as.integer(gsub("V", "", row_id)),
    intake = exp(log_intake)  
  ) %>%
  left_join(predict_data_global %>% mutate(row_id = row_number()), by = "row_id")

bins <- c(0, 25, 50, 100, 150, 200, Inf)


posterior_long_global <- posterior_long_global %>%
  mutate(
    intake_level = cut(intake, breaks = bins, right = FALSE, include.lowest = TRUE,
                       labels = paste0("level_", 1:(length(bins) - 1))),
    group_id = paste(gender, age, education, residence, sep = "_")
  )


pi_draws_global <- posterior_long_global %>%
  group_by(age, gender, education, residence, draw_id) %>%
  mutate(weight = 1 / n()) %>%
  group_by(age, gender, education, residence, draw_id, intake_level) %>%
  summarise(pi = sum(weight), .groups = "drop") %>%
  pivot_wider(names_from = intake_level, values_from = pi, values_fill = 0)

pi_draws_global <- pi_draws_global %>%
  dplyr::select(age, gender, education, residence, draw_id, 
         level_1, level_2, level_3, level_4, level_5, level_6)


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
RR_matrix[, 1] <- 1.0        
RR_matrix[, 2] <- 1.2         
RR_matrix[, 3] <- 1.4        
RR_matrix[, 4] <- 1.6         
RR_matrix[, 5] <- 1.8      
RR_matrix[, 6] <- RR_fixed   


calculate_paf <- function(pi_vec, rr_vec) {
  numerator <- sum(pi_vec * (rr_vec - 1))
  denominator <- sum(pi_vec * rr_vec)
  
  if(denominator <= 0) return(0)
  
  paf <- numerator / denominator
  return(max(0, min(1, paf)))
}


posterior_long_global_scaled <- posterior_long_global %>%
  group_by(draw_id) %>%
  mutate(
    
    intake_scaled = pmax(pmin(intake, quantile(intake, 0.9)), quantile(intake, 0.1)),
    intake = 10 + (intake_scaled - min(intake_scaled)) / (max(intake_scaled) - min(intake_scaled)) * (250 - 10)
  ) %>%
  ungroup() %>%
  
  mutate(
    intake_level = cut(intake, breaks = bins, right = FALSE, include.lowest = TRUE,
                       labels = paste0("level_", 1:(length(bins) - 1)))
  )


pi_draws_by_gender_global <- posterior_long_global_scaled %>%
  group_by(gender, draw_id) %>%
  mutate(weight = 1 / n()) %>%
  group_by(gender, draw_id, intake_level) %>%
  summarise(pi = sum(weight), .groups = "drop") %>%
  pivot_wider(names_from = intake_level, values_from = pi, values_fill = 0)


pi_draws_by_age_global <- posterior_long_global_scaled %>%
  group_by(age, draw_id) %>%
  mutate(weight = 1 / n()) %>%
  group_by(age, draw_id, intake_level) %>%
  summarise(pi = sum(weight), .groups = "drop") %>%
  pivot_wider(names_from = intake_level, values_from = pi, values_fill = 0)


pi_draws_by_education_global <- posterior_long_global_scaled %>%
  group_by(education, draw_id) %>%
  mutate(weight = 1 / n()) %>%
  group_by(education, draw_id, intake_level) %>%
  summarise(pi = sum(weight), .groups = "drop") %>%
  pivot_wider(names_from = intake_level, values_from = pi, values_fill = 0)


pi_draws_by_residence_global <- posterior_long_global_scaled %>%
  group_by(residence, draw_id) %>%
  mutate(weight = 1 / n()) %>%
  group_by(residence, draw_id, intake_level) %>%
  summarise(pi = sum(weight), .groups = "drop") %>%
  pivot_wider(names_from = intake_level, values_from = pi, values_fill = 0)

calculate_group_paf <- function(pi_data, group_var) {
  groups <- unique(pi_data[[group_var]])
  results <- list()
  
  for(group in groups) {
    group_data <- pi_data[pi_data[[group_var]] == group, ]
    level_cols <- paste0("level_", 1:6)
    
    paf_draws <- numeric(1000)
    for(j in 1:1000) {
      pi_vec <- as.numeric(group_data[j, level_cols])
      rr_vec <- RR_matrix[j, ]
      paf_draws[j] <- calculate_paf(pi_vec, rr_vec)
    }
    
    results[[as.character(group)]] <- data.frame(
      group_type = group_var,
      group_name = group,
      mean_paf = mean(paf_draws)*100 ,
      lower_paf = quantile(paf_draws, 0.025)*100,
      upper_paf = quantile(paf_draws, 0.975)*100
    )
  }
  
  return(do.call(rbind, results))
}


paf_gender <- calculate_group_paf(pi_draws_by_gender_global, "gender")
paf_age <- calculate_group_paf(pi_draws_by_age_global, "age")
paf_education <- calculate_group_paf(pi_draws_by_education_global, "education")
paf_residence <- calculate_group_paf(pi_draws_by_residence_global, "residence")


all_paf_results <- rbind(paf_gender, paf_age, paf_education, paf_residence)

plot_data <- all_paf_results %>%
  mutate(
    category = case_when(
      group_type == "gender" ~ "Sex",
      group_type == "education" ~ "Education level", 
      group_type == "residence" ~ "Area of residence",
      group_type == "age" ~ "Age category (years)"
    ),
    group_label = case_when(
      group_name == "Male" ~ "Male",
      group_name == "Female" ~ "Female",
      group_name == "0-6" ~ "Low",
      group_name == "6-12" ~ "Medium", 
      group_name == "≥12" ~ "High",
      group_name == "Rural" ~ "Rural",
      group_name == "Urban" ~ "Urban",
      TRUE ~ as.character(group_name)
    )
  ) %>%

  mutate(
    group_label = ifelse(category == "Age category (years)", 
                         ifelse(as.numeric(as.character(group_name)) >= 97.5,
                                "95+",  
                                paste0(as.numeric(as.character(group_name)) - 2.5, "-", 
                                       as.numeric(as.character(group_name)) + 2.5 - 1)),
                         group_label),
  
    age_order = ifelse(category == "Age category (years)",
                       as.numeric(as.character(group_name)),
                       NA)
  )

plot_data$category <- factor(plot_data$category,
                             levels = c("Sex", "Education level", "Area of residence",
                                        "Age category (years)"))

paf_by_dm <- ggplot(plot_data, aes(x = mean_paf, 
                                   y = case_when(
                                     category == "Education level" ~ factor(group_label, levels = c("High", "Medium", "Low")),
                                     category == "Age category (years)" ~ factor(group_label, levels = c("20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80-84", "85-89", "90-94", "95+")),
                                     TRUE ~ reorder(group_label, mean_paf)
                                   ))) +
  geom_col(aes(fill = category), alpha = 0.8, width = 0.7) +
  geom_errorbarh(aes(xmin = lower_paf, xmax = upper_paf),
                 height = 0.3, alpha = 0.7) +
  facet_wrap(~category, scales = "free_y", ncol = 1) +
  labs(
    x = "Gout incidence attributable to SSB intake (%)",
    y = NULL
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(hjust = 0, face = "bold", size = 11),
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 10)
  ) +
  scale_fill_viridis_d(option = "plasma", begin = 0.2, end = 0.8)

print(paf_by_dm)


ages <- seq(22.5, 92.5, by = 5)
age_ranges <- paste0(ages - 2.5, "-", ages + 2.5 - 1)
age_mapping <- setNames(age_ranges, as.character(ages))
age_mapping["97.5"] <- "95+"

get_group_asir <- function(group_type, group_name) {
  if(group_name == "Male") {
   
    asir_row <- inv_data %>% 
      filter(sex == "Male", age == "Age-standardized", metric == "Rate") 

  } else if(group_name == "Female") {
  
    asir_row <- inv_data %>% 
      filter(sex == "Female", age == "Age-standardized", metric == "Rate") 
     
  } else if(group_type == "age") {
    
    gbd_age <- paste0(age_mapping[group_name], " years")
    asir_row <- inv_data %>% 
      filter(sex == "Both", age == gbd_age, metric == "Rate") 
      
  } else {
 
    asir_row <- inv_data %>% 
      filter(sex == "Both", age == "Age-standardized", metric == "Rate") %>%
      filter(row_number() == 1)
  }
  
 
  if(nrow(asir_row) > 0) {
    return(asir_row$val[1])
  }
}


all_burden_results <- all_paf_results %>%
  rowwise() %>%
  mutate(
    asir = get_group_asir(group_type, group_name),
  
    attributable_incidence_mean = (mean_paf/100) * asir ,
    attributable_incidence_lower = (lower_paf/100) * asir ,
    attributable_incidence_upper = (upper_paf/100) * asir 
  ) %>%
  

  ungroup()

all_burden_results_for_csv <- all_burden_results %>%
  mutate(
    mean_paf = mean_paf / 100,
    lower_paf = lower_paf / 100,
    upper_paf = upper_paf / 100
  )

plot_data1 <- all_burden_results %>%
  mutate(
    category = case_when(
      group_type == "gender" ~ "Sex",
      group_type == "education" ~ "Education level",
      group_type == "residence" ~ "Area of residence",
      group_type == "age" ~ "Age category (years)"
    ),
    group_label = case_when(
      group_name == "Male" ~ "Male",
      group_name == "Female" ~ "Female",
      group_name == "0-6" ~ "Low",
      group_name == "6-12" ~ "Medium",
      group_name == "≥12" ~ "High",
      group_name == "Rural" ~ "Rural",
      group_name == "Urban" ~ "Urban",
      TRUE ~ as.character(group_name)
    )
  ) %>%
  
  mutate(
    group_label = ifelse(category == "Age category (years)",
                         ifelse(as.numeric(as.character(group_name)) >= 97.5,
                                "95+",
                                paste0(as.numeric(as.character(group_name)) - 2.5, "-",
                                       as.numeric(as.character(group_name)) + 2.5 - 1)),
                         group_label)
  )


plot_data1 <- plot_data1 %>%
  mutate(
    group_label = case_when(
      category == "Education level" & group_label == "Low" ~ "Low",
      category == "Education level" & group_label == "Medium" ~ "Medium", 
      category == "Education level" & group_label == "High" ~ "High",
      category == "Sex" & group_label == "Female" ~ "Female",
      category == "Sex" & group_label == "Male" ~ "Male",
      category == "Area of residence" & group_label == "Rural" ~ "Rural",
      category == "Area of residence" & group_label == "Urban" ~ "Urban",
      TRUE ~ group_label
    ),
   
    group_label = case_when(
      category == "Education level" ~ factor(group_label, levels = c("High", "Medium", "Low")),
      category == "Age category (years)" ~ factor(group_label, levels = c("20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80-84", "85-89", "90-94", "95+")),
      category == "Sex" ~ factor(group_label, levels = c("Female", "Male")),
      category == "Area of residence" ~ factor(group_label, levels = c("Rural", "Urban")),
      TRUE ~ as.factor(group_label)
    )
  )
plot_data1$category <- factor(plot_data1$category,
                             levels = c("Sex", "Education level", "Area of residence",
                                        "Age category (years)"))

attribute_by_dm <- ggplot(plot_data1, aes(x = attributable_incidence_mean, y = group_label)) +
  geom_col(aes(fill = category), alpha = 0.8, width = 0.7) +
  geom_errorbarh(aes(xmin = attributable_incidence_lower, xmax = attributable_incidence_upper),
                 height = 0.3, alpha = 0.7) +
  facet_wrap(~category, scales = "free_y", ncol = 1) +
  labs(
    x = "Attributable gout incidence (cases per 100,000 adults)",
    y = NULL
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(hjust = 0, face = "bold", size = 11),
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 10)
  ) +
  scale_fill_viridis_d(option = "plasma", begin = 0.2, end = 0.8)

print(attribute_by_dm)


write.csv(all_burden_results_for_csv, file.path(output, "01.Attributable_burden_results.csv"), row.names = FALSE)

ggsave(file.path(output, "02.PAF_by_demographics.png"), 
       plot = paf_by_dm, 
       width = 12, height = 10, dpi = 300)

ggsave(file.path(output, "03.Attributable_incidence_by_demographics.png"), 
       plot = attribute_by_dm, 
       width = 12, height = 10, dpi = 300)

ggsave(file.path(output, "02.PAF_by_demographics.pdf"), 
       plot = paf_by_dm, 
       width = 12, height = 10)

ggsave(file.path(output, "03.Attributable_incidence_by_demographics.pdf"), 
       plot = attribute_by_dm, 
       width = 12, height = 10)

