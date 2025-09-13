rm(list = ls()); gc()

ORIGINAL_DIR <- ""
output <- file.path(ORIGINAL_DIR, "01_burdenanalysis")

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

setwd(ORIGINAL_DIR)

library(dplyr)       
library(MASS)
library(knitr)
library(kableExtra)
library(tidyr)

set.seed(42)

gbd_data <- read.csv('00_rawdata/IHME-GBD_2021_DATA-ec7daba0-1.csv', header = T)
extract_gbd_data <- function(data, target_age, m_name, target_year, value_prefix) {
  
  result <- data %>%
    filter(age_name == target_age, metric_name == m_name, year == target_year) %>%
    dplyr::select(
      location_name, 
      sex_name, 
      measure_name, 
      !!value_prefix := val,                    
      !!paste0(value_prefix, "_lower") := lower,  
      !!paste0(value_prefix, "_upper") := upper   
    )
  
  return(result)
}

num_2021_data <- extract_gbd_data(gbd_data, "All ages", "Number", 2021, "n_2021") %>%
  distinct(location_name, sex_name, measure_name, .keep_all = TRUE) 
asr_2021_data <- extract_gbd_data(gbd_data, "Age-standardized", "Rate", 2021, "asr_2021") %>%
  distinct(location_name, sex_name, measure_name, .keep_all = TRUE)
asr_1990_data <- extract_gbd_data(gbd_data, "Age-standardized", "Rate", 1990, "asr_1990") %>%
  distinct(location_name, sex_name, measure_name, .keep_all = TRUE)

percentage_change <- function(data_2021, data_1990, 
                              correlation = 0.9, n_sim = 1000){ 
  
  merged_asr <- merge(data_2021, data_1990, 
                      by = c("location_name", "sex_name", "measure_name"), 
                      suffixes = c("_2021", "_1990"))

  merged_asr$se_2021 <- (merged_asr$asr_2021_upper - merged_asr$asr_2021_lower) / (2 * 1.96)
  merged_asr$se_1990 <- (merged_asr$asr_1990_upper - merged_asr$asr_1990_lower) / (2 * 1.96)
  
  merged_asr$mean_change <- numeric(nrow(merged_asr))    
  merged_asr$lower_change <- numeric(nrow(merged_asr))   
  merged_asr$upper_change <- numeric(nrow(merged_asr)) 
  
  for(i in 1:nrow(merged_asr)) {
    
    mean_1990 <- merged_asr$asr_1990[i]
    mean_2021 <- merged_asr$asr_2021[i]
    se_1990 <- merged_asr$se_1990[i]
    se_2021 <- merged_asr$se_2021[i]
 
    var_1990 <- se_1990^2
    var_2021 <- se_2021^2
    cov_12 <- correlation * se_1990 * se_2021
    sigma_matrix <- matrix(c(var_1990, cov_12, cov_12, var_2021), nrow = 2)
    mu_vector <- c(mean_1990, mean_2021)

    samples <- mvrnorm(n_sim, mu_vector, sigma_matrix)
    samples[samples <= 0] <- 0.01

    perc_change <- (samples[,2] - samples[,1]) / samples[,1] * 100
    merged_asr$mean_change[i] <- round(mean(perc_change), 1)
    merged_asr$lower_change[i] <- round(quantile(perc_change, 0.025), 1)
    merged_asr$upper_change[i] <- round(quantile(perc_change, 0.975),1)
  }
  return(merged_asr)
}

percentage_change_result <- percentage_change(asr_2021_data, asr_1990_data)

format_ui <- function(data, val_col, lower_col, upper_col, digit = 1, big_mark = ",") {

  val <- round(data[[val_col]], digit)
  lower <- round(data[[lower_col]], digit)
  upper <- round(data[[upper_col]], digit)
  val_formatted <- trimws(format(val, big.mark = big_mark, scientific = FALSE))
  lower_formatted <- trimws(format(lower, big.mark = big_mark, scientific = FALSE))
  upper_formatted <- trimws(format(upper, big.mark = big_mark, scientific = FALSE))

  data$formatted_result <- paste0(val_formatted,"<br>","(",lower_formatted,"-",upper_formatted,")")
  return(data)
}

format_percentage_change <- format_ui(percentage_change_result, "mean_change", "lower_change", "upper_change")
format_asr2021<- format_ui(asr_2021_data, "asr_2021", "asr_2021_lower", "asr_2021_upper", 0)
format_num2021 <- format_ui(num_2021_data, "n_2021", "n_2021_lower", "n_2021_upper", 0)

num_clean <- format_num2021 %>%
  dplyr::select(location_name, sex_name, measure_name, 
                n_2021_formatted = formatted_result)

asr_clean <- format_asr2021 %>%
  dplyr::select(location_name, sex_name, measure_name,
                asr_2021_formatted = formatted_result)

change_clean <- format_percentage_change %>%
  dplyr::select(location_name, sex_name, measure_name,
                change_formatted = formatted_result)

final_table <- merge(merge(num_clean, asr_clean, 
                           by = c("location_name", "sex_name", "measure_name")),
                     change_clean,
                     by = c("location_name", "sex_name", "measure_name"))

groups <- list(
  Global = filter(final_table, location_name == "Global", sex_name == "Both"),
  Male = filter(final_table, location_name == "Global", sex_name == "Male"),
  Female = filter(final_table, location_name == "Global", sex_name == "Female"),
  SDI_Regions = filter(final_table, location_name %in% c("High SDI", "High-middle SDI", "Middle SDI", "Low-middle SDI", "Low SDI"), sex_name == "Both") %>%
    mutate(location_name = factor(location_name, levels = c("High SDI", "High-middle SDI", "Middle SDI", "Low-middle SDI", "Low SDI"))) %>%
    arrange(location_name),
  GBD_Regions = filter(final_table, !(location_name %in% c("Global", "High SDI", "High-middle SDI", "Middle SDI", "Low-middle SDI", "Low SDI")), sex_name == "Both")
  
)

super_regions <- c(
  "High-income",
  "Central Europe, Eastern Europe, and Central Asia", 
  "Latin America and Caribbean",
  "North Africa and Middle East",
  "South Asia",
  "Southeast Asia, East Asia, and Oceania",
  "Sub-Saharan Africa"
)

super_regions_with_children <- c(
  "High-income",
  "Central Europe, Eastern Europe, and Central Asia",
  "Latin America and Caribbean", 
  "Southeast Asia, East Asia, and Oceania",
  "Sub-Saharan Africa"
)

duplicate_super_regions <- c(
  "North Africa and Middle East",
  "South Asia"
)

make_summary_table <- function(df) {
  
  incidence <- df %>%
    filter(measure_name == "Incidence") %>%
    dplyr::select(location_name, sex_name, n_2021_formatted, asr_2021_formatted, change_formatted) %>%
    rename_with(~paste0("incidence_", .), -c(location_name, sex_name))
  
  prevalence <- df %>%
    filter(measure_name == "Prevalence") %>%
    dplyr::select(location_name, sex_name, n_2021_formatted, asr_2021_formatted, change_formatted) %>%
    rename_with(~paste0("prevalence_", .), -c(location_name, sex_name))
  
  dalys <- df %>%
    filter(measure_name == "DALYs (Disability-Adjusted Life Years)") %>%
    dplyr::select(location_name, sex_name, n_2021_formatted, asr_2021_formatted, change_formatted) %>%
    rename_with(~paste0("dalys_", .), -c(location_name, sex_name))
  
  summary_table <- incidence %>%
    inner_join(prevalence, by = c("location_name", "sex_name")) %>%
    inner_join(dalys, by = c("location_name", "sex_name"))
  
  return(summary_table)
}

final_summary <- tibble()

for (group_name in names(groups)) {
  df <- groups[[group_name]]
  
  if (group_name == "GBD_Regions") {
  
    df_with_super <- df %>%
      mutate(
        super_region = case_when(
          location_name %in% super_regions ~ location_name,
          location_name == "High-income Asia Pacific" ~ "High-income",
          location_name == "High-income North America" ~ "High-income", 
          location_name == "Western Europe" ~ "High-income",
          location_name == "Australasia" ~ "High-income",
          location_name == "Southern Latin America" ~ "High-income",
          
          location_name == "Central Europe" ~ "Central Europe, Eastern Europe, and Central Asia",
          location_name == "Eastern Europe" ~ "Central Europe, Eastern Europe, and Central Asia",
          location_name == "Central Asia" ~ "Central Europe, Eastern Europe, and Central Asia",
          
          location_name == "Andean Latin America" ~ "Latin America and Caribbean",
          location_name == "Central Latin America" ~ "Latin America and Caribbean",
          location_name == "Tropical Latin America" ~ "Latin America and Caribbean",
          location_name == "Caribbean" ~ "Latin America and Caribbean",
          
          location_name == "North Africa and Middle East" ~ "North Africa and Middle East",
          
          location_name == "South Asia" ~ "South Asia",
          
          location_name == "East Asia" ~ "Southeast Asia, East Asia, and Oceania",
          location_name == "Southeast Asia" ~ "Southeast Asia, East Asia, and Oceania",
          location_name == "Oceania" ~ "Southeast Asia, East Asia, and Oceania",
          
          location_name == "Central Sub-Saharan Africa" ~ "Sub-Saharan Africa",
          location_name == "Eastern Sub-Saharan Africa" ~ "Sub-Saharan Africa",
          location_name == "Southern Sub-Saharan Africa" ~ "Sub-Saharan Africa",
          location_name == "Western Sub-Saharan Africa" ~ "Sub-Saharan Africa",
          TRUE ~ "Other"
        ),
        level = case_when(
          location_name %in% super_regions ~ "super_region",
          TRUE ~ "specific_region"
        ),
        show_row = case_when(
          TRUE ~ TRUE
        )
      ) %>%
      filter(show_row == TRUE)
    
    summary_table <- make_summary_table(df_with_super)
    summary_table <- summary_table %>%
      mutate(
        group = group_name,
        super_region = df_with_super$super_region[match(location_name, df_with_super$location_name)],
        level = df_with_super$level[match(location_name, df_with_super$location_name)]
      )
    
    summary_table <- summary_table %>%
      mutate(
        sub_group = super_region,
        display_group = case_when(
       
          level == "super_region" ~ location_name,      
       
          level == "specific_region" ~ paste0("\u00A0\u00A0\u00A0\u00A0", location_name), 
          TRUE ~ location_name
        )
      ) %>%
      arrange(
        super_region,  
        case_when(
          level == "super_region" ~ 1,    
          level == "specific_region" ~ 2   
        ),
        location_name 
      )
    
  } else {
    summary_table <- make_summary_table(df)
    summary_table <- summary_table %>%
      mutate(group = group_name)
    
    if (group_name %in% c("Male", "Female")) {
      summary_table$sub_group <- summary_table$sex_name
      summary_table$display_group <- group_name
    } else if (group_name == "Global") {
      summary_table$sub_group <- summary_table$location_name
      summary_table$display_group <- group_name
    } else {
      summary_table$sub_group <- summary_table$location_name
      summary_table$display_group <- paste0("  ", summary_table$location_name)
    }
  }
  
  final_summary <- bind_rows(final_summary, summary_table)
}

final_summary <- final_summary %>%
  dplyr::select(-any_of(c("sex_name", "location_name")))

display_table <- final_summary %>%
  dplyr::select(group, display_group, 
                starts_with("incidence_"), 
                starts_with("prevalence_"), 
                starts_with("dalys_"))

kable_data <- display_table %>%
  dplyr::select(-group)

current_row <- 1
pack_rows_info <- list()

for(group_name in unique(display_table$group)) {
  group_rows <- sum(display_table$group == group_name)
  pack_rows_info[[group_name]] <- c(current_row, current_row + group_rows - 1)
  current_row <- current_row + group_rows
}

result_table <- kable_data %>%
  kable(escape = FALSE, 
        format = "html",
        align = c("l", rep("c", 9)), 
        col.names = c(
          "Groups",
          "N,2021", "Age-standardized rate (per 100,000),2021", "Percentage change,1990-2021(%)",
          "N,2021", "Age-standardized rate (per 100,000),2021", "Percentage change,1990-2021(%)",
          "N,2021", "Age-standardized rate (per 100,000),2021", "Percentage change,1990-2021(%)"
        )) %>%
  add_header_above(c(" " = 1, "Incidence" = 3, "Prevalence" = 3, "DALYs (Disability-Adjusted Life Years)" = 3)) %>%
  kable_styling(full_width = TRUE, bootstrap_options = c("striped", "hover"))

for(group_name in names(pack_rows_info)) {
  start_row <- pack_rows_info[[group_name]][1]
  end_row <- pack_rows_info[[group_name]][2]

  if(group_name %in% c("Global", "Male", "Female") && (end_row - start_row + 1) == 1) {
    next
  }
  
  result_table <- result_table %>%
    pack_rows(group_name, start_row, end_row,
              label_row_css = "background-color: #f5f5f5; border-top: 2px solid #ddd; font-weight: normal !important;")
}

print(result_table)

display_table_clean <- display_table %>%
  mutate(across(where(is.character), ~ gsub("<br>", " ", .)))

write.csv(display_table_clean, file.path(output, file = "01.summary_table_burdenanalysis.csv"), row.names = FALSE)
writeLines(as.character(result_table), 
           file.path(output, "02.summary_table_burdenanalysis.html"))
