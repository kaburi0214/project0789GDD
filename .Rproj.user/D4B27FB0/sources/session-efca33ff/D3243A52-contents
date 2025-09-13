rm(list = ls()); gc()

ORIGINAL_DIR <- ""
output <- file.path(ORIGINAL_DIR, "04_altitude_gout")

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

setwd(ORIGINAL_DIR)

set.seed(42)


library(ggplot2)
library(dplyr)
library(geodata)
library(rnaturalearth)
library(exactextractr)
library(sf)
library(terra)


gbd_sample <- read.csv("00_rawdata/IHME-GBD_2021_DATA-7cf7e249-1.csv")


if (file.exists(file.path(output, "01.elevation_data.csv"))) {
  elevation_data <- read.csv(file.path(output, "01.elevation_data.csv"))
} else {

  world_countries <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  gbd_countries_names <- unique(gbd_sample$location_name)
  gbd_countries <- world_countries[world_countries$name %in% gbd_countries_names, ]
  
 
  gbd_countries$iso_a3[gbd_countries$name == "France"] <- "FRA"
  gbd_countries$iso_a3[gbd_countries$name == "Norway"] <- "NOR"
  
  
  global_pop <- geodata::population(year = 2020, res = 0.5, path = tempdir())
  country_elevation <- numeric(nrow(gbd_countries))
  
  for(i in 1:nrow(gbd_countries)) {
    country_iso3 <- gbd_countries$iso_a3[i] 
    
    if(is.na(country_iso3) || country_iso3 == "-99") {
      country_elevation[i] <- NA
      next
    }
    
    tryCatch({
      
      elevation_raster <- geodata::elevation_30s(country = country_iso3, path = tempdir())
      
    
      country_pop <- terra::crop(global_pop, gbd_countries[i,])
      country_pop <- terra::mask(country_pop, gbd_countries[i,])
      country_pop_resampled <- terra::resample(country_pop, elevation_raster, method = "bilinear")
      
      elev_vals <- terra::values(elevation_raster)
      pop_vals <- terra::values(country_pop_resampled)
      
      valid_idx <- !is.na(elev_vals) & !is.na(pop_vals) & pop_vals > 0
      
      if(sum(valid_idx) > 0) {
        mean_elev <- weighted.mean(elev_vals[valid_idx], pop_vals[valid_idx])
      } else {
        mean_elev <- exactextractr::exact_extract(elevation_raster, gbd_countries[i,], 'mean')
      }
      
      if(is.na(mean_elev) || mean_elev == 0) {
        country_elevation[i] <- NA
      } else {
        country_elevation[i] <- mean_elev
      }
      
    }, error = function(e) {
      country_elevation[i] <- NA
    })
  }
  
 
  elevation_data <- data.frame(
    Country = gbd_countries$name,
    ISO3 = gbd_countries$iso_a3,
    Elevation = country_elevation
  )
  elevation_data <- elevation_data[!is.na(elevation_data$Elevation), ]
  write.csv(elevation_data, file.path(output, "01.elevation_data.csv"), row.names = FALSE)
  }


process_gout_indicator <- function(indicator_type) {
  

  measure_filter <- switch(indicator_type,
                           "asir" = "Incidence",
                           "aspr" = "Prevalence", 
                           "asdr" = "DALYs"
  )
  

  gbd_processed <- gbd_sample %>%
    filter(sex_name == "Both", grepl(measure_filter, measure_name, ignore.case = TRUE)) %>%
    dplyr::select(location_name, val) %>%
    rename(Country = location_name) %>%
    rename(!!indicator_type := val)
  

  merged_data <- merge(elevation_data, gbd_processed, by = "Country", all = FALSE)
  merged_data <- merged_data[complete.cases(merged_data), ]
  

  cor_result <- cor.test(merged_data$Elevation, merged_data[[indicator_type]], method = "spearman")
  significant <- abs(cor_result$estimate) > 0.3 & cor_result$p.value < 0.05 
  
  
  plot <- ggplot(merged_data, aes(x = Elevation, y = .data[[indicator_type]])) +
    geom_point(color = "#8B0000", size = 2.5, alpha = 0.7) +
    geom_smooth(method = "loess", se = TRUE, color = "#8B0000", fill = "#FFB6C1", alpha = 0.3) +
    annotate("text", x = Inf, y = Inf, 
             label = paste0("r = ", round(cor_result$estimate, 3), 
                            ", p ", 
                            ifelse(cor_result$p.value < 0.001, "< 0.001", 
                                   paste("=", round(cor_result$p.value, 3)))),
             hjust = 1.1, vjust = 1.5, size = 4) +
    labs(title = paste(toupper(indicator_type), "under different altitudes"),
         x = "Elevation (meters)", y = paste0(toupper(indicator_type), " (per 100,000)")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
          panel.border = element_rect(color = "black", fill = NA))
  
  return(list(data = merged_data, correlation = cor_result$estimate, 
              p_value = cor_result$p.value, significant = significant, plot = plot))
}


asir_results <- process_gout_indicator("asir")
aspr_results <- process_gout_indicator("aspr")
asdr_results <- process_gout_indicator("asdr")


correlation_summary <- data.frame(
  Indicator = c("ASIR", "ASPR", "ASDR"),
  Correlation = c(asir_results$correlation, aspr_results$correlation, asdr_results$correlation),
  P_value = c(asir_results$p_value, aspr_results$p_value, asdr_results$p_value),
  Significant = c(asir_results$significant, aspr_results$significant, asdr_results$significant),
  Sample_Size = c(nrow(asir_results$data), nrow(aspr_results$data), nrow(asdr_results$data))
)

print(correlation_summary)


print(asir_results$plot)
print(aspr_results$plot)
print(asdr_results$plot)


write.csv(correlation_summary, file.path(output, "02.correlation_results.csv"), row.names = FALSE)
ggsave(file.path(output, "03.asir_plot.png"), plot = asir_results$plot, width = 10, height = 6, dpi = 300)
ggsave(file.path(output,"04.aspr_plot.png"), plot = aspr_results$plot, width = 10, height = 6, dpi = 300)
ggsave(file.path(output,"05.asdr_plot.png"), plot = asdr_results$plot, width = 10, height = 6, dpi = 300)
ggsave(file.path(output, "03.asir_plot.pdf"), plot = asir_results$plot, width = 10, height = 6)
ggsave(file.path(output,"04.aspr_plot.pdf"), plot = aspr_results$plot, width = 10, height = 6)
ggsave(file.path(output,"05.asdr_plot.pdf"), plot = asdr_results$plot, width = 10, height = 6)

