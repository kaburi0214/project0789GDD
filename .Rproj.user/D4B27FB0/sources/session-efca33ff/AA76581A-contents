rm(list = ls()); gc()

ORIGINAL_DIR <- ""
output <- file.path(ORIGINAL_DIR, "08_2050pred")

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

setwd(ORIGINAL_DIR)

library(segmented)
library(forecast)
library(tseries)
library(dplyr)
library(ggplot2)

gbd_32 <- read.csv("00_rawdata/IHME-GBD_2021_DATA-ec7daba0-1.csv")

global_data <- gbd_32 %>%
  filter(
    location_name == "Global",
    age_name == "Age-standardized", 
    sex_name == "Both",
    metric_name == "Rate"
  ) %>%
  arrange(year) %>%
  dplyr::select(measure_name, year, val, upper, lower)

measures <- unique(global_data$measure_name)

global_data <- global_data %>%
  group_by(measure_name) %>%
  mutate(trend = predict(lm(val ~ year))) %>%
  ungroup()

arima_results_joinpoint <- list()

for(measure in measures) {

  measure_data <- global_data %>%
    filter(measure_name == measure)
  
  measure_data$log_val <- log(measure_data$val)
  

  lm_model <- lm(log_val ~ year, data = measure_data)
  

  tryCatch({
 
    seg_model <- segmented(lm_model, seg.Z = ~year, npsi = 1)

    slopes <- slope(seg_model)
    last_slope <- slopes$year[nrow(slopes$year), 1]
    breakpoint <- seg_model$psi[,"Est."]
    
   
    recent_years <- which(measure_data$year >= breakpoint)
    if(length(recent_years) > 2) {
     
      recent_data <- measure_data[recent_years, ]
      recent_trend <- lm(val ~ year, data = recent_data[1:(length(recent_years)-2), ])
      
 
      predicted_vals <- predict(recent_trend, newdata = measure_data[recent_years, ])
      measure_data$val_adjusted <- measure_data$val
      measure_data$val_adjusted[recent_years[(length(recent_years)-1):length(recent_years)]] <- 
        predicted_vals[(length(predicted_vals)-1):length(predicted_vals)]
    } else {
      measure_data$val_adjusted <- measure_data$val
    }
    
  }, error = function(e) {
  
    measure_data$val_adjusted <- measure_data$val
  })
  

  ts_data <- ts(measure_data$val_adjusted, start = 1990, end = 2021, frequency = 1)
  

  arima_model <- auto.arima(ts_data)
  

  forecast_result <- forecast(arima_model, h = 29)
  

  arima_results_joinpoint[[measure]] <- list(
    model = arima_model,
    forecast = forecast_result,
    historical_data = measure_data,
    original_data = measure_data$val,
    adjusted_data = measure_data$val_adjusted
  )
}

measure_names <- c(
  "DALYs (Disability-Adjusted Life Years)" = "01.ASDR",
  "Incidence" = "02.ASIR", 
  "Prevalence" = "03.ASPR"
)


measure_names_png <- c(
  "DALYs (Disability-Adjusted Life Years)" = "04.ASDR",
  "Incidence" = "05.ASIR", 
  "Prevalence" = "06.ASPR"
)


for(measure in measures) {
  

  result <- arima_results_joinpoint[[measure]]
  hist_data <- result$historical_data
  forecast_data <- result$forecast
  

  arima_fitted <- data.frame(
    year = 1990:2021,
    value = as.numeric(fitted(result$model)),
    type = "ARIMA拟合"
  )
  

  arima_forecast <- data.frame(
    year = 2022:2050,
    value = as.numeric(forecast_data$mean),
    lower95 = as.numeric(forecast_data$lower[,2]),
    upper95 = as.numeric(forecast_data$upper[,2]),
    lower80 = as.numeric(forecast_data$lower[,1]),
    upper80 = as.numeric(forecast_data$upper[,1]),
    type = "ARIMA predict"
  )
  
  filename <- measure_names[measure]
  write.csv(arima_forecast, file.path(output, paste0(filename, "_forecast.csv")), row.names = FALSE)
  

  p <- ggplot() +

    geom_point(data = hist_data, aes(x = year, y = val), 
               color = "grey", size = 0.6, alpha = 0.5) +
    

    geom_line(data = hist_data, aes(x = year, y = val_adjusted), 
              color = "black", linewidth = 0.8) +
    

    geom_line(data = hist_data, aes(x = year, y = trend), 
              color = "green", linewidth = 0.8) +
    

    geom_line(data = arima_fitted, aes(x = year, y = value), 
              color = "blue", linewidth = 0.8) +
    

    geom_line(data = arima_forecast, aes(x = year, y = value), 
              color = "blue", linewidth = 0.8) +
    

    geom_ribbon(data = arima_forecast, 
                aes(x = year, ymin = lower95, ymax = upper95), 
                fill = "lightgray", alpha = 0.6) +
    

    geom_ribbon(data = arima_forecast, 
                aes(x = year, ymin = lower80, ymax = upper80), 
                fill = "blue", alpha = 0.3) +
    

    geom_vline(xintercept = 2021.5, linetype = "dotted", color = "gray") +
    
    labs(
      title = paste("Forecasts from ARIMA:", measure),
      x = "year",
      y = "ASR"
      
    ) +
    
    scale_x_continuous(breaks = seq(1990, 2050, 10)) +
    theme_classic()
  

  print(p)
  

  ggsave(file.path(output, paste0(measure_names_png[measure], "_forecast.png")), 
         plot = p, width = 10, height = 6, dpi = 300)
  ggsave(file.path(output, paste0(measure_names_png[measure], "_forecast.pdf")), 
         plot = p, width = 10, height = 6)
}
