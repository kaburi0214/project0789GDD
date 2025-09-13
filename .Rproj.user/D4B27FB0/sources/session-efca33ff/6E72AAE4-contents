rm(list = ls()); gc()

ORIGINAL_DIR <- ""
output <- file.path(ORIGINAL_DIR, "05_interaction")

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

setwd(ORIGINAL_DIR)

library(stats)
library(readxl)
library(countrycode)
library(dplyr)
library(tidyr)
library(ggplot2)

ssb_countries <- read.csv('00_rawdata/Country-level estimates/v15_cnty.csv')
gbd_data <- read.csv('00_rawdata/IHME-GBD_2021_DATA-7cf7e249-1.csv')
SDI_data <- read_excel('00_rawdata/IHME_GBD_SDI_2021_SDI_1990_2021_Y2024M05D16.XLSX')
elevation_data <- read.csv('04_altitude_gout/01.elevation_data.csv')

SDI_2018 <- SDI_data %>%
  dplyr::select(location = 1, sdi_2018 = 30) %>%
  filter(row_number() > 1, 
         !is.na(sdi_2018),
         location != "Location") %>%
  mutate(sdi_2018 = as.numeric(sdi_2018))

SDI_2018_countries <- SDI_2018%>%
  mutate(iso3 = countrycode(location, "country.name", "iso3c")) %>%
  filter(!is.na(iso3))

ssb_countries <- ssb_countries %>%
  filter(age != 999, female != 999, urban != 999, edu != 999, year == 2018) %>%
  arrange(iso3, female, age, edu, urban)

ssb_with_sdi <- merge(ssb_countries, SDI_2018_countries, 
                      by = "iso3", 
                      all.x = TRUE, sort = FALSE) %>%
  arrange(iso3, female, age, edu, urban)


ssb_with_elevation <- merge(ssb_with_sdi, elevation_data, 
                            by.x = "iso3", by.y = "ISO3", 
                            all.x = TRUE, sort = FALSE) %>%
  arrange(iso3, female, age, edu, urban)

gbd_gout_2021 <- gbd_data %>%
  filter(
         sex_name %in% c("Male", "Female")) %>%  
  mutate(female = ifelse(sex_name == "Female", 1, 0)) %>%  
  dplyr::select(location_name, female, measure_name, val, lower, upper) %>%
  pivot_wider(names_from = measure_name, 
              values_from = c(val, lower, upper)) %>%
  mutate(iso3 = countrycode(location_name, "country.name", "iso3c"))

gout_ssb <- merge(ssb_with_elevation, gbd_gout_2021, 
                    by = c("iso3", "female"), 
                    all.x = TRUE) %>%
  arrange(iso3, female, age, edu, urban)

gout_ssbz <- gout_ssb %>%
  dplyr::select(
  
    iso3, female, age, urban, edu,
    
    median,
     
    sdi_2018, Elevation,
 
    val_Incidence, val_Prevalence, `val_DALYs (Disability-Adjusted Life Years)`
  ) %>%
  filter(!is.na(median), !is.na(Elevation),
         age %in% seq(22.5, 97.5, by = 5)) %>%  
  mutate(
    SSB_z = scale(median)[,1], 
    Elevation_z = scale(Elevation)[,1]
  )


model_ASIR <- lm(val_Incidence ~ SSB_z + Elevation_z + SSB_z:Elevation_z + 
                   sdi_2018 + factor(urban) + factor(female) + factor(age) + factor(edu), 
                 data = gout_ssbz)

model_ASPR <- lm(val_Prevalence ~ SSB_z + Elevation_z + SSB_z:Elevation_z + 
                   sdi_2018 + factor(urban) + factor(female) + factor(age) + factor(edu), 
                 data = gout_ssbz)

model_ASDR <- lm(`val_DALYs (Disability-Adjusted Life Years)` ~ SSB_z + Elevation_z + SSB_z:Elevation_z + 
                   sdi_2018 + factor(urban) + factor(female) + factor(age) + factor(edu), 
                 data = gout_ssbz)


summary(model_ASIR)$coefficients["SSB_z:Elevation_z", ]
summary(model_ASPR)$coefficients["SSB_z:Elevation_z", ]
summary(model_ASDR)$coefficients["SSB_z:Elevation_z", ]


conf_asir <- confint(model_ASIR)["SSB_z:Elevation_z", ]
conf_aspr <- confint(model_ASPR)["SSB_z:Elevation_z", ]
conf_asdr <- confint(model_ASDR)["SSB_z:Elevation_z", ]


tertile_breaks <- quantile(gout_ssbz$median, probs = seq(0, 1, 1/3), na.rm = TRUE)

gout_ssbz <- gout_ssbz %>%
  mutate(SSB_group = cut(median, 
                         breaks = tertile_breaks,
                         labels = c("T1 (Low)", "T2 (Medium)", "T3 (High)"),
                         include.lowest = TRUE))


p_ASIR <- ggplot(gout_ssbz, aes(x = Elevation, y = val_Incidence)) +
  geom_point(alpha = 0.6, size = 0.8, color = "pink") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "red", fill = "lightpink") +
  facet_wrap(~SSB_group) +
  labs(x = "Elevation (m)", 
       y = "ASIR (per 100,000)",
       title = "Interaction Effects of SSBs Intake and Altitude on Gout ASIR") +
  theme_bw() +
  theme(plot.title = element_text(size = 12, hjust = 0.5))


p_ASPR <- ggplot(gout_ssbz, aes(x = Elevation, y = val_Prevalence)) +
  geom_point(alpha = 0.6, size = 0.8, color = "pink") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "red", fill = "lightpink") +
  facet_wrap(~SSB_group) +
  labs(x = "Elevation (m)", 
       y = "ASPR (per 100,000)",
       title = "Interaction Effects of SSBs Intake and Altitude on Gout ASPR") +
  theme_bw() +
  theme(plot.title = element_text(size = 12, hjust = 0.5))

 
p_ASDR <- ggplot(gout_ssbz, aes(x = Elevation, y = `val_DALYs (Disability-Adjusted Life Years)`)) +
  geom_point(alpha = 0.6, size = 0.8, color = "pink") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "red", fill = "lightpink") +
  facet_wrap(~SSB_group) +
  labs(x = "Elevation (m)", 
       y = "ASDR (DALYs per 100,000)",
       title = "Interaction Effects of SSBs Intake and Altitude on Gout ASDR") +
  theme_bw() +
  theme(plot.title = element_text(size = 12, hjust = 0.5))


interaction_summary <- data.frame(
  indicator = c("ASIR", "ASPR", "ASDR"),
  coefficient = c(
    summary(model_ASIR)$coefficients["SSB_z:Elevation_z", "Estimate"],
    summary(model_ASPR)$coefficients["SSB_z:Elevation_z", "Estimate"], 
    summary(model_ASDR)$coefficients["SSB_z:Elevation_z", "Estimate"]
  ),
  p_value = c(
    summary(model_ASIR)$coefficients["SSB_z:Elevation_z", "Pr(>|t|)"],
    summary(model_ASPR)$coefficients["SSB_z:Elevation_z", "Pr(>|t|)"],
    summary(model_ASDR)$coefficients["SSB_z:Elevation_z", "Pr(>|t|)"]
  )
)

write.csv(interaction_summary, file.path(output, "01.interaction_results.csv"), row.names = FALSE)
ggsave(file.path(output, "02.ASIR_interaction.png"), p_ASIR, width = 12, height = 6, dpi = 300)
ggsave(file.path(output, "03.ASPR_interaction.png"), p_ASPR, width = 12, height = 6, dpi = 300)
ggsave(file.path(output, "04.ASDR_interaction.png"), p_ASDR, width = 12, height = 6, dpi = 300)
ggsave(file.path(output, "02.ASIR_interaction.pdf"), p_ASIR, width = 12, height = 6)
ggsave(file.path(output, "03.ASPR_interaction.pdf"), p_ASPR, width = 12, height = 6)
ggsave(file.path(output, "04.ASDR_interaction.pdf"), p_ASDR, width = 12, height = 6)