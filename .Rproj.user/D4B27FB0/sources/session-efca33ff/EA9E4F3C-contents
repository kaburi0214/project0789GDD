rm(list = ls()); gc()

ORIGINAL_DIR <- ""
output <- file.path(ORIGINAL_DIR, "00_rawdata")

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

setwd(ORIGINAL_DIR)

required_files <- c(
  "00_rawdata/IHME-GBD_2021_DATA-ec7daba0-1.csv",
  "00_rawdata/Global estimates/v15_global.csv",
  "00_rawdata/Regional estimates/v15_superregion2.csv",
  "00_rawdata/Country-level estimates/v15_cnty.csv",
  "00_rawdata/IHME-GBD_2021_DATA-ec7daba0-1.csv",
  "00_rawdata/IHME-GBD_2021_DATA-7cf7e249-1.csv",
  "00_rawdata/IHME_GBD_SDI_2021_SDI_1990_2021_Y2024M05D16.XLSX",
  "00_rawdata/IHME-GBD_2021_DATA-7010f9b0-1.csv",
  "00_rawdata/IHME-GBD_2021_DATA-d92a7383-1.csv",
  "00_rawdata/GDD 2018 Codebook_Jan 10 2022.xlsx"
)

cat("inspecting...\n")
missing_files <- c()

for(file in required_files) {
  if(file.exists(file)) {
    cat(sprintf("✓ %s\n", file))
  } else {
    cat(sprintf("✗ %s\n", file))
    missing_files <- c(missing_files, file)
  }
}

if(length(missing_files) > 0) {
  stop(sprintf("\nlack files:\n%s\n", 
               paste(missing_files, collapse = "\n")))
}

cat("\ninspection completed！\n\n")