# 当前脚本应该在01_project0789目录下，而不是在分析点目录00_rawdata里面
# 清除当前环境所有变量并返回内存使用情况
rm(list = ls()); gc()

# 定义两个变量,output和ORIGINAL_DIR
ORIGINAL_DIR <- "/data/nas1/zhangtongrui_OD/project/01_project0789"
output <- file.path(ORIGINAL_DIR, "00_rawdata")

# 确保输出目录存在
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

# 切换到项目目录下
setwd(ORIGINAL_DIR)

# 检查数据集存在
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