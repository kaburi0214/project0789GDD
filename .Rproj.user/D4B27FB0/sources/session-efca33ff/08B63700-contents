# 当前脚本应该在01_project0789目录下，而不是在分析点目录02_ssbintk里面
# 清除当前环境所有变量并返回内存使用情况
rm(list = ls()); gc()

# 定义两个变量,output和ORIGINAL_DIR
ORIGINAL_DIR <- "/data/nas1/zhangtongrui_OD/project/01_project0789"
output <- file.path(ORIGINAL_DIR, "02_ssbintk")

# 确保输出目录存在
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

# 切换到项目目录下
setwd(ORIGINAL_DIR)

# 加载包
library(dplyr)
library(brms)
library(tidyr)
library(purrr)

# 设置全局随机种子以确保结果可重复
set.seed(42)

# 加载数据集
ssb_global <- read.csv('00_rawdata/Global estimates/v15_global.csv')
ssb_countries <- read.csv('00_rawdata/Country-level estimates/v15_cnty.csv')
# 定义按人口数排名前30的国家ISO3代码
top30 <- c("CHN", "IND", "USA", "IDN", "BRA","PAK",
           "RUS", "JPN", "BGD", "NGA", "MEX","DEU",
           "VNM", "PHL", "EGY", "IRN", "TUR","THA",
           "ETH", "GBR", "ITA", "FRA", "KOR", "ESP",
           "COD", "ZAF", "UKR", "MMR", "COL","ARG")

# 筛选并清理数据，"999"代表all levels
ssb_global_18 <- ssb_global %>%
  filter(
    year == 2018,
    age != 999,
    female != 999,
    urban != 999,
    edu != 999,
    !is.na(median),
    !is.na(lowerci_95),
    !is.na(upperci_95)
  )

ssb_countries_18 <- ssb_countries %>%
  filter(
    year == 2018,
    age!= 999,
    female != 999,
    urban != 999,
    edu != 999,
    !is.na(median),
    !is.na(lowerci_95),
    !is.na(upperci_95)
  )

#检查数据分布,SSB摄入量取对数后更接近正态分布
medians <- ssb_global_18$median
hist(medians, breaks = 40,
     main = "Histogram of SSB Intake Median",
     xlab = "Median Intake (g/day)",
     col = "skyblue", border = "white")
log_medians <- log(medians)
hist(log_medians, breaks = 40,
     main = "Histogram of log(Median SSB Intake)",
     xlab = "Log(Median Intake)", col = "lightgreen")
qqnorm(log_medians); qqline(log_medians, col = "red")

# 将分类变量转换为因子
ssb_countries_18 <- ssb_countries_18 %>%
  mutate(
    age = factor(age, levels = c(0.5, 1.5, 3.5, seq(7.5, 97.5, by = 5))),
    gender = factor(female, levels = c(0, 1), labels = c("Male", "Female")),
    education = factor(edu, levels = c(1, 2, 3), labels = c("0-6", "6-12", "≥12")),
    residence = factor(urban, levels = c(0, 1), labels = c("Rural", "Urban")),
    superregion = factor(superregion2, levels = c("Asia", "FSU", "HIC", "LAC", "MENA", "SAARC", "SSA")),
    country = factor(iso3)
  )

# 拟合贝叶斯模型
if (file.exists(file.path(output, "ssb_all_country_bm.rds"))) {
  model_cnty <- readRDS(file.path(output, "ssb_all_country_bm.rds"))
} else {
  model_cnty <- brm(
    formula = log(median) ~ age + gender + education + residence + (1 | superregion / country), # 贝叶斯分层模型：log变换的SSB摄入量 ~ 固定效应(人口学因素) + 随机效应(地理层级)
    data = ssb_countries_18,
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
  saveRDS(model_cnty, file.path(output, "ssb_all_country_bm.rds"))
}

model_cnty <- readRDS(file.path(output, "ssb_all_country_bm.rds"))

# 创建分层预测用的组合数据框
age_levels <- levels(ssb_countries_18$age)
gender_levels <- levels(ssb_countries_18$gender)
education_levels <- levels(ssb_countries_18$education)
residence_levels <- levels(ssb_countries_18$residence)

strata_grid <- expand.grid(
  age = age_levels,
  gender = gender_levels,
  education = education_levels,
  residence = residence_levels
)

pred_data_global <- strata_grid %>%
  mutate(superregion = NA, country = NA)

pred_data_region <- strata_grid %>%
  crossing(superregion = levels(ssb_countries_18$superregion)) %>%
  mutate(country = NA)

# 获取人口数前30country-superregion对应关系表
country_region_lookup <- ssb_countries_18 %>%
  distinct(country, superregion) %>%
  filter(as.character(country) %in% top30) %>%  # 转为character比较
  mutate(country = as.character(country))       # 转为character

pred_data_country <- strata_grid %>%
  crossing(country = top30) %>%  
  left_join(country_region_lookup, by = "country") %>%
  mutate(country = factor(country, levels = top30))

# 获取后验预测样本
get_posterior_predict <- function(model, newdata, re_formula = NULL) {
  pred <- posterior_epred(
    model,
    newdata = newdata,
    re_formula = re_formula,
    allow_new_levels = TRUE
  )
  exp(pred) # 返回指数变换后的预测值
}

# 生成预测矩阵
pred_global <- get_posterior_predict(model_cnty, pred_data_global, re_formula = NA)
pred_region <- get_posterior_predict(model_cnty, pred_data_region, re_formula = ~ (1 | superregion))
pred_country <- get_posterior_predict(model_cnty, pred_data_country, re_formula = ~ (1 | superregion / country))

# 计算中位数和95%区间的函数
calculate_summary <- function(pred_matrix) {
  # pred_matrix：samples x observations
  median_val <- apply(pred_matrix, 2, median)
  lower_val <- apply(pred_matrix, 2, quantile, 0.025)
  upper_val <- apply(pred_matrix, 2, quantile, 0.975)

  tibble(
    median = median_val,
    lower = lower_val,
    upper = upper_val,
    formatted = sprintf("%.1f (%.1f–%.1f)", median_val, lower_val, upper_val)
  )
}

# 按组计算overall及各分组层级的表
# 修正后的分层表格生成函数
generate_stratified_tables <- function(pred_matrix, pred_data, level_var, level_value) {
  
  # 参数说明：
  # pred_matrix: 后验预测矩阵
  # pred_data: 预测用的协变量数据框
  # level_var: 层级变量名
  # level_value: 层级值
  # 确保数据对齐
  complete_idx <- complete.cases(pred_data[, c("age", "gender", "education", "residence")])
  pred_data <- pred_data[complete_idx, ]
  pred_matrix <- pred_matrix[, complete_idx]

  # 添加层级信息的辅助函数
  append_level_info <- function(df) {
    df[[level_var]] <- level_value
    df
  }

  # 1. Overall（总体）
  overall_samples <- rowMeans(pred_matrix)  # 每个posterior draw的平均值
  overall <- tibble(
    median = median(overall_samples),
    lower = quantile(overall_samples, 0.025),
    upper = quantile(overall_samples, 0.975),
    formatted = sprintf("%.1f (%.1f–%.1f)", median, lower, upper),
    group = "Overall",
    subgroup = NA
  ) %>% append_level_info()

  # 2. 按性别分组
  gender_table <- map_dfr(unique(pred_data$gender), function(g) {
    idx <- which(pred_data$gender == g)
    group_samples <- rowMeans(pred_matrix[, idx, drop = FALSE])
    tibble(
      median = median(group_samples),
      lower = quantile(group_samples, 0.025),
      upper = quantile(group_samples, 0.975),
      formatted = sprintf("%.1f (%.1f–%.1f)", median, lower, upper),
      group = "Sex",
      subgroup = as.character(g)
    ) %>% append_level_info()
  })

  # 3. 按教育水平分组
  edu_table <- map_dfr(unique(pred_data$education), function(e) {
    idx <- which(pred_data$education == e)
    group_samples <- rowMeans(pred_matrix[, idx, drop = FALSE])
    tibble(
      median = median(group_samples),
      lower = quantile(group_samples, 0.025),
      upper = quantile(group_samples, 0.975),
      formatted = sprintf("%.1f (%.1f–%.1f)", median, lower, upper),
      group = "Education",
      subgroup = as.character(e)
    ) %>% append_level_info()
  })

  # 4. 按居住地分组
  res_table <- map_dfr(unique(pred_data$residence), function(r) {
    idx <- which(pred_data$residence == r)
    group_samples <- rowMeans(pred_matrix[, idx, drop = FALSE])
    tibble(
      median = median(group_samples),
      lower = quantile(group_samples, 0.025),
      upper = quantile(group_samples, 0.975),
      formatted = sprintf("%.1f (%.1f–%.1f)", median, lower, upper),
      group = "Residence",
      subgroup = as.character(r)
    ) %>% append_level_info()
  })

  # 5. 按年龄分组（≥20岁）
  age_keep <- as.numeric(as.character(pred_data$age)) >= 22.5
  age_levels <- unique(pred_data$age[age_keep])
  age_table <- map_dfr(age_levels, function(a) {
    idx <- which(pred_data$age == a)
    group_samples <- rowMeans(pred_matrix[, idx, drop = FALSE])
    age_num <- as.numeric(as.character(a))
    if(age_num >= 97.5) {
      age_label <- "95+ years"
    } else {
      age_label <- paste0(age_num - 2.5, "-", age_num + 1.5, " years")
    }
    tibble(
      median = median(group_samples),
      lower = quantile(group_samples, 0.025),
      upper = quantile(group_samples, 0.975),
      formatted = sprintf("%.1f (%.1f–%.1f)", median, lower, upper),
      group = "Age",
      subgroup = age_label
    ) %>% append_level_info()
  })

  # 合并所有结果
  bind_rows(overall, gender_table, edu_table, res_table, age_table)
}
# 生成各层级表格
table_global <- generate_stratified_tables(pred_global, pred_data_global, "superregion", "World")

table_region <- map_dfr(levels(ssb_countries_18$superregion), function(region) {
  idx <- which(pred_data_region$superregion == region)
  generate_stratified_tables(pred_region[, idx], pred_data_region[idx, ], "superregion", region)
})

table_country <- map_dfr(top30, function(cty) {
  if(cty %in% pred_data_country$country) {  
    idx <- which(pred_data_country$country == cty)
    generate_stratified_tables(pred_country[, idx], pred_data_country[idx, ], "country", cty)
  }
})

# 合并所有结果
final_table_stratified <- bind_rows(table_global, table_region, table_country)

# 查看结果
print(final_table_stratified)

# 从final_table_stratified提取数据生成两张表格
#表1: 2018年全球和七大世界区域按年龄、性别、教育程度和居住地区分列的成人含糖饮料摄入量的全球和地区平均值(95%UI)
table_1_pivot <- final_table_stratified %>%
  # 筛选全球和区域层级
  filter(superregion %in% c("World", "Asia", "FSU", "HIC", "LAC", "MENA", "SAARC", "SSA")) %>%
  # 选择需要的分组
  filter(group %in% c("Overall", "Sex", "Age", "Education", "Residence")) %>%
  # 创建行标识符
  mutate(
    row_label = case_when(
      group == "Overall" ~ "Overall",
      group == "Sex" & subgroup == "Female" ~ "Female",
      group == "Sex" & subgroup == "Male" ~ "Male",
      group == "Age" ~ paste("Age", subgroup),
      group == "Education" & subgroup == "0-6" ~ "Education 0-6 years",
      group == "Education" & subgroup == "6-12" ~ "Education 6-12 years",
      group == "Education" & subgroup == "≥12" ~ "Education ≥12 years",
      group == "Residence" & subgroup == "Rural" ~ "Rural",
      group == "Residence" & subgroup == "Urban" ~ "Urban"
    )
  ) %>%
  # 选择需要的列
  dplyr::select(superregion, row_label, formatted) %>%
  # 透视表：行变列
  pivot_wider(
    names_from = superregion,
    values_from = formatted,
    names_sort = TRUE
  ) %>%
  # 重新排列列顺序
  dplyr::select(row_label, World, Asia, FSU, HIC, LAC, MENA, SAARC, SSA)

# 表2: 2018年人口最多的30个国家按性别、教育程度和居住地区分列的成年人含糖饮料摄入量的全国平均值(95%UI)
table_2_pivot <- final_table_stratified %>%
  # 筛选国家层级
  filter(!is.na(country)) %>%
  # 选择需要的分组（不包含年龄）
  filter(group %in% c("Overall", "Sex", "Education", "Residence")) %>%
  # 创建列标识符
  mutate(
    col_label = case_when(
      group == "Overall" ~ "Overall",
      group == "Sex" & subgroup == "Female" ~ "Female",
      group == "Sex" & subgroup == "Male" ~ "Male",
      group == "Education" & subgroup == "0-6" ~ "Education 0-6 years",
      group == "Education" & subgroup == "6-12" ~ "Education 6-12 years",
      group == "Education" & subgroup == "≥12" ~ "Education ≥12 years",
      group == "Residence" & subgroup == "Rural" ~ "Rural",
      group == "Residence" & subgroup == "Urban" ~ "Urban"
    )
  ) %>%
  # 选择需要的列
  dplyr::select(country, col_label, formatted) %>%
  # 透视表：列变行
  pivot_wider(
    names_from = col_label,
    values_from = formatted
  ) %>%
  # 重新排列列顺序
  dplyr::select(country, Overall, Female, Male, `Education 0-6 years`, `Education 6-12 years`, `Education ≥12 years`, Rural, Urban)

# 保存结果文件
write.csv(table_1_pivot,
          file.path(output, "01.ssbsintake_global_region.csv"),
          row.names = FALSE)

write.csv(table_2_pivot,
          file.path(output, "02.ssbsintake_top30countries.csv"),
          row.names = FALSE)
