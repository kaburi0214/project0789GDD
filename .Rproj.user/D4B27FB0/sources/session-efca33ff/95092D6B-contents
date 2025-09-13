# 当前脚本应该在01_project0789目录下，而不是在分析点目录05_interaction里面
# 清除当前环境所有变量并返回内存使用情况
rm(list = ls()); gc()

# 定义两个变量,output和ORIGINAL_DIR
ORIGINAL_DIR <- "/data/nas1/zhangtongrui_OD/project/01_project0789"
output <- file.path(ORIGINAL_DIR, "05_interaction")

# 确保输出目录存在
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

# 切换到项目目录下
setwd(ORIGINAL_DIR)

# 加载必要的包
library(stats)
library(readxl)
library(countrycode)
library(dplyr)
library(tidyr)
library(ggplot2)

# 读取并整理数据
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

# 提取按性别分层的痛风数据
gbd_gout_2021 <- gbd_data %>%
  filter(
         sex_name %in% c("Male", "Female")) %>%  # 按性别分开
  mutate(female = ifelse(sex_name == "Female", 1, 0)) %>%  # 转换为0/1编码
  dplyr::select(location_name, female, measure_name, val, lower, upper) %>%
  pivot_wider(names_from = measure_name, 
              values_from = c(val, lower, upper)) %>%
  mutate(iso3 = countrycode(location_name, "country.name", "iso3c"))

# 合并痛风数据（需要同时匹配iso3和female）
gout_ssb <- merge(ssb_with_elevation, gbd_gout_2021, 
                    by = c("iso3", "female"), 
                    all.x = TRUE) %>%
  arrange(iso3, female, age, edu, urban)

# 在gout_ssbz数据处理步骤中添加年龄筛选
gout_ssbz <- gout_ssb %>%
  dplyr::select(
    # 分组变量
    iso3, female, age, urban, edu,
    # 暴露变量
    median,
    # 协变量  
    sdi_2018, Elevation,
    # 响应变量
    val_Incidence, val_Prevalence, `val_DALYs (Disability-Adjusted Life Years)`
  ) %>%
  filter(!is.na(median), !is.na(Elevation),
         age %in% seq(22.5, 97.5, by = 5)) %>%  # 22.5, 27.5, 32.5, ..., 97.5
  mutate(
    SSB_z = scale(median)[,1], # 标准化SSB摄入量
    Elevation_z = scale(Elevation)[,1] # 标准化海拔高度
  )

# 构建包含主效应和交互项的多元线性回归模型
# 对三个痛风负担指标分别建模
# 因变量：痛风疾病负担指标
# 自变量：SSB摄入量、海拔高度及其交互项
# 控制变量：SDI、居住地区、性别、年龄、教育水平
model_ASIR <- lm(val_Incidence ~ SSB_z + Elevation_z + SSB_z:Elevation_z + 
                   sdi_2018 + factor(urban) + factor(female) + factor(age) + factor(edu), 
                 data = gout_ssbz)

model_ASPR <- lm(val_Prevalence ~ SSB_z + Elevation_z + SSB_z:Elevation_z + 
                   sdi_2018 + factor(urban) + factor(female) + factor(age) + factor(edu), 
                 data = gout_ssbz)

model_ASDR <- lm(`val_DALYs (Disability-Adjusted Life Years)` ~ SSB_z + Elevation_z + SSB_z:Elevation_z + 
                   sdi_2018 + factor(urban) + factor(female) + factor(age) + factor(edu), 
                 data = gout_ssbz)

# 查看交互项系数以及显著性
summary(model_ASIR)$coefficients["SSB_z:Elevation_z", ]
summary(model_ASPR)$coefficients["SSB_z:Elevation_z", ]
summary(model_ASDR)$coefficients["SSB_z:Elevation_z", ]

# 提取交互项95%置信区间
conf_asir <- confint(model_ASIR)["SSB_z:Elevation_z", ]
conf_aspr <- confint(model_ASPR)["SSB_z:Elevation_z", ]
conf_asdr <- confint(model_ASDR)["SSB_z:Elevation_z", ]

# 计算三分位点（包含边界）
tertile_breaks <- quantile(gout_ssbz$median, probs = seq(0, 1, 1/3), na.rm = TRUE)

gout_ssbz <- gout_ssbz %>%
  mutate(SSB_group = cut(median, 
                         breaks = tertile_breaks,
                         labels = c("T1 (Low)", "T2 (Medium)", "T3 (High)"),
                         include.lowest = TRUE))

# ASIR图
p_ASIR <- ggplot(gout_ssbz, aes(x = Elevation, y = val_Incidence)) +
  geom_point(alpha = 0.6, size = 0.8, color = "pink") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "red", fill = "lightpink") +
  facet_wrap(~SSB_group) +
  labs(x = "Elevation (m)", 
       y = "ASIR (per 100,000)",
       title = "Interaction Effects of SSBs Intake and Altitude on Gout ASIR") +
  theme_bw() +
  theme(plot.title = element_text(size = 12, hjust = 0.5))

# ASPR图
p_ASPR <- ggplot(gout_ssbz, aes(x = Elevation, y = val_Prevalence)) +
  geom_point(alpha = 0.6, size = 0.8, color = "pink") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "red", fill = "lightpink") +
  facet_wrap(~SSB_group) +
  labs(x = "Elevation (m)", 
       y = "ASPR (per 100,000)",
       title = "Interaction Effects of SSBs Intake and Altitude on Gout ASPR") +
  theme_bw() +
  theme(plot.title = element_text(size = 12, hjust = 0.5))

# ASDR图  
p_ASDR <- ggplot(gout_ssbz, aes(x = Elevation, y = `val_DALYs (Disability-Adjusted Life Years)`)) +
  geom_point(alpha = 0.6, size = 0.8, color = "pink") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "red", fill = "lightpink") +
  facet_wrap(~SSB_group) +
  labs(x = "Elevation (m)", 
       y = "ASDR (DALYs per 100,000)",
       title = "Interaction Effects of SSBs Intake and Altitude on Gout ASDR") +
  theme_bw() +
  theme(plot.title = element_text(size = 12, hjust = 0.5))

# 保存分析结果和图片
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