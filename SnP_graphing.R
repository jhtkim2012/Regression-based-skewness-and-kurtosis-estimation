# 필요한 패키지 로드
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# 데이터 불러오기
#df_100 <- read_csv("C:\\Data\\skew_kurt_result_100.csv")
df_100 <- read_csv("skew_kurt_result_100.csv")


### 왜도
# 1. 데이터 형식 변환 및 skew 데이터만 추출
skew_data <- df_100 %>%
  select(Date, contains("skew")) %>%
  pivot_longer(-Date, names_to = "Method", values_to = "Skewness") %>%
  mutate(
    Date = as.Date(Date) - 50,  # x축 왼쪽으로 50일 이동
    Method = case_when(
      Method == "New5_smallresid_skew" ~ "New",
      Method == "Quantile2_skew" ~ "Groeneveld",
      Method == "lmom_skew" ~ "Lmom",
      Method == "modal_skew" ~ "Modal",
      Method == "scott_skew" ~ "KDE",
      Method == "Empirical_skew" ~ "Empirical",
      TRUE ~ Method
    )
  ) %>%
  filter(Method %in% c("Empirical", "New", "Lmom", "Groeneveld", "Modal", "KDE"))

skew_data$Method <- factor(skew_data$Method, levels = c("Empirical", "New", "Lmom", "Groeneveld", "Modal", "KDE"))


color_mapping <- c(
  "New" = "black",
  "Empirical" = "blue",
  "Lmom" = "green",
  "Groeneveld" = "gray45",
  "Modal" = "gray70",
  "KDE" = "red"
)

linetype_mapping <- c(
  "New" = "solid",
  "Empirical" = "solid",
  "Lmom" = "solid",
  "Groeneveld" = "dashed",
  "Modal" = "dotted",
  "KDE" = "solid"
)


ggplot(skew_data, aes(x = Date, y = Skewness, color = Method, linetype = Method)) +
  geom_line(size = 0.5) +
  scale_color_manual(values = color_mapping) +
  scale_linetype_manual(values = linetype_mapping) +
  scale_x_date(
    date_breaks = "6 month",
    date_labels = "%Y.%m",
    limits = c(as.Date("2008-05-27") , as.Date("2023-12-29") ),
    expand = c(0, 0)
  ) +
  labs(title = "",
       x = "", y = "Skewness") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    text = element_text(size = 14)
  )



### 첨도
# 1. 데이터 형식 변환 및 kurt 데이터만 추출
kurt_data <- df_100 %>%
  select(Date, contains("kurt")) %>%
  pivot_longer(-Date, names_to = "Method", values_to = "Kurtosis") %>%
  mutate(
    Date = as.Date(Date) - 50,  # x축 왼쪽으로 50일 이동
    Method = case_when(
      Method == "New5_smallresid_kurt" ~ "New",
      Method == "Quantile1_kurt" ~ "Moors",
      Method == "Quantile3_kurt" ~ "Crow",
      Method == "lmom_kurt" ~ "Lmom",
      Method == "scott_kurt" ~ "KDE",
      Method == "Empirical_kurt" ~ "Empirical",
      TRUE ~ Method
    )
  ) %>%
  filter(Method %in% c("Empirical", "New", "Lmom", "Moors", "Crow", "KDE"))

# Method 순서 고정
kurt_data$Method <- factor(kurt_data$Method, levels = c("Empirical", "New", "Lmom", "Moors", "Crow", "KDE"))

color_mapping <- c(
  "Empirical" = "blue",
  "New" = "black",
  "Lmom" = "green",
  "Moors" = "gray45",
  "Crow" = "gray70",
  "KDE" = "red"
)

linetype_mapping <- c(
  "Empirical" = "solid",
  "New" = "solid",
  "Lmom" = "solid",
  "Moors" = "dashed",
  "Crow" = "dotted",
  "KDE" = "solid"
)

ggplot(kurt_data, aes(x = Date, y = Kurtosis, color = Method, linetype = Method)) +
  geom_line(size = 0.5) +
  scale_color_manual(values = color_mapping) +
  scale_linetype_manual(values = linetype_mapping) +
  scale_x_date(
    date_breaks = "6 month",
    date_labels = "%Y.%m",
    limits = c(as.Date("2008-05-27") , as.Date("2023-12-29") ),
    expand = c(0, 0)
  ) +
  labs(title = "",
       x = "", y = "Kurtosis") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    text = element_text(size = 14)
  )
