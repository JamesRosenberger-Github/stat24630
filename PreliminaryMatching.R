# Following Example 8

library("MatchIt")
library(tidyverse)

fish <- read_csv("nhanes_fish.csv")

fish_clean <- fish %>%
  mutate(high_fish = ifelse(fish.level == "high", 1, 0)) %>%
  mutate(race = as.factor(race)) %>%
  select(high_fish, gender, age, income, income.missing, race, education,
         smoking.ever, smoking.now)

model <- glm(high_fish ~ ., data = fish_clean, family = "binomial")
eps <- predict(model, type = "response")
lps <- predict(model)


m.out0 <- matchit(high_fish ~ gender + age + income + income.missing + race +
                    education + smoking.ever + smoking.now, data = fish_clean,
                  method = NULL, distance = "glm")

summary(m.out0)
plot(summary(m.out0), abs = F)


m.out1 <- matchit(high_fish ~ gender + age + income + income.missing + race +
                    education + smoking.ever + smoking.now, data = fish_clean,
                  estimand = "ATT", 
                  method = "nearest", distance = lps)
summary(m.out1)

