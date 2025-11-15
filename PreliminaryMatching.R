
# Following Example 8

library(MatchIt)
library(tidyverse)
library(ggplot2)

fish <- read_csv("nhanes_fish.csv")

fish_clean <- fish %>%
  mutate(high_fish = ifelse(fish.level == "high", 1, 0)) %>%
  mutate(race = as.factor(race)) %>%
  select(high_fish, gender, age, income, income.missing, race, education,
         smoking.ever, smoking.now)

# Propensity Score Calculation
model <- glm(high_fish ~ ., data = fish_clean, family = "binomial")
eps <- predict(model, type = "response")
lps <- predict(model)

temp.data <- data.frame(lps = lps, treated = as.factor(fish_clean$high_fish))
ggplot(temp.data, aes(x = lps, fill = treated, color = treated)) + 
  geom_histogram(alpha = 0.5, position = "identity") + 
  xlab("Linearized propensity score") 



# Original Covariate Balancing

m.out0 <- matchit(high_fish ~ gender + age + income + income.missing + race +
                    education + smoking.ever + smoking.now, data = fish_clean,
                  method = NULL, distance = "glm")

summary(m.out0)
plot0 = plot(summary(m.out0), abs = F)
plot0

# Matched Covariate Balancing
m.out1 <- matchit(high_fish ~ gender + age + income + income.missing + race +
                    education + smoking.ever + smoking.now, data = fish_clean,
                  estimand = "ATT", 
                  method = "nearest", distance = lps)
summary(m.out1)
plot1 = plot(summary(m.out1), abs = F)
plot1

matched_fish <- match_data(m.out1)

matched_model <- glm(high_fish ~ ., data = matched_fish, family = "binomial")
matched_eps <- predict(matched_model, type = "response")
matched_lps <- predict(matched_model)
matched_temp.data <- data.frame(lps = matched_lps, treated = as.factor(matched_fish$high_fish))

ggplot(matched_temp.data, aes(x = matched_lps, fill = treated, color = treated)) + 
  geom_histogram(alpha = 0.5, position = "identity") + 
  xlab("Linearized propensity score") 

