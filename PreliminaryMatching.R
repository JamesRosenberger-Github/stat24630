# Following Example 8

library(MatchIt)
library(tidyverse)
library(ggplot2)

fish <- read_csv("nhanes_fish.csv")

# Cleaning data
fish_clean <- fish %>%
  mutate(high_fish = ifelse(fish.level == "high", 1, 0)) %>%
  mutate(race = as.factor(race)) %>%
  rename(mercury = o.LBXTHG) %>%
  select(mercury, high_fish, gender, age, income, income.missing, race, education,
         smoking.ever, smoking.now)

# Propensity Score Calculation
model <- glm(high_fish ~ ., data = fish_clean[-1], family = "binomial")
eps <- predict(model, type = "response")
lps <- predict(model)

# LPS Histogram
temp.data <- data.frame(lps = lps, treated = as.factor(fish_clean$high_fish))
ggplot(temp.data, aes(x = lps, fill = treated, color = treated)) + 
  geom_histogram(alpha = 0.5, position = "identity") + 
  xlab("Linearized propensity score")

# EPS Histogram
temp.data_eps <- data.frame(eps = eps, treated = as.factor(fish_clean$high_fish))
ggplot(temp.data_eps, aes(x = eps, fill = treated, color = treated)) + 
  geom_histogram(alpha = 0.5, position = "identity") + 
  xlab("Linearized propensity score")


# Covariate Balancing Before matching
m.out0 <- matchit(high_fish ~ gender + age + income + income.missing + race +
                    education + smoking.ever + smoking.now, data = fish_clean,
                  method = NULL, distance = "glm")

summary(m.out0)
plot0 = plot(summary(m.out0), abs = F)
plot0

# Covariate Balancing after Greedy Matching
m.out1 <- matchit(high_fish ~ gender + age + income + income.missing + race +
                    education + smoking.ever + smoking.now, data = fish_clean,
                  estimand = "ATT", 
                  method = "nearest", distance = lps)
summary(m.out1)
plot1 = plot(summary(m.out1), abs = F)
plot1

matched_fish_g <- match_data(m.out1)


#Optimal Matched Covariate Balancing (in this case the result is equivalent to greedy)
m.out2 <- matchit(high_fish ~ gender + age + income + income.missing + race +
                    education + smoking.ever + smoking.now, data = fish_clean,
                  estimand = "ATT", 
                  method = "optimal", distance = lps)
summary(m.out2)
plot2 = plot(summary(m.out2), abs = F)
plot2

matched_fish_opt <- match_data(m.out2)




#regression adjsutment on matching

