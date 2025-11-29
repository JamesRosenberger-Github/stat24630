library(MatchIt)
library(tidyverse)
library(ggplot2)
library(sandwich)

fish <- read_csv("nhanes_fish.csv")
source("utils.R")

# Cleaning data
fish_clean <- fish %>%
  mutate(high_fish = ifelse(fish.level == "high", 1, 0)) %>%
  mutate(race = as.factor(race)) %>%
  rename(mercury = o.LBXTHG) %>%
  select(mercury, high_fish, gender, age, income, income.missing, race, education,
         smoking.ever, smoking.now)

model <- glm(high_fish ~ . , data = fish_clean[-1], family = "binomial")
eps <- predict(model, type = "response")

n.treated <- sum(fish_clean$high_fish == 1)
n.control <- sum(fish_clean$high_fish == 0)
weights <- ifelse(fish_clean$high_fish == 1, 1/eps, 1/(1 - eps))

#Histogram of weights, log transformation on x-axis
temp.data <- data.frame(weights = weights, treated = as.factor(fish_clean$high_fish))
ggplot(temp.data, aes(x = weights, fill = treated, color = treated)) + 
  geom_histogram(alpha = 0.5, position = "identity") + 
  scale_x_log10() +
  xlab("Weights") +
  theme_bw()

#Histogram of eps
temp.data.eps <- data.frame(eps = eps, treated = as.factor(fish_clean$high_fish))
ggplot(temp.data.eps, aes(x = eps, fill = treated, color = treated)) + 
  geom_histogram(alpha = 0.5, position = "identity") +
  ggtitle("Histogram of eps before trimming") +
  theme_bw()

#chancing race into indicator variables
#note that there is no race equal to 5 because of NHANES recoding
fish_indic <- fish_clean %>%
  mutate(mexican_american = ifelse(race == 1, 1, 0),
         other_hispanic = ifelse(race == 2, 1, 0),
         white = ifelse(race == 3, 1, 0),
         black = ifelse(race == 4, 1, 0),
         asian = ifelse(race == 6, 1, 0),
         other = ifelse(race == 7, 1, 0)) %>%
  select(-race)

#Love plot before trimming
raw.smd <- love.plot(fish_indic[, -c(1:2)], fish_indic$high_fish)
weighted.smd <- love.plot(fish_indic[, -c(1:2)], fish_indic$high_fish, weights = weights)

plot.data <- data.frame(smd = c(raw.smd, weighted.smd), 
                        covariates = c(names(raw.smd), names(weighted.smd)),
                        category = c(rep("Original", length(raw.smd)), rep("IPW", length(weighted.smd))))
range <- max(abs(plot.data$smd))

ggplot(plot.data) + geom_point(aes(x = as.numeric(smd), y = covariates, color = category)) +
  geom_vline(xintercept = c(-0.1, -0.05, 0, 0.05, 0.1),
             linetype = c("solid", "dashed", "solid", "dashed", "solid")) + 
  xlim(-range, range) +
  labs(x = 'Standardized Difference in Means') +
  theme_bw()



###-Trimming-----------------------------------(always makes it worse)
rm.idx <- which(eps < 0.05 | eps > 0.95)
temp.data.trimmed <- temp.data[-rm.idx, ]
fish_indic.trimmed <- fish_indic[-rm.idx, ]

# temp.data.trimmed <- temp.data
# fish_indic.trimmed <- fish_indic

## Refitting the propensity score model
model.trimmed <- glm(high_fish ~ . , data = fish_indic.trimmed[-1], family = "binomial")
eps.trimmed <- predict(model.trimmed, type = "response")
weights.trimmed <- ifelse(fish_indic.trimmed$high_fish == 1, 1/eps.trimmed, 1/(1 - eps.trimmed))
temp.data.trimmed$weights <- weights.trimmed

#love plot after trimming is worse, so definitely stick with untrimmed
raw.smd <- love.plot(fish_indic.trimmed[, -c(1:2)], fish_indic.trimmed$high_fish)
weighted.smd <- love.plot(fish_indic.trimmed[, -c(1:2)], fish_indic.trimmed$high_fish, weights = weights.trimmed)

plot.data <- data.frame(smd = c(raw.smd, weighted.smd), 
                        covariates = c(names(raw.smd), names(weighted.smd)),
                        category = c(rep("Original", length(raw.smd)), rep("IPW", length(weighted.smd))))
range <- max(abs(plot.data$smd))

ggplot(plot.data) + geom_point(aes(x = as.numeric(smd), y = covariates, color = category)) +
  geom_vline(xintercept = c(-0.1, -0.05, 0, 0.05, 0.1),
             linetype = c("solid", "dashed", "solid", "dashed", "solid")) + 
  xlim(-range, range) +
  labs(x = 'Standardized Difference in Means') +
  theme_bw()


#------Because trimming made it worse, we estimate IPW using the untrimmed data

lm.result <- lm(mercury ~ high_fish, weights = weights, data = fish_indic)
summary(lm.result)

#sandwich variance estimate
tau_hat <- lm.result$coefficients[2]
SE <- sqrt(diag(vcovHC(lm.result, type = "HC2")))[2]

## get the 95% CI
result <- c(tau_hat, SE, c(tau_hat- 1.96 * SE, tau_hat + 1.96 * SE))
names(result) <- c("est", "sd", "CI_lower", "CI_upper")
result


# using bootstrap
X <- model.matrix(model)
Y <- fish_indic$mercury
W <- fish_indic$high_fish

SE_boostrap <- IPW_bootstrap(W, Y, X, 5000)[2]
result_bootstrap <- c(tau_hat, SE_boostrap, 
                      c(tau_hat- 1.96 * SE_boostrap, tau_hat + 1.96 * SE_boostrap))
names(result_bootstrap) <- c("est", "sd (bootstrap)", "CI_lower (bootstrap)", "CI_upper (bootstrap)")
result_bootstrap
