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
  xlab("Estimated propensity score")

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

# Optimal Matched Covariate Balancing (in this case the result is 
# functionally equivalent to greedy)
m.out2 <- matchit(high_fish ~ gender + age + income + income.missing + race +
                    education + smoking.ever + smoking.now, data = fish_clean,
                  estimand = "ATT", 
                  method = "optimal", distance = lps)
summary(m.out2)
plot2 = plot(summary(m.out2), abs = F)
plot2

matched_fish_opt <- match_data(m.out2)





#Neyman's Approach without Regression Adjustment (Pairwise Random Exp.)
m.data <- match.data(m.out2)

tau_hat_vec <- sapply(levels(m.data$subclass), function(sc) {
  treated <- m.data$mercury[m.data$subclass == sc & m.data$high_fish == 1]
  control <- m.data$mercury[m.data$subclass == sc & m.data$high_fish == 0]
  return(treated - control)
})

tau_hat <- mean(tau_hat_vec, na.rm = T)
sd_tau_hat <- sd(tau_hat_vec)/sqrt(length(tau_hat_vec))
print(c(tau_hat, sd_tau_hat))

## 95% CI
print(c(tau_hat- 1.96 * sd_tau_hat, tau_hat + 1.96 * sd_tau_hat))


#Neyman's Approach with Regression Adjustment

## First run a linear regression on the control (remember to remove the treatment assignment variable!)

fit0 <- lm(mercury ~ gender + age + income +
             income.missing + race + education + 
             smoking.ever + smoking.now, 
           data = m.data[m.data$high_fish == 0, -c(2)])

m.data$predicted_y0 <- predict(fit0, m.data)

## Neyman's approach with regression adjustment
## First run a linear regression on the control (remember to remove the treatment assignment variable!)

tau_hat_vec_adjusted <- sapply(levels(m.data$subclass), function(sc) {
  treated <- m.data$mercury[m.data$subclass == sc & m.data$high_fish == 1]
  control <- m.data$mercury[m.data$subclass == sc & m.data$high_fish == 0]
  
  bias <- m.data$predicted_y0[m.data$subclass == sc & m.data$high_fish == 1] - 
    m.data$predicted_y0[m.data$subclass == sc & m.data$high_fish == 0]
  control.adjusted <- control + bias
  
  return(treated - control.adjusted)
})

tau_hat_adjusted <- mean(tau_hat_vec_adjusted)
sd_tau_hat_adjusted <- sd(tau_hat_vec_adjusted)/sqrt(length(tau_hat_vec_adjusted))
print(c(tau_hat_adjusted, sd_tau_hat_adjusted))

## 95% CI
print(c(tau_hat_adjusted - 1.96 * sd_tau_hat_adjusted, tau_hat_adjusted + 1.96 * sd_tau_hat_adjusted))


## Histograms of covariate distributions before and after matching

# Notice that matching greatly improves overlap
# creating a new column with appropriate labels
help <- function(df, label) {
  df2 <- df[, c("age", "income", "high_fish")]
  df2$sample <- label
  df2
}

# defining a new larger object for graphing
hist_dat <- rbind(help(fish_clean, "Original"), help(m.data, "Matched"))
hist_dat$sample <- factor(hist_dat$sample, levels = c("Original", "Matched"))
hist_dat$high_fish  <- factor(hist_dat$high_fish,
                              levels = c(0, 1),
                              labels = c("Control", "Treated"))

ggplot(hist_dat, aes(x = age, fill = high_fish)) +
  geom_histogram(position = "identity", alpha = 0.6) +
  facet_grid(sample ~ .) +
  labs(title = "Age Distributions by Different Procedures (Matching)",
       x = "age", y = "Count")

ggplot(hist_dat, aes(x = income, fill = high_fish)) +
  geom_histogram(position = "identity", alpha = 0.6) +
  facet_grid(sample ~ .) +
  labs(title = "Income Distributions by Different Procedures (Matching)",
       x = "Income", y = "Count")
