library(MatchIt)
library(tidyverse)
library(ggplot2)

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
pscore <- model$fitted.values

# EPS histogram (before trimming)
temp.data <- data.frame(eps = pscore, treated = as.factor(fish_clean$high_fish))
ggplot(temp.data, 
       aes(x = eps, fill = treated, color = treated)) + 
  geom_histogram(alpha = 0.5, position = "identity") + 
  xlab("Estimated propensity score") 

# ------ To be consistent with IPW, we trim in a similar fashion prior to stratify
rm.idx <- which(pscore < 0.05 | pscore > 0.95)

fish_clean.trimmed <- fish_clean[-rm.idx, ]
pscore.trim <- pscore[-rm.idx]


## Refitting the propensity score model
model.trimmed <- glm(high_fish ~ . , data = fish_clean.trimmed[-1], family = "binomial")
eps.trimmed <- predict(model.trimmed, type = "response")
lps.trimmed <- predict(model.trimmed)

# Histogram of eps after trimming
temp.data.eps.trimmed <- data.frame(eps = eps.trimmed, treated = as.factor(fish_clean.trimmed$high_fish))
ggplot(temp.data.eps.trimmed, aes(x = eps, fill = treated, color = treated)) + 
  geom_histogram(alpha = 0.5, position = "identity") +
  ggtitle("Histogram of eps after trimming") +
  theme_bw()

# ---- end of changed code

# Stratifying after trimming procedure (quantiles)
x <- model.matrix(model.trimmed)[, -1]
n <- table(fish_clean.trimmed$high_fish)
z <- fish_clean.trimmed$high_fish

nn <- 5
q.pscore <- quantile(pscore.trim, (1:( nn -1)) / nn)
ps.strata <- cut(pscore.trim, breaks = c(0 , q.pscore, 1), labels = 1:nn)
balance_check <- apply(x, 2, function(v) NeymanSRE(fish_clean.trimmed$high_fish, v, ps.strata)) 

temp <- data.frame(x = colnames(x), y = balance_check[3, ])
ggplot(temp, aes(x, y)) + geom_point() + geom_hline(yintercept = 0) +
  geom_hline(yintercept = -1.96, color = "red")+ 
  geom_hline(yintercept = 1.96, color = "red") +
  xlab("") + ylab("t statistics") + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

# Looking at groups that were made
tble <- table(ps.strata, fish_clean.trimmed$high_fish)
colnames(tble) <- c("Control", "Treated")

tble


#Sequential splitting
## First, create a data frame to store the current grouping information
temp <- data.frame(e = lps.trimmed, treat = fish_clean.trimmed$high_fish, b = 1)

t.max <-  1.96
# check whether t. stat is above t.max for first iteration
t.stat <- t.test(x = temp$e[temp$treat == 1], y =  temp$e[temp$treat == 0], var.equal=T)$statistic
condition = t.stat > t.max
# minimum n.t or n.c in each block
size <- 3
size.new <- 20

# continue until all t statistics are below t.max
set.seed(0)
while(condition)
{
  # calculate size in each group
  # we want to not split a block if it is too small
  b <- max(temp$b)
  ignore <- sapply(1:b, function(j) {
    n.t <- sum(temp$treat == 1 & temp$b == j)
    n.c <- sum(temp$treat == 0 & temp$b == j)
    return(n.t < size | n.c < size | (n.t + n.c < size.new * 2))
  })
  
  
  # split unbalanced blocks into more blocks
  split <- which((abs(t.stat) > t.max) & (!ignore))
  
  if(length(split) == 0)
    break
  
  
  ## we need to keep a current copy of the block information and which block to ignore as block assignments are going to change later
  b.current <- temp$b
  
  for(j in split)
  {
    
    cutoff <- median(temp$e[b.current == j])
    ## We split units into two new blocks
    ## extract the index of units belonging to each new stratum
    idx.s <- which(b.current == j & temp$e < cutoff)
    idx.l <- which(b.current == j & temp$e > cutoff)
    ## randomly put half of the ties into one category
    idx.e <- which(temp$e == cutoff & b.current == j)
    n.tie <- length(idx.e)
    if (n.tie >= 1) {
      if (n.tie > 1) {
        idx.e <- sample(idx.e)
        idx.s <- c(idx.s, idx.e[1:round(n.tie/2)])
      }
      idx.l <- c(idx.l, idx.e[(round(n.tie/2)+ 1):n.tie])
    }
    
    ## we split only when new stratum has at least size number of control/treated units
    if (sum(temp$treat[idx.s]==1) > size && sum(temp$treat[idx.s]==0) > size && 
        sum(temp$treat[idx.l]==1) > size && sum(temp$treat[idx.l]==0) > size) {
      # anything above the current will have to be moved up 1
      temp$b[b.current > j] <- temp$b[b.current > j] + 1
      ## the upper new stratum will also have the block idx added by 1
      temp$b[idx.l] <- temp$b[idx.l] + 1
    }
    ## We don't do anything if we do not want to split
  }
  
  # calculate t statistic for each block
  b <- max(temp$b)
  t.stat <- sapply(1:b, function(j) {
    t.test(x = temp$e[temp$treat == 1 & temp$b == j], 
           y = temp$e[temp$treat == 0 & temp$b == j], var.equal=T)$statistic
  })
  
  ## Update condition
  # check whether ANY blocks are above t.max
  # AND are not too small
  condition <- any(abs(t.stat) > t.max)
}

fish_clean.trimmed$blocks <- temp$b
print("number of individuals per strata")

table(fish_clean.trimmed[, c("high_fish", "blocks")])

## check the range of estimated propensity scores
lps_blocks <- sapply(1:max(fish_clean.trimmed$blocks), 
                     function(j) range(temp$e[temp$b == j]))
eps_blocks <- exp(lps_blocks)/(1 + exp(lps_blocks))
colnames(eps_blocks) <- paste("strata", 1:8)
rownames(eps_blocks) <- c("min eps", "max eps")
eps_blocks



# Balance check post sequential splitting (can see it is slightly improved)
balance_check <- apply(x, 2, function(v) NeymanSRE(fish_clean.trimmed$high_fish, v, fish_clean.trimmed$blocks)) 


temp <- data.frame(x = colnames(x), y = balance_check[3, ])
ggplot(temp, aes(x, y)) + geom_point() + geom_hline(yintercept = 0) +
  geom_hline(yintercept = -1.96, color = "red")+ 
  geom_hline(yintercept = 1.96, color = "red") +
  xlab("") + ylab("t statistics") + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) 


#Neyman's ATE Estimate (data as Stratified Random Experiment)
result <- NeymanSRE(fish_clean.trimmed$high_fish, fish_clean.trimmed$mercury, fish_clean.trimmed$blocks)
result <- c(result[1:2], c(result[1] - 1.96 * result[2], result[1] + 1.96 * result[2]))
names(result) <- c("est", "sd", "CI_lower", "CI_upper")
result

## Histograms of covariate distributions before and after sequential splitting
## Note that original graphs are contained in PreliminaryMatching.R

# creating a new df with needed labels
hist_dat <- fish_clean.trimmed[, c("age", "income", "high_fish", "blocks")]

hist_dat$high_fish <- factor(
  hist_dat$high_fish,
  levels = c(0, 1),
  labels = c("Control", "Treated")
)

# age blocks
ggplot(hist_dat, aes(x = age, fill = high_fish, color = high_fish)) +
  geom_histogram(position = "identity", alpha = 0.6) +
  facet_wrap(~ blocks) +
  labs(title = "Age Distributions by Block After Sequential Splitting",
       x = "Age", y = "Count") +
  theme_bw()

# income blocks
ggplot(hist_dat, aes(x = income, fill = high_fish, color = high_fish)) +
  geom_histogram(position = "identity", alpha = 0.6) +
  facet_wrap(~ blocks) +
  labs(title = "Income Distributions by Block After Sequential Splitting",
       x = "Income", y = "Count") +
  theme_bw()