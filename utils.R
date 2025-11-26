#Script for functions
NeymanSRE <- function(W, Y, strata.labels) {
  n <- length(W)
  groups <- unique(strata.labels)
  ests <- sapply(groups, function(gg) {
    wt <- sum(strata.labels == gg)/n
    est <- mean(Y[strata.labels == gg & W == 1]) - mean(Y[strata.labels == gg & W == 0])
    est.var <- var(Y[strata.labels == gg & W == 1])/sum(strata.labels == gg & W == 1) + 
      var(Y[strata.labels == gg & W == 0])/sum(strata.labels == gg & W == 0)
    return(c(wt, est, est.var))
  })
  # print(ests)
  neyman.est <- sum(ests[1, ] * ests[2, ])
  neyman.var <- sum(ests[1, ]^2 * ests[3, ])
  return(c(est = neyman.est, se = sqrt(neyman.var), 
           tstat = neyman.est/sqrt(neyman.var)))
}

IPW_estimator <- function(W, Y, X) {
  ## Estimate propensity score
  model <- glm(W ~ X , family = "binomial")
  eps <- predict(model, type = "response")
  
  ## Calculate the weights
  weights <- ifelse(W == 1, 1/eps, 1/(1 - eps))
  
  ## Calculate weighted mean difference between treated and control group
  est <- lm(Y ~ W, weights = weights)$coef[2]
  return(est)
}

IPW_bootstrap <- function(W, Y, X, n.boot = 200){
  est <- IPW_estimator(W, Y, X)
  IPWboot <- sapply(1:n.boot, function(i) {
    id.boot <- sample(1:length(W), replace = T)
    IPW_estimator(W[id.boot], Y[id.boot], X[id.boot, ])
  })
  return(c(est, sd(IPWboot)))
}


## Draw love plot
love.plot = function(cov, treat,  ## cov is the matrix of covariates and treat is a vector of treatment assignment
                     weights = rep(1, length(treat)),
                     plot = F) 
{
  
  ## mean with normalized weights \sum w_i x_i / (\sum w_i)
  treat.means <- colSums(cov[treat == 1,] * weights[treat == 1])/sum(weights[treat == 1])
  treat.var <- colSums(t(t(cov[treat == 1,]) - treat.means)^2 *
                         weights[treat == 1])/sum(weights[treat == 1])
  
  control.means <- colSums(cov[treat == 0,] * weights[treat == 0])/sum(weights[treat == 0])
  control.var <- colSums(t(t(cov[treat == 0,]) - control.means)^2 *
                           weights[treat == 0])/sum(weights[treat == 0])
  
  ## the standardized mean differences for every covariate
  smd <- (treat.means - control.means)/sqrt((treat.var + control.var)/2)
  names(smd) <- colnames(cov)
  
  if (plot == T) {
    plot.data <- data.frame(smd = smd, covariates = names(smd))
    range <- max(abs(smd))
    ggplot(plot.data) + geom_point(aes(x = as.numeric(smd), y = covariates)) +
      geom_vline(xintercept = 0) + xlim(-range, range) +
      labs(x = 'Standardized Difference in Means')
  }
  return(smd)
}