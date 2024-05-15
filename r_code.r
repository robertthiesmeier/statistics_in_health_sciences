##########################################################################
# Example I: Theoretical vs. empirical distribution
##########################################################################

rm(list = ls())

# function to simulate study
sim_exp1 <- function() {
  n <- 1000
  py_1 <- 0.3
  py_0 <- 0.1
  x_p <- 0.3
  
  x <- rbinom(n, 1, x_p)
  y <- ifelse(x == 1, rbinom(n, 1, py_1), rbinom(n, 1, py_0))
  
  df <- data.frame(y = y, x = x)
  model1 <- glm(y ~ x, data = df, family = binomial(link = "logit"))
  
  est_b1 <- coef(model1)[2]
  vcov <- vcov(model1)
  est_se_b1 <- sqrt(vcov[2, 2])
  return(c(est_b1 = est_b1, est_se_b1 = est_se_b1))
}

# simulate the once
sim_exp1()

# simulate the study many times
reps <- 10000
est_b1 <- rep(NA, reps)
est_se_b1 <- rep(NA, reps)

for (i in seq_len(reps)) {
  cat(".")
  if (i %% 10 == 0 || i == reps) {
    cat(i)
  }
  if (i %% 50 == 0 && i != reps) {
    cat("\n")
  }
  est <- sim_exp1()
  est_b1[i] <- est[1]
  est_se_b1[i] <- est[2]
}

# summary statistics
summary(est_b1)
exp(mean(est_b1))
summary(est_se_b1)

# Expected mean and std deviation
exp_mean <- (60/140)/(80/720)
exp_sd <- sqrt(1/60+1/140+1/80+1/720)
exp_mean
exp_sd

##########################################################################
# Example II: Confounding ##########################################################################

rm(list = ls())
library(broom)

# define function for logit and inverse of logit (invlogit)
logit <- function(p) {log(p/(1 - p))}
invlogit <- function(p){1/(1+exp(-p))}

# define the parameters of the study
n <- 1000
p_y <- 0.1
p_x_z0 <- 0.3
p_w <- 0.5
p_z <- 0.4

# w is an outcome predictor but not a confounder 
w <- rbinom(n, 1, p_w)

# z is a confounding variable (common cause of exposure and outcome)
z <- rbinom(n, 1, p_z)

# x is the main predictor 
x <- rbinom(n, 1, invlogit(logit(p_x_z0)+log(3)*z))

# y is the outcome
y <- rbinom(n, 1, invlogit(logit(p_y)+log(1.5)*x+log(2)*w+log(6)*z))

df <- data.frame( y = y, x = x, z = z, w = w)

# Model 0. Omit w and z 
model0 <- glm(y ~ x, data = df, family = binomial(link = "logit"))
sum_model0 <- tidy(model0)
print(sum_model0)

# Model 1. Include an outcome predictor w but not the confounder z 
model1 <- glm(y ~ x + w, data = df, family = binomial(link = "logit"))
sum_model1 <- tidy(model1)
print(sum_model1)

# Model 2. Include only the confounder z but not the outcome predictor w
model2 <- glm(y ~ x + z, data = df, family = binomial(link = "logit"))
sum_model2 <- tidy(model2)
print(sum_model2)

# Model 3. Include both the confounder z and the outcome predictor w
model3 <- glm(y ~ x + z + w, data = df, family = binomial(link = "logit"))
sum_model3 <- tidy(model3)
print(sum_model3)

# simulating the unadjusted and adjusted effect for z a many times
sim_exp2 <- function(){
  n <- 1000
  p_y <- 0.1
  p_x_z0 <- 0.3
  p_w <- 0.5
  p_z <- 0.4
  
  w <- rbinom(n, 1, p_w)
  z <- rbinom(n, 1, p_z)
  x <- rbinom(n, 1, invlogit(logit(p_x_z0)+log(3)*z))
  y <- rbinom(n, 1, invlogit(logit(p_y)+log(1.5)*x+log(2)*w+log(6)*z))
  
  df <- data.frame( y = y, x = x, z = z, w = w)
  
  model0 <- glm(y ~ x, data = df, family = binomial(link = "logit"))
    est_b1_model0 <- coef(model0)[2]
    vcov <- vcov(model0)
    est_se_b1_model0 <- sqrt(vcov[2, 2])
    
  model1 <- glm(y ~ x + w, data = df, family = binomial(link = "logit"))
    est_b1_model1 <- coef(model1)[2]
    vcov <- vcov(model1)
    est_se_b1_model1 <- sqrt(vcov[2, 2])
    
  model2 <- glm(y ~ x + z, data = df, family = binomial(link = "logit"))
    est_b1_model2 <- coef(model2)[2]
    vcov <- vcov(model2)
    est_se_b1_model2 <- sqrt(vcov[2, 2])
    
  model3 <- glm(y ~ x + z + w, data = df, family = binomial(link = "logit"))
    est_b1_model3 <- coef(model3)[2]
    vcov <- vcov(model3)
    est_se_b1_model3 <- sqrt(vcov[2, 2])
    p_value_model3 <- coef(summary(model3))["x", "Pr(>|z|)"]
    
    return(c(est_b1_model0 = est_b1_model0, est_se_b1_model0 = est_se_b1_model0,
             est_b1_model1 = est_b1_model1, est_se_b1_model1 = est_se_b1_model1,
             est_b1_model2 = est_b1_model2, est_se_b1_model2 = est_se_b1_model2,
             est_b1_model3 = est_b1_model3, est_se_b1_model3 = est_se_b1_model3,
             p_value_model3 = p_value_model3))
}

# simulate the study once
sim_exp2()

# simulate the study many times
reps <- 1000
est_b1_model0 <- rep(NA, reps)
est_se_b1_model0 <- rep(NA, reps)
est_b1_model1 <- rep(NA, reps)
est_se_b1_model1 <- rep(NA, reps)
est_b1_model2 <- rep(NA, reps)
est_se_b1_model2 <- rep(NA, reps)
est_b1_model3 <- rep(NA, reps)
est_se_b1_model3 <- rep(NA, reps)
p_value_model3 <- rep(NA, reps)

for (i in seq_len(reps)) {
  cat(".")
      if (i %% 10 == 0 || i == reps) {
           cat(i)
      }
      if (i %% 50 == 0 && i != reps) {
          cat("\n")
      }
  est <- sim_exp2()
  est_b1_model0[i] <- est[1]
  est_se_b1_model0[i] <- est[2]
  est_b1_model1[i] <- est[3]
  est_se_b1_model1[i] <- est[4]
  est_b1_model2[i] <- est[5]
  est_se_b1_model2[i] <- est[6]
  est_b1_model3[i] <- est[7]
  est_se_b1_model3[i] <- est[8]
  p_value_model3[i] <- est[9]
}

# summary statistics
summary(est_b1_model0)
exp(mean(est_b1_model0))
summary(est_b1_model1)
exp(mean(est_b1_model1))
summary(est_b1_model2)
exp(mean(est_b1_model2))
summary(est_b1_model3)
exp(mean(est_b1_model3))

# power
power <- sum(p_value_model3 < 0.05)
perc_power <- power / reps * 100
result_power <- paste("The simulated power is =", round(perc_power, 2), "%")
print(result_power)

##########################################################################
# Example II: Missing data ##########################################################################

rm(list = ls())
library(broom)
library(mice)

# define function for logit and inverse of logit (invlogit)
logit <- function(p) {log(p/(1 - p))}
invlogit <- function(p){1/(1+exp(-p))}

# define the parameters
n <- 1000
p_y <-  0.1
p_x <-  0.3
p_z <-  0.4

z <- rbinom(n, 1, p_z)
x <- rbinom(n, 1, invlogit(logit(p_x)+log(3)*z))
y <- rbinom(n, 1, invlogit(logit(p_y)+log(1.5)*x+log(6)*z))

df <- data.frame(y = y, x = x, z = z)

# model without any missing data
model0 <- glm(y ~ x + z, data = df, family = binomial(link = "logit"))
sum_model0 <- tidy(model0)
print(sum_model0)
estimates_full <- summary(model0)

# MCAR for the confounder 70% 
df$z[runif(nrow(df)) < 0.7] <- NA

# complete case analysis
model1 <- glm(y ~ x + z, data = df, na.action = na.exclude, family = binomial(link = "logit"))
sum_model1 <- tidy(model1)
print(sum_model1)
estimates_cc <- summary(model1)

# imputing missing data
df$z <- factor(df$z)
imputed_data <- mice(df, m = 70, maxit = 5, method = "logreg")
model2 <- with(imputed_data, glm(y ~ x + z, family = binomial(link = "logit")))
summary(pool(model2)) 
estimates_mi <- summary(model2)
estimates_full

# simulate the example many times
sim_exp3 <- function() {
  n <- 1000
  p_y <- 0.1
  p_x <- 0.3
  p_z <- 0.4
  
  z <- rbinom(n, 1, p_z)
  x <- rbinom(n, 1, invlogit(logit(p_x) + log(3) * z))
  y <- rbinom(n, 1, invlogit(logit(p_y) + log(1.5) * x + log(6) * z))
  
  df <- data.frame(y = y, x = x, z = z)
  
  # Model 0
  model0 <- glm(y ~ x + z, data = df, family = binomial(link = "logit"))
  est_b1_model0 <- coef(model0)[2]
  est_se_b1_model0 <- sqrt(vcov(model0)[2, 2])
  
  # Introduce missingness in z
  df$z[runif(nrow(df)) < 0.7] <- NA
  
  # Model 1
  model1 <- glm(y ~ x + z, data = df, na.action = na.exclude, family = binomial(link = "logit"))
  est_b1_model1 <- coef(model1)[2]
  est_se_b1_model1 <- sqrt(vcov(model1)[2, 2])
  
  # Imputation with MICE
  df$z <- factor(df$z) # z is recognized as categorical variable
  imputed_data <- mice(df, m = 10, maxit = 5, method = "logreg", printFlag = FALSE)
  model2 <- with(imputed_data, glm(y ~ x + z, family = binomial(link = "logit")))
  
  # Combine results from multiple imputations
  est_b1_model2 <- summary(pool(model2))$estimate[2]
  est_se_b1_model2 <- summary(pool(model2))$std.error[2]
  
  return(c(
    est_b1_model0 = est_b1_model0,
    est_se_b1_model0 = est_se_b1_model0,
    est_b1_model1 = est_b1_model1,
    est_se_b1_model1 = est_se_b1_model1,
    est_b1_model2 = est_b1_model2,
    est_se_b1_model2 = est_se_b1_model2
  ))
}

# simulate the study once
sim_exp3()

# simulate the study many times
reps <- 1000
est_b1_model0 <- rep(NA, reps)
est_se_b1_model0 <- rep(NA, reps)
est_b1_model1 <- rep(NA, reps)
est_se_b1_model1 <- rep(NA, reps)
est_b1_model2 <- rep(NA, reps)
est_se_b1_model2 <- rep(NA, reps)

for (i in seq_len(reps)) {
  cat(".")
  if (i %% 10 == 0 || i == reps) {
    cat(i)
  }
  if (i %% 50 == 0 && i != reps) {
    cat("\n")
  }
  est <- sim_exp3()
  est_b1_model0[i] <- est[1]
  est_se_b1_model0[i] <- est[2]
  est_b1_model1[i] <- est[3]
  est_se_b1_model1[i] <- est[4]
  est_b1_model2[i] <- est[5]
  est_se_b1_model2[i] <- est[6]
}

# summary statistics
summary(est_b1_model0)
exp(mean(est_b1_model0)) # no missing data
summary(est_b1_model1)
exp(mean(est_b1_model1)) # complete case
summary(est_b1_model2)
exp(mean(est_b1_model2
