#' Example code for Quantifying Within- and Between-Cluster Variance in Logistic Multilevel Models
#' @author:  Sean Devine
#' @contact: seandamiandevine@gmail.com
#' @github:  https://github.com/seandamiandevine/logisticicc
#' @paper:   TODO

# Load data ---------------------------------------------------------------
# 1. For simplicity, recode effort choices so that 1 = high effort and 0 = low effort
# 2. subset the data such that we only look at participants in the control condition 

data               <- read.csv('Bogdanovetal2021/DST_data_osf2.csv')
data$effort_choice <- abs(data$effort_choice-1) 
data.Ctl           <- data[data$condition=='control', ] # focus just on control condition for now

# Logistic multilevel model for choice -------------------------------------------------------------------
# load the lme4 package (assuming the package is installed)
# if the lme4 is not installed, run this command first:
# install.packages('lme4')

library(lme4)

# Fit an intercept only model of choice
logistic_MLM0 <- glmer(effort_choice ~ 1 + (1|PID), data=data.Ctl, family='binomial')

# print summary statistics of this model
summary(logistic_MLM0)

# convert the intercept effect from log odds to probability scale and print
pChooseHighEffort <- exp(fixef(logistic_MLM0))/(1+exp(fixef(logistic_MLM0)))
cat('estimated group average probability of choosing the high-effort cue = ', pChooseHighEffort, '\n')


# # Reproduce Figure 2 from paper -----------------------------------------

# extract random intercepts from model and add them to the fixed effect (gamma00) 
eb <- coef(logistic_MLM0)$PID
head(eb)

# plot how individual differences in EB manifest in the behavioural data

# get P(High effort) from each participant
phigheff = tapply(data.Ctl$effort_choice, data.Ctl$PID, mean)

# plot correlation between EB
plot(eb$`(Intercept)`, phigheff, 
     ylim=c(0,1), xlim=c(-3,3),
     pch=21, cex=2, 
     col='black', bg='lightblue', 
     xlab=expression(gamma[`00`] + U[`0j`]), 
     ylab='P(Choose High Effort)')
abline(h=.5, lty=2)

# extract the estimated random intercept variance (tau^2_0)
tau <- VarCorr(logistic_MLM0)$PID[1]

# plot the histogram of the random intercepts and overlay an estimated normal density curve with 
# mu = gamma00 and sd = sqrt(tau^2_0)
hist(eb[[1]], main='', xlab=expression(gamma[`00`] + U[`0j`]), freq = F, xlim=c(-3,3), ylim=c(0,1))
curve(dnorm(x, fixef(logistic_MLM0), sqrt(tau)), add=T, lwd=2)

# finally, shade in areas of demand-avoidance (b0 < 0) and demand-seeking respectively (b0 > 0)
abline(v=0, lty='dashed', lwd=3)
rect(0, 0, 10, 100, col = scales::alpha('red', .1))
rect(0, 0, -10, 100, col = scales::alpha('blue', .1))
text(-.75, .95, labels = 'Demand\navoidant')
text(.75, .95, labels = 'Demand\nseeking')

# ICC: Latent variable approach -------------------------------------------

tau20     <- VarCorr(logistic_MLM0)$PID[1]
sigma2    <- pi^2/3
threshICC <- tau20/(tau20+sigma2)

# print result
cat('ICC using the Latent Threshold Approach = ',threshICC,'\n')


# ICC: Simulation method ---------------------------------------------------------------------

# 0. Extract intercept-variance and fixed effects, specify number of simulations, and set random seed

set.seed(2022) # for reproducibility 
tau20       <- VarCorr(logistic_MLM0)$PID[1]
gamma00     <- fixef(logistic_MLM0)[[1]]
M           <- 1e5

# 1. Simulate random effects 
# (note, we take the square root of tau20 here, because we want the standard dev., not the variance)
U0j         <- rnorm(M, 0, sqrt(tau20))

# 2. Predict probabilities 
logit_p_hat <- gamma00+U0j
p_hat       <- exp(logit_p_hat)/(1+exp(logit_p_hat))

# 3. Compute level-1 variance (Bernouilli variance)
var_L1      <- p_hat*(1-p_hat) 

# 4. Compute ICC
sigma2      <- mean(var_L1)
simICC      <- var(p_hat)/(var(p_hat) + sigma2)

# print result
cat('ICC using the Simulation Model = ',simICC,'\n')

# ICC: Linearization -----------------------------------------------------------

# 0. Extract relevant parameters (tau0 and gamma00)
tau20   <- VarCorr(logistic_MLM0)$PID[1]
gamma00 <- fixef(logistic_MLM0)[[1]]

# 1. Evaluate the probability of success at the mean of random effects (i.e., the fixed effect)
p       <- exp(gamma00)/(1+exp(gamma00))

# 2. Compute the Bernouilli variance of this fixed estimate first
var1    <- p*(1-p)

# 3. Compute the variance in the level-1 outcome
var2    <- tau20*p^2*(1 + exp(gamma00))^(-2)

# 4. Compute ICC
linICC  <- var2/(var1+var2)

# print result
cat('ICC using the Linearization method = ',linICC,'\n')

# Median Odds Ratio ----------------------------------------------------------------

tau20 <- VarCorr(logistic_MLM0)$PID[1]
phi75 <- qnorm(.75) # 75th percentile of normal CDF
MOR   <- exp(sqrt(2*tau20)*phi75)
 
cat('Median Odds Ratio = ',MOR,'\n')

# Bootstrapping ------------------------------------------------------

source('fx.R') # import custom functions

bootstrap_thr <- bootstrap_icc(logistic_MLM0, gr='PID', method='icc_thr', B = 100, seed = 2022)
bootstrap_sim <- bootstrap_icc(logistic_MLM0, gr='PID', method='icc_sim', B = 100, seed = 2022)
bootstrap_lin <- bootstrap_icc(logistic_MLM0, gr='PID', method='icc_lin', B = 100, seed = 2022)
bootstrap_MOR <- bootstrap_icc(logistic_MLM0, gr='PID', method='MOR', B = 100, seed = 2022)

quantile(bootstrap_thr, c(.025, .975))
quantile(bootstrap_sim, c(.025, .975))
quantile(bootstrap_lin, c(.025, .975))
quantile(bootstrap_MOR, c(.025, .975))



# Bootstrap by hand, for those who are interested.

# # 0. Set constants for bootstrapping procedure
# B       <- 1000                    # number of bootstraps
# ids     <- logistic_MLM0@frame$PID # extract id vector from model data
# K       <- length(unique(ids))     # number of clusters (subjects)
# nTrials <- table(ids)              # number of trials per subject  
# tau20   <- VarCorr(logistic_MLM0)$PID[1]
# g00     <- fixef(logistic_MLM0)[[1]]
# output  <- matrix(NA, B, 7, dimnames = list(NULL, c('iteration', 'threshICC', 'simICC', 'linICC', 'MOR', 'g00', 'tau20')))
# 
# # 1. Cycle through iterations (i) of bootstrapped samples
# for(i in 1:B) {
#   if(i%%100==0) cat('bootstrapping is', round(i/B,2)*100, '% complete.')
#   # 1.1 Simulate data
#   U0j           <- rnorm(K, 0, tau20)        # random deviations per subject
#   LOR_j         <- g00 + U0j                 # log odds of response==1 per subject
#   p_j           <- exp(LOR_j)/(1+exp(LOR_j)) # convert LOR to probability
#   y_ij          <- sapply(1:K, function(k) rbinom(nTrials[k], 1, p_j[k]))
#   y_ij          <- unlist(y_ij)              # break out of a list format
#   
#   # 1.2 Fit new model and compute values of interest
#   thisMLM       <- glmer(y_ij ~ 1 + (1|ids), family='binomial')
#   thisG00       <- fixef(thisMLM)[[1]]
#   thistau20     <- VarCorr(thisMLM)$ids[1]
#   
#   # 1.2.1. Latent Threshold ICC
#   thisthreshICC <- thistau20/(thistau20+pi^2/3)
# 
#   # 1.2.2. Simulation ICC
#   thisU0j       <- rnorm(1e5, 0, sqrt(thistau20))
#   logit_p_hat   <- thisG00+U0j
#   p_hat         <- exp(logit_p_hat)/(1+exp(logit_p_hat))
#   var_L1        <- p_hat*(1-p_hat) 
#   sigma2        <- mean(var_L1)
#   thissimICC    <- var(p_hat)/(var(p_hat) + sigma2)
#   
#   # 1.2.3. Linearlization 
#   p             <- exp(thisG00)/(1+exp(thisG00))
#   var1          <- p*(1-p)
#   var2          <- thistau20*p^2*(1 + exp(gamma00))^(-2)
#   thislinICC    <- var2/(var1+var2)
#   
#   # 1.2.4. MOR
#   thisMOR       <- exp(sqrt(2*thistau20)*qnorm(.75))
#   
#   # 1.3. Save output
#   output[i,] = c(i, thisthreshICC, thissimICC, thislinICC, thisMOR, thisG00, thistau20)
#   
# }


# ICC with one within-person predictor and random intercepts -----------------------------------

# Compute trial per block to use as a predictor
data.Ctl <- data.Ctl[order(c(data.Ctl$PID, data.Ctl$block)), ]
data.Ctl <- data.Ctl[!is.na(data.Ctl$PID),]

trial <- c()
for(i in unique(data.Ctl$PID)) {
  for(j in unique(data.Ctl$block)) {
    trial <- c(trial, 1:nrow(data.Ctl[data.Ctl$PID==i & data.Ctl$block==j, ]))
  }
}

data.Ctl$trial <- trial

# Model effort choice trial-by-trial
data.Ctl$trial_c <- data.Ctl$trial - median(data.Ctl$trial)
data.Ctl$trial0  <- data.Ctl$trial_c/max(data.Ctl$trial_c)
  
logistic_MLM1 <- glmer(effort_choice ~ trial0 + (1|PID), 
                       data=data.Ctl, family='binomial')
summary(logistic_MLM1)

icc_thr(logistic_MLM0, 'PID')
icc_thr(logistic_MLM1, 'PID')

icc_sim(logistic_MLM0, 'PID')
icc_sim(logistic_MLM1, 'PID')

icc_lin(logistic_MLM0, 'PID')
icc_lin(logistic_MLM1, 'PID')

MOR(logistic_MLM0, 'PID')
MOR(logistic_MLM1, 'PID')

# Make figure 3A
# compute the ICC at each different steps along `trial` to see how it shifts
# we will visualize this with the linearization
iccs_across_trials <- icc_lin(logistic_MLM1, 'PID', avg_preds = F)

plot(1:75, iccs_across_trials,
     pch=21, cex=2,
     col='black', bg='lightblue', 
     xlab="Trial Number", 
     ylab='ICC (Linearization Method)')
abline(h=icc_lin(logistic_MLM1, 'PID'), lty=2, col='red')

# ICC with one between-person predictor and random intercepts -----------------------------------
# here we model effort choice by condition (control and stress)

# dummy code condition
data$condition <- as.factor(data$condition)
contrasts(data$condition) <- contr.sum(2)

# Fit model
logistic_MLM2 <- glmer(effort_choice ~ condition + (1|PID), 
                       data=data, family='binomial')
summary(logistic_MLM2)

# Compute ICC for each gender
icc_thr(logistic_MLM2, 'PID') # 0.1507089
icc_sim(logistic_MLM2, 'PID', avg_preds = F, seed=2023) # 0.1127125 0.1106655
icc_lin(logistic_MLM2, 'PID', avg_preds = F) # 0.1245555 0.1189923
MOR(logistic_MLM2, 'PID') # 2.072635

# Compute averaged ICC
icc_thr(logistic_MLM2, 'PID') # 0.1507089
icc_sim(logistic_MLM2, 'PID',seed=2023) # 0.1115791
icc_lin(logistic_MLM2, 'PID') # 0.1221064
MOR(logistic_MLM2, 'PID') # 2.072635

# Make Figure 3B
barplot(icc_lin(logistic_MLM2, 'PID', avg_preds = F), 
        ylim=c(0.10,0.14), xpd=F,
        names.arg = c('Control', 'Stress'), 
        ylab='ICC (Linearization Method)', 
        xlab='Experimental Condition')
abline(h=icc_lin(logistic_MLM2, 'PID'), lty=2, col='red')


# Supplemental ------------------------------------------------------------

# Linear regression of RT --

# effects code the "Switch" variable
linear_regression <- lm(TS_RT ~ Switch, data=data.Ctl)
print(summary(linear_regression))
cat('The residual variance (sigma^2) of the linear regession equals', sigma(linear_regression)^2, '\n')

# Linear MLM of RT  --
# Note: this requires lme4 to be loaded

linear_MLM <- lmer(TS_RT ~ Switch + (1|PID), data=data.Ctl)
tau20      <- VarCorr(linear_MLM)$PID[1] # between-subject variance 
sigma2     <- sigma(linear_MLM)^2        # within (residual) subject variance
ICC        <- tau20/(tau20+sigma2)

print(summary(linear_MLM))
cat('The ICC for the linear MLM is', ICC, '\n')

# Reproduce Figure S1
# If this looks confusing, see the comments in section "Reproduce Figure 2" above. 
# The steps are basically the same. 

layout(1)
U0js <- coef(linear_MLM)$PID[,1]
g00  <- fixef(linear_MLM)[1]
hist(U0js, xlab = expression(gamma['00']+U['0j']), ylab='Number of Particpants', main='')
abline(v=g00, col='red', lty=2, lwd=2)
legend('topright', bty='n', col='red', lty=2, lwd=2, legend=expression(gamma['00']))
