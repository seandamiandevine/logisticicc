#' Example code for Quantifying Within- and Between-Cluster Variance in Logistic Multilevel Models
#' @author:  Sean Devine
#' @contact: seandamiandevine@gmail.com
#' @github:  https://github.com/seandamiandevine/logisticicc
#' @paper:   TODO

# Load data ---------------------------------------------------------------
# Load the data from .csv and do the following things: 
  # 1. subset the data such that we only look at participants in the control condition 
  # 2. Add participant demographic information into the dataframe of participant choice
  # 3. For simplicity, recode effort choices so that 1 = high effort and 0 = low effort

data                   <- read.csv('Bogdanovetal2021/DST_data_osf2.csv')
data.Ctl               <- data[data$condition=='control', ] # focus just on control condition
demo                   <- read.csv('Bogdanovetal2021/DST_demographics_osf.csv')
colnames(demo)[which(colnames(demo)=='ID')] <- 'PID'
data.Ctl               <- merge(data.Ctl, demo[,c('PID', 'Gender')], by='PID')
data.Ctl$effort_choice <- abs(data.Ctl$effort_choice-1) 

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
eb                 <- coef(logistic_MLM0)$PID

# extract the estimated random intercept variance (tau^2_0)
tau = VarCorr(logistic_MLM0)$PID[1]

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

# ICC: Simulation method (Conditional ICC) ---------------------------------------------------------------------

# Convert gender to a factor in the data-framre and effects code it 
data.Ctl$Gender            <- factor(data.Ctl$Gender)
contrasts(data.Ctl$Gender) <- contr.sum(2)
logisticMLM_gender         <- glmer(effort_choice ~ Gender + (1|PID), data=data.Ctl, family='binomial')

# 0. Extract intercept-variance and fixed effects and specify number of simulations
set.seed(2022) # for reproducibility 
tau20         <- VarCorr(logisticMLM_gender)$PID[1]
gamma         <- fixef(logisticMLM_gender)
M             <- 1e5

# 1. Simulate random effects
U0j           <- rnorm(M, 0, sqrt(tau20))

# 2. Predict probabilities for repeat and switch trials
logit_p_hat_m <- gamma[1]+U0j + gamma[2]*1
logit_p_hat_f <- gamma[1]+U0j + gamma[2]*-1

p_hat_m       <- exp(logit_p_hat_m)/(1+exp(logit_p_hat_m))
p_hat_f       <- exp(logit_p_hat_f)/(1+exp(logit_p_hat_f))

# 3. Compute level-1 variance (Bernouilli variance)
var_L1_m      <- p_hat_m*(1-p_hat_m) 
var_L1_f      <- p_hat_f*(1-p_hat_f) 

# 4. Compute ICC for repeat and switch trials
sigma2_m      <- mean(var_L1_m)
sigma2_f      <- mean(var_L1_f)

simICC_m      <- var(p_hat_m)/(var(p_hat_m) + sigma2_m)
simICC_f      <- var(p_hat_f)/(var(p_hat_f) + sigma2_f)

# print result
cat('ICC for men using the Simulation Model = ',simICC_m,'\n')
cat('ICC for women using the Simulation Model = ',simICC_f,'\n')

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
MOR   <- exp(sqrt(2*tau0)*phi75)
 
cat('Median Odds Ratio = ',MOR,'\n')

# Bootstrapping ------------------------------------------------------

# 0. Set constants for bootstrapping procedure
B       <- 1000                    # number of bootstraps
ids     <- logistic_MLM0@frame$PID # extract id vector from model data
K       <- length(unique(ids))     # number of clusters (subjects)
nTrials <- table(ids)              # number of trials per subject  
tau20   <- VarCorr(logistic_MLM0)$PID[1]
g00     <- fixef(logistic_MLM0)[[1]]
output  <- matrix(NA, B, 7, dimnames = list(NULL, c('iteration', 'threshICC', 'simICC', 'linICC', 'MOR', 'g00', 'tau20')))

# 1. Cycle through iterations (i) of bootstrapped samples
for(i in 1:B) {
  if(i%%100==0) cat('bootstrapping is', round(i/B,2)*100, '% complete.')
  # 1.1 Simulate data
  U0j           <- rnorm(K, 0, tau20)        # random deviations per subject
  LOR_j         <- g00 + U0j                 # log odds of response==1 per subject
  p_j           <- exp(LOR_j)/(1+exp(LOR_j)) # convert LOR to probability
  y_ij          <- sapply(1:K, function(k) rbinom(nTrials[k], 1, p_j[k]))
  y_ij          <- unlist(y_ij)              # break out of a list format
  
  # 1.2 Fit new model and compute values of interest
  thisMLM       <- glmer(y_ij ~ 1 + (1|ids), family='binomial')
  thisG00       <- fixef(thisMLM)[[1]]
  thistau20     <- VarCorr(thisMLM)$ids[1]
  
  # 1.2.1. Latent Threshold ICC
  thisthreshICC <- thistau20/(thistau20+pi^2/3)

  # 1.2.2. Simulation ICC
  thisU0j       <- rnorm(1e5, 0, sqrt(thistau20))
  logit_p_hat   <- thisG00+U0j
  p_hat         <- exp(logit_p_hat)/(1+exp(logit_p_hat))
  var_L1        <- p_hat*(1-p_hat) 
  sigma2        <- mean(var_L1)
  thissimICC    <- var(p_hat)/(var(p_hat) + sigma2)
  
  # 1.2.3. Linearlization 
  p             <- exp(thisG00)/(1+exp(thisG00))
  var1          <- p*(1-p)
  var2          <- thistau20*p^2*(1 + exp(gamma00))^(-2)
  thislinICC    <- var2/(var1+var2)
  
  # 1.2.4. MOR
  thisMOR       <- exp(sqrt(2*thistau20)*qnorm(.75))
  
  # 1.3. Save output
  output[i,] = c(i, thisthreshICC, thissimICC, thislinICC, thisMOR, thisG00, thistau20)
  
}


# Reproduce Figure 3 from paper -------------------------------------------

layout(matrix(1:6,2,3,byrow=T))
pretty_titles <- list('threshICC'='Latent Threshold', 'simICC'='Simulation', 'linICC' = 'Linearization', 
                      'MOR'='MOR', 'g00'=expression(bold(gamma['00'])), 'tau20'=expression(bold(tau[0]^2)))
for(var in c('threshICC', 'simICC', 'linICC', 'MOR', 'g00', 'tau20')) {
  
  hist(output[,var], main=pretty_titles[[var]], col='white', xlab='')
  
}


# Supplemental ------------------------------------------------------------

# Linear regression of RT --

# effects code the "Switch" variable
contrasts(data$Switch) <- contr.sum(2) 
linear_regression      <- lm(TS_RT ~ Switch, data=data.Ctl)
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



