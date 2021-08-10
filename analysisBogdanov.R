
# Load data ---------------------------------------------------------------

data <- read.csv('Bogdanovetal2021/DST_data_osf2.csv')
data.Ctl <- data[data$condition=='control', ] # focus just on control condition
demo <- read.csv('Bogdanovetal2021/DST_demographics_osf.csv')
colnames(demo)[which(colnames(demo)=='ID')] <- 'PID'
data.Ctl <- merge(data.Ctl, demo[,c('PID', 'Gender')], by='PID')

# Descriptives ------------------------------------------------------------

N <- length(unique(data$PID)) # number of subjects
nTrials <- table(data$PID)    # number of trials per subject

# Linear regression of RT -------------------------------------------------

contrasts(data$Switch) <- contr.sum(2)

linear_regression <- lm(TS_RT ~ Switch, data=data.Ctl)

summary(linear_regression)

sigma(linear_regression)^2


# Linear MLM of RT --------------------------------------------------------

library(lme4)

linear_MLM <- lmer(TS_RT ~ Switch + (1|PID), data=data.Ctl)

summary(linear_MLM)

tau20 <- attr(VarCorr(linear_MLM)[[1]], 'stddev')^2
sigma2 <- sigma(linear_MLM)^2

vpc_linear_MLM <- tau20/(tau20+sigma2)

linear_MLM2 <- lmer(TS_RT ~ Switch + block + (1|PID), data=data.Ctl)

sigma2_M2 <- sigma(linear_MLM2)^2  

pseudoR2 <- (sigma2-sigma2_M2)/sigma2
  
# Logistic model for choice -------------------------------------------------------------------

glm0 <- glm(effort_choice ~ 1, data=data.Ctl, family='binomial')
logistic_MLM0 <- glmer(effort_choice ~ 1 + (1|PID), data=data.Ctl, family='binomial')

summary(logistic_MLM0)

effAvoidance <- exp(fixef(logistic_MLM0))/(1+exp(fixef(logistic_MLM0)))

eb <- coef(logistic_MLM0)$PID
eb$mean_eff_choice <- tapply(data.Ctl$effort_choice, data.Ctl$PID, mean)

hist(eb[[1]], main='', xlab=expression(U[`0j`]+gamma[`00`]), ylim=c(0, 25))
#plot(density(eb[[1]]))
abline(v=0, lty='dashed', lwd=3)
text(-.75, 20, labels = 'Demand\nseeking')
text(.75, 20, labels = 'Demand\navoidant')


# VPC: Latent variable approach -------------------------------------------

tau20 <- attr(VarCorr(logistic_MLM0)[[1]], 'stddev')^2
sigma2 <- pi^2/3
threshVPC <- tau20/(tau20+sigma2)

threshVPC


# VPC: Simulation method ---------------------------------------------------------------------

# 0. Extract intercept-variance and fixed effects and specify number of simulations
tau0 <- attr(VarCorr(logistic_MLM0)[[1]], 'stddev') 
gamma00 <- fixef(logistic_MLM0)[[1]]
M <- 1e5
set.seed(10408) # for reproducibility 

# 1. Simulate random effects
U0j <- rnorm(M, 0, tau0)

# 2. Predict probabilities 
logit_p_hat <- gamma00+U0j
p_hat <- exp(logit_p_hat)/(1+exp(logit_p_hat))

# 3. Compute level-1 variance (Bernouilli variance)
var_L1 <- p_hat*(1-p_hat) 

# 4. Compute VPC
sigma2 <- mean(var_L1)
simVPC <- var(p_hat)/(var(p_hat) + sigma2)

simVPC

# For a more complex model! 

data.Ctl$Gender <- factor(data.Ctl$Gender)
contrasts(data.Ctl$Gender) <- contr.sum(2)
logisticMLM_gender <- glmer(effort_choice ~ Gender + (1|PID), data=data.Ctl, family='binomial')

# 0. Extract intercept-variance and fixed effects and specify number of simulations
tau0 <- attr(VarCorr(logisticMLM_gender)[[1]], 'stddev') 
gamma <- fixef(logisticMLM_gender)
M <- 1e5
set.seed(10408) # for reproducibility 

# 1. Simulate random effects
U0j <- rnorm(M, 0, tau0)

# 2. Predict probabilities for repeat and switch trials
logit_p_hat_m <- gamma[1]+U0j + gamma[2]*1
logit_p_hat_f <- gamma[1]+U0j + gamma[2]*-1

p_hat_m <- exp(logit_p_hat_m)/(1+exp(logit_p_hat_m))
p_hat_f <- exp(logit_p_hat_f)/(1+exp(logit_p_hat_f))

# 3. Compute level-1 variance (Bernouilli variance)
var_L1_m <- p_hat_m*(1-p_hat_m) 
var_L1_f <- p_hat_f*(1-p_hat_f) 

# 4. Compute VPC for repeat and switch trials
sigma2_m <- mean(var_L1_m)
sigma2_f <- mean(var_L1_f)

simVPC_m <- var(p_hat_m)/(var(p_hat_m) + sigma2_m)
simVPC_f <- var(p_hat_f)/(var(p_hat_f) + sigma2_f)

c('male'=simVPC_m, 'female'=simVPC_f)


# Linearization -----------------------------------------------------------



# VPC: MOR ----------------------------------------------------------------

tau0 <- attr(VarCorr(logistic_MLM0)[[1]], 'stddev') 
phi75 <- qnorm(.75) # 75th percentile of normal CDF
MOR <- exp(sqrt(2*tau0)*phi75)

as.numeric(MOR)

# For a more complex model! 

tau0 <- attr(VarCorr(logisticMLM_gender)[[1]], 'stddev') 
phi75 <- qnorm(.75) # 75th percentile of normal CDF
MOR <- exp(sqrt(2*tau0)*phi75)

MOR

# Bootstrapping ------------------------------------------------------

B <- 1000
K <- length(unique(logistic_MLM0@frame$PID))
tau0 <- attr(VarCorr(logistic_MLM0)[[1]], 'stddev') 

for(iteration in B) {
  U0j <- rnorm(K, 0, tau0)
  
}




# misc --------------------------------------------------------------------


source('ICC/LogitICC.R')
LogitICC(m0, method = 'MOR')
LogitICC(logistic_MLM0, method = 'simGold')
LogitICC(m0, method = 'threshold')

source('R2Dicho.R')
R2dicho(m0, 'threshold')
R2dicho(m0, 'simGold')



