R2dicho <- function(model, sigmaR='', M=5e4, nCov=5, simReturnFunction=mean, seed=NULL) {
  tau2 = as.numeric(as.data.frame(VarCorr(model))['vcov'])
  sigmaF = var(predict(model))
  switch (sigmaR,
    'threshold' = {
      sigmaR = pi^2/3
      }, 
    'simGold' = {
      oldSeed <- .Random.seed # store old seed so you don't mess with anything outside the function
      
      if(is.null(seed)) {
        set.seed(10408) # change seed to my custom seed for reproducibility
      } else {
        set.seed(seed)
      }
      
      # 1. Extract relevant information from model
      tau2 <- as.numeric(as.data.frame(VarCorr(model))['vcov'])  # extract tau from model
      Gxj <- matrix(fixef(model))                                # fixed effects (gammas)
      if(length(Gxj)==1) {                                       # if null model
        X <-  matrix(0)          
      } else {
        X <- model@frame[,2:(ncol(model@frame)-1), drop=F] 
      }
      Y <- model@frame[,1, drop=F]                               # (binary) outcome variable
      J <- model@frame[,ncol(model@frame), drop=F]               # grouping factor (clusters)
      tauHat <- rnorm(M, 0, sqrt(tau2))                          # sample between-cluster variance
      
      # 2. Specify a covariate pattern upon which to base the simulations
      covRange <- list()                                      # sample predictor values
      for(pred in 1:ncol(X)) {                                   
        if(is.factor(X[,pred])) {                             # deal with categorical predictors
          codings <- contrasts(X[,pred])
          X[,pred] <- codings[X[,pred]]
          covRange[[pred]] <- range(X[,pred])
        } else {
          covRange[[pred]] <- seq(min(X[,pred]), max(X[,pred]), length.out=nCov) 
        }
      }
      covRange <- append(covRange, 0, 0)                     # append X values for intercept
      covPattern <- expand.grid(covRange)                    # pattern of covariates for simulation
      
      # 3. Simulate!
      vars <- c()
      for(row in 1:nrow(covPattern)) {
        p <- exp(sum(t(Gxj)*covPattern[row,]) + tauHat)/(1+exp(sum(t(Gxj)*covPattern[row,]) + tauHat))
        varP <- p*(1-p)
        vars[row] = mean(varP)
      }
      sigmaR = simReturnFunction(vars)
    }, 
    {sigmaR=pi^2/3}
  )
  
  sigmaF/sum(sigmaF, tau2, sigmaR)
}