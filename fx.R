

icc_thr <- function(fit, gr) {
  #' Computes the ICC using the Latent Threshold Approach
  #'
  #' @param fit A fit glmer model 
  #' @param gr  A string of the group identifier (e.g., "PID")
  #'
  #' @return A singular ICC value
  #'
  #' @examples icc_thr(model, 'PID')

  tau2 <- lme4::VarCorr(fit)[[gr]][1]
  tau2/(tau2+pi^2/3)
  
}


icc_sim <- function (fit, gr, avg_preds = T, seed = 12345, nsims = 5e3) {
  
  #' Computes the ICC using a modified version of Goldstein's (2002) Simulation Aproach
  #' Specifically, it computes the simulated ICC for the average predicted fixed effect
  #'
  #' @param fit A fit glmer model 
  #' @param gr A string of the group identifier (e.g., "PID")
  #' @param avg_preds Should the average of predictor values be used to compute the ICC?
  #' @param seed Random seed value for reproducibility
  #' @param nsims Number of simulations (should be a large number)
  #'
  #' @return A singular ICC value
  #'
  #' @examples icc_sim(model, 'PID', avg_preds = TRUE, seed=2022, nsims=1e4)
  
  set.seed(seed)
  tau2          <- lme4::VarCorr(fit)[[gr]][1]
  y_lgt_pred    <- predict(fit, re.form = NA)
  if(avg_preds) {
    # use the average of the predicted values to compute ICC
    mn_y_lgt_pred <- mean(y_lgt_pred)
    p_hat         <- plogis(rnorm(nsims, mn_y_lgt_pred, sqrt(tau2)))
    v_p_hat       <- var(p_hat)
    m_p_hat       <- mean(p_hat)
    return(v_p_hat / (m_p_hat - m_p_hat ^ 2))
  }
  # otherwise, compute a unique ICC for each unique value in `y_lgt_pred`
  y_lgt_pred_un <- unique(y_lgt_pred)
  sapply(y_lgt_pred_un, function(x) {
    p_hat         <- plogis(rnorm(nsims, x, sqrt(tau2)))
    v_p_hat       <- var(p_hat)
    m_p_hat       <- mean(p_hat)
    v_p_hat / (m_p_hat - m_p_hat ^ 2)
  })

}


icc_num <- function (fit, gr, avg_preds = T) {
  #' Computes the ICC by estimating the moments of the logit-normal distribution via numerical integration
  #'
  #' @param fit A fit glmer model 
  #' @param gr  A string of the group identifier (e.g., "PID")
  #' @param avg_preds Should the average of predictor values be used to compute the ICC?

  #' @return A singular ICC value
  #'
  #' @examples icc_num(model, 'PID', avg_preds=TRUE)
  
  tau2 <- lme4::VarCorr(fit)[[gr]][1]
  y_lgt_pred <- predict(fit, re.form = NA)
  
  if(avg_preds) {
    #  use the average of the predicted values to compute ICC
    mn_y_lgt_pred <- mean(y_lgt_pred)
    mu_p <- integrate(function (x) {
      plogis(x) * dnorm(x, mn_y_lgt_pred, sqrt(tau2))
    }, -Inf, Inf)$value
    tau2_p <- integrate(function (x) {
      (plogis(x) - mu_p) ^ 2 * dnorm(x, mn_y_lgt_pred, sqrt(tau2))
    }, -Inf, Inf)$value
    return(tau2_p / (mu_p - mu_p ^ 2))
  }
  # otherwise, compute a unique ICC for each unique value in `y_lgt_pred`
  y_lgt_pred_un <- unique(y_lgt_pred)
  sapply(y_lgt_pred_un, function(y_pred) {
    mn_y_lgt_pred <- y_pred
    mu_p <- integrate(function (x) {
      plogis(x) * dnorm(x, mn_y_lgt_pred, sqrt(tau2))
    }, -Inf, Inf)$value
    tau2_p <- integrate(function (x) {
      (plogis(x) - mu_p) ^ 2 * dnorm(x, mn_y_lgt_pred, sqrt(tau2))
    }, -Inf, Inf)$value
    tau2_p / (mu_p - mu_p ^ 2)
  })

}


icc_lin <- function (fit, gr, avg_preds=T) {
  #' Computes the ICC using linearization
  #'
  #' @param fit A fit glmer model 
  #' @param gr  A string of the group identifier (e.g., "PID")
  #' @param avg_preds Should the average of predictor values be used to compute the ICC?
  #'
  #' @return A singular ICC value
  #'
  #' @examples icc_lin(model, 'PID')
  
  tau2 <- lme4::VarCorr(fit)[[gr]][1]
  y_lgt_pred <- predict(fit, re.form = NA)
  
  if(avg_preds) {
    #  use the average of the predicted values to compute ICC
    mn_y_lgt_pred <- mean(y_lgt_pred)
    pi_ <- plogis(mn_y_lgt_pred)
    v_2 <- tau2 * pi_ ^ 2 / (1 + exp(mn_y_lgt_pred)) ^ 2
    return(v_2 / (v_2 + pi_ - pi_ ^ 2))
  }
  # otherwise, compute a unique ICC for each unique value in `y_lgt_pred`
  y_lgt_pred_un <- unique(y_lgt_pred)
  sapply(y_lgt_pred_un, function(y_pred) {
    mn_y_lgt_pred <- y_pred
    pi_ <- plogis(mn_y_lgt_pred)
    v_2 <- tau2 * pi_ ^ 2 / (1 + exp(mn_y_lgt_pred)) ^ 2
    v_2 / (v_2 + pi_ - pi_ ^ 2)
  })
}

MOR <- function(fit, gr) {
  #' Computes Merlo's et al. (2006) Median Odds Ratio
  #'
  #' @param fit A fit glmer model 
  #' @param gr  A string of the group identifier (e.g., "PID")
  #'
  #' @return A singular MOR value
  #'
  #' @examples MOR(model, 'PID')
  tau2  <- lme4::VarCorr(fit)[[gr]][1]
  exp(sqrt(2*tau2)*qnorm(.75))
  
}

bootstrap_icc <- function(fit, gr, method, B=5e3, seed=1234) {
  #' Computes bootstrapped values for methods used above
  #'
  #' @param fit A fit glmer model
  #' @param gr A string of the group identifier (e.g., "PID")
  #' @param method A string of the method to compute (e.g., "icc_lin")
  #' @param B  The number of bootstrapped samples
  #' @param seed Random seed value for reproducibility
  #'
  #' @return A vector (1xB) of bootstrapped estimates
  #'
  #' @examples bootstrap_icc(model,"PID","icc_lin",B=1e3,seed=2022)

  sims <- simulate(fit, B, seed)
  pb <- txtProgressBar(0,B)
  out <- rep(NA, B)
  for(i in 1:B) {
    setTxtProgressBar(pb,i)
    mod <- lme4::refit(fit, newresp = sims[[i]])
    out[i] <- eval(parse(text=paste0(method,'(mod,"',gr,'")')))
  }
  out
}
