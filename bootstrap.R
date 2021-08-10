B = 3
BSModels = list()
for(bs in 1:B) {
  
  cat('Bootstrap #', bs, '/', B, '\n')
  clusters = model@frame[,ncol(model@frame)]
  clusterName = tail(names(model@frame), 1)
  K = length(unique(clusters))
  
  tau2 = as.numeric(as.data.frame(VarCorr(model))['vcov'])
  tau2Hat = rnorm(K, 0 , tau2)
  gammas = fixef(model)
  Bhat = matrix(gammas[-1])
  
  Yij = c()
  for(j in 1:K) {
    DesignMatrix = model.matrix(model)[model@frame[,clusterName] == unique(clusters)[j], -1]
    LORj = gammas[1] + tau2Hat[j] + DesignMatrix %*% Bhat
    pj = exp(LORj)/(1+exp(LORj))
    Yij = c(Yij, rbinom(n=length(pj), size=1, prob = pj))
    if(sum(is.na(Yij)) > 0) break
  }
  
  DesignMatrix = model.matrix(model)[, -1]
  BSModels[[bs]] = lme4::glmer(Yij ~  DesignMatrix + (1|clusters), family='binomial')
}


