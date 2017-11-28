# Purpose: test differences between pairwise genetic distances from two different samples.
# The program compares two samples at the time, but it has the capacity of doing "batches of pairs,"
# as for example when one wishes to test two different regions from a set of patients to see whether
# or not they evolved at different evolutionary rates.
# Written by: EEGiorgi, egiorgi@lanl.gov -- Latest update: 2/23/2011
# Authors: EE Giorgi and T Bhattacharya
# LA-CC-11-024 (2/24/2011)
# 
#
# adapted by Joy Pai for integration into HIV pairwise genetic distance analysis (2017)
# Nussenzweig Lab, Rockefeller University

OurVarEst <- function(s11,s12,s21,s22,n1,n2){
  
  temp1 <- 4*((n1*(n1-1)*(n1-2)*(n1-3))^(-1))
  temp2 <- 4*((n2*(n2-1)*(n2-2)*(n2-3))^(-1))
  
  ourv <- temp1*(2*s11+s12)+temp2*(2*s21+s22)
  return(ourv)
  
}

myvarmu <- function(n1, sigm1, sigm2){
  
  temp1 <- 4*((n1*(n1-1)*(n1-2)*(n1-3))^(-1))
  return(temp1*(2*sigm1+sigm2))
  
}

mynu <- function(ns, sigm1, sigm2){
  
  num <- ((2/((ns-2)*(ns-3)))*(sigm2+2*sigm1))^2
  den1 <- (4*(ns-1)/((ns-2)^2))*((2/((ns-1)*(ns-2)))^2)*((sigm2+sigm1)^2)
  den2 <- ((2*ns)/((ns-3)*(ns-2)^2))*((2/(ns*(ns-2)*(ns-3)))*((ns-4)*sigm2-2*sigm1))^2
  return(num/(den1+den2))
  
}

two_sample_test <- function(dist_dfs) {
  nbases <- rep(0,2)
  nseq <- rep(0,2)
  mult <- rep(0,2)
  dvec <- vector("list", length=2)
  vecone <- vector("list", length=2)
  vectwo <- vector("list", length=2)
  
  for(i in 1:2){	
    df <- dist_dfs[[i]]
    vecone[[i]] <- df$row
    vectwo[[i]] <- df$col
    nbases[i] <- df$length[1]
    dvec[[i]] <- df$value
  
    nseq[i] <- 1+max(df$row)
    mult[i] <- nseq[i]*(nseq[i]-1)*(0.5)
    
  }
  
  if(nbases[1] != nbases[2]){ 
    dvec[[1]] <- dvec[[1]]/nbases[1]
    dvec[[2]] <- dvec[[2]]/nbases[2]
  }
  
  sample1 <- dvec[[1]]
  sample2 <- dvec[[2]] ## sample2 is the larger (for computational convenience)
  
  N1 <- nseq[1]
  N2 <- nseq[2]
  
  mu1 <- sum(sample1)/length(sample1)
  mu2 <- sum(sample2)/length(sample2)
  
  sigma12 <- sum((sample1-mu1)^2)
  sigma22 <- sum((sample2-mu2)^2)
  
  #### calculate the two means and the two statistics and then store the p-value
  #### need to compute the pairwise distances
  
  sigma11 <-  0 #sample1
  sigma21 <-  0 #sample2
  
  for(i in 1:(N1-2)){
    for(j in (i+1):(N1-1)){
      
      d1ij <- sample1[which((vecone[[1]]==i)&(vectwo[[1]]==j))]
      
      for(l in (j+1):N1){
        d1il <- sample1[which((vecone[[1]]==i)&(vectwo[[1]]==l))]
        d1jl <- sample1[which((vecone[[1]]==j)&(vectwo[[1]]==l))]
        sigma11 <- sigma11+(d1ij-mu1)*(d1il-mu1)+(d1ij-mu1)*(d1jl-mu1)+(d1il-mu1)*(d1jl-mu1)
      }
    }
  }
  
  for(i in 1:(N2-2)){
    for(j in (i+1):(N2-1)){
      
      d2ij <- sample2[which((vecone[[2]]==i)&(vectwo[[2]]==j))]
      
      for(l in (j+1):N2){
        d2il <- sample2[which((vecone[[2]]==i)&(vectwo[[2]]==l))]
        d2jl <- sample2[which((vecone[[2]]==j)&(vectwo[[2]]==l))]
        sigma21 <- sigma21+(d2ij-mu2)*(d2il-mu2)+(d2ij-mu2)*(d2jl-mu2)+(d2il-mu2)*(d2jl-mu2)
      }
    }
    
  }
  
  #### T-stats
  
  OurT <- (mu1-mu2)/sqrt(OurVarEst(sigma11,sigma12,sigma21,sigma22,N1,N2))
  
  npooled <- (1/N1 + 1/N2)^(-1)
  varmu1 <- myvarmu(N1, sigma11, sigma12)
  varmu2 <- myvarmu(N2, sigma21, sigma22)
  coeff1 <- varmu1/(varmu1+varmu2)
  coeff2 <- varmu2/(varmu1+varmu2)
  nu1 <- mynu(N1, sigma11, sigma12)
  nu2 <- mynu(N2, sigma21, sigma22)
  nu <- ((coeff1^2)*nu1^(-1) + (coeff1^2)*nu2^(-1))^(-1)
  
  OurPval_norm <- 2*(1-pnorm(abs(OurT)))
  OurPval <- 2*(1-pt(abs(OurT), df=nu, lower.tail = TRUE))
  
  print(paste("T=", OurT, "df=", nu, sep=" "))
  print(paste("Z-test P=", OurPval_norm, sep=" "))
  print(paste("T-test P=", OurPval, sep=" "))
  
  results <- list()
  results$t.statistic <- OurT
  results$p.ztest <- OurPval_norm
  results$p.ttest <- OurPval
  
  return(results)

}
