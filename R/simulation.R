library(readr)
library(dplyr)
library(tidyr)
library(tidytext)
library(Matrix)
library(e1071)
library(mvtnorm)
library(sampling)
library(mase)
library(spdep)
library(maptools)
library(ggmap)
library(maps)
library(spatialreg)
library(tm)
library(ggthemes)
source('R/models.R')


### Construct spatial adjacency matrix
m1 <- map(database='state', fill=T, plot=F)
df1 <- data.frame(polyname=m1$names) %>% left_join(state.fips)
m2 <- map2SpatialPolygons(m1, IDs=df1$fips)
W <-as.matrix(as_dgRMatrix_listw(nb2listw(poly2nb(m2), style='B')))
W <- W[-8,-8]
W <- W/rowSums(W)
eig <- eigen(W)
bf <- eig$vectors[,1:25]

### ANES survey data
anes <- read_delim('Data/anes.txt','|')
anes <- anes %>% filter(dem_hisp %in% c(" 1"," 2"))
anes$caseid <- as.numeric(anes$caseid)

### Free response text data
anesText <- read_csv('Data/anesText.csv')
probs1 <- anesText %>% select(caseid, mip_prob1)
probs1 <- probs1 %>% unnest_tokens(word, mip_prob1) %>% anti_join(stop_words) 
topWords <- probs1 %>% count(word, sort=T)
probs1 <- probs1 %>% filter(word %in% topWords$word[1:1000])
probs1$ones <- 1
probs1 <- probs1 %>% group_by(caseid, word) %>% summarize(n=sum(ones))
probs1 <- probs1 %>% spread(word, n, fill=0)




anes2 <- anes %>% select(caseid, fedspend_welfare, dem_hisp, weight_full, gender_respondent_x, sample_state, sample_stfips, 
                         presapp_job)
wordDF <- data.frame(caseid=probs1$caseid, words=I(as.matrix(probs1[,-1])))
anes2 <- anes2 %>% left_join(wordDF, by='caseid')
anes2$words[is.na(anes2$words)] <- 0

anes2 <- anes2 %>% mutate(dem_hisp=factor(dem_hisp), gender=factor(gender_respondent_x))
anes2$R <- anes2$presapp_job==" 1"
anes2 <- anes2 %>% filter(!(sample_state %in% c("DC", "AK", "HI")))
anes2$words2 <- anes2$words>0
anes2$Psi <- I(bf[as.numeric(as.factor(anes2$sample_stfips)),])
fullX1 <-  model.matrix(~  dem_hisp + gender , data=anes2)


### Set up simulation storage
nRep <- 100
grid1 <- expand.grid(unique(anes2$sample_state))
mod1Est <- mod2Est <- mod3Est <- HTEst <-  UWEst <- matrix(NA, nrow=nrow(grid1), ncol=nRep)
mod1Sd <- mod2Sd <- mod3Sd <- mod1L <- mod2L <- mod3L <- mod1H <- mod2H <- mod3H <- matrix(NA, nrow=nrow(grid1), ncol=nRep)
sampSizes <- matrix(NA, nrow=nrow(grid1), ncol=nRep)
HTsd <- matrix(NA, nrow=nrow(grid1), ncol=nRep)

### True population values
truth <- grid1
truth$R <- NA
for(j in 1:nrow(truth)){
  subD <- anes2$R[anes2$sample_state==truth$Var1[j]]
  truth$R[j] <- mean(subD)
}


### Run empirical simulation
ss <- 1000
set.seed(1)
class(anes2$Psi) <- 'matrix'
for(i in 1:nRep){
  repeat{ ## ensure all states have a sample so that we can get HT estimate
    prob <- inclusionprobabilities(anes2$weight_full + 0.7*(anes2$R), ss)
    Ind <- UPpoisson(prob)
    sample <- anes2[as.logical(Ind),]
    sample$P <- prob[as.logical(Ind)]
    sample$W <- 1/sample$P
    sample$scaledWGT <- (sample$W)*nrow(sample)/sum(sample$W)
    if(length(unique(sample$sample_state))==48) break
  }
  sampSizes[,i] <- (grid1 %>% left_join(sample %>% group_by(sample_state) %>% summarize(ss=n()), by=c("Var1"="sample_state")))$ss
  
  
  fun1 <- function(st, n=10){ ## imputes text data
    poss <- which(sample$sample_state==st)
    wordSet <- poss[sample.int(length(poss),n, replace=T, prob=sample$P[poss])]
    return(wordSet)
  }
  
  ## Generate random weights
  elm <- matrix(rnorm(240*ncol(cbind(1,sample$words2)))*rbinom(240*ncol(cbind(1,sample$words2)), 1, 0.9), nrow=ncol(sample$words2)+1) 
  
  ## Fit baseline model
  modX1 <- model.matrix(~  dem_hisp + gender , data=sample)
  mod1 <- pgVBbase(X=cbind(modX1, sample$Psi), Y=sample$R, wgt=sample$scaledWGT) ## Base model
  
  ## Fit proposed
  elmo1 <- cbind(1,sample$words2)%*%elm
  elmo1 <- sigmoid(elmo1)
  mod2 <- pgVB(X=cbind(modX1, sample$Psi), Z=elmo1,  
               Y=sample$R, wgt=sample$scaledWGT, Au=0.5, Bu=0.5) ## ELM Model
  
  
  ## Get predictions
  beta1 <- rmvnorm(1000, mean=mod1$BU$mean, sigma = as.matrix(mod1$BU$Cov))
  beta2 <- rmvnorm(1000, mean=mod2$BU$mean, sigma = as.matrix(mod2$BU$Cov))
  
  preds1 <- sigmoid(cbind(fullX1, anes2$Psi) %*% t(beta1))
  wordInd <- t(sapply(anes2$sample_state, fun1, n=1000))
  
  preds2 <- matrix(NA, nrow=nrow(anes2), ncol=1000)
  pb <- txtProgressBar(min=0, max=1000, style=3)
  for(ind in 1:1000){
    elmo2 <- Matrix(cbind(1,sample$words2[wordInd[,ind], ]))%*%elm   
    elmo2 <- sigmoid(elmo2)
    preds2[,ind] <- as.numeric(sigmoid(cbind(fullX1, as.matrix(anes2$Psi), elmo2) %*% t(beta2)[,ind]) )
    setTxtProgressBar(pb, ind)
  }
  

  elmo3 <- cbind(1,anes2$words2)%*%elm
  elmo3 <- sigmoid(elmo3)
  preds3 <- sigmoid(cbind(fullX1, as.matrix(anes2$Psi), elmo3) %*% t(beta2))
  
  
  
  preds1 <- matrix(rbinom(length(preds1), 1, c(preds1)), ncol=1000)
  preds2 <- matrix(rbinom(length(preds2), 1, c(preds2)), ncol=1000)
  preds3 <- matrix(rbinom(length(preds2), 1, c(preds3)), ncol=1000)
  
  ## Calculate domain estimates and other quantities
  for(j in 1:nrow(grid1)){
    subD <- preds1[anes2$sample_state==grid1$Var1[j],]
    subD <- apply(subD,2,mean)
    mod1Est[j,i] <- mean(subD)
    mod1Sd[j,i] <- sd(subD)
    mod1L[j,i] <- quantile(subD, prob=0.025)
    mod1H[j,i] <- quantile(subD, prob=0.975)
    
    
    subD <- preds2[anes2$sample_state==grid1$Var1[j],]
    subD <- apply(subD,2,mean)
    mod2Est[j,i] <- mean(subD)
    mod2Sd[j,i] <- sd(subD)
    mod2L[j,i] <- quantile(subD, prob=0.025)
    mod2H[j,i] <- quantile(subD, prob=0.975)
    
    subD <- preds3[anes2$sample_state==grid1$Var1[j],]
    subD <- apply(subD,2,mean)
    mod3Est[j,i] <- mean(subD)
    mod3Sd[j,i] <- sd(subD)
    mod3L[j,i] <- quantile(subD, prob=0.025)
    mod3H[j,i] <- quantile(subD, prob=0.975)
    
    
    subD <- sample %>% filter(sample_state==grid1$Var1[j])
    HTEst[j,i] <- horvitzThompson(as.numeric(subD$R), pi=1/subD$W)$pop_mean
    HTsd[j,i] <- sqrt(horvitzThompson(as.numeric(subD$R), pi=1/subD$W, var_est=T, var_method="LinHB")$pop_mean_var)
    UWEst[j,i] <- mean(subD$R)
  }
  
  print(paste0("Done with simulation ",i))
}
##### Notes #####
## mod1 is PLLR
## mod2 is BURN-IT with imputation
## mod3 is BURN-KT without imputation (known text covariates)


### MSE
mean((mod1Est - truth$R)^2)
mean((mod2Est - truth$R)^2)
mean((mod3Est - truth$R)^2)
mean((HTEst - truth$R)^2)
mean((UWEst - truth$R)^2)
mean((mod2Est - truth$R)^2)/mean((mod1Est - truth$R)^2)

### Squared Bias
mean((rowMeans(mod1Est) - truth$R)^2)
mean((rowMeans(mod2Est) - truth$R)^2)
mean((rowMeans(mod3Est) - truth$R)^2)
mean((rowMeans(HTEst) - truth$R)^2)
mean((rowMeans(UWEst) - truth$R)^2)



### Generate output
df2 <- data.frame(region=tolower(state.name), abb=state.abb)
results <- data.frame(abb=grid1$Var1,PLLR=rowMeans(mod1Sd), 
                      BURN_IT=rowMeans(mod2Sd),
                      BURN_KT=rowMeans(mod3Sd),
                      Direct=rowMeans(sqrt(HTsd), na.rm=T))
sts <- map_data('state') %>% left_join(df2, by=c("region")) %>% left_join(results, by="abb")
sts <- sts %>% pivot_longer(8:11, names_to="Model", values_to="SE")
sts$Model <- if_else(sts$Model=="BURN_IT", "BURN-IT", if_else(sts$Model=="BURN_KT", "BURN-KT", sts$Model))
sts$Model <- factor(sts$Model, levels=c("Direct", "PLLR", "BURN-IT", "BURN-KT"))

# Standard error plot
ggplot(sts)+
  geom_polygon(aes(x=long, y=lat, group=group, fill=SE))+
  facet_wrap(~Model, nrow=2)+
  scale_fill_viridis_c()+
  theme_map()


# MSE plot
df2 <- data.frame(BURN=rowMeans((mod2Est - truth$R)^2), Dir=rowMeans((HTEst - truth$R)^2), ss=rowMeans(sampSizes))
ggplot(df2, aes(x=Dir, y=BURN))+
  geom_point()+
  geom_abline(slope=1, color="red", linetype='dashed')+
  theme_classic()+
  ylab("BURN-IT")+
  xlab("Direct")


# Ratio of estimates
df2 <- data.frame(Ratio=rowMeans(HTEst/mod2Est), ss=rowMeans(sampSizes))
ggplot(df2, aes(x=ss, y=Ratio))+
  geom_point()+
  geom_abline(slope=0, intercept=1, color='red', linetype='dashed')+
  theme_classic()+
  xlab("Sample Size")+
  ylab("Ratio of Estimates")



