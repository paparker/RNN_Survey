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
library(cowplot)
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
probs1 <- probs1 %>% unnest_tokens(word, mip_prob1) %>% anti_join(stop_words) #
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
#anes2$R <- anes2$fedspend_welfare==" 2"
anes2$R <- anes2$presapp_job==" 1"
anes2 <- anes2 %>% filter(!(sample_state %in% c("DC", "AK", "HI")))
anes2$words2 <- anes2$words>0
anes2$Psi <- I(bf[as.numeric(as.factor(anes2$sample_stfips)),])
anes2$bfID <- as.numeric(as.factor(anes2$sample_stfips))


### poststrat cells
pop <- read_csv('Data/popest.csv')
pop <- pop %>% filter(AGE >= 18 & STATE %in% unique(as.numeric(anes2$sample_stfips)) &
                        SEX != 0 & ORIGIN !=0)
pop <- pop %>% group_by(STATE, SEX, ORIGIN) %>% summarize(N=sum(POPESTIMATE2012))
pop <- pop %>% mutate(dem_hisp=factor(if_else(ORIGIN==2,1,2)), gender=factor(SEX))

fullX1 <-  model.matrix(~  dem_hisp + gender , data=anes2)
sample <- anes2
sample$scaledWGT <- nrow(sample)*sample$weight_full/sum(sample$weight_full)
modX1 <- model.matrix(~  dem_hisp + gender , data=sample)

### Random weights
elm <- matrix(rnorm(240*ncol(cbind(1,sample$words2)))*rbinom(240*ncol(cbind(1,sample$words2)), 1, 0.9), nrow=ncol(sample$words2)+1) 



### Fit proposed model
elmo1 <- cbind(1,sample$words2)%*%elm
elmo1 <- sigmoid(elmo1)
mod2 <- pgVB(X=cbind(modX1, sample$Psi), Z=elmo1,  
             Y=sample$R, wgt=sample$scaledWGT, Au=0.5, Bu=0.5) ## Text Model

### Get population estimates
popX <- model.matrix(~  dem_hisp + gender , data=pop)
fun1 <- function(st, n=10){
  poss <- which(sample$STATE==st)
  wordSet <- poss[sample.int(length(poss),n, replace=T, prob=1/sample$weight_full[poss])]
  return(wordSet)
}
sample$STATE <- as.numeric(sample$sample_stfips)
pop <- pop %>% left_join(unique(sample %>% select(STATE, bfID)))
betas <- rmvnorm(1000, mean=mod2$BU$mean, sigma = as.matrix(mod2$BU$Cov))
EstDF <- data.frame(STATE=unique(pop$STATE))
EstMat <- matrix(NA, nrow=48, ncol=1000)
for(i in 1:1000){
  pop2 <- pop
  pop2$est <- NA
  for(j in 1:nrow(pop)){
    tt <- fun1(pop$STATE[j], n=pop$N[j])
    temp <- (cbind(as.numeric(names(table(tt))), as.numeric(table(tt))))
    tt2 <- cbind(matrix(rep(c(popX[j,], bf[pop$bfID[j], ]), nrow(temp)), nrow=nrow(temp), byrow = T), elmo1[temp[,1], ])
    tt2 <- plogis(as.numeric(tt2%*%betas[i,]))
    pop2$est[j] <- sum(rbinom(nrow(temp), size=temp[,2], prob=tt2))
  }
  pop3 <- pop2 %>% group_by(STATE) %>% summarize(Est=sum(est)/sum(N))
  EstMat[,i] <- pop3$Est
  print(i)
}
EstDF$Est <- rowMeans(EstMat)
EstDF$SE <- apply(EstMat,1,sd)


### Get direct estimates
EstDF$DirEst <- NA
EstDF$HTVar <- NA
sample$weight_full <- sum(pop$N)*sample$weight_full/sum(sample$weight_full) ## HT weights to sum to pop total
for(st in EstDF$STATE){
  subD <- sample %>% filter(STATE==st)
  EstDF$DirEst[EstDF$STATE==st] <- horvitzThompson(as.numeric(subD$R), pi=1/subD$weight_full, var_est = T)$pop_mean
  EstDF$HTVar[EstDF$STATE==st] <- sqrt(horvitzThompson(as.numeric(subD$R), pi=1/subD$weight_full, var_est = T)$pop_mean_var)
}





#### Plot results
EstDF <- EstDF %>% left_join(unique(data.frame(fips=state.fips$fips, abb=state.fips$abb)), by=c("STATE"="fips"))
dfMean <- EstDF %>% select(-SE, -HTVar, Model=Est, Direct=DirEst) %>% pivot_longer(2:3, names_to="Estimator", values_to="Proportion")
dfMean <- map_data('state') %>% left_join(data.frame(region=tolower(state.name), abb=state.abb), by=c("region")) %>% 
  left_join(dfMean, by="abb") %>% filter(!is.na(Proportion))
ggplot(dfMean, aes(x=long, y=lat, group=group))+
  geom_polygon(color='black', size=.1, aes(fill=Proportion))+
  facet_wrap(~Estimator, nrow=2)+
  scale_fill_viridis_c(name="Estimated Proportion")+
  theme_map()+
  coord_map()+
  theme(legend.position = c(-0.3,0))



dfSE <- EstDF %>% select(-Est, -DirEst, Model=SE, Direct=HTVar) %>% pivot_longer(2:3, names_to="Estimator", values_to="SE")
dfSE <- map_data('state') %>% left_join(data.frame(region=tolower(state.name), abb=state.abb), by=c("region")) %>% 
  left_join(dfSE, by="abb") %>% filter(!is.na(SE))
p1 <- ggplot(dfSE, aes(x=long, y=lat, group=group))+
  geom_polygon(color='white', size=.1, aes(fill=SE))+
  facet_wrap(~Estimator)+
  scale_fill_viridis_c(name="Standard Error")+
  theme_map()+
  coord_map()+
  theme(legend.position = c(0,-0.3))
p2 <- ggplot(dfSE)+
  geom_density(alpha=0.3, bw=.01, fill="green", aes(x=SE))+
  facet_wrap(~Estimator)+
  theme_classic()

plot_grid(p1, p2, nrow=2)



dfSE <- EstDF %>% select(-Est, -DirEst, Model=SE, Direct=HTVar)
dfSE <- map_data('state') %>% left_join(data.frame(region=tolower(state.name), abb=state.abb), by=c("region")) %>% 
  left_join(dfSE, by="abb") %>% filter(!is.na(Model))
ggplot(dfSE, aes(x=long, y=lat, group=group))+
  geom_polygon(color='white', size=.1, aes(fill=Model/Direct))+
  scale_fill_viridis_c(name="Standard Error Ratio (Model/Direct)")+
  theme_map()+
  coord_map()


ggplot(dfSE)+
  geom_boxplot(aes(x=Model/Direct))+
  xlab("Ratio of Standard Errors (Model/Direct)")+
  theme_classic()






