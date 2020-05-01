library(dplyr)

####Generate the population####

#Set Seed for Repeatable Results
set.seed(20831748)

N<-10000 #Population size
K<-4 #Number of strata
alpha<-0.025 #95% Confidence

#Matrix to House the Results
pop<-matrix(ncol=2, nrow=N)

#Code the Strata
Nh<-c(1000, 2000, 3000, 4000)
pop[1:1000,1]<-1
pop[1001:3000,1]<-2
pop[3001:6000,1]<-3
pop[6001:N,1]<-4

#Simulate Data
for(i in 1:N){pop[i,2]<-rnorm(1, mean = sqrt(pop[i,1]) )}

#Calculate population mean
y.mu<-mean(pop[,2])

#Stratum samples sizes; 400 observations total
nh<-c(48, 81, 122, 149)

#Mean Variance
V<-aggregate(pop[,2]~pop[,1], FUN=var)[,2] %>%
  "*" ((1/nh)*(1-nh/Nh)*(Nh/N)^2) %>%
  sum()

#### Simulation ####
n.sim<-100 #100 Simulations

#This will house the simulation  results
Sim.Results<-list()
Sim.Results$Analytical<-data.frame(matrix(nrow=0, ncol=3))
Sim.Results$Jackknife<-data.frame(matrix(nrow=0, ncol=3))
Sim.Results$`Rao Wu`<-data.frame(matrix(nrow=0, ncol=3))
Sim.Results$`Mirror Match`<-data.frame(matrix(nrow=0, ncol=3))

for(l in 1:n.sim){
  
  #### Sample Observations ####
  
  #Sample observations
  Sample<-list() #It's easier to work with these algorithms if the stratum data are in a list.
  for(i in 1:4){
    df<-subset(pop, pop[,1]==i)
    Sample[[i]]<-df[sample(1:nrow(df), size=nh[i]), ]
  }
  
  #Estimate for the mean (y.bar)
  y.bar<-lapply(Sample, function(x) mean(x[,2])) %>%
    unlist() %>%
    "*"(Nh/sum(Nh)) %>% 
    sum()
  
  #### Analytical Estimator ####
  
  #Stratum Variances
  sd.yst<-lapply(Sample, function(x) var(x[,2]))%>%
    unlist()
  
  #Variance Estimator
  v.yst<-((Nh/N)^2*(1-nh/Nh)*sd.yst/nh) %>%
    sum()
  
  #C.I. Lower Bound
  CI.yst.LB<-y.bar-qnorm(p = 1-alpha)*sqrt(v.yst)
  
  #C.I. Upper Bound
  CI.yst.UB<-y.bar+qnorm(p = 1-alpha)*sqrt(v.yst)
  
  #Save Results
  Sim.Results$Analytical<-rbind(Sim.Results$Analytical,
                                c(signif(v.yst, 3), signif(CI.yst.LB,4), signif(CI.yst.UB,4)))
  
  #### Jackknife (Jones) ####
  
  #Matrix for stratum variance
  theta.h<-matrix(nrow=K, ncol=2)
  
  #Jackknife
  for(i in 1:K){
    
    #Vector to house the estimations of theta_{hi}
    theta.hi<-vector(length=nh[i])
    
    #Sample size after deleting one observation
    Nh.j<-Nh
    Nh.j[i]<-Nh.j[i]-1
    
    for(j in 1:nrow(Sample[[i]])){
      
      #Sample the jth observation deleted
      Sample.j<-Sample
      Sample.j[[i]]<-Sample.j[[i]][-j,]
      
      #Jackknife Estimate
      theta.hi[j]<-lapply(Sample.j, function(x) mean(x[,2])) %>%
        unlist() %>%
        "*"(Nh.j/sum(Nh.j)) %>% 
        sum()
    }
    
    #Stratum Mean
    theta.h[i, 1]<-mean(theta.hi)
    
    #Stratum variance
    theta.h[i, 2]<-var(theta.hi)
  }
  
  #Mean Estimator
  y.j<-mean(theta.h[,1])
  
  #Variance Estimator
  v.j<-sum((nh-1)^2*((Nh-nh)/(Nh*nh))*theta.h[,2])
  
  #C.I. Lower Bound
  CI.j.LB<-y.bar-qnorm(p = 1-alpha)*sqrt(v.j)
  
  #C.I. Upper Bound
  CI.j.UB<-y.bar+qnorm(p = 1-alpha)*sqrt(v.j)
  
  #Save Results
  Sim.Results$Jackknife<-rbind(Sim.Results$Jackknife, c(signif(v.j, 3), signif(CI.j.LB,4), signif(CI.j.UB,4)))
  
  
  #### Bootstrap (Rao and Wu) ####
  
  #Calculate the weights w_{hi}
  w <- Nh/(N*nh)
  
  #Bootstrap iterations
  B<-1000
  
  #This vector will house the bootstrap results
  theta.rw<-vector(length=B)
  
  #Boostrap
  for(b in 1:B){
    
    #Make a copy of the sample
    Sample.RW<-Sample
    
    #Boostrap sample
    for(i in 1:K){
      
      #Import sample
      Sample.rw<-Sample[[i]]
      
      #Sample with replacement
      Sample.mh<-sample(1:nh[i], size = nh[i]-3, replace = T)
      
      #Compile m*, the number of times each observation of S_h appears in the boostrap Sample.mh
      m<-vector(length = nh[i])
      for(j in 1:nh[i]){
        m[j]<-length(which(Sample.mh==j))
      } 
      
      #Bootstrap weights
      w.star<-(1-sqrt((nh[i]-3)/(nh[i]-1))+sqrt((nh[i]-3)/(nh[i]-1))*(nh[i]/(nh[i]-3))*m)*w[i]
      Sample.rw[,2]<-Sample.rw[,2]*w.star
      
      #Replace in population
      Sample.RW[[i]]<-Sample.rw
    }
    
    #Return estimators
    Sample.RW<-do.call(rbind, Sample.RW)
    theta.rw[b]<-sum(Sample.RW[,2])
  }
  
  #Variance
  v.rw<-var(theta.rw)*(B-1)/B
  
  #Jackknife variance
  v.j.rw<-0
  for(i in 1:K){
    
    #Sample size after deleting one observation
    Nh.j<-Nh
    Nh.j[i]<-Nh.j[i]-1
    
    
    theta.j.rw<-vector(length=nh[i])
    for(j in 1:nh[i]){
      
      #Copy sample
      Sample.j.rw<-Sample
      
      #Delete one observations
      Sample.j.rw[[i]]<-Sample.j.rw[[i]][-j,]
      
      #Bootstrap
      theta.j.rw[j]<-lapply(Sample.j.rw, function(x) mean(x[,2])) %>%
        unlist() %>%
        "*"(Nh.j/sum(Nh.j)) %>% 
        sum()}
    
    #Tally jackknife variance
    v.j.rw<-v.j.rw+((nh[i]-1)/nh[i])*sum((theta.j.rw-mean(theta.j.rw))^2)
  }
  
  #Calculate t-statistic using the boostrap
  t.CI<-(theta.rw-y.bar)/sqrt(v.j.rw)
  
  #C.I. Lower Bound
  CI.rw.LB<-y.bar-quantile(t.CI, 1-alpha)*sqrt(v.j.rw)
  
  #C.I. Upper Bound
  CI.rw.UB<-y.bar-quantile(t.CI, alpha)*sqrt(v.j.rw)
  
  #Save Results
  Sim.Results$`Rao Wu`<-rbind(Sim.Results$`Rao Wu`,
                              c(signif(v.rw, 3), signif(CI.rw.LB,4), signif(CI.rw.UB,4)))
  
  #### Mirror Match (Sitter) ####
  
  #Number of subsamples
  k<-table(pop[,1])/nh %>% as.numeric()
  
  #Randomization procedure, since the k are likely not integers
  k<-sapply(k, function(x) ifelse(
    ((1/x-1/ceiling(x))/(1/floor(x)-1/ceiling(x)))<runif(1), 
    floor(x), ceiling(x)))
  
  #Size of each subsample
  n.mm<-nh/k %>% as.numeric()
  
  #Bootstrap Iterations
  B<-1000
  
  #This vector will house the bootstrap results
  theta.mm<-vector(length=B)
  
  #This vector will be used to build confidence intervals
  theta.mm.ci<-vector(length=B)
  
  for(b in 1:B){
    
    #Make a copy of the sample
    Sample.MM<-Sample
    
    for(i in 1:K){
      
      #Randomization procedure, since n.mm are not integers
      n.mm[i]<-ifelse((ceiling(n.mm[i])-n.mm[i])<runif(1), 
                      floor(n.mm[i]), ceiling(n.mm[i]))
      
      #This dataframe will house the observations
      Sample.mm<-matrix(nrow=0, ncol=ncol(Sample[[i]]))
      
      #Create mirror sample by sampling with replacement
      for(j in 1:as.numeric(k[i])){
        Sample.mm<-rbind(Sample.mm, Sample[[i]][sample(1:nrow(Sample[[i]]), size = n.mm[i]),])
      }
      
      #Replace stratum sample
      Sample.MM[[i]]<-Sample.mm
    }
    
    
    #Save bootstrap estimate
    theta.mm[b]<-sapply(Sample.MM, function(x) mean(x[,2])) %>%
      "*"(lengths(Sample.MM)/2) %>%
      "/"(sum(k*n.mm)) %>%
      sum()
    
    
    #Confidence interval -  a second bootstrap, with the difference being the sample
    
    for(i in 1:K){
      
      #This dataframe will house the observations
      Sample.mm<-matrix(nrow=0, ncol=ncol(Sample.MM[[i]]))
      
      #Create mirror sample by sampling with replacement
      for(j in 1:as.numeric(k[i])){
        Sample.mm<-rbind(Sample.mm, Sample.MM[[i]][sample(1:nrow(Sample.MM[[i]]), size = n.mm[i]),])}
      
      #Replace stratum sample
      Sample.MM[[i]]<-Sample.mm}
    
    #Save confidence intervals bootstrap estimate
    theta.mm.ci[b]<-sapply(Sample.MM, function(x) mean(x[,2])) %>%
      "*"(lengths(Sample.MM)/2) %>%
      "/"(sum(k*n.mm)) %>%
      sum()
  }
  
  #Variance estimate
  v.mm<-(B-1)*var(theta.mm)/B
  
  #Variance of of the confidence interval bootstrap estimates
  v.mm.CI<-(B-1)*var(theta.mm.ci)/B
  
  #Statistic for C.I. calculation
  t.CI<-(theta.mm-y.bar)/sqrt(v.mm.CI)
  
  #C.I. Lower Bound
  CI.mm.LB<-y.bar-quantile(t.CI, 1-alpha)*sqrt(v.mm)
  
  #C.I. Upper Bound
  CI.mm.UB<-y.bar-quantile(t.CI, alpha)*sqrt(v.mm)
  
  #Save Results
  Sim.Results$`Mirror Match`<-rbind(Sim.Results$`Mirror Match`,
                                    c( signif(v.mm, 3), signif(CI.mm.LB,4), signif(CI.mm.UB,4)))
}

#### Interpretation of Simulation Results ####

#Assign names to the columns
Sim.Results<-lapply(Sim.Results, setNames, c("Variance", "95% CI (Lower Bound)", "95% CI (Upper Bound)"))

#Matrix for the Results
Results<-data.frame(matrix(nrow=4, ncol=0))
row.names(Results)<-names(Sim.Results)

#Relative bias
Results$`Relative Bias`<-sapply(Sim.Results, function(x) mean((x$Variance-V)/V) %>% signif(digits = 3))

#MSE
Results$MSE<-sapply(Sim.Results, function(x) mean((x$Variance-V)^2) %>% signif(digits = 3))

#Average Length
Results$`Average Length`<-sapply(Sim.Results, function(x) mean(x$`95% CI (Upper Bound)`-x$`95% CI (Lower Bound)`))

#Lower Tail Error Rate
Results$`Lower Tail Error`<-sapply(Sim.Results, function(x) mean(ifelse(x$`95% CI (Lower Bound)`>y.mu, 1, 0)))

#Upper Tail Error Rate
Results$`Upper Tail Error`<-sapply(Sim.Results, function(x) mean(ifelse(x$`95% CI (Upper Bound)`<y.mu, 1, 0)))

#Coverage Probability Rate
Results$`Coverage Probability`<-sapply(Sim.Results, 
                                       function(x) mean(ifelse(x$`95% CI (Upper Bound)`>y.mu &x$`95% CI (Lower Bound)`<y.mu , 1, 0)))
