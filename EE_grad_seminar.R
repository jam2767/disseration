library(mcmcplots)
####plot of all individuals####
denplot(test.SD[,2:7],collapse=TRUE,xlim=c(0,1.25),main=c("Residual",
                                                          "Individual",
                                                          "Landscape",
                                                          "Region",
                                                          "Microsite",
                                                          "Year"))

##summary stats for making barplot of median variance explained
summary.all <- summary(test.SD) 
var.all<-summary.all$quantiles[2:7,3]
#proportion of variance explained
proportion.var.all<- var.all/sum(var.all)

##graduate student seminar analysis##
# Analysis of species turnover - Bayes
#Acer rubrum
acru.qc <- which(combined.dmg.clean.qc$SPP == "ACRU")
acru.dat.qc.bayes <- combined.dmg.clean.qc[acru.qc,]


##########ACRU model
########## model just fitting the mean with year, site, plot, sub random effect
s.acru<-as.numeric(as.factor(acru.dat.qc.bayes$SITE))
### HACK to fix indexing by plots ###
s.acru<-replace(s.acru,which(s.acru==5),4)
s.acru<-replace(s.acru,which(s.acru==6),5)
s.acru<-replace(s.acru,which(s.acru==8),6)
###
p.acru<-as.numeric(as.factor(paste0(as.character(acru.dat.qc.bayes$SITE),
                                    as.character(acru.dat.qc.bayes$PLOT))))
sp.acru <- as.numeric(as.factor(paste0(as.character(acru.dat.qc.bayes$SITE),
                                       as.character(acru.dat.qc.bayes$PLOT),
                                       as.character(acru.dat.qc.bayes$SUB))))

x.acru = 1/var(as.vector(as.matrix(acru.dat.qc.bayes[,33:62])),na.rm=TRUE)
y.acru = mean(as.matrix(acru.dat.qc.bayes[,33:62]),na.rm=TRUE)
mean30yr.qc.acru = rowMeans(acru.dat.qc.bayes[33:62],na.rm=TRUE)
column.means.acru = colMeans(acru.dat.qc.bayes[33:62],na.rm=TRUE)
us.acru = tapply(X=mean30yr.qc.acru,INDEX=list(s.acru),FUN=mean,na.rm=TRUE)
up.acru = tapply(X=mean30yr.qc.acru,INDEX=list(p.acru),FUN=mean,na.rm=TRUE)
usp.acru = tapply(X=mean30yr.qc.acru,INDEX=list(sp.acru),FUN=mean,na.rm=TRUE)
ty.acru = 1/var(column.means.acru)
ts.acru = 1/var(us.acru)
#ts =  1/mean(tapply(X=mean30yr.qc,INDEX=list(combined.dmg.clean.qc$SITE),FUN=mean,na.rm=TRUE))
tp.acru= 1/var(up.acru)
#tp = 1/mean(tapply(X=mean30yr.qc,INDEX=list(combined.dmg.clean.qc$PLOT),FUN=mean,na.rm=TRUE))
tsp.acru= 1/var(usp.acru)
as.acru = us.acru - mean(us.acru)
ap.acru = up.acru - mean(up.acru) 
asp.acru =usp.acru - mean(usp.acru)
ai.acru = mean30yr.qc.acru - mean(y.acru)#as.vector(as.matrix(combined.dmg.clean.qc[,33:62])) - mean(y)
ay.acru = column.means.acru - mean(column.means.acru)

data2.acru = list(x=as.matrix(acru.dat.qc.bayes[,33:62]), nyear=30,
             nrep = 153,s=s.acru,p=p.acru,sp=sp.acru,ns=length(unique(s.acru)),np = length(unique(p.acru)),nsp=length(unique(sp.acru)))
init.acru = list(mu=y.acru,prec=x.acru,tau.y=ts.acru,tau.s=ts.acru,tau.p=tp.acru,tau.sp=tsp.acru,tau.ind=x.acru,
            alpha.site=as.acru,alpha.plot=ap.acru,alpha.sub=asp.acru,alpha.ind=ai.acru,
            alpha.y=ay.acru)

FitRand.acru = "
model {
mu ~ dnorm(0,0.001)
prec ~ dgamma(0.1,0.1)
tau.y ~ dgamma(0.1,0.1)  #prior year-effect
tau.s ~ dgamma(0.1,0.1)  #prior site
tau.p ~ dgamma(0.1,0.1)  #prior plot
tau.sp ~ dgamma(0.1,0.1) #prior subplot
tau.ind ~ dgamma(0.1,0.1) #prior individual

for(s in 1:ns){
alpha.site[s] ~ dnorm(0,tau.s)
}
for(p in 1:np){
alpha.plot[p] ~ dnorm(0,tau.p)
}
for(sp in 1:nsp){
alpha.sub[sp] ~ dnorm(0,tau.sp)
}
for(i in 1:nrep){
alpha.ind[i] ~ dnorm(0,tau.ind)
}


for(y in 1:nyear){
alpha.y[y] ~dnorm(0,tau.y)
for(i in 1:nrep){
Ex[i,y] <- mu + alpha.y[y] + alpha.site[s[i]] + alpha.plot[p[i]] + alpha.sub[sp[i]] + alpha.ind[i]
#Ex[i,y] <- mu + alpha.y[y] + alpha.site[s] + alpha.plot[p] + alpha.sub[sp] + alpha.ind[i]

x[i,y] ~dnorm(Ex[i,y],prec)
}
}
}"

## compile JAGS model
j.model.acru   <- jags.model (file = textConnection(FitRand.acru),
                              data = data2.acru,
                              inits = init.acru,
                              n.chains = 3)
## burn-in
b.acru   <- coda.samples (model = j.model.acru,
                          variable.names = c("mu","prec","tau.y","tau.s","tau.p","tau.sp","tau.ind"),
                          n.iter = 5000,
                          thin = 2)

# bt   <- coda.samples (model = j.model,
#                      variable.names = c("mu","prec","tau.y","tau.s","tau.p","tau.sp","tau.ind"),
#                      n.iter = 5000)

update.jags<-function(model,previous.coda,variable.names,iter){
  samp <- coda.samples (model = model,
                        variable.names = variable.names,
                        n.iter = iter,
                        thin=2)
  new.list <- list()
  for(i in 1:length(samp)){
    new.list[[i]] = as.mcmc(rbind(previous.coda[[i]],samp[[i]]))
  }
  bmcmc.acru <<- mcmc.list(new.list)
  #return(bmcmc)
}

update.jags(j.model.acru,b.acru,c("mu","prec","tau.y","tau.s","tau.p","tau.sp","tau.ind"),10000)
# b1<- coda.samples (model = j.model,
#                    variable.names = c("mu","prec","tau.y","tau.s","tau.p","tau.sp","tau.ind"),#,"alpha.y","alpha.site","alpha.plot","alpha.sub"),#,"alpha.ind"),
#                    n.iter = 100)#,

b.out.acru<-b.acru
bmcmc1.acru<-bmcmc.acru
for(i in 1:7){
  plot(bmcmc1.acru[,i]) ## mcmc trace and density plots
}

update.jags(j.model.acru,bmcmc1.acru,c("mu","prec","tau.y","tau.s","tau.p","tau.sp","tau.ind"),20000)
bmcmc2.acru<-bmcmc.acru
for(i in 1:7){
  plot(bmcmc2.acru[,i]) ## mcmc trace and density plots
}

update.jags(j.model.acru,bmcmc2.acru,c("mu","prec","tau.y","tau.s","tau.p","tau.sp","tau.ind"),40000)
bmcmc3.acru<-bmcmc.acru
for(i in 1:7){
  plot(bmcmc3.acru[,i]) ## mcmc trace and density plots
}

update.jags(j.model.acru,bmcmc3.acru,c("mu","prec","tau.y","tau.s","tau.p","tau.sp","tau.ind"),40000)
bmcmc4.acru<-bmcmc.acru
for(i in 1:7){
  plot(bmcmc4.acru[,i]) ## mcmc trace and density plots
}
update.jags(j.model.acru,bmcmc4.acru,c("mu","prec","tau.y","tau.s","tau.p","tau.sp","tau.ind"),50000)
bmcmc5.acru<-bmcmc.acru
# for(i in 1:7){
#   plot(bmcmc5.acru[,i]) ## mcmc trace and density plots
# }

# update.jags(j.model.acru,bmcmc5.acru,c("mu","prec","tau.y","tau.s","tau.p","tau.sp","tau.ind"),50000)
# bmcmc6.acru<-bmcmc.acru
# for(i in 1:7){
#   plot(bmcmc6.acru[,i]) ## mcmc trace and density plots
# }
#transform posteriors to standard deviations
test.sd2 <- list()
test.sd2 <- for(i in 1:3){
  test.sd2[[i]]<-as.mcmc(sqrt(1/bmcmc5.acru[[i]]))
  
  test.SD.acru <- mcmc.list(test.sd2)
}

par(mfrow=c(1,1))
for(i in 1:7){
  plot(test.SD.acru[,i],cex.lab=1.8,xlim=c(0,1.5))
}

#plots of posteriors for ACRU
denplot(test.SD.acru[,2:7],collapse=TRUE,xlim=c(0,1.25),main=c("Residual",
                                                          "Individual",
                                                          "Landscape",
                                                          "Region",
                                                          "Microsite",
                                                          "Year"))

##summary stats for making barplot of median variance explained
summary.acru <- summary(test.SD.acru) 
var.acru<-summary.acru$quantiles[2:7,3]
#proportion of variance explained
proportion.var.acru<- var.acru/sum(var.acru)
########### end model ##############

#Quercus Rubra
quru.qc <- which(combined.dmg.clean.qc$SPP == "QURU")
quru.dat.qc.bayes <- combined.dmg.clean.qc[quru.qc,]
##########QURU model
########## model just fitting the mean with year, site, plot, sub random effect
s.quru<-as.numeric(as.factor(quru.dat.qc.bayes$SITE))
### HACK to fix indexing by plots ###
s.quru<-replace(s.quru,which(s.quru==3),1)
s.quru<-replace(s.quru,which(s.quru==4),2)
s.quru<-replace(s.quru,which(s.quru==5),3)
s.quru<-replace(s.quru,which(s.quru==6),4)
s.quru<-replace(s.quru,which(s.quru==7),5)

###
p.quru<-as.numeric(as.factor(paste0(as.character(quru.dat.qc.bayes$SITE),
                                    as.character(quru.dat.qc.bayes$PLOT))))
sp.quru <- as.numeric(as.factor(paste0(as.character(quru.dat.qc.bayes$SITE),
                                       as.character(quru.dat.qc.bayes$PLOT),
                                       as.character(quru.dat.qc.bayes$SUB))))

x.quru = 1/var(as.vector(as.matrix(quru.dat.qc.bayes[,33:62])),na.rm=TRUE)
y.quru = mean(as.matrix(quru.dat.qc.bayes[,33:62]),na.rm=TRUE)
mean30yr.qc.quru = rowMeans(quru.dat.qc.bayes[33:62],na.rm=TRUE)
column.means.quru = colMeans(quru.dat.qc.bayes[33:62],na.rm=TRUE)
us.quru = tapply(X=mean30yr.qc.quru,INDEX=list(s.quru),FUN=mean,na.rm=TRUE)
up.quru = tapply(X=mean30yr.qc.quru,INDEX=list(p.quru),FUN=mean,na.rm=TRUE)
usp.quru = tapply(X=mean30yr.qc.quru,INDEX=list(sp.quru),FUN=mean,na.rm=TRUE)
ty.quru = 1/var(column.means.quru)
ts.quru = 1/var(us.quru)
#ts =  1/mean(tapply(X=mean30yr.qc,INDEX=list(combined.dmg.clean.qc$SITE),FUN=mean,na.rm=TRUE))
tp.quru= 1/var(up.quru)
#tp = 1/mean(tapply(X=mean30yr.qc,INDEX=list(combined.dmg.clean.qc$PLOT),FUN=mean,na.rm=TRUE))
tsp.quru= 1/var(usp.quru)
as.quru = us.quru - mean(us.quru)
ap.quru = up.quru - mean(up.quru) 
asp.quru =usp.quru - mean(usp.quru)
ai.quru = mean30yr.qc.quru - mean(y.quru)#as.vector(as.matrix(combined.dmg.clean.qc[,33:62])) - mean(y)
ay.quru = column.means.quru - mean(column.means.quru)

data2.quru = list(x=as.matrix(quru.dat.qc.bayes[,33:62]), nyear=30,
                  nrep = 96,s=s.quru,p=p.quru,sp=sp.quru,ns=length(unique(s.quru)),np = length(unique(p.quru)),nsp=length(unique(sp.quru)))
init.quru = list(mu=y.quru,prec=x.quru,tau.y=ts.quru,tau.s=ts.quru,tau.p=tp.quru,tau.sp=tsp.quru,tau.ind=x.quru,
                 alpha.site=as.quru,alpha.plot=ap.quru,alpha.sub=asp.quru,alpha.ind=ai.quru,
                 alpha.y=ay.quru)

FitRand.quru = "
model {
mu ~ dnorm(0,0.001)
prec ~ dgamma(0.1,0.1)
tau.y ~ dgamma(0.1,0.1)  #prior year-effect
tau.s ~ dgamma(0.1,0.1)  #prior site
tau.p ~ dgamma(0.1,0.1)  #prior plot
tau.sp ~ dgamma(0.1,0.1) #prior subplot
tau.ind ~ dgamma(0.1,0.1) #prior individual

for(s in 1:ns){
alpha.site[s] ~ dnorm(0,tau.s)
}
for(p in 1:np){
alpha.plot[p] ~ dnorm(0,tau.p)
}
for(sp in 1:nsp){
alpha.sub[sp] ~ dnorm(0,tau.sp)
}
for(i in 1:nrep){
alpha.ind[i] ~ dnorm(0,tau.ind)
}


for(y in 1:nyear){
alpha.y[y] ~dnorm(0,tau.y)
for(i in 1:nrep){
Ex[i,y] <- mu + alpha.y[y] + alpha.site[s[i]] + alpha.plot[p[i]] + alpha.sub[sp[i]] + alpha.ind[i]
#Ex[i,y] <- mu + alpha.y[y] + alpha.site[s] + alpha.plot[p] + alpha.sub[sp] + alpha.ind[i]

x[i,y] ~dnorm(Ex[i,y],prec)
}
}
}"

## compile JAGS model
j.model.quru   <- jags.model (file = textConnection(FitRand.quru),
                              data = data2.quru,
                              inits = init.quru,
                              n.chains = 3)
## burn-in
b.quru   <- coda.samples (model = j.model.quru,
                          variable.names = c("mu","prec","tau.y","tau.s","tau.p","tau.sp","tau.ind"),
                          n.iter = 5000,
                          thin = 2)

# bt   <- coda.samples (model = j.model,
#                      variable.names = c("mu","prec","tau.y","tau.s","tau.p","tau.sp","tau.ind"),
#                      n.iter = 5000)

update.jags<-function(model,previous.coda,variable.names,iter){
  samp <- coda.samples (model = model,
                        variable.names = variable.names,
                        n.iter = iter,
                        thin=2)
  new.list <- list()
  for(i in 1:length(samp)){
    new.list[[i]] = as.mcmc(rbind(previous.coda[[i]],samp[[i]]))
  }
  bmcmc.quru <<- mcmc.list(new.list)
  #return(bmcmc)
}

update.jags(j.model.quru,b.quru,c("mu","prec","tau.y","tau.s","tau.p","tau.sp","tau.ind"),10000)
# b1<- coda.samples (model = j.model,
#                    variable.names = c("mu","prec","tau.y","tau.s","tau.p","tau.sp","tau.ind"),#,"alpha.y","alpha.site","alpha.plot","alpha.sub"),#,"alpha.ind"),
#                    n.iter = 100)#,

b.out.quru<-b.quru
bmcmc1.quru<-bmcmc.quru
# for(i in 1:7){
#   plot(bmcmc1.quru[,i]) ## mcmc trace and density plots
# }

update.jags(j.model.quru,bmcmc1.quru,c("mu","prec","tau.y","tau.s","tau.p","tau.sp","tau.ind"),20000)
bmcmc2.quru<-bmcmc.quru
# for(i in 1:7){
#   plot(bmcmc2.quru[,i]) ## mcmc trace and density plots
# }

update.jags(j.model.quru,bmcmc2.quru,c("mu","prec","tau.y","tau.s","tau.p","tau.sp","tau.ind"),40000)
bmcmc3.quru<-bmcmc.quru
# for(i in 1:7){
#   plot(bmcmc3.quru[,i]) ## mcmc trace and density plots
# }

update.jags(j.model.quru,bmcmc3.quru,c("mu","prec","tau.y","tau.s","tau.p","tau.sp","tau.ind"),40000)
bmcmc4.quru<-bmcmc.quru
# for(i in 1:7){
#   plot(bmcmc4.quru[,i]) ## mcmc trace and density plots
# }
update.jags(j.model.quru,bmcmc4.quru,c("mu","prec","tau.y","tau.s","tau.p","tau.sp","tau.ind"),50000)
bmcmc5.quru<-bmcmc.quru
# for(i in 1:7){
#   plot(bmcmc5.quru[,i]) ## mcmc trace and density plots
# }

# update.jags(j.model.quru,bmcmc5.quru,c("mu","prec","tau.y","tau.s","tau.p","tau.sp","tau.ind"),50000)
# bmcmc6.quru<-bmcmc.quru
# for(i in 1:7){
#   plot(bmcmc6.quru[,i]) ## mcmc trace and density plots
# }
#transform posteriors to standard deviations
test.sd3 <- list()
test.sd3 <- for(i in 1:3){
  test.sd3[[i]]<-as.mcmc(sqrt(1/bmcmc5.quru[[i]]))
  
  test.SD.quru <- mcmc.list(test.sd3)
}

# par(mfrow=c(1,1))
# for(i in 1:7){
#   plot(test.SD.quru[,i],cex.lab=1.8,xlim=c(0,1.5))
# }

#plots of posteriors for quru
denplot(test.SD.quru[,2:7],collapse=TRUE,xlim=c(0,1.25),main=c("Residual",
                                                               "Individual",
                                                               "Landscape",
                                                               "Region",
                                                               "Microsite",
                                                               "Year"))
##summary stats for making barplot of median variance explained
summary.quru <- summary(test.SD.quru) 
var.quru<-summary.quru$quantiles[2:7,3]
#proportion of variance explained
proportion.var.quru<- var.quru/sum(var.quru)
                                                               
################ End of QURU model ##########                                                               
                                                               
#Acer saccharum
acsa.qc <- which(combined.dmg.clean.qc$SPP == "ACSA3")
acsa.dat.qc.bayes <- combined.dmg.clean.qc[acsa.qc,]

#Fagus grandifolia n = 75
fagr.qc <- which(combined.dmg.clean.qc$SPP == "FAGR")
fagr.dat.qc.bayes <- combined.dmg.clean.qc[fagr.qc,]

#Tsuga canadensis n = 156
tsca.qc <- which(combined.dmg.clean.qc$SPP == "TSCA")
tsca.dat.qc.bayes <- combined.dmg.clean.qc[tsca.qc,]



###what explains individual effects###
############ combined.dmg.clean.qc used for BAYES ##############################################
################## MOVING TO BAYES!!!! #########################################################
################################################################################################
########## model just fitting the mean with year, site, plot, sub random effect
########## simplest model to be used for poster
#############This is the correct model that works USE THIS!!!!!!!!!!! ##########################
s<-as.numeric(as.factor(combined.dmg.clean.qc$SITE))
p<-as.numeric(as.factor(paste0(as.character(combined.dmg.clean.qc$SITE),
                               as.character(combined.dmg.clean.qc$PLOT))))
sp <- as.numeric(as.factor(paste0(as.character(combined.dmg.clean.qc$SITE),
                                  as.character(combined.dmg.clean.qc$PLOT),
                                  as.character(combined.dmg.clean.qc$SUB))))

x = 1/var(as.vector(as.matrix(combined.dmg.clean.qc[,33:62])),na.rm=TRUE)
y = mean(as.matrix(combined.dmg.clean.qc[,33:62]),na.rm=TRUE)
mean30yr.qc = rowMeans(combined.dmg.clean.qc[33:62],na.rm=TRUE)
column.means = colMeans(combined.dmg.clean.qc[33:62],na.rm=TRUE)
us = tapply(X=mean30yr.qc,INDEX=list(s),FUN=mean,na.rm=TRUE)
up = tapply(X=mean30yr.qc,INDEX=list(p),FUN=mean,na.rm=TRUE)
usp = tapply(X=mean30yr.qc,INDEX=list(sp),FUN=mean,na.rm=TRUE)
ty = 1/var(column.means)
ts = 1/var(us)
#ts =  1/mean(tapply(X=mean30yr.qc,INDEX=list(combined.dmg.clean.qc$SITE),FUN=mean,na.rm=TRUE))
tp= 1/var(up)
#tp = 1/mean(tapply(X=mean30yr.qc,INDEX=list(combined.dmg.clean.qc$PLOT),FUN=mean,na.rm=TRUE))
tsp= 1/var(usp)
as = us - mean(us)
ap = up - mean(up) 
asp =usp - mean(usp)
ai = mean30yr.qc - mean(y)#as.vector(as.matrix(combined.dmg.clean.qc[,33:62])) - mean(y)
ay = column.means - mean(column.means)

data2 = list(x=as.matrix(combined.dmg.clean.qc[,33:62]), nyear=30, 
             nrep = 1380,s=s,p=p,sp=sp,ns=max(s),np = max(p),nsp=max(sp))
init = list(mu=y,prec=x,tau.y=ts,tau.s=ts,tau.p=tp,tau.sp=tsp,tau.ind=x,
            alpha.site=as,alpha.plot=ap,alpha.sub=asp,alpha.ind=ai,
            alpha.y=ay)

FitRand = "
model {
mu ~ dnorm(0,0.001)
prec ~ dgamma(0.1,0.1)
tau.y ~ dgamma(0.1,0.1)  #prior year-effect
tau.s ~ dgamma(0.1,0.1)  #prior site
tau.p ~ dgamma(0.1,0.1)  #prior plot
tau.sp ~ dgamma(0.1,0.1) #prior subplot
tau.ind ~ dgamma(0.1,0.1) #prior individual

for(s in 1:ns){
alpha.site[s] ~ dnorm(0,tau.s)
}
for(p in 1:np){
alpha.plot[p] ~ dnorm(0,tau.p)
}
for(sp in 1:nsp){
alpha.sub[sp] ~ dnorm(0,tau.sp)
}
for(i in 1:nrep){
alpha.ind[i] ~ dnorm(0,tau.ind)
}


for(y in 1:nyear){
alpha.y[y] ~dnorm(0,tau.y)
for(i in 1:nrep){
Ex[i,y] <- mu + alpha.y[y] + alpha.site[s[i]] + alpha.plot[p[i]] + alpha.sub[sp[i]] + alpha.ind[i]
#Ex[i,y] <- mu + alpha.y[y] + alpha.site[s] + alpha.plot[p] + alpha.sub[sp] + alpha.ind[i]

x[i,y] ~dnorm(Ex[i,y],prec)
}
}
}"

## compile JAGS model
j.model   <- jags.model (file = textConnection(FitRand),
                         data = data2,
                         inits = init,
                         n.chains = 3)
## burn-in
b   <- coda.samples (model = j.model,
                     variable.names = c("mu","prec","tau.y","tau.s","tau.p","tau.sp","tau.ind"),
                     n.iter = 5000,
                     thin = 2)

# bt   <- coda.samples (model = j.model,
#                      variable.names = c("mu","prec","tau.y","tau.s","tau.p","tau.sp","tau.ind"),
#                      n.iter = 5000)

update.jags<-function(model,previous.coda,variable.names,iter){
  samp <- coda.samples (model = model,
                        variable.names = variable.names,
                        n.iter = iter,
                        thin=2)
  new.list <- list()
  for(i in 1:length(samp)){
    new.list[[i]] = as.mcmc(rbind(previous.coda[[i]],samp[[i]]))
  }
  bmcmc <<- mcmc.list(new.list)
  #return(bmcmc)
}

update.jags(j.model,b,c("mu","prec","tau.y","tau.s","tau.p","tau.sp","tau.ind"),5000)
# b1<- coda.samples (model = j.model,
#                    variable.names = c("mu","prec","tau.y","tau.s","tau.p","tau.sp","tau.ind"),#,"alpha.y","alpha.site","alpha.plot","alpha.sub"),#,"alpha.ind"),
#                    n.iter = 100)#,

b.out<-b
bmcmc1<-bmcmc
update.jags(j.model,bmcmc1,c("mu","prec","tau.y","tau.s","tau.p","tau.sp","tau.ind"),1000)
bmcmc2<-bmcmc