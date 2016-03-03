library(rjags)
trees.dmg <- read.csv("/Users/josh/Dropbox/Dietze_Lab_Undergrads/JAM - Xsite/Growth/plot data/Combined_2014_Adult_Field_Data_DMG.csv")
rings <- Read_Tuscon("/Users/josh/Desktop/Windendro_All_Sites/")

combined.dmg <- matchInventoryRings(trees=trees.dmg,rings)
library(data.table)
setnames(combined.dmg,old='X.',new='STEM')
combined.dmg.clean <- cbind(combined.dmg[,1:15],combined.dmg[,27:76])

#remove cores with suspiciously large growth increments
#did by hand because couldn't figure out the indexing - fml
#write.csv(combined.dmg.clean,file="/Users/josh/Desktop/combined.dmg.clean")
combined.dmg.clean.qc <- read.csv(file="/home/jam2767/combined.dmg.clean.qc.csv")#/Users/josh/Documents/meetings/AGU_2014/combined.dmg.clean.qc.csv")
combined.dmg.clean.qc <- combined.dmg.clean.qc[,2:66]

colnames<-c("SITE","PLOT","SUB","TAG","STEM","SPP","X","Y","DBH09","DBH11","DBH12",
            "DBH13","DBH14","DMG14","DMG14_BOOLEAN","1964","1965","1966","1967","1968",
            "1969","1970","1971","1972","1973","1974","1975","1976","1977","1978","1979",
            "1980","1981","1982","1983","1984","1985","1986","1987","1988","1989","1990",
            "1991","1992","1993","1994","1995","1996","1997","1998","1999","2000","2001",
            "2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012",
            "2013")
colnames(combined.dmg.clean.qc) <- colnames#colnames(combined.dmg.clean) #get rid of stupid "x" before year

#recreate data frame to be compliant with structure needed to include time as variable (50 years)
ids.rep.new <- data.frame(matrix(NA,nrow=1,ncol=15))                    
colnames(ids.rep.new) <- colnames(combined.dmg.clean.qc[1:15])
for(i in 1:nrow(combined.dmg.clean.qc)){
  ids.rep <- combined.dmg.clean.qc[rep(i,times=50),1:15]
  ids.rep.new <- rbind(ids.rep.new,ids.rep)
}

#remove row of NAs
combined.ids <- ids.rep.new[2:nrow(ids.rep.new),]

#transpose and concatenate increments for all trees
inc.new <- vector("numeric",length=0)
for(i in 1:nrow(combined.dmg.clean.qc)){
  inc <- as.matrix(combined.dmg.clean.qc[i,16:65])
  inc.new <- as.numeric(c(inc.new,inc))
}

#add ring increments
combined.ids.incs <- cbind(combined.ids,inc.new)

#rename inc variable
setnames(combined.ids.incs,old='inc.new',new='INCREMENT')

YEAR <- seq(1964,2013,1)
YEAR <- rep(YEAR,1380)

#final version with 50yrs of growth
combined.all <- cbind(combined.ids.incs,YEAR)

#final version with 30yrs of growth
dat <- data.frame(matrix(NA,nrow=1,ncol=17)) 
colnames(dat) <- colnames(combined.all)

for(i in 1981:2010){
  combined.all.30 <- which(combined.all$YEAR == i)
  dat.new<-combined.all[combined.all.30,]
  dat <- rbind(dat.new,dat)
}

############# code for barplots of common spp ###################

mean30yr.qc = rowMeans(combined.dmg.clean.qc[33:62]) #2012-2003 average growth
combined.dmg.clean.qc.mean30yr = cbind(combined.dmg.clean.qc,mean30yr.qc)

fit.30.qc1 = aov(mean30yr.qc ~ SITE/PLOT/SUB, 
                 data = combined.dmg.clean.qc.mean30yr)
var.30.qc1 = anova(fit.30.qc1)[,2]
var.30.qc1 = var.30.qc1/sum(var.30.qc1)

color.scheme <- c(rgb(255, 255, 204, maxColorValue = 255),
                  rgb(161, 218, 180, maxColorValue = 255),
                  rgb(65, 182, 196, maxColorValue = 255),
                  rgb(44, 127, 184, maxColorValue = 255),
                  rgb(37, 52, 148, maxColorValue = 255))

#site, plot, subplot
barplot(matrix(rev(var.30.qc1),4,1),col=rev(color.scheme[c(1:3,5)]),
        ylab="",legend.text=c("Individual","Microsite","Landscape","Region"), 
        args.legend = list(x=.72,y=.31,cex = 1.5, pt.cex = 1.5))

################
# Analysis of species turnover
#Acer rubrum
acru.qc <- which(combined.dmg.clean.qc.mean30yr$SPP == "ACRU")
acru.dat.qc <- combined.dmg.clean.qc.mean30yr[acru.qc,]

fit.acru.qc = aov(mean30yr.qc ~ SITE/PLOT/SUB, data = acru.dat.qc)
var.acru.qc = anova(fit.acru.qc)[,2] #sum of squares
var.acru.qc = var.acru.qc/sum(var.acru.qc)  #total sum of squares

#site, plot, subplot
#Acer Rubrum n = 153
barplot(matrix(rev(var.acru.qc),4,1),col=rev(color.scheme[c(1:3,5)]),yaxt='n')#,
#ylab="Proportion of Variance Explained", main = "ACRU  n=154")

#Quercus Rubra
quru.qc <- which(combined.dmg.clean.qc.mean30yr$SPP == "QURU")
quru.dat.qc <- combined.dmg.clean.qc.mean30yr[quru.qc,]

fit.quru.qc = aov(mean30yr.qc ~ SITE/PLOT/SUB, data = quru.dat.qc)
var.quru.qc = anova(fit.quru.qc)[,2] #sum of squares
var.quru.qc = var.quru.qc/sum(var.quru.qc)  #total sum of squares
#site, plot, subplot
#Quercus Rubra n = 96
barplot(matrix(rev(var.quru.qc),4,1),col=rev(color.scheme[c(1:3,5)]))#,
#ylab="Proportion of Variance Explained", main = "ACRU  n=154")

#Acer saccharum
acsa.qc <- which(combined.dmg.clean.qc.mean30yr$SPP == "ACSA3")
acsa.dat.qc <- combined.dmg.clean.qc.mean30yr[acsa.qc,]

fit.acsa.qc = aov(mean30yr.qc ~ SITE/PLOT/SUB, data = acsa.dat.qc)
var.acsa.qc = anova(fit.acsa.qc)[,2] #sum of squares
var.acsa.qc = var.acsa.qc/sum(var.acsa.qc)  #total sum of squares

#site, plot, subplot
#Acer saccharum n = 115
barplot(matrix(rev(var.acsa.qc),4,1),col=rev(color.scheme[c(1:3,5)]),yaxt='n')#,
#ylab="Proportion of Variance Explained", main = "ACRU  n=154")

#Fagus grandifolia n = 75
fagr.qc <- which(combined.dmg.clean.qc.mean30yr$SPP == "FAGR")
fagr.dat.qc <- combined.dmg.clean.qc.mean30yr[fagr.qc,]

fit.fagr.qc = aov(mean30yr.qc ~ SITE/PLOT/SUB, data = fagr.dat.qc)
var.fagr.qc = anova(fit.fagr.qc)[,2] #sum of squares
var.fagr.qc = var.fagr.qc/sum(var.fagr.qc)  #total sum of squares

#site, plot, subplot
#Fagus grandifolia n = 75
barplot(matrix(rev(var.fagr.qc),4,1),col=rev(color.scheme[c(1:3,5)]),yaxt='n')#,
#ylab="Proportion of Variance Explained", main = "ACRU  n=154")

#Tsuga canadensis n = 156
tsca.qc <- which(combined.dmg.clean.qc.mean30yr$SPP == "TSCA")
tsca.dat.qc <- combined.dmg.clean.qc.mean30yr[tsca.qc,]

fit.tsca.qc = aov(mean30yr.qc ~ SITE/PLOT/SUB, data = tsca.dat.qc)
var.tsca.qc = anova(fit.tsca.qc)[,2] #sum of squares
var.tsca.qc = var.tsca.qc/sum(var.tsca.qc)  #total sum of squares

#site, plot, subplot
#Tsuga canadensis n = 156
barplot(matrix(rev(var.tsca.qc),4,1),col=rev(color.scheme[c(1:3,5)]),yaxt='n')#,
#ylab="Proportion of Variance Explained", main = "ACRU  n=154")



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
update.jags(j.model,bmcmc2,c("mu","prec","tau.y","tau.s","tau.p","tau.sp","tau.ind"),2000)
bmcmc3<-bmcmc
update.jags(j.model,bmcmc3,c("mu","prec","tau.y","tau.s","tau.p","tau.sp","tau.ind"),1000)
bmcmc4<-bmcmc
update.jags(j.model,bmcmc4,c("mu","prec","tau.y","tau.s","tau.p","tau.sp","tau.ind"),10000)
bmcmc5<-bmcmc
update.jags(j.model,bmcmc5,c("mu","prec","tau.y","tau.s","tau.p","tau.sp","tau.ind"),20000)
bmcmc6<-bmcmc
update.jags(j.model,bmcmc6,c("mu","prec","tau.y","tau.s","tau.p","tau.sp","tau.ind"),20000)
bmcmc7<-bmcmc
update.jags(j.model,bmcmc7,c("mu","prec","tau.y","tau.s","tau.p","tau.sp","tau.ind"),10000)
bmcmc8<-bmcmc
update.jags(j.model,bmcmc8,c("mu","prec","tau.y","tau.s","tau.p","tau.sp","tau.ind"),20000)
bmcmc9<-bmcmc

#transform posteriors to standard deviations
test.sd2 <- list()
test.sd2 <- for(i in 1:3){
  test.sd2[[i]]<-as.mcmc(sqrt(1/bmcmc9[[i]]))

  test.SD <- mcmc.list(test.sd2)
}

for(i in 1:7){
  plot(test.SD[,i],cex.lab=1.8,xlim=c(0,1.25))
}

str(summary(bmcmc9))
#for(i in 1:7){
  summary(bmcmc9)#[,i])    ## summary table
#}

for(i in 1:7){
  plot(bmcmc9[,i])
  #plot(sqrt(1/(bmcmc9[,i]))) ## mcmc trace and density plots
}

for(i in 1:7){
  plot(bmcmc7[,i]) ## mcmc trace and density plots
}

for(i in 1:7){
  plot(bmcmc6[,i]) ## mcmc trace and density plots
}

for(i in 1:7){
  plot(bmcmc5[,i]) ## mcmc trace and density plots
}
for(i in 1:7){
  plot(bmcmc3[,i]) ## mcmc trace and density plots
}
#bmcmc is new output
for(i in 1:7){
  plot(bmcmc2[,i]) ## mcmc trace and density plots
}
for(i in 1:7){
  plot(b.out[,i]) ## mcmc trace and density plots
}
for(i in 1:7){
autocorr.plot(b.out[,i]) ## autocorrelation
}
for(i in 1:7){
cumuplot(b.out[,i]) ## quantile plot
}
for(i in 1:7){
gelman.plot(b.out[,i])      ## GRB statistic
}
for(i in 1:7){
summary(b.out[,i])    ## summary table
}

##### done visualizing outputs #####

update.jags<-function(model,previous.coda,variable.names,iter){
  samp <- coda.samples (model = model,
                variable.names = variable.names,
                n.iter = iter)
  bmcmc <<- mcmc(do.call(rbind,c(previous.coda,samp)))
  #return(bmcmc)
}

update.jags(j.model,b,c("mu","prec","tau.y","tau.s","tau.p","tau.sp","tau.ind"),100)
# b3 <- coda.samples (model = j.model,
#                     variable.names = c("mu","prec"),#,"alpha.y","alpha.site","alpha.plot","alpha.sub"),#,"alpha.ind"),
#                     n.iter = 500)#,
# update.jags(model=j.model,previous.coda=b3,variable.names=c("mu","prec"),iter=500)

# bmcmc3 = b3
# plot(bmcmc3)



mu2 = as.data.frame(as.matrix(bmcmc2))$mu
quantile(mu2,c(0.025,0.5,0.975))



###############################################################################################
###############################################################################################
############### The JAGS models below here are wrong and don't work! Don't use!################
###############################################################################################
###############################################################################################


s<-as.numeric(as.factor(combined.dmg.clean.qc$SITE))
p<-as.numeric(as.factor(paste0(as.character(combined.dmg.clean.qc$SITE),
                               as.character(combined.dmg.clean.qc$PLOT))))
sp <- as.numeric(as.factor(paste0(as.character(combined.dmg.clean.qc$SITE),
                                  as.character(combined.dmg.clean.qc$PLOT),
                                  as.character(combined.dmg.clean.qc$SUB))))
data2 = list(x=as.matrix(combined.dmg.clean.qc[,33:62]), nyear=30, 
             nrep = 1380,s=s,p=p,sp=sp)
init = NULL

FitRand = "
model {
mu ~ dnorm(0,0.001)
prec ~ dgamma(0.1,0.1)
tau.y ~ dgamma(0.1,0.1)  #prior year-effect
tau.s ~ dgamma(0.1,0.1)  #prior site
tau.p ~ dgamma(0.1,0.1)  #prior plot
tau.sp ~ dgamma(0.1,0.1) #prior subplot
tau.ind ~ dgamma(0.1,0.1) #prior individual

for(s in 1:length(s)){
alpha.site[s] ~ dnorm(0,tau.s)
}
for(p in 1:length(p)){
alpha.plot[p] ~ dnorm(0,tau.p)
}
for(sp in 1:length(sp)){
alpha.sub[sp] ~ dnorm(0,tau.sp)
}
for(i in 1:nrep){
alpha.ind[i] ~ dnorm(0,tau.ind)
}


for(y in 1:nyear){
alpha.y[y] ~dnorm(0,tau.y)
for(i in 1:nrep){
Ex[i,y] <- mu + alpha.y[y] + alpha.site[s[i]] + alpha.plot[p[i]] + alpha.sub[sp[i]] + alpha.ind[i]

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
b2   <- coda.samples (model = j.model,
                      variable.names = c("mu","prec","alpha.y","alpha.site","alpha.plot","alpha.sub","alpha.ind"),
                      n.iter = 20000,
                      thin=100)

#init.cond2 <- function(){  ## starts each chain from the same initial conditions
#  list(mu=5,prec=2) }
#init.cond3  <- function(){  ## generates random initial conditions
#  list(mu=rnorm(1,5,2),prec=runif(1,0.5,2.0))  }

bmcmc2 = b2

plot(bmcmc2)         ## mcmc history and density plot
autocorr.plot(bmcmc2)        ## autocorrelation
cumuplot(bmcmc2)     ## quantile plot
gelman.plot(bmcmc2)      ## GRB statistic
summary(bmcmc2)      ## summary table

mu2 = as.data.frame(as.matrix(bmcmc2))$mu
quantile(mu2,c(0.025,0.5,0.975))


# #### toy examples to get stuff runnning earlier:
# library(rjags)
# #example of fitting a normal distribution
# FitNorm = "
# model {
# mu ~ dnorm(0,0.001)
# prec ~ dgamma(0.1,0.1)
# for(i in 1:n){
# x[i] ~ dnorm(mu,prec)
# }
# }"
# 
# # x = rnorm(10,3,0.5)               ## pseudo-data
# # data = list(x=x,n=10)
# # init = NULL
# 
# #my data
# data = list(x=as.matrix(combined.dmg.clean.qc[,33:62]),nrep = 1380,nyear=30 )
# #data = list(x=combined.all.30$INCREMENT)
# init = NULL
# #model just fitting the mean
# FitNorm = "
# model {
# mu ~ dnorm(0,0.001)
# prec ~ dgamma(0.1,0.1)
# for(i in 1:nrep){
# for(y in 1:nyear){
# x[i] ~ dnorm(mu,prec)
# }
# #}
# }"
# 
# ## compile JAGS model
# j.model   <- jags.model (file = textConnection(FitNorm),
#                          data = data,
#                          inits = init,
#                          n.chains = 3)
# ## burn-in
# b1   <- coda.samples (model = j.model,
#                       variable.names = c("mu","prec"),
#                       n.iter = 2000)
# 
# #init.cond2 <- function(){  ## starts each chain from the same initial conditions
# #  list(mu=5,prec=2) }
# #init.cond3  <- function(){  ## generates random initial conditions
# #  list(mu=rnorm(1,5,2),prec=runif(1,0.5,2.0))  }
# 
# bmcmc = b1
# 
# plot(bmcmc)         ## mcmc history and density plot
# autocorr.plot(bmcmc)        ## autocorrelation
# cumuplot(bmcmc)     ## quantile plot
# gelman.plot(bmcmc)      ## GRB statistic
# summary(bmcmc)      ## summary table
# 
# mu = as.data.frame(as.matrix(bmcmc))$mu
# quantile(mu,c(0.025,0.5,0.975))
# 
# 
# ########## model just fitting the mean with year random effect
# data2 = list(x=as.matrix(combined.dmg.clean.qc[,33:62]), nyear=30, nrep = 1380)
# init = NULL
# 
# FitYear = "
# model {
# mu ~ dnorm(0,0.001)
# prec ~ dgamma(0.1,0.1)
# tau.y ~ dgamma(0.1,0.1) #prior year-effect
# for(y in 1:nyear){
# alpha.y[y] ~dnorm(0,tau.y)
# Ex[y] <- mu + alpha.y[y]
# for(i in 1:nrep){
# x[i,y] ~dnorm(Ex[y],prec)
# }
# }
# }"
# 
# ## compile JAGS model
# j.model   <- jags.model (file = textConnection(FitYear),
#                          data = data2,
#                          inits = init,
#                          n.chains = 3)
# ## burn-in
# b2   <- coda.samples (model = j.model,
#                       variable.names = c("mu","prec","alpha.y"),
#                       n.iter = 5000,
#                       thin = 2)
# 
# #init.cond2 <- function(){  ## starts each chain from the same initial conditions
# #  list(mu=5,prec=2) }
# #init.cond3  <- function(){  ## generates random initial conditions
# #  list(mu=rnorm(1,5,2),prec=runif(1,0.5,2.0))  }
# 
# bmcmc2 = b2
# 
# plot(bmcmc2)         ## mcmc history and density plot
# autocorr.plot(bmcmc2)        ## autocorrelation
# cumuplot(bmcmc2)     ## quantile plot
# gelman.plot(bmcmc2)      ## GRB statistic
# summary(bmcmc2)      ## summary table
# 
# mu2 = as.data.frame(as.matrix(bmcmc2))$mu
# quantile(mu2,c(0.025,0.5,0.975))






############ combined.all.30 used for lmer (ditched) #############
# ####### ANOVAS #########
# library(optimx)
# #only year and site/plot/sub effects
# simple <- lmer(INCREMENT ~ YEAR + (1|SITE/PLOT/SUB), data=combined.all.30)
# summary(simple)
# plot((simple)) #residuals plot
# #residual var: .69032
# 
# #year and site/individual random effects only
# #model failed to converge!!
# year.test <- lmer(INCREMENT ~ YEAR + (1|SITE/PLOT/SUB/TAG), data=combined.all.30)
# summary(year.test)
# plot((year.test))
# #residual var: .38174
# 
# #year, site/individual effects, spp
# year.test1 <- lmer(INCREMENT ~ YEAR + (1|SITE/PLOT/SUB/TAG) + (1|SPP), data=combined.all.30)
# summary(year.test1)
# plot((year.test1))
# #residual var: .38174
# 
# #year,site/individual,damage
# #model failed to converge
# year.test.dam <- lmer(INCREMENT ~ YEAR + (1|SITE/PLOT/SUB/TAG) + (1|DMG14_BOOLEAN), data=combined.all.30)
# summary(year.test.dam)
# #residual var: .38174
# 
# #year, site,individ,spp,year:site interaction
# #fails to converge!!
# year.test2 <- lmer(INCREMENT ~ YEAR + (1|SITE/PLOT/SUB/TAG) + (1|SPP) + (1|YEAR:SITE), 
#                    data=combined.all.30)
# #year.test2.LBFGSB <- update(year.test2, control=lmerControl(optimizer="optimx",
# # optCtrl=list(method="L-BFGS-B")))
# summary(year.test2)
# plot((year.test2)) #plot of residuals
# #residual var: .37362
# 
# #year fixed, dbh fixed, individual, spp
# year.test3 <- lmer(INCREMENT ~ YEAR + DBH12 + (1|SITE/PLOT/SUB/TAG) + (1|SPP),data=combined.all.30)
# summary(year.test3)
# #residual var: .36997
# 
# #year/dbh fixed, site, spp rand
# year.test4 <- lmer(INCREMENT ~ YEAR + DBH12 + (1|SITE/PLOT/SUB/TAG) + (1|SPP),data=combined.all.30)
# summary(year.test4)
# plot((year.test4))
# #residual var: .36997
# 
# #year/dbh fixed, site, spp, year:site rand
# #fails to converge!!
# year.test5 <- lmer(INCREMENT ~ YEAR + DBH12 + (1|SITE/PLOT/SUB/TAG) + (1|SPP) + (1|YEAR:SITE),data=combined.all.30)
# summary(year.test5)
# plot((year.test5))
# #residual var: .36194
# 
# #year/dbh fixed, site, spp, year:site, dbh:site rand
# #fails to converge!!
# year.test6 <- lmer(INCREMENT ~ YEAR + DBH12 + (1|SITE/PLOT/SUB/TAG) + (1|SPP) + (1|YEAR:SITE) + (1:DBH12:SITE),data=combined.all.30)
# summary(year.test6)
# #residual var: .75884
# 
# #year fixed, dbh fixed, individual included
# year.test7 <- lmer(INCREMENT ~ DBH12 + YEAR + (1|SITE/PLOT/SUB/TAG) + (1|SPP),data=combined.all.30)
# summary(year.test7)
# #residual var: .36997
# 
# #year random, spp random, no year:site
# #model failed to converge!!
# year.test8 <- lmer(INCREMENT ~ DBH12 + (1|YEAR) + (1|SITE/PLOT/SUB/TAG) + (1|SPP) + (1|YEAR:SITE),data=combined.all.30)
# summary(year.test8)
# #residual var: .36187
# 
# #year fixed, dbh fixed, spp fixed
# #model failed to converge!!
# year.test9 <- lmer(INCREMENT ~ DBH12 + YEAR + SPP + (1|SITE/PLOT/SUB/TAG) + (1|YEAR:SITE),data=combined.all.30)
# summary(year.test9)
# #residual var: .36194

