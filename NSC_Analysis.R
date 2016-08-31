#### code to read and match 96-well plate data (sugar and starch concentrations)
#### and adult census data and do basic exploratory analyses

#load libraries
library(lme4)
library(ggplot2)

#helper functions
#function to find last value of vector
last <- function(x) { return( x[length(x)] ) }

#give sample size on boxplot for ggplot
give.n <- function(x){
  return(c(y = mean(x), label = length(x)))
}

#color-blind-friendly palettes
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#read data
#read Adult Census Data
adult <- read.csv("/Users/Josh/Dropbox/Dietze_Lab_Undergrads/JAM - Xsite/Field Data Entry/Data Sheets/Adult_Field_Data_JAM_MCD.csv")

#read NSC tracking data - contains sample mass used for assays
mass <- read.csv("/Users/Josh/Dropbox/NSC_Runs/All/Tracking/NSC_Mass_All_8172016.csv")

#all files starch files in one directory
data.dir.starch = "/Users/Josh/Dropbox/NSC_Runs/All/Starch"
data.dir.sugar = "/Users/Josh/Dropbox/NSC_Runs/All/Sugar"

#starch files
files.starch <- dir(data.dir.starch,"Xsite")
baseline.files.starch <- files.starch[grep(pattern = "Baseline",files.starch)]
HK.files.starch <- files.starch[grep(pattern = "HK",files.starch)]

baseline.starch = NULL
for(i in seq_along(baseline.files.starch)){
  tmp <- read.csv(file.path(data.dir.starch,baseline.files.starch[i]),sep="",colClasses="character")
  baseline.starch = rbind(baseline.starch,tmp)
}
mean.baseline.starch <- tapply(X = as.numeric(baseline.starch$Meas.), INDEX = baseline.starch$Sample, FUN = mean)
baseline.starch.id <- unique(baseline.starch$Sample)

#these files contain a run after Hexokinase is added
#measures starch content via breakdown of glucose & production of NADPH
HK.starch = NULL
for(i in seq_along(HK.files.starch)){
  tmp.HK.starch <- read.csv(file.path(data.dir.starch,HK.files.starch[i]),sep="",colClasses="character")
  HK.starch = rbind(HK.starch,tmp.HK.starch)
}
max.HK.starch <- tapply(X = as.numeric(HK.starch$Meas.), INDEX = HK.starch$Sample, FUN = max)
HK.starch.id <- unique(HK.starch$Sample)

# Check that samples finished converting all Starch. If not, flag for reanalysis
dif.HK.starch <- tapply(X = as.numeric(HK.starch$Meas.), INDEX = HK.starch$Sample, FUN = diff)

#difference between baseline and HK
total.starch <- max.HK.starch-mean.baseline.starch

##TODO use each standard curve for each individual plate
##currently only using one standard curve
#create standard curve equation

standard.abs.starch <- c(mean(total.starch[3:4]),mean(total.starch[1:2]),
                         mean(total.starch[5:6]),mean(total.starch[7:8]))
#standard.conc.starch <- c(0,0.5,1.0,1.5) #standard glucose concentrations in mmol for plate 2
standard.conc.starch <- c(0,20,40,60) #standard concentrations * 40 uL
standard.lm.starch <- lm(standard.abs.starch ~ standard.conc.starch)
summary(standard.lm.starch) #slope, intercept, R^2

#visualize standard curve
plot(standard.conc.starch,standard.abs.starch)
abline(standard.lm.starch)

#standard curve equation - calculate mmol glucose equivalents of samples
sample.conc.starch <- (total.starch - standard.lm.starch$coefficients[1])/standard.lm.starch$coefficients[2]
sample.conc.starch <- as.matrix(sample.conc.starch)
sample.ids.starch <- strsplit(row.names(sample.conc.starch),"_")

###create ID to match Core_ID in adult census data
core.id.starch <- array()
for(i in 1:length(sample.conc.starch)){
  core.id.starch[i] <- paste0(sample.ids.starch[[i]][1],
                              sample.ids.starch[[i]][2],
                              sample.ids.starch[[i]][3])
}

##bind core.id to sample matrix
sample.starch <- data.frame(core.id.starch,sample.conc.starch)
sample.agg.starch <- aggregate(sample.conc.starch~core.id.starch,sample.starch,mean)
matches.starch <- match(as.character(sample.agg.starch$core.id),as.character(adult$Core_ID14))
row.names(sample.agg.starch) <- NULL
nsc.adult.dat.starch <- data.frame(sample.agg.starch,adult[matches.starch,])

#match sample masses to absorbance
matches.mass.starch <- match(as.character(nsc.adult.dat.starch$Core_ID14),as.character(mass$Core_ID))
nsc.adult.starch.mass <- data.frame(nsc.adult.dat.starch,mass[matches.mass.starch,])

#mass column is Mass_mg
#calculate concentrations based on mass of sample used
nsc.adult.starch.mass$ug.sample <- nsc.adult.starch.mass$sample.conc.starch / 40 * 840 #micrograms/sample 40 ul drawn from 840 ul total volume of Starch Extraction
nsc.adult.starch.mass$Starch.mg.g <- nsc.adult.starch.mass$ug.sample / nsc.adult.starch.mass$Mass_mg #mg/gram of sample

############ Sugar Files ############
files.sugar <- dir(data.dir.sugar,"Xsite")
baseline.files.sugar <- files.sugar[grep(pattern = "Baseline", files.sugar)]
glucose.files.sugar <- files.sugar[grep(pattern = "Glucose", files.sugar)]
fructose.files.sugar <- files.sugar[grep(pattern = "Fructose", files.sugar)]
sucrose.files.sugar <- files.sugar[grep(pattern = "Sucrose", files.sugar)]

baseline.sugar = NULL
for(i in seq_along(baseline.files.sugar)){
  tmp <- read.csv(file.path(data.dir.sugar,baseline.files.sugar[i]),sep="",colClasses="character")
  baseline.sugar = rbind(baseline.sugar,tmp)
}
mean.baseline.sugar <- tapply(X = as.numeric(baseline.sugar$Meas.), INDEX = baseline.sugar$Sample, FUN = mean) #consider replacing mean with last
baseline.sugar.id <- unique(baseline.sugar$Sample)

#these files contain a run after Hexokinase is added
#measures sugar content via breakdown of glucose & production of NADPH
glucose.sugar = NULL
for(i in seq_along(glucose.files.sugar)){
  tmp.glucose <- read.csv(file.path(data.dir.sugar,glucose.files.sugar[i]),sep="",colClasses="character")
  glucose.sugar = rbind(glucose.sugar,tmp.glucose)
}
max.glucose.sugar <- tapply(X = as.numeric(glucose.sugar$Meas.), INDEX = glucose.sugar$Sample, FUN = max)
glucose.sugar.id <- unique(glucose.sugar$Sample)

#difference between baseline and HK (glucose content)
total.glucose <- max.glucose.sugar-mean.baseline.sugar

#these files contain  a run after hexokinase and PGI have been added
#measures sugar content via breakdown of fructose & production of NADPH
fructose.sugar = NULL
for(i in seq_along(fructose.files.sugar)){
  tmp.fructose <- read.csv(file.path(data.dir.sugar,fructose.files.sugar[i]),sep="",colClasses="character")
  fructose.sugar = rbind(fructose.sugar,tmp.fructose)
}
max.fructose.sugar <- tapply(X = as.numeric(fructose.sugar$Meas.), INDEX = fructose.sugar$Sample, FUN = max)
fructose.sugar.id <- unique(fructose.sugar$Sample)

#difference between HK and PGI (Fructose content)
last.glucose.sugar <- tapply(X = as.numeric(glucose.sugar$Meas.), INDEX = glucose.sugar$Sample, FUN = last)
total.fructose <- max.fructose.sugar-last.glucose.sugar

#these files contain  a run after hexokinase, PGI and Invertase have been added
#measures sugar content via breakdown of sucrose & production of NADPH
sucrose.sugar = NULL
for(i in seq_along(sucrose.files.sugar)){
  tmp.sucrose <- read.csv(file.path(data.dir.sugar,sucrose.files.sugar[i]),sep="",colClasses="character")
  sucrose.sugar = rbind(sucrose.sugar,tmp.sucrose)
}
max.sucrose.sugar <- tapply(X = as.numeric(sucrose.sugar$Meas.), INDEX = sucrose.sugar$Sample, FUN = max)
sucrose.sugar.id <- unique(sucrose.sugar$Sample)

#difference between HK and PGI (Fructose content)
last.fructose.sugar <- tapply(X = as.numeric(fructose.sugar$Meas.), INDEX = fructose.sugar$Sample, FUN = last)
total.sucrose <- max.sucrose.sugar-last.fructose.sugar

##TODO use each standard curve for each individual plate
##currenlty only using one standard curve
#create standard curve equation

standard.abs.sugar <- c(mean(total.glucose[3:4]),mean(total.glucose[1:2]),
                         mean(total.glucose[5:6]),mean(total.glucose[7:8]))
#standard.conc.starch <- c(0,0.5,1.0,1.5) #standard glucose concentrations in mmol for plate 2
standard.conc.sugar <- c(0,20,40,60) #standard concentrations * 40 uL
standard.lm.sugar <- lm(standard.abs.sugar ~ standard.conc.sugar)
summary(standard.lm.sugar) #slope, intercept, R^2

#visualize standard curve
plot(standard.conc.sugar,standard.abs.sugar)
abline(standard.lm.sugar)

#standard curve equation - calculate mmol glucose equivalents of samples
sample.conc.glucose <- ((total.glucose) - standard.lm.sugar$coefficients[1])/standard.lm.sugar$coefficients[2] #only 20ul of sample was used, compared to 40ul of standard
sample.conc.glucose <- as.matrix(sample.conc.glucose)
sample.ids.glucose <- strsplit(row.names(sample.conc.glucose),"_")

sample.conc.fructose <- ((total.fructose) - standard.lm.sugar$coefficients[1])/standard.lm.sugar$coefficients[2] #only 20ul of sample was used, compared to 40ul of standard
sample.conc.fructose <- as.matrix(sample.conc.fructose)
sample.ids.fructose <- strsplit(row.names(sample.conc.fructose),"_")

sample.conc.sucrose <- ((total.sucrose) - standard.lm.sugar$coefficients[1])/standard.lm.sugar$coefficients[2] #only 20ul of sample was used, compared to 40ul of standard
sample.conc.sucrose <- as.matrix(sample.conc.sucrose)
sample.ids.sucrose <- strsplit(row.names(sample.conc.sucrose),"_")

###create ID to match Core_ID in adult census data
core.id.sugar <- array()
for(i in 1:length(sample.conc.glucose)){
  core.id.sugar[i] <- paste0(sample.ids.glucose[[i]][1],
                              sample.ids.glucose[[i]][2],
                              sample.ids.glucose[[i]][3])
}

##bind core.id to sample matrix
sample.sugar <- data.frame(core.id.sugar,sample.conc.glucose,sample.conc.fructose,sample.conc.sucrose)
sample.agg.glucose <- aggregate(sample.conc.glucose ~ core.id.sugar, sample.sugar, mean)
sample.agg.fructose <- aggregate(sample.conc.fructose ~ core.id.sugar, sample.sugar, mean)
sample.agg.sucrose <- aggregate(sample.conc.sucrose ~ core.id.sugar, sample.sugar, mean)
sample.agg.sugar <- data.frame(sample.agg.glucose$core.id.sugar, sample.agg.glucose$sample.conc.glucose,
                               sample.agg.fructose$sample.conc.fructose, sample.agg.sucrose$sample.conc.sucrose)

matches.sugar <- match(as.character(sample.agg.sugar$sample.agg.glucose.core.id),as.character(adult$Core_ID14))
row.names(sample.agg.sugar) <- NULL
nsc.adult.dat.sugar <- data.frame(sample.agg.sugar,adult[matches.sugar,])

#match sample masses to absorbance
matches.mass.sugar <- match(as.character(nsc.adult.dat.sugar$Core_ID14),as.character(mass$Core_ID))
nsc.adult.sugar.mass <- data.frame(nsc.adult.dat.sugar,mass[matches.mass.sugar,])

#mass column is Mass_mg
#calculate concentrations of glucose based on mass of sample used
nsc.adult.sugar.mass$ug.glucose <- nsc.adult.sugar.mass$sample.agg.glucose.sample.conc.glucose / 20 * 845 #micrograms/sample 20 ul drawn from 845 ul total volume of Sugar Extraction
nsc.adult.sugar.mass$glucose.mg.g <- nsc.adult.sugar.mass$ug.glucose / nsc.adult.sugar.mass$Mass_mg #mg/gram of sample

#calculate concentrations of fructose based on mass of sample used
nsc.adult.sugar.mass$ug.fructose <- nsc.adult.sugar.mass$sample.agg.fructose.sample.conc.fructose / 20 * 845 #micrograms/sample 20 ul drawn from 845 ul total volume of Sugar Extraction
nsc.adult.sugar.mass$fructose.mg.g <- nsc.adult.sugar.mass$ug.fructose / nsc.adult.sugar.mass$Mass_mg #mg/gram of sample

#calculate concentrations of fructose based on mass of sample used
nsc.adult.sugar.mass$ug.sucrose <- nsc.adult.sugar.mass$sample.agg.sucrose.sample.conc.sucrose / 20 * 845 #micrograms/sample 20 ul drawn from 845 ul total volume of Sugar Extraction
nsc.adult.sugar.mass$sucrose.mg.g <- nsc.adult.sugar.mass$ug.sucrose / nsc.adult.sugar.mass$Mass_mg #mg/gram of sample

###########Combine Sugar and Starch for total NSC, far fewer samples##############
starch.only <- data.frame(nsc.adult.starch.mass$core.id.starch,nsc.adult.starch.mass$Starch.mg.g)
matches.sugar.starch <- match(as.character(nsc.adult.sugar.mass$sample.agg.glucose.core.id.sugar),
                              as.character(starch.only$nsc.adult.starch.mass.core.id.starch))
row.names(nsc.adult.sugar.mass) <- NULL
row.names(starch.only) <- NULL
nsc.adult.total.mass <- data.frame(nsc.adult.sugar.mass,starch.only[matches.sugar.starch,])
nsc.adult.total.mass$total.nsc <- mapply(nsc.adult.total.mass$nsc.adult.starch.mass.Starch.mg.g,
                                      nsc.adult.total.mass$glucose.mg.g,
                                      nsc.adult.total.mass$fructose.mg.g,
                                      nsc.adult.total.mass$sucrose.mg.g, FUN=sum, na.rm=FALSE)                                      
#rename starch column
colnames(nsc.adult.total.mass)[115] <- "core.id.starch"
colnames(nsc.adult.total.mass)[116] <- "starch.mg.g"

#remove NAs for total.nsc
nsc.adult.total.mass <- nsc.adult.total.mass[which(nsc.adult.total.mass$total.nsc != "NA"),]
#write.csv(nsc.adult.total.mass, file = "/Users/josh/Dropbox/Dissertation/CH2_NSC_Empirical/Data/nsc.adult.total.mass.csv")
#nsc.adult.total.mass$Spp <- as.factor(nsc.adult.total.mass$Spp)
########## Total NSC Exploratory Analysis ###########
#####################################################
##adding PFTs to dataframe
#functional group
###Reference: http://plants.usda.gov/java/
### Also USFS silvics manual
Early.Hardwood <- c(
              which(nsc.adult.total.mass$Spp == "ACPE"),#
              which(nsc.adult.total.mass$Spp == "AIAL"),#
              which(nsc.adult.total.mass$Spp == "BEAL2"),#
              which(nsc.adult.total.mass$Spp == "BELE"),#
              which(nsc.adult.total.mass$Spp == "BEPA"),#
              which(nsc.adult.total.mass$Spp == "BETUL"),#
              which(nsc.adult.total.mass$Spp == "COFL2"),#
              which(nsc.adult.total.mass$Spp == "JUCI"),#
              which(nsc.adult.total.mass$Spp == "JUNI"),#
              which(nsc.adult.total.mass$Spp == "LIST2"),#
              which(nsc.adult.total.mass$Spp == "LITU"),#
              which(nsc.adult.total.mass$Spp == "POGR4"),#
              which(nsc.adult.total.mass$Spp == "POTR5"),#
              which(nsc.adult.total.mass$Spp == "PRSE2"),#
              which(nsc.adult.total.mass$Spp == "ROPS"),#
              which(nsc.adult.total.mass$Spp == "SAAL5")#
              ) #14 Spp

Mid.Hardwood <- c(
                  which(nsc.adult.total.mass$Spp == "FRAM2"),#
                  which(nsc.adult.total.mass$Spp == "FRNI"),#
                  which(nsc.adult.total.mass$Spp == "FRQU"),#
                  which(nsc.adult.total.mass$Spp == "ACNE2"),#
                  which(nsc.adult.total.mass$Spp == "ACRU"),#
                  which(nsc.adult.total.mass$Spp == "CACA18"),#
                  which(nsc.adult.total.mass$Spp == "CACO15"),#
                  which(nsc.adult.total.mass$Spp == "CAGL8"),#
                  which(nsc.adult.total.mass$Spp == "CALA21"),#
                  which(nsc.adult.total.mass$Spp == "CAOV2"),#
                  which(nsc.adult.total.mass$Spp == "CATO6"),#
                  which(nsc.adult.total.mass$Spp == "CECA4"),#
                  which(nsc.adult.total.mass$Spp == "MAAC"),#
                  which(nsc.adult.total.mass$Spp == "MAFR"),#
                  which(nsc.adult.total.mass$Spp == "MAVI2"),#
                  which(nsc.adult.total.mass$Spp == "PLOC"),#
                  which(nsc.adult.total.mass$Spp == "QUAL"),#
                  which(nsc.adult.total.mass$Spp == "QUCO2"),#
                  which(nsc.adult.total.mass$Spp == "QUERC"),#
                  which(nsc.adult.total.mass$Spp == "QUIM"),#
                  which(nsc.adult.total.mass$Spp == "QULA2"),#
                  which(nsc.adult.total.mass$Spp == "QULA3"),#
                  which(nsc.adult.total.mass$Spp == "QUMO4"),#
                  which(nsc.adult.total.mass$Spp == "QUMU"),#
                  which(nsc.adult.total.mass$Spp == "QUNI"),#
                  which(nsc.adult.total.mass$Spp == "QUPH"),#
                  which(nsc.adult.total.mass$Spp == "QURU"),#
                  which(nsc.adult.total.mass$Spp == "QUST"),#
                  which(nsc.adult.total.mass$Spp == "QUVE"),#
                  which(nsc.adult.total.mass$Spp == "QUVI"),#
                  which(nsc.adult.total.mass$Spp == "ULAM")#
                  ) #23 Spp

Late.Hardwood <- c(which(nsc.adult.total.mass$Spp == "ACFL"),#
                   which(nsc.adult.total.mass$Spp == "ACSA3"),#
                   which(nsc.adult.total.mass$Spp == "AEGL"),#
                   which(nsc.adult.total.mass$Spp == "DIVI5"),#
                   which(nsc.adult.total.mass$Spp == "FAGR"),#
                   which(nsc.adult.total.mass$Spp == "MAGR4"),#
                   which(nsc.adult.total.mass$Spp == "NYSY"),#
                   which(nsc.adult.total.mass$Spp == "OXAR"),#
                   which(nsc.adult.total.mass$Spp == "TIAM")#
                   ) #9 Spp


#set PFTs
nsc.adult.total.mass$PFT <- NA
nsc.adult.total.mass$PFT[Early.Hardwood] <- "Early.Hardwood"
nsc.adult.total.mass$PFT[Mid.Hardwood] <- "Mid.Hardwood"
nsc.adult.total.mass$PFT[Late.Hardwood] <- "Late.Hardwood"


Early.Conifer <- c(
              which(nsc.adult.total.mass$Spp == "ABBA"),#
              which(nsc.adult.total.mass$Spp == "JUVI"),#
              which(nsc.adult.total.mass$Spp == "PIPA2"),#
              which(nsc.adult.total.mass$Spp == "PIST"),#
              which(nsc.adult.total.mass$Spp == "PIVI"),##should be PIVI2 for Virginia Pine
              which(nsc.adult.total.mass$Spp == "PIRE"),#
              which(nsc.adult.total.mass$Spp == "PIRI"),#
              which(nsc.adult.total.mass$Spp == "PITA")#
              ) #6 Spp

#Mid.Conifer <- c()

Late.Conifer <- c(
              which(nsc.adult.total.mass$Spp == "LALA"),#
              which(nsc.adult.total.mass$Spp == "PIGL"),#
              which(nsc.adult.total.mass$Spp == "PIMA"),#
              which(nsc.adult.total.mass$Spp == "PIRU"),#
              which(nsc.adult.total.mass$Spp == "THOC2"),#
              which(nsc.adult.total.mass$Spp == "TSCA")#
              ) #6 Spp

nsc.adult.total.mass$PFT[Early.Conifer] <- "Early.Conifer"
nsc.adult.total.mass$PFT[Late.Conifer] <- "Late.Conifer"

#adding to data frame
###Wood anatomy
###Reference http://www.wood-database.com/northern-white-cedar/
###Disagree with reference for hickory classification - put all as semi-ring porous. conferred with Christy

large.resin <- c(which(nsc.adult.total.mass$Spp == "PIPA2"),
                 which(nsc.adult.total.mass$Spp == "PIST"),
                 which(nsc.adult.total.mass$Spp == "PIVI"), #should be PIVI2 for Virginia Pine
                 which(nsc.adult.total.mass$Spp == "PIRE"),
                 which(nsc.adult.total.mass$Spp == "PIRI"),
                 which(nsc.adult.total.mass$Spp == "PITA")
                 )

small.resin <- c(which(nsc.adult.total.mass$Spp == "LALA"),
                 which(nsc.adult.total.mass$Spp == "PIGL"),
                 which(nsc.adult.total.mass$Spp == "PIMA"),
                 which(nsc.adult.total.mass$Spp == "PIRU")
                 )

no.resin <- c(which(nsc.adult.total.mass$Spp == "ABBA"),
              which(nsc.adult.total.mass$Spp == "JUVI"),
              which(nsc.adult.total.mass$Spp == "THOC2"),
              which(nsc.adult.total.mass$Spp == "TSCA")
              )

ring.porous <- c(which(nsc.adult.total.mass$Spp == "AIAL"),              
                 which(nsc.adult.total.mass$Spp == "FRAM2"),
                 which(nsc.adult.total.mass$Spp == "FRNI"),
                 which(nsc.adult.total.mass$Spp == "FRQU"),
                 which(nsc.adult.total.mass$Spp == "QUAL"),
                 which(nsc.adult.total.mass$Spp == "QUCO2"),
                 which(nsc.adult.total.mass$Spp == "QUERC"),
                 which(nsc.adult.total.mass$Spp == "QUIM"),
                 which(nsc.adult.total.mass$Spp == "QULA2"),
                 which(nsc.adult.total.mass$Spp == "QULA3"),
                 which(nsc.adult.total.mass$Spp == "QUMO4"),
                 which(nsc.adult.total.mass$Spp == "QUMU"),
                 which(nsc.adult.total.mass$Spp == "QUNI"),
                 which(nsc.adult.total.mass$Spp == "QUPH"),
                 which(nsc.adult.total.mass$Spp == "QURU"),
                 which(nsc.adult.total.mass$Spp == "QUST"),
                 which(nsc.adult.total.mass$Spp == "QUVE"),
                 which(nsc.adult.total.mass$Spp == "QUVI"),
                 which(nsc.adult.total.mass$Spp == "ROPS"),
                 which(nsc.adult.total.mass$Spp == "SAAL5")
                 )

diffuse.porous <- c(which(nsc.adult.total.mass$Spp == "ACNE2"),
                    which(nsc.adult.total.mass$Spp == "ACPE"),
                    which(nsc.adult.total.mass$Spp == "ACRU"),
                    which(nsc.adult.total.mass$Spp == "ACSA3"),
                    which(nsc.adult.total.mass$Spp == "AEGL"),
                    which(nsc.adult.total.mass$Spp == "BEAL2"),
                    which(nsc.adult.total.mass$Spp == "BELE"),
                    which(nsc.adult.total.mass$Spp == "BEPA"),
                    which(nsc.adult.total.mass$Spp == "BETUL"),
                    which(nsc.adult.total.mass$Spp == "CACA18"),
                    which(nsc.adult.total.mass$Spp == "COFL2"),
                    which(nsc.adult.total.mass$Spp == "LIST2"),
                    which(nsc.adult.total.mass$Spp == "LITU"),
                    which(nsc.adult.total.mass$Spp == "MAAC"),
                    which(nsc.adult.total.mass$Spp == "MAFR"),
                    which(nsc.adult.total.mass$Spp == "MAVI2"),                  
                    which(nsc.adult.total.mass$Spp == "PLOC"),
                    which(nsc.adult.total.mass$Spp == "FAGR"),
                    which(nsc.adult.total.mass$Spp == "MAGR4"),
                    which(nsc.adult.total.mass$Spp == "NYSY"),
                    which(nsc.adult.total.mass$Spp == "OXAR"),
                    which(nsc.adult.total.mass$Spp == "TIAM"),
                    which(nsc.adult.total.mass$Spp == "POGR4"),
                    which(nsc.adult.total.mass$Spp == "POTR5"),
                    which(nsc.adult.total.mass$Spp == "PRSE2")               
                    )

semiring.porous <- c(which(nsc.adult.total.mass$Spp == "CACO15"),
                     which(nsc.adult.total.mass$Spp == "CAGL8"),
                     which(nsc.adult.total.mass$Spp == "CALA21"),
                     which(nsc.adult.total.mass$Spp == "CAOV2"),
                     which(nsc.adult.total.mass$Spp == "CATO6"),
                     which(nsc.adult.total.mass$Spp == "CECA4"),
                     which(nsc.adult.total.mass$Spp == "DIVI5"), #bordering diffuse porous
                     which(nsc.adult.total.mass$Spp == "JUCI"),
                     which(nsc.adult.total.mass$Spp == "JUNI"),
                     which(nsc.adult.total.mass$Spp == "ULAM")                                   
                     )

#set anatomy
nsc.adult.total.mass$wood <- NA
nsc.adult.total.mass$wood[large.resin] <- "large.resin"
nsc.adult.total.mass$wood[small.resin] <- "small.resin"
nsc.adult.total.mass$wood[no.resin] <- "no.resin"
nsc.adult.total.mass$wood[ring.porous] <- "ring.porous"
nsc.adult.total.mass$wood[diffuse.porous] <- "diffuse.porous"
nsc.adult.total.mass$wood[semiring.porous] <- "semiring.porous"



#make PFT factor
#nsc.adult.total.mass$PFT <- as.factor(nsc.adult.total.mass$PFT)

#only hardwood dataframe
early.h <- which(nsc.adult.total.mass$PFT == "Early.Hardwood")
mid.h <- which(nsc.adult.total.mass$PFT == "Mid.Hardwood")
late.h <- which(nsc.adult.total.mass$PFT == "Late.Hardwood")

nsc.adult.total.mass.hard <- rbind(
                             nsc.adult.total.mass[early.h,],
                             nsc.adult.total.mass[mid.h,],
                             nsc.adult.total.mass[late.h,]
                             )
#drop unused hardwood factor levels
# nsc.adult.total.mass.hard$Spp <- droplevels(nsc.adult.total.mass.hard$Spp)
# nsc.adult.total.mass.hard$PFT <- droplevels(nsc.adult.total.mass.hard$PFT)

nsc.adult.total.mass.hard$Spp <- as.character(nsc.adult.total.mass.hard$Spp)

nsc.adult.total.mass.hard$Spp <- factor(nsc.adult.total.mass.hard$Spp, 
                  levels = c("ACPE","AIAL","BEAL2","BELE",
                  "BEPA","BETUL","JUCI","JUNI","LIST2","LITU","POGR4","POTR5",
                  "PRSE2","ROPS","SAAL5","FRAM2","FRNI","FRQU","ACNE2",
                  "ACRU","CACA18","CACO15","CAGL8","CALA21","CAOV2",
                  "CATO6","CECA4","MAAC","MAFR","QUAL","QUCO2","QUIM",
                  "QUMO4","QUMU","QUPH","QURU","QUVE","QUVI","ULAM","ACFL",
                  "ACSA3","AEGL","FAGR","NYSY","OXAR","TIAM"), ordered = TRUE)

#levels(nsc.adult.total.mass.hard$Spp)
#levels(nsc.adult.total.mass.hard$PFT) <- c("Early.Hardwood","Mid.Hardwood","Late.Hardwood")
#levels(nsc.adult.total.mass.hard$PFT)

#nsc.adult.total.mass.hard <- nsc.adult.total.mass.hard[order(nsc.adult.total.mass.hard$Spp),]
#nsc.adult.total.mass.hard <- nsc.adult.total.mass.hard[order(nsc.adult.total.mass.hard$PFT),]

#only conifer dataframe
early.c <- which(nsc.adult.total.mass$PFT == "Early.Conifer")
late.c <- which(nsc.adult.total.mass$PFT == "Late.Conifer")

nsc.adult.total.mass.conif <- rbind(
  nsc.adult.total.mass[early.c,],
  nsc.adult.total.mass[late.c,]
)

#drop unused conifer factors
# nsc.adult.total.mass.conif$Spp <- droplevels(nsc.adult.total.mass.conif$Spp)
# nsc.adult.total.mass.conif$PFT <- droplevels(nsc.adult.total.mass.conif$PFT)

nsc.adult.total.mass.conif$Spp <- as.character(nsc.adult.total.mass.conif$Spp)

#relevel factors
nsc.adult.total.mass.conif$Spp <- factor(nsc.adult.total.mass.conif$Spp, 
                              levels = c("ABBA","JUVI","PIPA2","PIST","PIVI","PIRE",
                                         "PIRI","PITA","LALA","PIGL","PIMA","PIRU",
                                         "THOC2","TSCA"), ordered = TRUE)
#########################
#########################
#Wood anatomy ggplot boxplot
ggplot(nsc.adult.total.mass, aes(x = wood, y = total.nsc, fill = wood)) +
  geom_boxplot() + theme(axis.text.x = element_text(angle=90,hjust=1,size=16)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.title.x = element_text(size=20)) +
  theme(axis.title.y = element_text(size=20)) +
  theme(legend.text = element_text(size=16)) +
  stat_summary(fun.data = give.n, geom= "text", aes(y= 155)) + 
  ylab(expression(paste("Total NSC mg g"^"-1"))) + xlab("Wood Anatomy") + ggtitle("All Wood Types") +
#   geom_segment(aes(x=0,y=med.early.conif,xend=7.5,
#                    yend=med.early.conif),linetype="longdash",color="indianred4") +
#   geom_segment(aes(x=7.5,y=med.late.conif,xend=12,
#                    yend=med.late.conif),linetype="longdash",color="darkgreen") +
  scale_fill_manual(values=rev(cbPalette))

#wood anatomy anova
#explains ~4.5% more variance than PFT
lm.wood.total <- lm(total.nsc ~ wood, data = nsc.adult.total.mass)
anova(lm.wood.total)
summary(lm.wood.total)

# ALL Conifer ggplot boxplot
med.early.conif <- median(nsc.adult.total.mass$total.nsc[early.c])
med.late.conif <- median(nsc.adult.total.mass$total.nsc[late.c])
#med.total <- median(nsc.adult.total.mass.hard$total.nsc)

ggplot(nsc.adult.total.mass.conif, aes(x = Spp, y=total.nsc, fill=PFT)) + 
  geom_boxplot() + theme(axis.text.x = element_text(angle=90,hjust=1,size=16)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.title.x = element_text(size=20)) +
  theme(axis.title.y = element_text(size=20)) +
  theme(legend.text = element_text(size=16)) +
  stat_summary(fun.data = give.n, geom= "text", aes(y= 65)) + 
  #geom_hline(aes(yintercept=med.total),linetype="longdash") +
  ylab(expression(paste("Total NSC mg g"^"-1"))) + xlab("Species") + ggtitle("All") +
  geom_segment(aes(x=0,y=med.early.conif,xend=7.5,
                   yend=med.early.conif),linetype="longdash",color="indianred4") +
  geom_segment(aes(x=7.5,y=med.late.conif,xend=12,
                   yend=med.late.conif),linetype="longdash",color="darkgreen") +
  scale_fill_manual(values=rev(cbPalette))

#species anova within all conifers
spp.conif.alive.lm <- lm(total.nsc ~ Spp, data = nsc.adult.total.mass.conif.alive)
anova(spp.conif.alive.lm)
summary(spp.conif.alive.lm)

#pft anova within all hardwoods
pft.conif.alive.lm <- lm(total.nsc ~ PFT, data = nsc.adult.total.mass.conif.alive)
anova(pft.conif.alive.lm)
summary(pft.conif.alive.lm)

#ALIVE Conifer ggplot boxplot
nsc.adult.total.mass.alive <- which(nsc.adult.total.mass$Mortality14 == "A")
conif.alive <- which(nsc.adult.total.mass.conif$Mortality14 == "A")
nsc.adult.total.mass.conif.alive <- nsc.adult.total.mass.conif[conif.alive,]

early.c.alive <- which(nsc.adult.total.mass.conif.alive$PFT == "Early.Conifer")
late.c.alive <- which(nsc.adult.total.mass.conif.alive$PFT == "Late.Conifer")

med.early.conif.alive <- median(nsc.adult.total.mass.conif.alive$total.nsc[early.c.alive])
med.late.conif.alive <- median(nsc.adult.total.mass.conif.alive$total.nsc[late.c.alive])

#Plot
ggplot(nsc.adult.total.mass.conif.alive, aes(x = Spp, y=total.nsc, fill=PFT)) + 
  geom_boxplot() + theme(axis.text.x = element_text(angle=90,hjust=1,size=16)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.title.x = element_text(size=20)) +
  theme(axis.title.y = element_text(size=20)) +
  theme(legend.text = element_text(size=16)) +
  stat_summary(fun.data = give.n, geom= "text", aes(y= 85)) + 
  #geom_hline(aes(yintercept=med.total),linetype="longdash") +
  ylab(expression(paste("Total NSC mg g"^"-1"))) + xlab("Species") + ggtitle("Alive") +
  geom_segment(aes(x=0.5,y=med.early.conif.alive,xend=7.5,
                   yend=med.early.conif.alive),linetype="longdash",color="indianred4") +
  geom_segment(aes(x=7.5,y=med.late.conif.alive,xend=10.5,
                   yend=med.late.conif.alive),linetype="longdash",color="darkgreen") +
  scale_fill_manual(values=rev(cbPalette))

#species anova within alive conifers
spp.conif.alive.lm <- lm(total.nsc ~ Spp, data = nsc.adult.total.mass.conif.alive)
anova(spp.conif.alive.lm)
summary(spp.conif.alive.lm)

#pft anova within alive hardwoods
pft.conif.alive.lm <- lm(total.nsc ~ PFT, data = nsc.adult.total.mass.conif.alive)
anova(pft.conif.alive.lm)
summary(pft.conif.alive.lm)


#STRESSED Conifer ggplot boxplot
conif.stress <- which(nsc.adult.total.mass.conif$Mortality14 == "IBS")
nsc.adult.total.mass.conif.stress <- nsc.adult.total.mass.conif[conif.stress,]

early.c.stress <- which(nsc.adult.total.mass.conif.stress$PFT == "Early.Conifer")
late.c.stress <- which(nsc.adult.total.mass.conif.stress$PFT == "Late.Conifer")

med.early.conif.stress <- median(nsc.adult.total.mass.conif.stress$total.nsc[early.c.stress])
med.late.conif.stress <- median(nsc.adult.total.mass.conif.stress$total.nsc[late.c.stress])

#Plot
ggplot(nsc.adult.total.mass.conif.stress, aes(x = Spp, y=total.nsc, fill=PFT)) + 
  geom_boxplot() + theme(axis.text.x = element_text(angle=90,hjust=1,size=16)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.title.x = element_text(size=20)) +
  theme(axis.title.y = element_text(size=20)) +
  theme(legend.text = element_text(size=16)) +
  stat_summary(fun.data = give.n, geom= "text", aes(y= 42.5)) + 
  #geom_hline(aes(yintercept=med.total),linetype="longdash") +
  ylab(expression(paste("Total NSC mg g"^"-1"))) + xlab("Species") + ggtitle("Stressed") +
  geom_segment(aes(x=0.5,y=med.early.conif.stress,xend=3.5,
                   yend=med.early.conif.stress),linetype="longdash",color="indianred4") +
  geom_segment(aes(x=3.5,y=med.late.conif.stress,xend=5.5,
                   yend=med.late.conif.stress),linetype="longdash",color="darkgreen") +
  scale_fill_manual(values=rev(cbPalette))

#species anova within stressed conifers
spp.conif.stress.lm <- lm(total.nsc ~ Spp, data = nsc.adult.total.mass.conif.stress)
anova(spp.conif.stress.lm)
summary(spp.conif.stress.lm)

#pft anova within stressed conifers
pft.conif.stress.lm <- lm(total.nsc ~ PFT, data = nsc.adult.total.mass.conif.stress)
anova(pft.conif.stress.lm)
summary(pft.conif.stress.lm)


##DEAD Conifer ggplot boxplot
##5/16 dead conifers A, M, NA in 2013
##11/16 dead conifers Dead in 2013
#16/16 dead conifers still dead in 2015
conif.dead <- which(nsc.adult.total.mass.conif$Mortality14 == "D")
nsc.adult.total.mass.conif.dead <- nsc.adult.total.mass.conif[conif.dead,]

early.c.dead <- which(nsc.adult.total.mass.conif.dead$PFT == "Early.Conifer")
late.c.dead <- which(nsc.adult.total.mass.conif.dead$PFT == "Late.Conifer")

med.early.conif.dead <- median(nsc.adult.total.mass.conif.dead$total.nsc[early.c.dead])
med.late.conif.dead <- median(nsc.adult.total.mass.conif.dead$total.nsc[late.c.dead])

#Plot
ggplot(nsc.adult.total.mass.conif.dead, aes(x = Spp, y=total.nsc, fill=PFT)) + 
  geom_boxplot() + theme(axis.text.x = element_text(angle=90,hjust=1,size=16)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.title.x = element_text(size=20)) +
  theme(axis.title.y = element_text(size=20)) +
  theme(legend.text = element_text(size=16)) +
  stat_summary(fun.data = give.n, geom= "text", aes(y= 42.5)) + 
  #geom_hline(aes(yintercept=med.total),linetype="longdash") +
  ylab(expression(paste("Total NSC mg g"^"-1"))) + xlab("Species") + ggtitle("Dead") +
  geom_segment(aes(x=0.5,y=med.early.conif.dead,xend=3.5,
                   yend=med.early.conif.dead),linetype="longdash",color="indianred4") +
  geom_segment(aes(x=3.5,y=med.late.conif.dead,xend=6.75,
                   yend=med.late.conif.dead),linetype="longdash",color="darkgreen") +
  scale_fill_manual(values=rev(cbPalette))

#species anova within dead conifers
spp.conif.dead.lm <- lm(total.nsc ~ Spp, data = nsc.adult.total.mass.conif.dead)
anova(spp.conif.dead.lm)
summary(spp.conif.dead.lm)

#pft anova within dead hardwoods
pft.conif.dead.lm <- lm(total.nsc ~ PFT, data = nsc.adult.total.mass.conif.dead)
anova(pft.conif.dead.lm)
summary(pft.conif.dead.lm)


#########################


#Species
nsc.adult.total.mass$Spp <- as.factor(nsc.adult.total.mass$Spp)
#relevel factors, changed Spp to factor 
#and set BEAL2 as ref group because close to median of all species
#nsc.adult.total.mass$Spp <- relevel(nsc.adult.total.mass$Spp, ref = "ACNE2")

#boxplot n = 757
b <- boxplot(total.nsc ~ droplevels(Spp), data = nsc.adult.total.mass,
             las=2, ylab = "Concentration mg/g") # Glucose Equivalency Units")
boxplot(total.nsc ~ droplevels(Spp), data = nsc.adult.total.mass,
        las=2, ylab = expression(paste("Total NSC mg g"^"-1")),ylim=c(0,160))
text(seq_along(droplevels(nsc.adult.total.mass$Spp)), 160, b$n, cex=.9) #b$stats[4,], b$n, cex = .75)
abline(h=median(nsc.adult.total.mass$total.nsc, na.rm = TRUE),col=2)

#Species violin plot
ggplot(nsc.adult.total.mass, aes(x=Spp, y=total.nsc)) + 
  geom_violin() + stat_summary(fun.y=median,geom="point")

lm.spp.total <- lm(total.nsc ~ Spp, data = nsc.adult.total.mass)
anova(lm.spp.total)
summary(lm.spp.total)
#plot(lm.spp.suc) #assumptions of anova not met: heteroskedasticity, 
#non-normal residuals, residuals increase with fitted values

#PFTs violin plot
nsc.adult.total.mass$PFT <- as.factor(nsc.adult.total.mass$PFT)
ggplot(nsc.adult.total.mass, aes(x=PFT, y=total.nsc)) + 
  geom_violin() + stat_summary(fun.y=median,geom="point")
lm.pft.total <- lm(total.nsc ~ PFT, data = nsc.adult.total.mass)
anova(lm.pft.total)
summary(lm.pft.total)

#Hardwood Spp plot
#boxplot
boxplot(total.nsc ~ droplevels(Spp), data = nsc.adult.total.mass.hard,
        las=2, ylab = expression(paste("Total NSC mg g"^"-1")),ylim=c(0,160))

#violin plot
# too thin, can't use
# ggplot(nsc.adult.total.mass.hard, aes(x=Spp, y=total.nsc, fill=PFT)) + 
#   geom_violin() + stat_summary(fun.y=median,geom="point") + 
#   theme(axis.text.x = element_text(angle=90,hjust=1))

##########################
# hardwood ggplot boxplot for ALL individuals
#median for total.nsc
med.early.hard <- median(nsc.adult.total.mass$total.nsc[early.h])
med.mid.hard <- median(nsc.adult.total.mass$total.nsc[mid.h])
med.late.hard <- median(nsc.adult.total.mass$total.nsc[late.h])
med.total <- median(nsc.adult.total.mass.hard$total.nsc)

ggplot(nsc.adult.total.mass.hard, aes(x = Spp, y=total.nsc, fill=PFT)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle=90,hjust=1,size=16)) + 
  #geom_point(aes(x = Spp, y=total.nsc, colour = PFT)) + #works to color points, uncomment for mean + 2SE
  #geom_point(stat = mean, mapping = aes(x = Spp, y=total.nsc, colour = PFT), na.rm = TRUE) +
  #geom_errorbar(ymax = (mean(nsc.adult.total.mass.hard$total.nsc) + sd(nsc.adult.total.mass.hard$total.nsc)), ymin = (mean(nsc.adult.total.mass.hard$total.nsc) - sd(nsc.adult.total.mass.hard$total.nsc))) + theme(axis.text.x = element_text(angle=90,hjust=1,size=16)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.title.x = element_text(size=20)) +
  theme(axis.title.y = element_text(size=20)) +
  theme(legend.text = element_text(size=16)) +
  #stat_summary(mapping = aes(x = Spp, y=total.nsc, colour = PFT), fun.y = mean, geom = "point") + #can't get mean to plot!
  #stat_summary(fun.data = "mean_cl_boot", colour = "red", size = 2)
  #stat_summary(fun.data = mean_se, mult = 2, geom = "errorbar", aes(x = Spp, y=total.nsc, colour = PFT)) + #uncomment to get mean + 2SE
  stat_summary(fun.data = give.n, geom= "text", aes(y= 160)) + 
  #geom_hline(aes(yintercept=med.total),linetype="longdash") +
  ylab(expression(paste("Total NSC mg g"^"-1"," Sapwood"))) + xlab("Species") + ggtitle("All") +
  geom_segment(aes(x=0,y=med.early.hard,xend=13.5,
                   yend=med.early.hard),linetype="longdash",color="indianred4") +
  geom_segment(aes(x=13.5,y=med.mid.hard,xend=35.5,
                   yend=med.mid.hard),linetype="longdash",color="dodgerblue4") +
  geom_segment(aes(x=35.5,y=med.late.hard,xend=40.5,
                   yend=med.late.hard),linetype="longdash",color="darkgreen") +
  scale_fill_manual(values=cbPalette)

#species anova within all hardwoods
spp.hard.lm <- lm(total.nsc ~ Spp, data = nsc.adult.total.mass.hard)
anova(spp.hard.lm)
summary(spp.hard.lm)

#pft anova within all hardwoods
pft.hard.lm <- lm(total.nsc ~ PFT, data = nsc.adult.total.mass.hard)
anova(pft.hard.lm)
summary(pft.hard.lm)

##species boxplot for ALIVE individuals only
nsc.adult.total.mass.alive <- which(nsc.adult.total.mass$Mortality14 == "A")
hard.alive <- which(nsc.adult.total.mass.hard$Mortality14 == "A")
nsc.adult.total.mass.hard.alive <- nsc.adult.total.mass.hard[hard.alive,]

early.h.alive <- which(nsc.adult.total.mass.hard.alive$PFT == "Early.Hardwood")
mid.h.alive <- which(nsc.adult.total.mass.hard.alive$PFT == "Mid.Hardwood")
late.h.alive <- which(nsc.adult.total.mass.hard.alive$PFT == "Late.Hardwood")

med.early.hard.alive <- median(nsc.adult.total.mass.hard.alive$total.nsc[early.h.alive])
med.mid.hard.alive <- median(nsc.adult.total.mass.hard.alive$total.nsc[mid.h.alive])
med.late.hard.alive <- median(nsc.adult.total.mass.hard.alive$total.nsc[late.h.alive])
med.total.alive <- median(nsc.adult.total.mass.hard$total.nsc[hard.alive])

#plot
ggplot(nsc.adult.total.mass.hard.alive, aes(x = Spp, y=total.nsc, fill=PFT)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle=90,hjust=1,size=16)) + 
  #geom_point(aes(x = Spp, y=total.nsc, colour = PFT)) + #works to color points, uncomment for mean + 2SE
  #geom_point(stat = mean, mapping = aes(x = Spp, y=total.nsc, colour = PFT), na.rm = TRUE) +
  #geom_errorbar(ymax = (mean(nsc.adult.total.mass.hard$total.nsc) + sd(nsc.adult.total.mass.hard$total.nsc)), ymin = (mean(nsc.adult.total.mass.hard$total.nsc) - sd(nsc.adult.total.mass.hard$total.nsc))) + theme(axis.text.x = element_text(angle=90,hjust=1,size=16)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.title.x = element_text(size=20)) +
  theme(axis.title.y = element_text(size=20)) +
  theme(legend.text = element_text(size=16)) +
  #stat_summary(mapping = aes(x = Spp, y=total.nsc, colour = PFT), fun.y = mean, geom = "point") + #can't get mean to plot!
  #stat_summary(fun.data = "mean_cl_boot", colour = "red", size = 2)
  #stat_summary(fun.data = mean_se, mult = 2, geom = "errorbar", aes(x = Spp, y=total.nsc, colour = PFT)) + #uncomment to get mean + 2SE
  stat_summary(fun.data = give.n, geom= "text", aes(y= 160)) + 
  #geom_hline(aes(yintercept=med.total),linetype="longdash") +
  ylab(expression(paste("Total NSC mg g"^"-1"," Sapwood"))) + xlab("Species") + ggtitle("Alive") +
  geom_segment(aes(x=0,y=med.early.hard.alive,xend=13.5,
                   yend=med.early.hard.alive),linetype="longdash",color="indianred4") +
  geom_segment(aes(x=13.5,y=med.mid.hard.alive,xend=34.5,
                   yend=med.mid.hard.alive),linetype="longdash",color="dodgerblue4") +
  geom_segment(aes(x=34.5,y=med.late.hard.alive,xend=40.5,
                   yend=med.late.hard.alive),linetype="longdash",color="darkgreen") +
  scale_fill_manual(values=cbPalette)
  
#species anova within alive hardwoods
spp.hard.alive.lm <- lm(total.nsc ~ Spp, data = nsc.adult.total.mass.hard.alive)
anova(spp.hard.alive.lm)
summary(spp.hard.alive.lm)

#pft anova within alive hardwoods
pft.hard.alive.lm <- lm(total.nsc ~ PFT, data = nsc.adult.total.mass.hard.alive)
anova(pft.hard.alive.lm)
summary(pft.hard.alive.lm)

##species boxplot for STRESSED individuals only
nsc.adult.total.mass.stress <- which(nsc.adult.total.mass$Mortality14 == "IBS")
hard.stress <- which(nsc.adult.total.mass.hard$Mortality14 == "IBS")
nsc.adult.total.mass.hard.stress <- nsc.adult.total.mass.hard[hard.stress,]

early.h.stress <- which(nsc.adult.total.mass.hard.stress$PFT == "Early.Hardwood")
mid.h.stress <- which(nsc.adult.total.mass.hard.stress$PFT == "Mid.Hardwood")
late.h.stress <- which(nsc.adult.total.mass.hard.stress$PFT == "Late.Hardwood")

med.early.hard.stress <- median(nsc.adult.total.mass.hard.stress$total.nsc[early.h.stress])
med.mid.hard.stress <- median(nsc.adult.total.mass.hard.stress$total.nsc[mid.h.stress])
med.late.hard.stress <- median(nsc.adult.total.mass.hard.stress$total.nsc[late.h.stress])
med.total.stress <- median(nsc.adult.total.mass.hard$total.nsc[hard.stress])

#Plot
ggplot(nsc.adult.total.mass.hard.stress, aes(x = Spp, y=total.nsc, fill=PFT)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle=90,hjust=1,size=16)) + 
  #geom_point(aes(x = Spp, y=total.nsc, colour = PFT)) + #works to color points, uncomment for mean + 2SE
  #geom_point(stat = mean, mapping = aes(x = Spp, y=total.nsc, colour = PFT), na.rm = TRUE) +
  #geom_errorbar(ymax = (mean(nsc.adult.total.mass.hard$total.nsc) + sd(nsc.adult.total.mass.hard$total.nsc)), ymin = (mean(nsc.adult.total.mass.hard$total.nsc) - sd(nsc.adult.total.mass.hard$total.nsc))) + theme(axis.text.x = element_text(angle=90,hjust=1,size=16)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.title.x = element_text(size=20)) +
  theme(axis.title.y = element_text(size=20)) +
  theme(legend.text = element_text(size=16)) +
  #stat_summary(mapping = aes(x = Spp, y=total.nsc, colour = PFT), fun.y = mean, geom = "point") + #can't get mean to plot!
  #stat_summary(fun.data = "mean_cl_boot", colour = "red", size = 2)
  #stat_summary(fun.data = mean_se, mult = 2, geom = "errorbar", aes(x = Spp, y=total.nsc, colour = PFT)) + #uncomment to get mean + 2SE
  stat_summary(fun.data = give.n, geom= "text", aes(y= 160)) + 
  #geom_hline(aes(yintercept=med.total),linetype="longdash") +
  ylab(expression(paste("Total NSC mg g"^"-1"," Sapwood"))) + xlab("Species") + ggtitle("Stressed") +
  geom_segment(aes(x=0,y=med.early.hard.stress,xend=9.5,
                   yend=med.early.hard.stress),linetype="longdash",color="indianred4") +
  geom_segment(aes(x=9.5,y=med.mid.hard.stress,xend=21.5,
                   yend=med.mid.hard.stress),linetype="longdash",color="dodgerblue4") +
  geom_segment(aes(x=21.5,y=med.late.hard.stress,xend=24.5,
                   yend=med.late.hard.stress),linetype="longdash",color="darkgreen") +
  scale_fill_manual(values=cbPalette)

#species anova within stressed hardwoods
spp.hard.stress.lm <- lm(total.nsc ~ Spp, data = nsc.adult.total.mass.hard.stress)
anova(spp.hard.stress.lm)
summary(spp.hard.stress.lm)

#pft anova within stressed hardwoods
pft.hard.stress.lm <- lm(total.nsc ~ PFT, data = nsc.adult.total.mass.hard.stress)
anova(pft.hard.stress.lm)
summary(pft.hard.stress.lm)

##species boxplot for DEAD individuals only
##17/39 dead hardwoods were also dead in 2013
##22/39 dead hardwoods were alive or NA in 2013
##1/39 dead hardwoods was IBS in 2015 - zombie tree!!!
nsc.adult.total.mass.dead <- which(nsc.adult.total.mass$Mortality14 == "D")
hard.dead <- which(nsc.adult.total.mass.hard$Mortality14 == "D")
nsc.adult.total.mass.hard.dead <- nsc.adult.total.mass.hard[hard.dead,]

early.h.dead <- which(nsc.adult.total.mass.hard.dead$PFT == "Early.Hardwood")
mid.h.dead <- which(nsc.adult.total.mass.hard.dead$PFT == "Mid.Hardwood")
late.h.dead <- which(nsc.adult.total.mass.hard.dead$PFT == "Late.Hardwood")

med.early.hard.dead <- median(nsc.adult.total.mass.hard.dead$total.nsc[early.h.dead])
med.mid.hard.dead <- median(nsc.adult.total.mass.hard.dead$total.nsc[mid.h.dead])
med.late.hard.dead <- median(nsc.adult.total.mass.hard.dead$total.nsc[late.h.dead])
med.total.dead <- median(nsc.adult.total.mass.hard$total.nsc[hard.dead])

#Plot
ggplot(nsc.adult.total.mass.hard.dead, aes(x = Spp, y=total.nsc, fill=PFT)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle=90,hjust=1,size=16)) + 
  #geom_point(aes(x = Spp, y=total.nsc, colour = PFT)) + #works to color points, uncomment for mean + 2SE
  #geom_point(stat = mean, mapping = aes(x = Spp, y=total.nsc, colour = PFT), na.rm = TRUE) +
  #geom_errorbar(ymax = (mean(nsc.adult.total.mass.hard$total.nsc) + sd(nsc.adult.total.mass.hard$total.nsc)), ymin = (mean(nsc.adult.total.mass.hard$total.nsc) - sd(nsc.adult.total.mass.hard$total.nsc))) + theme(axis.text.x = element_text(angle=90,hjust=1,size=16)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.title.x = element_text(size=20)) +
  theme(axis.title.y = element_text(size=20)) +
  theme(legend.text = element_text(size=16)) +
  #stat_summary(mapping = aes(x = Spp, y=total.nsc, colour = PFT), fun.y = mean, geom = "point") + #can't get mean to plot!
  #stat_summary(fun.data = "mean_cl_boot", colour = "red", size = 2)
  #stat_summary(fun.data = mean_se, mult = 2, geom = "errorbar", aes(x = Spp, y=total.nsc, colour = PFT)) + #uncomment to get mean + 2SE
  stat_summary(fun.data = give.n, geom= "text", aes(y= 70)) + 
  #geom_hline(aes(yintercept=med.total),linetype="longdash") +
  ylab(expression(paste("Total NSC mg g"^"-1"," Sapwood"))) + xlab("Species") + ggtitle("Dead") +
  geom_segment(aes(x=0,y=med.early.hard.dead,xend=7.5,
                   yend=med.early.hard.dead),linetype="longdash",color="indianred4") +
  geom_segment(aes(x=7.5,y=med.mid.hard.dead,xend=17.5,
                   yend=med.mid.hard.dead),linetype="longdash",color="dodgerblue4") +
  geom_segment(aes(x=17.5,y=med.late.hard.dead,xend=20,
                   yend=med.late.hard.dead),linetype="longdash",color="darkgreen") +
  scale_fill_manual(values=cbPalette)

#species anova within hardwoods
spp.hard.dead.lm <- lm(total.nsc ~ Spp, data = nsc.adult.total.mass.hard.dead)
anova(spp.hard.dead.lm)
summary(spp.hard.dead.lm)

#pft anova within hardwoods
pft.hard.dead.lm <- lm(total.nsc ~ PFT, data = nsc.adult.total.mass.hard.dead)
anova(pft.hard.dead.lm)
summary(pft.hard.dead.lm)

library(pwr)
#power analysis for hardwoods
pwr.f2.test( u = 39, v = 764, f2 = 0.25, sig.level = 0.05, power = ) #regression
pwr.anova.test(k = 39, n = 30, f = 0.2, sig.level = 0.05, power = ) #anova

##hardwood mean, SE plot, ggplot
tot.nsc.summary <- stat_summary(data = nsc.adult.total.mass.hard, aes(x = Spp, y=total.nsc))

#species anova within hardwoods
spp.hard.lm <- lm(total.nsc ~ Spp, data = nsc.adult.total.mass.hard)
anova(spp.hard.lm)
summary(spp.hard.lm)

#pft anova within hardwoods
pft.hard.lm <- lm(total.nsc ~ PFT, data = nsc.adult.total.mass.hard)
anova(pft.hard.lm)
summary(pft.hard.lm)

#species anova all individuals
spp.lm <- lm(total.nsc ~ Spp, data = nsc.adult.total.mass)
anova(spp.lm)
summary(spp.lm)

#pft anova all species
pft.lm <- lm(total.nsc ~ PFT, data = nsc.adult.total.mass)
anova(pft.lm)
summary(pft.lm)

#species anova conifers
spp.conif.lm <- lm(total.nsc ~ Spp, data = nsc.adult.total.mass.conif)
anova(spp.conif.lm)
summary(spp.conif.lm)

#pft anova conifers
pft.conif.lm <- lm(total.nsc ~ PFT, data = nsc.adult.total.mass.conif)
anova(pft.conif.lm)
summary(pft.conif.lm)
###########################  

#Conifer Spp plot
ggplot(nsc.adult.total.mass.conif, aes(x=Spp, y=total.nsc)) + 
  geom_violin() + stat_summary(fun.y=median,geom="point")

#Hardwood PFTs plot
ggplot(nsc.adult.total.mass.hard, aes(x=PFT, y=total.nsc)) + 
  geom_violin() + stat_summary(fun.y=median,geom="point")

#Conifer PFTs plot
ggplot(nsc.adult.total.mass.conif, aes(x=PFT, y=total.nsc)) + 
  geom_violin() + stat_summary(fun.y=median,geom="point")

#mortality
nsc.adult.total.mass$Mortality14 <- as.factor(nsc.adult.total.mass$Mortality14)
#nsc.adult.sugar.mass$Mortality14 <- as.character(nsc.adult.sugar.mass$Mortality14)
b <- boxplot(total.nsc ~ droplevels(Mortality14), data = nsc.adult.total.mass)
boxplot(total.nsc ~ droplevels(Mortality14) + PFT, data = nsc.adult.total.mass,
        ylab = expression(paste("Total NSC mg g"^"-1")),ylim=c(0,160), las=2)
text(seq_along(droplevels(nsc.adult.total.mass$Mortality14)), 160, b$n, cex=.7)
lm.mort.total <- lm(total.nsc ~ Mortality14, data = nsc.adult.total.mass)
anova(lm.mort.total)
summary(lm.mort.total)
#plot(lm.mort.total)

#Mortality Violin plot
ggplot(nsc.adult.total.mass, aes(x=Mortality14, y=total.nsc)) + 
  geom_violin() + stat_summary(fun.y=median,geom="point")


#mortality boxplot
ggplot(nsc.adult.total.mass, aes(x=Mortality14, y=total.nsc)) + 
  geom_boxplot() +
  scale_x_discrete(limits=c("A","IBS","D")) + #ignores "M" and "NA" for now and reorders plot
  stat_summary(fun.data = give.n, geom= "text", aes(y= 165)) + 
  ylab(expression(paste("Total NSC mg g"^"-1"," Sapwood"))) + xlab("Mortality Status")
  
  
  
# Spp distribution of dead
#sort dead
mort.total.sort<-which(nsc.adult.total.mass$Mortality14 == "D")
count.mort.spp <- count(nsc.adult.total.mass$Spp[mort.total.sort])
barplot(count.mort.spp[,2],names.arg=count.mort.spp[,1],las=2, ylab = "Number of Trees",
        main="Total NSC Dead")

#size distribution of alive, living, dead
#alive
nsc.adult.total.mass.alive <- which(nsc.adult.total.mass$Mortality14 == "A")
hist(nsc.adult.total.mass$DBH14[nsc.adult.total.mass.alive])
mean(nsc.adult.total.mass$DBH14[nsc.adult.total.mass.alive])
median(nsc.adult.total.mass$DBH14[nsc.adult.total.mass.alive])

#IBS
nsc.adult.total.mass.ibs <- which(nsc.adult.total.mass$Mortality14 == "IBS")
hist(nsc.adult.total.mass$DBH14[nsc.adult.total.mass.ibs])
mean(nsc.adult.total.mass$DBH14[nsc.adult.total.mass.ibs])
median(nsc.adult.total.mass$DBH14[nsc.adult.total.mass.ibs])

#Dead
nsc.adult.total.mass.dead <- which(nsc.adult.total.mass$Mortality14 == "D")
hist(nsc.adult.total.mass$DBH14[nsc.adult.total.mass.dead])
mean(nsc.adult.total.mass$DBH14[nsc.adult.total.mass.dead],na.rm=TRUE)
median(nsc.adult.total.mass$DBH14[nsc.adult.total.mass.dead],na.rm=TRUE)

#anova for size and mortality
lm.mort.size <- lm(DBH14 ~ Mortality14, data = nsc.adult.total.mass)
anova(lm.mort.size)
summary(lm.mort.size)
#Dead and IBS trees are significantly smaller, 6.7 and 3.8 cm respectively

#plot total NSC versus size
plot(nsc.adult.total.mass$DBH14, nsc.adult.total.mass$total.nsc, 
     ylab = expression(paste("Total NSC mg g"^"-1")),xlab = "DBH cm")
lm.dbh.total <- lm(total.nsc ~ DBH14, data = nsc.adult.total.mass)
abline(lm.dbh.total)
summary(lm.dbh.total) #significant but explains very little

##Canopy status and total NSC ALIVE individuals
nsc.adult.total.mass$Canopy15 <- as.factor(nsc.adult.total.mass$Canopy15)
b <- boxplot(total.nsc[nsc.adult.total.mass.alive] ~ Canopy15[nsc.adult.total.mass.alive], data = nsc.adult.total.mass)
boxplot(total.nsc[nsc.adult.total.mass.alive] ~ Canopy15[nsc.adult.total.mass.alive], data = nsc.adult.total.mass,
        ylab = expression(paste("Total NSC mg g"^"-1")),ylim=c(0,170), 
        names = c("Understory","Co-Dominant","Emergent"), main = "Alive")
text(seq_along(droplevels(nsc.adult.total.mass$Canopy15[nsc.adult.total.mass.alive])), 165, b$n, cex=.7)
lm.canopy.total.alive <- lm(total.nsc[nsc.adult.total.mass.alive] ~ Canopy15[nsc.adult.total.mass.alive], data = nsc.adult.total.mass)
anova(lm.canopy.total.alive)
summary(lm.canopy.total.alive)

##Canopy status and total NSC STRESSED individuals
nsc.adult.total.mass$Canopy15 <- as.factor(nsc.adult.total.mass$Canopy15)
b <- boxplot(total.nsc[nsc.adult.total.mass.ibs] ~ Canopy15[nsc.adult.total.mass.ibs], data = nsc.adult.total.mass)
boxplot(total.nsc[nsc.adult.total.mass.ibs] ~ Canopy15[nsc.adult.total.mass.ibs], data = nsc.adult.total.mass,
        ylab = expression(paste("Total NSC mg g"^"-1")),ylim=c(0,170), 
        names = c("Understory","Co-Dominant","Emergent"), main = "Stressed")
text(seq_along(droplevels(nsc.adult.total.mass$Canopy15[nsc.adult.total.mass.ibs])), 165, b$n, cex=.7)
lm.canopy.total.ibs <- lm(total.nsc[nsc.adult.total.mass.ibs] ~ Canopy15[nsc.adult.total.mass.ibs], data = nsc.adult.total.mass)
anova(lm.canopy.total.ibs)
summary(lm.canopy.total.ibs)

##Canopy status and total NSC DEAD individuals
nsc.adult.total.mass$Canopy15 <- as.factor(nsc.adult.total.mass$Canopy15)
b <- boxplot(total.nsc[nsc.adult.total.mass.dead] ~ Canopy15[nsc.adult.total.mass.dead], data = nsc.adult.total.mass)
boxplot(total.nsc[nsc.adult.total.mass.dead] ~ Canopy15[nsc.adult.total.mass.dead], data = nsc.adult.total.mass,
        ylab = expression(paste("Total NSC mg g"^"-1")),ylim=c(0,170), 
        names = c("Understory","Co-Dominant","Emergent"), main = "Dead")
text(seq_along(droplevels(nsc.adult.total.mass$Canopy15[nsc.adult.total.mass.dead])), 165, b$n, cex=.7)
lm.canopy.total.dead <- lm(total.nsc[nsc.adult.total.mass.dead] ~ Canopy15[nsc.adult.total.mass.dead], data = nsc.adult.total.mass)
anova(lm.canopy.total.dead)
summary(lm.canopy.total.dead)

#Site violin plot
ggplot(nsc.adult.total.mass, aes(x=Site, y=total.nsc)) + 
  geom_violin() + stat_summary(fun.y=median,geom="point")

#Canopy violin plot
nsc.adult.total.mass$Canopy15 <- as.factor(nsc.adult.total.mass$Canopy15)
ggplot(nsc.adult.total.mass, aes(x=Canopy15, y=total.nsc)) + 
  geom_violin() + stat_summary(fun.y=median,geom="point")

#pairs plots
#glucose:starch
plot(nsc.adult.total.mass$glucose.mg.g, nsc.adult.total.mass$starch.mg.g)
lm.glucose.starch.tot <- lm(glucose.mg.g ~ starch.mg.g, data = nsc.adult.total.mass)
anova(lm.glucose.starch.tot)
summary(lm.glucose.starch.tot) #R2 = .03
abline(lm.glucose.starch.tot,col=2,lwd=2)

#sucrose:starch
plot(nsc.adult.total.mass$starch.mg.g, nsc.adult.total.mass$sucrose.mg.g,
     xlab = expression(paste("Starch mg g"^"-1")),
     ylab = expression(paste("Sucrose mg g"^"-1")))
lm.sucrose.starch.tot <- lm(sucrose.mg.g ~ starch.mg.g, data = nsc.adult.total.mass)
anova(lm.sucrose.starch.tot)
summary(lm.sucrose.starch.tot) #R2 = .30
abline(lm.sucrose.starch.tot,col=2,lwd=2)

#fructose:starch
plot(nsc.adult.total.mass$fructose.mg.g, nsc.adult.total.mass$starch.mg.g)
lm.fructose.starch.tot <- lm(fructose.mg.g ~ starch.mg.g, data = nsc.adult.total.mass)
anova(lm.fructose.starch.tot)
summary(lm.fructose.starch.tot) #R2 = .08
abline(lm.fructose.starch.tot,col=2,lwd=2)

#starch:total
plot(nsc.adult.total.mass$starch.mg.g, nsc.adult.total.mass$total.nsc,
     xlab = expression(paste("Starch mg g"^"-1")),
     ylab = expression(paste("Total NSC mg g"^"-1")))
lm.starch.tot <- lm(total.nsc ~ starch.mg.g, data = nsc.adult.total.mass)
anova(lm.starch.tot)
summary(lm.starch.tot)
abline(lm.starch.tot,col=2,lwd=2)
text(x=15,y=140,labels=expression(paste("R"^"2"," = 0.89")))

#starch mortality
lm.mort.starch <- lm(starch.mg.g ~ Mortality14, data = nsc.adult.total.mass)
anova(lm.mort.starch)
summary(lm.mort.starch)

#glucose:total
plot(nsc.adult.total.mass$glucose.mg.g, nsc.adult.total.mass$total.nsc,
     xlab = expression(paste("Glucose mg g"^"-1")),
     ylab = expression(paste("Total NSC mg g"^"-1")))
lm.gluc.tot <- lm(total.nsc ~ glucose.mg.g, data = nsc.adult.total.mass)
anova(lm.gluc.tot)
summary(lm.gluc.tot) #R2 = .13
abline(lm.gluc.tot,col=2,lwd=2)

#sucrose:total
plot(nsc.adult.total.mass$sucrose.mg.g, nsc.adult.total.mass$total.nsc,
     xlab = expression(paste("Sucrose mg g"^"-1")),
     ylab = expression(paste("Total NSC mg g"^"-1")))
lm.suc.tot <- lm(total.nsc ~ sucrose.mg.g, data = nsc.adult.total.mass)
anova(lm.suc.tot)
summary(lm.suc.tot) #R2 = .46
abline(lm.suc.tot,col=2,lwd=2)

#fructose:total
plot(nsc.adult.total.mass$fructose.mg.g, nsc.adult.total.mass$total.nsc,
     xlab = expression(paste("Fructose mg g"^"-1")),
     ylab = expression(paste("Total NSC mg g"^"-1")))
lm.fruc.tot <- lm(total.nsc ~ fructose.mg.g, data = nsc.adult.total.mass)
anova(lm.fruc.tot)
summary(lm.fruc.tot) #R2 = .23
abline(lm.fruc.tot,col=2,lwd=2)

#sucrose:glucose
plot(nsc.adult.total.mass$sucrose.mg.g, nsc.adult.total.mass$glucose.mg.g,
     xlab = expression(paste("Sucrose mg g"^"-1")),
     ylab = expression(paste("Glucose mg g"^"-1")))
lm.sucrose.glucose.tot <- lm(glucose.mg.g ~ sucrose.mg.g, data = nsc.adult.total.mass)
anova(lm.sucrose.glucose.tot)
summary(lm.sucrose.glucose.tot) #R2 = .02
abline(lm.sucrose.glucose.tot,col=2,lwd=2)

#sucrose:fructose
plot(nsc.adult.total.mass$sucrose.mg.g, nsc.adult.total.mass$fructose.mg.g,
     xlab = expression(paste("Sucrose mg g"^"-1")),
     ylab = expression(paste("Fructose mg g"^"-1")))
lm.sucrose.fructose.tot <- lm(fructose.mg.g ~ sucrose.mg.g, data = nsc.adult.total.mass)
anova(lm.sucrose.fructose.tot)
summary(lm.sucrose.fructose.tot) #R2 = 0.002
abline(lm.sucrose.fructose.tot,col=2,lwd=2)

#glucose:fructose
plot(nsc.adult.total.mass$glucose.mg.g, nsc.adult.total.mass$fructose.mg.g,
     xlab = expression(paste("Glucose mg g"^"-1")),
     ylab = expression(paste("Fructose mg g"^"-1")))
lm.fructose.glucose.tot <- lm(fructose.mg.g ~ glucose.mg.g, data = nsc.adult.total.mass)
anova(lm.fructose.glucose.tot)
summary(lm.fructose.glucose.tot) #R2 = .635
abline(lm.fructose.glucose.tot,col=2,lwd=2)
#Strong positive correlation between fructose and glucose, but not the others.

############Starch ONLY exploratory Analysis##############
##########################################################
nsc.adult.starch.mass <- nsc.adult.starch.mass[which(nsc.adult.starch.mass$Starch.mg.g != "NA"),]

##adding PFTs to dataframe
#functional group
Early.Hardwood <- c(
  which(nsc.adult.starch.mass$Spp == "ACPE"),
  which(nsc.adult.starch.mass$Spp == "AIAL"),
  which(nsc.adult.starch.mass$Spp == "BEAL2"),
  which(nsc.adult.starch.mass$Spp == "BELE"),
  which(nsc.adult.starch.mass$Spp == "BEPA"),
  which(nsc.adult.starch.mass$Spp == "JUCI"),
  which(nsc.adult.starch.mass$Spp == "JUNI"),
  which(nsc.adult.starch.mass$Spp == "LIST2"),
  which(nsc.adult.starch.mass$Spp == "LITU"),
  which(nsc.adult.starch.mass$Spp == "POGR4"),
  which(nsc.adult.starch.mass$Spp == "POTR5"),
  which(nsc.adult.starch.mass$Spp == "PRSE2"),
  which(nsc.adult.starch.mass$Spp == "ROPS"),
  which(nsc.adult.starch.mass$Spp == "SAAL5")
)

Mid.Hardwood <- c(
  which(nsc.adult.starch.mass$Spp == "FRAM2"),
  which(nsc.adult.starch.mass$Spp == "FRNI"),
  which(nsc.adult.starch.mass$Spp == "FRQU"),
  which(nsc.adult.starch.mass$Spp == "ACNE2"),
  which(nsc.adult.starch.mass$Spp == "ACRU"),
  which(nsc.adult.starch.mass$Spp == "CACA18"),
  which(nsc.adult.starch.mass$Spp == "CACO15"),
  which(nsc.adult.starch.mass$Spp == "CAGL8"),
  which(nsc.adult.starch.mass$Spp == "CALA21"),
  which(nsc.adult.starch.mass$Spp == "CAOV2"),
  which(nsc.adult.starch.mass$Spp == "CATO6"),
  which(nsc.adult.starch.mass$Spp == "CECA4"),
  which(nsc.adult.starch.mass$Spp == "MAAC"),
  which(nsc.adult.starch.mass$Spp == "MAFR"),
  which(nsc.adult.starch.mass$Spp == "QUAL"),
  which(nsc.adult.starch.mass$Spp == "QUCO2"),
  which(nsc.adult.starch.mass$Spp == "QUIM"),
  which(nsc.adult.starch.mass$Spp == "QUMO4"),
  which(nsc.adult.starch.mass$Spp == "QUMU"),
  which(nsc.adult.starch.mass$Spp == "QUPH"),
  which(nsc.adult.starch.mass$Spp == "QURU"),
  which(nsc.adult.starch.mass$Spp == "QUVE"),
  which(nsc.adult.starch.mass$Spp == "ULAM")
)

Late.Hardwood <- c(which(nsc.adult.starch.mass$Spp == "ACFL"),
                   which(nsc.adult.starch.mass$Spp == "ACSA3"),
                   which(nsc.adult.starch.mass$Spp == "AEGL"),
                   which(nsc.adult.starch.mass$Spp == "FAGR"),
                   which(nsc.adult.starch.mass$Spp == "NYSY"),
                   which(nsc.adult.starch.mass$Spp == "OXAR"),
                   which(nsc.adult.starch.mass$Spp == "TIAM")
)


#set PFTs
nsc.adult.starch.mass$PFT <- NA
nsc.adult.starch.mass$PFT[Early.Hardwood] <- "Early.Hardwood"
nsc.adult.starch.mass$PFT[Mid.Hardwood] <- "Mid.Hardwood"
nsc.adult.starch.mass$PFT[Late.Hardwood] <- "Late.Hardwood"


Early.Conifer <- c(
  which(nsc.adult.starch.mass$Spp == "ABBA"),
  which(nsc.adult.starch.mass$Spp == "JUVI"),
  which(nsc.adult.starch.mass$Spp == "PIST"),
  which(nsc.adult.starch.mass$Spp == "PIVI"), #should be PIVI2 for Virginia Pine
  which(nsc.adult.starch.mass$Spp == "PIRE"),
  which(nsc.adult.starch.mass$Spp == "PIRI")
)

#Mid.Conifer <- c()

Late.Conifer <- c(
  which(nsc.adult.starch.mass$Spp == "LALA"),
  which(nsc.adult.starch.mass$Spp == "PIGL"),
  which(nsc.adult.starch.mass$Spp == "PIMA"),
  which(nsc.adult.starch.mass$Spp == "PIRU"),
  which(nsc.adult.starch.mass$Spp == "THOC2"),
  which(nsc.adult.starch.mass$Spp == "TSCA")
)

nsc.adult.starch.mass$PFT[Early.Conifer] <- "Early.Conifer"
nsc.adult.starch.mass$PFT[Late.Conifer] <- "Late.Conifer"

#make PFT factor
#nsc.adult.total.mass$PFT <- as.factor(nsc.adult.total.mass$PFT)

#only hardwood dataframe
early.h <- which(nsc.adult.starch.mass$PFT == "Early.Hardwood")
mid.h <- which(nsc.adult.starch.mass$PFT == "Mid.Hardwood")
late.h <- which(nsc.adult.starch.mass$PFT == "Late.Hardwood")

nsc.adult.starch.mass.hard <- rbind(
  nsc.adult.starch.mass[early.h,],
  nsc.adult.starch.mass[mid.h,],
  nsc.adult.starch.mass[late.h,]
)


nsc.adult.starch.mass.hard$Spp <- as.character(nsc.adult.starch.mass.hard$Spp)

nsc.adult.starch.mass.hard$Spp <- factor(nsc.adult.starch.mass.hard$Spp, 
                                  levels = c("ACPE","AIAL","BEAL2","BELE",
                                             "BEPA","JUCI","JUNI","LIST2","LITU","POGR4","POTR5",
                                             "PRSE2","ROPS","SAAL5","FRAM2","FRNI","FRQU","ACNE2",
                                             "ACRU","CACA18","CACO15","CAGL8","CALA21","CAOV2",
                                             "CATO6","CECA4","MAAC","MAFR","QUAL","QUCO2","QUIM",
                                             "QUMO4","QUMU","QUPH","QURU","QUVE","ULAM","ACSA3",
                                             "FAGR","NYSY","OXAR","TIAM"), ordered = TRUE)

# #drop unused hardwood factor levels
# nsc.adult.starch.mass.hard$Spp <- droplevels(nsc.adult.starch.mass.hard$Spp)
# nsc.adult.starch.mass.hard$PFT <- droplevels(nsc.adult.starch.mass.hard$PFT)

#only conifer dataframe
early.c <- which(nsc.adult.starch.mass$PFT == "Early.Conifer")
late.c <- which(nsc.adult.starch.mass$PFT == "Late.Conifer")

nsc.adult.starch.mass.conif <- rbind(
  nsc.adult.starch.mass[early.c,],
  nsc.adult.starch.mass[late.c,]
)



nsc.adult.starch.mass.conif$Spp <- as.character(nsc.adult.starch.mass.conif$Spp)

#relevel factors
nsc.adult.starch.mass.conif$Spp <- factor(nsc.adult.starch.mass.conif$Spp, 
                                         levels = c("ABBA","JUVI","PIST","PIVI","PIRE",
                                                    "PIRI","LALA","PIGL","PIMA","PIRU",
                                                    "THOC2","TSCA"), ordered = TRUE)

# #drop unused conifer factors
# nsc.adult.starch.mass.conif$Spp <- droplevels(nsc.adult.starch.mass.conif$Spp)
# nsc.adult.starch.mass.conif$PFT <- droplevels(nsc.adult.starch.mass.conif$PFT)

# hardwood ggplot boxplot
#median for total.nsc
med.early.hard.starch <- median(nsc.adult.starch.mass$Starch.mg.g[early.h])
med.mid.hard.starch <- median(nsc.adult.starch.mass$Starch.mg.g[mid.h])
med.late.hard.starch <- median(nsc.adult.starch.mass$Starch.mg.g[late.h])
med.total.starch <- median(nsc.adult.starch.mass.hard$Starch.mg.g)

ggplot(nsc.adult.starch.mass.hard, aes(x = Spp, y=Starch.mg.g, fill=PFT)) + 
  geom_boxplot() + theme(axis.text.x = element_text(angle=90,hjust=1,size=16)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.title.x = element_text(size=20)) +
  theme(axis.title.y = element_text(size=20)) +
  theme(legend.text = element_text(size=16)) +
  stat_summary(fun.data = give.n, geom= "text", aes(y= 100)) + 
  #geom_hline(aes(yintercept=med.total),linetype="longdash") +
  ylab(expression(paste("Starch mg g"^"-1"," Sapwood"))) + xlab("Species") +
  geom_segment(aes(x=0,y=med.early.hard.starch,xend=14.5,
                   yend=med.early.hard.starch),linetype="longdash",color="grey30") +
  geom_segment(aes(x=14.5,y=med.mid.hard.starch,xend=36.3,
                   yend=med.mid.hard.starch),linetype="longdash",color="dodgerblue1") +
  geom_segment(aes(x=36.5,y=med.late.hard.starch,xend=41.5,
                   yend=med.late.hard.starch),linetype="longdash",color="darkgoldenrod4") +
  scale_fill_manual(values=cbPalette)


#ggplot boxplots of starch species
#Conifer ggplot boxplot
med.early.conif.starch <- median(nsc.adult.starch.mass$Starch.mg.g[early.c])
med.late.conif.starch <- median(nsc.adult.starch.mass$Starch.mg.g[late.c])
#med.total <- median(nsc.adult.total.mass.hard$total.nsc)

ggplot(nsc.adult.starch.mass.conif, aes(x = Spp, y=Starch.mg.g, fill=PFT)) + 
  geom_boxplot() + theme(axis.text.x = element_text(angle=90,hjust=1,size=16)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.title.x = element_text(size=20)) +
  theme(axis.title.y = element_text(size=20)) +
  theme(legend.text = element_text(size=16)) +
  stat_summary(fun.data = give.n, geom= "text", aes(y= 36)) + 
  #geom_hline(aes(yintercept=med.total),linetype="longdash") +
  ylab(expression(paste("Starch NSC mg g"^"-1"))) + xlab("Species") +
  geom_segment(aes(x=0,y=med.early.conif.starch,xend=6.5,
                   yend=med.early.conif.starch),linetype="longdash",color="indianred4") +
  geom_segment(aes(x=6.6,y=med.late.conif.starch,xend=12.5,
                   yend=med.late.conif.starch),linetype="longdash",color="darkgreen") +
  scale_fill_manual(values=rev(cbPalette))

#species anova within hardwoods
spp.hard.starch.lm <- lm(Starch.mg.g ~ Spp, data = nsc.adult.starch.mass.hard)
anova(spp.hard.starch.lm)
summary(spp.hard.starch.lm)

#pft anova within hardwoods
pft.hard.starch.lm <- lm(Starch.mg.g ~ PFT, data = nsc.adult.starch.mass.hard)
anova(pft.hard.starch.lm)
summary(pft.hard.starch.lm)

#species anova all individuals
spp.starch.lm <- lm(Starch.mg.g ~ Spp, data = nsc.adult.starch.mass)
anova(spp.starch.lm)
summary(spp.starch.lm)

#pft anova all species
pft.starch.lm <- lm(Starch.mg.g ~ PFT, data = nsc.adult.starch.mass)
anova(pft.starch.lm)
summary(pft.starch.lm)

#species anova conifers
spp.conif.starch.lm <- lm(Starch.mg.g ~ Spp, data = nsc.adult.starch.mass.conif)
anova(spp.conif.lm)
summary(spp.conif.lm)

#pft anova conifers
pft.conif.starch.lm <- lm(Starch.mg.g ~ PFT, data = nsc.adult.starch.mass.conif)
anova(pft.conif.lm)
summary(pft.conif.lm)

##End of starch only spp plots and anovas
############################################


#relevel factors, changed Spp to factor 
#and set BEAL2 as ref group because close to median of all species
nsc.adult.starch.mass$Spp <- relevel(nsc.adult.starch.mass$Spp, ref = "BEAL2")

#boxplot
b <- boxplot(Starch.mg.g ~ droplevels(Spp), data = nsc.adult.starch.mass,
        las=2)
boxplot(Starch.mg.g ~ droplevels(Spp), data = nsc.adult.starch.mass,
        las=2, ylab = expression(paste("Starch mg g"^"-1")),ylim=c(0,100))
text(seq_along(droplevels(nsc.adult.sugar.mass$Spp)), 99, b$n, cex=.9)
abline(h=median(nsc.adult.starch.mass$Starch.mg.g, na.rm = TRUE),col=2)

#species violin plot 
ggplot(nsc.adult.starch.mass, aes(x=Spp, y=Starch.mg.g)) + 
  geom_violin() + stat_summary(fun.y=median,geom="point")

lm.spp.starch <- lm(Starch.mg.g ~ Spp, data = nsc.adult.starch.mass)
anova(lm.spp.starch)
summary(lm.spp.starch)
#plot(lm.spp) #assumptions of anova not met: heteroskedasticity, 
#non-normal residuals, residuals increase with fitted values

# Site n = 1401
nsc.adult.starch.mass$Site <- relevel(nsc.adult.starch.mass$Site, ref = "D")
b <- boxplot(Starch.mg.g ~ droplevels(Site), data = nsc.adult.starch.mass)
boxplot(Starch.mg.g ~ droplevels(Site), data = nsc.adult.starch.mass,
        ylab = expression(paste("Starch mg g"^"-1")),ylim=c(-10,100))
text(seq_along(droplevels(nsc.adult.sugar.mass$Site)), 99, b$n, cex=.7)
abline(h=median(nsc.adult.starch.mass$Starch.mg.g, na.rm = TRUE))
lm.site <- lm(Starch.mg.g ~ Site, data = nsc.adult.starch.mass)
anova(lm.site)
summary(lm.site)
#plot(lm.site)

#site violin plot
ggplot(nsc.adult.starch.mass, aes(x=Site, y=Starch.mg.g)) + 
  geom_violin() + stat_summary(fun.y=median,geom="point")

#NSCs versus DBH in 2014 n = 1403
plot(nsc.adult.starch.mass$DBH14, nsc.adult.starch.mass$Starch.mg.g, 
     ylab = expression(paste("Starch mg g"^"-1")),xlab = "DBH cm")
lm.dbh <- lm(Starch.mg.g ~ DBH14, data = nsc.adult.starch.mass)
abline(lm.dbh)
summary(lm.dbh) #significant but explains very little

#NSCs and Mortality Status n = 1403
nsc.adult.starch.mass$Mortality14 <- as.factor(nsc.adult.starch.mass$Mortality14)
b <- boxplot(Starch.mg.g ~ Mortality14, data = nsc.adult.starch.mass)
boxplot(Starch.mg.g ~ Mortality14, data = nsc.adult.starch.mass,
        ylab = expression(paste("Starch mg g"^"-1")),ylim=c(-10,100))
text(seq_along(droplevels(nsc.adult.sugar.mass$Mortality14)), 100, b$n, cex=.9)
lm.mort <- lm(Starch.mg.g ~ Mortality14, data = nsc.adult.starch.mass)
anova(lm.mort)
summary(lm.mort)
#plot(lm.mort)

#mortality violin plot
ggplot(nsc.adult.starch.mass, aes(x=Mortality14, y=Starch.mg.g)) + 
  geom_violin() + stat_summary(fun.y=median,geom="point")

#sample sizes
table(nsc.adult.starch.mass$Mortality14)

#NSCs and average growth
nsc.adult.starch.mass$DBH09 <- as.numeric(nsc.adult.starch.mass$DBH09)
nsc.adult.starch.mass$DBH11 <- as.numeric(nsc.adult.starch.mass$DBH11)
nsc.adult.starch.mass$DBH12 <- as.numeric(nsc.adult.starch.mass$DBH12)
nsc.adult.starch.mass$DBH13 <- as.numeric(nsc.adult.starch.mass$DBH13)
nsc.adult.starch.mass$DBH14 <- as.numeric(nsc.adult.starch.mass$DBH14)

#avg growth #2 or 3 years is basically the same and sample size is very small starting in 09. n = 1353
nsc.adult.starch.mass$avg.growth <- (nsc.adult.starch.mass$DBH14 - nsc.adult.starch.mass$DBH12)/2

plot(nsc.adult.starch.mass$avg.growth,nsc.adult.starch.mass$Starch.mg.g,
     xlab = "Avg. Radial growth 2011-2014", ylab = expression(paste("Starch mg g"^"-1")),ylim=c(0,100))
lm.growth <- lm(nsc.adult.starch.mass$Starch.mg.g ~ nsc.adult.starch.mass$avg.growth)
abline(lm.growth)
summary(lm.growth)
#plot(lm.growth)

#NSCs and canopy status n = 1404
nsc.adult.starch.mass$Canopy15 <- as.factor(nsc.adult.starch.mass$Canopy15)
b <- boxplot(Starch.mg.g ~ Canopy15, data = nsc.adult.starch.mass)
boxplot(Starch.mg.g ~ Canopy15, data = nsc.adult.starch.mass,
        ylab = expression(paste("Starch mg g"^"-1")),ylim=c(0,100), 
        names = c("Understory","Co-Dominant","Emergent"))
text(seq_along(droplevels(nsc.adult.starch.mass$Canopy15)), 100, b$n, cex=.7)
lm.canopy <- lm(Starch.mg.g ~ Canopy15, data = nsc.adult.starch.mass)
anova(lm.canopy)
summary(lm.canopy)
#plot(lm.canopy)

#canopy violin plot
ggplot(nsc.adult.starch.mass, aes(x=Canopy15, y=Starch.mg.g)) + 
  geom_violin() + stat_summary(fun.y=median,geom="point")

#canopy status sample size
table(nsc.adult.starch.mass$Canopy15)

################Sugar ONLY Exploratory Analysis##################
#################################################################
########Sucrose only for now#########
nsc.adult.sucrose.mass <- nsc.adult.sugar.mass[which(nsc.adult.sugar.mass$sucrose.mg.g != "NA"),]

##adding PFTs to dataframe
#functional group
Early.Hardwood <- c(
  which(nsc.adult.sucrose.mass$Spp == "ACPE"),
  which(nsc.adult.sucrose.mass$Spp == "AIAL"),
  which(nsc.adult.sucrose.mass$Spp == "BEAL2"),
  which(nsc.adult.sucrose.mass$Spp == "BELE"),
  which(nsc.adult.sucrose.mass$Spp == "BEPA"),
  which(nsc.adult.sucrose.mass$Spp == "JUCI"),
  which(nsc.adult.sucrose.mass$Spp == "JUNI"),
  which(nsc.adult.sucrose.mass$Spp == "LIST2"),
  which(nsc.adult.sucrose.mass$Spp == "LITU"),
  which(nsc.adult.sucrose.mass$Spp == "POGR4"),
  which(nsc.adult.sucrose.mass$Spp == "POTR5"),
  which(nsc.adult.sucrose.mass$Spp == "PRSE2"),
  which(nsc.adult.sucrose.mass$Spp == "ROPS"),
  which(nsc.adult.sucrose.mass$Spp == "SAAL5")
)

Mid.Hardwood <- c(
  which(nsc.adult.sucrose.mass$Spp == "FRAM2"),
  which(nsc.adult.sucrose.mass$Spp == "FRNI"),
  which(nsc.adult.sucrose.mass$Spp == "FRQU"),
  which(nsc.adult.sucrose.mass$Spp == "ACNE2"),
  which(nsc.adult.sucrose.mass$Spp == "ACRU"),
  which(nsc.adult.sucrose.mass$Spp == "CACA18"),
  which(nsc.adult.sucrose.mass$Spp == "CACO15"),
  which(nsc.adult.sucrose.mass$Spp == "CAGL8"),
  which(nsc.adult.sucrose.mass$Spp == "CALA21"),
  which(nsc.adult.sucrose.mass$Spp == "CAOV2"),
  which(nsc.adult.sucrose.mass$Spp == "CATO6"),
  which(nsc.adult.sucrose.mass$Spp == "CECA4"),
  which(nsc.adult.sucrose.mass$Spp == "MAAC"),
  which(nsc.adult.sucrose.mass$Spp == "MAFR"),
  which(nsc.adult.sucrose.mass$Spp == "QUAL"),
  which(nsc.adult.sucrose.mass$Spp == "QUCO2"),
  which(nsc.adult.sucrose.mass$Spp == "QUIM"),
  which(nsc.adult.sucrose.mass$Spp == "QUMO4"),
  which(nsc.adult.sucrose.mass$Spp == "QUMU"),
  which(nsc.adult.sucrose.mass$Spp == "QUPH"),
  which(nsc.adult.sucrose.mass$Spp == "QURU"),
  which(nsc.adult.sucrose.mass$Spp == "QUVE"),
  which(nsc.adult.sucrose.mass$Spp == "ULAM")
)

Late.Hardwood <- c(which(nsc.adult.sucrose.mass$Spp == "ACSA3"),
                   which(nsc.adult.sucrose.mass$Spp == "FAGR"),
                   which(nsc.adult.sucrose.mass$Spp == "NYSY"),
                   which(nsc.adult.sucrose.mass$Spp == "OXAR"),
                   which(nsc.adult.sucrose.mass$Spp == "TIAM")
)


#set PFTs
nsc.adult.sucrose.mass$PFT <- NA
nsc.adult.sucrose.mass$PFT[Early.Hardwood] <- "Early.Hardwood"
nsc.adult.sucrose.mass$PFT[Mid.Hardwood] <- "Mid.Hardwood"
nsc.adult.sucrose.mass$PFT[Late.Hardwood] <- "Late.Hardwood"


Early.Conifer <- c(
  which(nsc.adult.sucrose.mass$Spp == "ABBA"),
  which(nsc.adult.sucrose.mass$Spp == "JUVI"),
  which(nsc.adult.sucrose.mass$Spp == "PIST"),
  which(nsc.adult.sucrose.mass$Spp == "PIVI"), #should be PIVI2 for Virginia Pine
  which(nsc.adult.sucrose.mass$Spp == "PIRE"),
  which(nsc.adult.sucrose.mass$Spp == "PIRI")
)

#Mid.Conifer <- c()

Late.Conifer <- c(
  which(nsc.adult.sucrose.mass$Spp == "LALA"),
  which(nsc.adult.sucrose.mass$Spp == "PIGL"),
  which(nsc.adult.sucrose.mass$Spp == "PIMA"),
  which(nsc.adult.sucrose.mass$Spp == "PIRU"),
  which(nsc.adult.sucrose.mass$Spp == "THOC2"),
  which(nsc.adult.sucrose.mass$Spp == "TSCA")
)

nsc.adult.sucrose.mass$PFT[Early.Conifer] <- "Early.Conifer"
nsc.adult.sucrose.mass$PFT[Late.Conifer] <- "Late.Conifer"

#make PFT factor
#nsc.adult.total.mass$PFT <- as.factor(nsc.adult.total.mass$PFT)

#only hardwood dataframe
early.h <- which(nsc.adult.sucrose.mass$PFT == "Early.Hardwood")
mid.h <- which(nsc.adult.sucrose.mass$PFT == "Mid.Hardwood")
late.h <- which(nsc.adult.sucrose.mass$PFT == "Late.Hardwood")

nsc.adult.sucrose.mass.hard <- rbind(
  nsc.adult.sucrose.mass[early.h,],
  nsc.adult.sucrose.mass[mid.h,],
  nsc.adult.sucrose.mass[late.h,]
)


nsc.adult.sucrose.mass.hard$Spp <- as.character(nsc.adult.sucrose.mass.hard$Spp)

nsc.adult.sucrose.mass.hard$Spp <- factor(nsc.adult.sucrose.mass.hard$Spp, 
                                         levels = c("ACPE","AIAL","BEAL2","BELE",
                                                    "BEPA","JUCI","JUNI","LIST2","LITU","POGR4","POTR5",
                                                    "PRSE2","ROPS","SAAL5","FRAM2","FRNI","FRQU","ACNE2",
                                                    "ACRU","CACA18","CACO15","CAGL8","CALA21","CAOV2",
                                                    "CATO6","CECA4","MAAC","MAFR","QUAL","QUCO2","QUIM",
                                                    "QUMO4","QUMU","QUPH","QURU","QUVE","ULAM","ACSA3",
                                                    "FAGR","NYSY","OXAR","TIAM"), ordered = TRUE)

# #drop unused hardwood factor levels
# nsc.adult.starch.mass.hard$Spp <- droplevels(nsc.adult.starch.mass.hard$Spp)
# nsc.adult.starch.mass.hard$PFT <- droplevels(nsc.adult.starch.mass.hard$PFT)

#only conifer dataframe
early.c <- which(nsc.adult.sucrose.mass$PFT == "Early.Conifer")
late.c <- which(nsc.adult.sucrose.mass$PFT == "Late.Conifer")

nsc.adult.sucrose.mass.conif <- rbind(
  nsc.adult.sucrose.mass[early.c,],
  nsc.adult.sucrose.mass[late.c,]
)



nsc.adult.sucrose.mass.conif$Spp <- as.character(nsc.adult.sucrose.mass.conif$Spp)

#relevel factors
nsc.adult.sucrose.mass.conif$Spp <- factor(nsc.adult.sucrose.mass.conif$Spp, 
                                          levels = c("ABBA","JUVI","PIST","PIVI","PIRE",
                                                     "PIRI","LALA","PIGL","PIMA","PIRU",
                                                     "THOC2","TSCA"), ordered = TRUE)

# #drop unused conifer factors
# nsc.adult.starch.mass.conif$Spp <- droplevels(nsc.adult.starch.mass.conif$Spp)
# nsc.adult.starch.mass.conif$PFT <- droplevels(nsc.adult.starch.mass.conif$PFT)

# hardwood ggplot boxplot
#median for total.nsc
med.early.hard.sucrose <- median(nsc.adult.sucrose.mass$sucrose.mg.g[early.h])
med.mid.hard.sucrose <- median(nsc.adult.sucrose.mass$sucrose.mg.g[mid.h])
med.late.hard.sucrose <- median(nsc.adult.sucrose.mass$sucrose.mg.g[late.h])
med.total.sucrose <- median(nsc.adult.sucrose.mass.hard$sucrose.mg.g)

ggplot(nsc.adult.sucrose.mass.hard, aes(x = Spp, y=sucrose.mg.g, fill=PFT)) + 
  geom_boxplot() + theme(axis.text.x = element_text(angle=90,hjust=1,size=16)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.title.x = element_text(size=20)) +
  theme(axis.title.y = element_text(size=20)) +
  theme(legend.text = element_text(size=16)) +
  stat_summary(fun.data = give.n, geom= "text", aes(y= 57)) + 
  #geom_hline(aes(yintercept=med.total),linetype="longdash") +
  ylab(expression(paste("Sucrose mg g"^"-1"," Sapwood"))) + xlab("Species") +
  geom_segment(aes(x=0,y=med.early.hard.sucrose,xend=13.5,
                   yend=med.early.hard.sucrose),linetype="longdash",color="indianred4") +
  geom_segment(aes(x=13.5,y=med.mid.hard.sucrose,xend=33.3,
                   yend=med.mid.hard.sucrose),linetype="longdash",color="dodgerblue4") +
  geom_segment(aes(x=33.5,y=med.late.hard.sucrose,xend=38.8,
                   yend=med.late.hard.sucrose),linetype="longdash",color="darkgreen") +
  scale_fill_manual(values=cbPalette)



#ggplot boxplots of sucrose species
#Conifer ggplot boxplot
med.early.conif.sucrose <- median(nsc.adult.sucrose.mass$sucrose.mg.g[early.c])
med.late.conif.sucrose <- median(nsc.adult.sucrose.mass$sucrose.mg.g[late.c])
#med.total <- median(nsc.adult.total.mass.hard$total.nsc)

ggplot(nsc.adult.sucrose.mass.conif, aes(x = Spp, y=sucrose.mg.g, fill=PFT)) + 
  geom_boxplot() + theme(axis.text.x = element_text(angle=90,hjust=1,size=16)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.title.x = element_text(size=20)) +
  theme(axis.title.y = element_text(size=20)) +
  theme(legend.text = element_text(size=16)) +
  stat_summary(fun.data = give.n, geom= "text", aes(y= 33)) + 
  #geom_hline(aes(yintercept=med.total),linetype="longdash") +
  ylab(expression(paste("Sucrose NSC mg g"^"-1"))) + xlab("Species") +
  geom_segment(aes(x=0,y=med.early.conif.sucrose,xend=6.5,
                   yend=med.early.conif.sucrose),linetype="longdash",color="indianred4") +
  geom_segment(aes(x=6.6,y=med.late.conif.sucrose,xend=9.5,
                   yend=med.late.conif.sucrose),linetype="longdash",color="darkgreen") +
  scale_fill_manual(values=rev(cbPalette))


#species anova within hardwoods
spp.hard.sucrose.lm <- lm(sucrose.mg.g ~ Spp, data = nsc.adult.sucrose.mass.hard)
anova(spp.hard.sucrose.lm)
summary(spp.hard.sucrose.lm)

#pft anova within hardwoods
pft.hard.sucrose.lm <- lm(sucrose.mg.g ~ PFT, data = nsc.adult.sucrose.mass.hard)
anova(pft.hard.sucrose.lm)
summary(pft.hard.sucrose.lm)

#species anova all individuals
spp.sucrose.lm <- lm(sucrose.mg.g ~ Spp, data = nsc.adult.sucrose.mass)
anova(spp.sucrose.lm)
summary(spp.sucrose.lm)

#pft anova all species
pft.sucrose.lm <- lm(sucrose.mg.g ~ PFT, data = nsc.adult.sucrose.mass)
anova(pft.sucrose.lm)
summary(pft.sucrose.lm)

#species anova conifers
spp.conif.sucrose.lm <- lm(sucrose.mg.g ~ Spp, data = nsc.adult.sucrose.mass.conif)
anova(spp.conif.sucrose.lm)
summary(spp.conif.sucrose.lm)

#pft anova conifers
pft.conif.sucrose.lm <- lm(sucrose.mg.g ~ PFT, data = nsc.adult.sucrose.mass.conif)
anova(pft.conif.sucrose.lm)
summary(pft.conif.sucrose.lm)


####Sucrose####
#Species, n = 1355
nsc.adult.sugar.mass$Spp <- as.factor(nsc.adult.sugar.mass$Spp)
#relevel factors, changed Spp to factor 
#and set BEAL2 as ref group because close to median of all species
nsc.adult.sugar.mass$Spp <- relevel(nsc.adult.sugar.mass$Spp, ref = "ACNE2")

#boxplot n = 757
b <- boxplot(sucrose.mg.g ~ droplevels(Spp), data = nsc.adult.sugar.mass,
        las=2, ylab = "Concentration mg/g") # Glucose Equivalency Units")
boxplot(sucrose.mg.g ~ droplevels(Spp), data = nsc.adult.sugar.mass,
               las=2, ylab = expression(paste("Sucrose mg g"^"-1")),ylim=c(0,55))
text(seq_along(droplevels(nsc.adult.sugar.mass$Spp)), 53, b$n, cex=.9) #b$stats[4,], b$n, cex = .75)
abline(h=median(nsc.adult.sugar.mass$sucrose.mg.g, na.rm = TRUE),col=2)

#species sucrose violin plot 
ggplot(nsc.adult.sugar.mass, aes(x=Spp, y=sucrose.mg.g)) + 
  geom_violin() + stat_summary(fun.y=median,geom="point")

lm.spp.suc <- lm(sucrose.mg.g ~Spp, data = nsc.adult.sugar.mass)
anova(lm.spp.suc)
summary(lm.spp.suc)
#plot(lm.spp.suc) #assumptions of anova not met: heteroskedasticity, 
#non-normal residuals, residuals increase with fitted values

####Fructose####
#Species boxplot n = 757
b <- boxplot(fructose.mg.g ~ droplevels(Spp), data = nsc.adult.sugar.mass,
             las=2, ylab = expression(paste("Fructose mg g"^"-1"))) # Glucose Equivalency Units")
boxplot(fructose.mg.g ~ droplevels(Spp), data = nsc.adult.sugar.mass,
        las=2, ylab = expression(paste("Fructose mg g"^"-1")),ylim=c(0,28))
text(seq_along(droplevels(nsc.adult.sugar.mass$Spp)), 28, b$n, cex=.7) #b$stats[4,], b$n, cex = .75)
abline(h=median(nsc.adult.sugar.mass$fructose.mg.g, na.rm = TRUE),col=2)

#species fructose violin plot
ggplot(nsc.adult.sugar.mass, aes(x=Spp, y=fructose.mg.g)) + 
  geom_violin() + stat_summary(fun.y=median,geom="point")


lm.spp.fruc <- lm(fructose.mg.g ~Spp, data = nsc.adult.sugar.mass)
anova(lm.spp.fruc)
summary(lm.spp.fruc)
#plot(lm.spp.fruc) #assumptions of anova not met: heteroskedasticity, 
#non-normal residuals, residuals increase with fitted values

####Glucose####
#Species boxplot n = 757
b <- boxplot(glucose.mg.g ~ droplevels(Spp), data = nsc.adult.sugar.mass,
             las=2, ylab = expression(paste("Glucose mg g"^"-1"))) # Glucose Equivalency Units")
boxplot(glucose.mg.g ~ droplevels(Spp), data = nsc.adult.sugar.mass,
        las=2, ylab = expression(paste("Glucose mg g"^"-1")),ylim=c(0,55))
text(seq_along(droplevels(nsc.adult.sugar.mass$Spp)), 50, b$n, cex=.7) #b$stats[4,], b$n, cex = .75)
abline(h=median(nsc.adult.sugar.mass$glucose.mg.g, na.rm = TRUE),col=2)

#violin glucose plot 
ggplot(nsc.adult.sugar.mass, aes(x=Spp, y=glucose.mg.g)) + 
  geom_violin() + stat_summary(fun.y=median,geom="point")


lm.spp.gluc <- lm(glucose.mg.g ~Spp, data = nsc.adult.sugar.mass)
anova(lm.spp.gluc)
summary(lm.spp.gluc)
#plot(lm.spp.gluc) #assumptions of anova not met: heteroskedasticity, 
#non-normal residuals, residuals increase with fitted values


####Sucrose####
# Site n = 798
#nsc.adult.sugar.mass$Site <- relevel(nsc.adult.sugar.mass$Site, ref = "D")
b <- boxplot(sucrose.mg.g ~ droplevels(Site), data = nsc.adult.sugar.mass, 
        ylab = expression(paste("Sucrose mg g"^"-1")))
boxplot(sucrose.mg.g ~ droplevels(Site), data = nsc.adult.sugar.mass, 
        ylab = expression(paste("Sucrose mg g"^"-1")),ylim=c(0,55))
text(seq_along(droplevels(nsc.adult.sugar.mass$Site)), 53, b$n, cex=.7)
abline(h=median(nsc.adult.sugar.mass$sucrose.mg.g, na.rm = TRUE))
lm.site.suc <- lm(sucrose.mg.g ~ Site, data = nsc.adult.sugar.mass)
anova(lm.site.suc)
summary(lm.site.suc)
#plot(lm.site.suc)

#Site sucrose violin plot
ggplot(nsc.adult.sugar.mass, aes(x=Site, y=sucrose.mg.g)) + 
  geom_violin() + stat_summary(fun.y=median,geom="point")


####Fructose####
#nsc.adult.sugar.mass$Site <- relevel(nsc.adult.sugar.mass$Site, ref = "D")
b <- boxplot(fructose.mg.g ~ droplevels(Site), data = nsc.adult.sugar.mass)
boxplot(fructose.mg.g ~ droplevels(Site), data = nsc.adult.sugar.mass,
        ylab = expression(paste("Fructose mg g"^"-1")),ylim=c(0,28))
text(seq_along(droplevels(nsc.adult.sugar.mass$Site)), 28, b$n, cex=.7)
abline(h=median(nsc.adult.sugar.mass$fructose.mg.g, na.rm = TRUE))
lm.site.fruc <- lm(fructose.mg.g ~ Site, data = nsc.adult.sugar.mass)
anova(lm.site.fruc)
summary(lm.site.fruc)
#plot(lm.site.fruc)

#Site fructose violin plot
ggplot(nsc.adult.sugar.mass, aes(x=Site, y=fructose.mg.g)) + 
  geom_violin() + stat_summary(fun.y=median,geom="point")


####Glucose####
#nsc.adult.sugar.mass$Site <- relevel(nsc.adult.sugar.mass$Site, ref = "D")
b <- boxplot(glucose.mg.g ~ droplevels(Site), data = nsc.adult.sugar.mass)
boxplot(glucose.mg.g ~ droplevels(Site), data = nsc.adult.sugar.mass,
        ylab = expression(paste("Glucose mg g"^"-1")),ylim=c(0,55))
text(seq_along(droplevels(nsc.adult.sugar.mass$Site)), 50, b$n, cex=.7)
abline(h=median(nsc.adult.sugar.mass$glucose.mg.g, na.rm = TRUE))
lm.site.gluc <- lm(glucose.mg.g ~ Site, data = nsc.adult.sugar.mass)
anova(lm.site.gluc)
summary(lm.site.gluc)
#plot(lm.site.gluc)

#Glucose Site violin plot
ggplot(nsc.adult.sugar.mass, aes(x=Site, y=glucose.mg.g)) + 
  geom_violin() + stat_summary(fun.y=median,geom="point")


#tried log transforming but didn't do much
####Sucrose####
#DBH14 n = 1403
plot(nsc.adult.sugar.mass$DBH14, nsc.adult.sugar.mass$sucrose.mg.g)
lm.dbh.suc <- lm(sucrose.mg.g ~ DBH14, data = nsc.adult.sugar.mass)
abline(lm.dbh.suc)
summary(lm.dbh.suc) #significant but explains very little

####Fructose####
plot(nsc.adult.sugar.mass$DBH14, nsc.adult.sugar.mass$fructose.mg.g)
lm.dbh.fruc <- lm(fructose.mg.g ~ DBH14, data = nsc.adult.sugar.mass)
abline(lm.dbh.fruc)
summary(lm.dbh.fruc) #significant but explains very little

####Glucose####
plot(nsc.adult.sugar.mass$DBH14, nsc.adult.sugar.mass$glucose.mg.g)
lm.dbh.gluc <- lm(glucose.mg.g ~ DBH14, data = nsc.adult.sugar.mass)
abline(lm.dbh.gluc)
summary(lm.dbh.gluc) #significant but explains very little


####Sucrose####
#Mortality Status n = 799
nsc.adult.sugar.mass$Mortality14 <- as.factor(nsc.adult.sugar.mass$Mortality14)
#nsc.adult.sugar.mass$Mortality14 <- as.character(nsc.adult.sugar.mass$Mortality14)
b <- boxplot(sucrose.mg.g ~ droplevels(Mortality14), data = nsc.adult.sugar.mass)
boxplot(sucrose.mg.g ~ droplevels(Mortality14), data = nsc.adult.sugar.mass,
        ylab = expression(paste("Sucrose mg g"^"-1")),ylim=c(0,55))
text(seq_along(droplevels(nsc.adult.sugar.mass$Mortality14)), 53, b$n, cex=.9)
lm.mort.suc <- lm(sucrose.mg.g ~ Mortality14, data = nsc.adult.sugar.mass)
anova(lm.mort.suc)
summary(lm.mort.suc)
#plot(lm.mort.suc)

####Fructose####
b <- boxplot(fructose.mg.g ~ droplevels(Mortality14), data = nsc.adult.sugar.mass)
boxplot(fructose.mg.g ~ droplevels(Mortality14), data = nsc.adult.sugar.mass,
        ylab = expression(paste("Fructose mg g"^"-1")),ylim=c(0,28))
#text(seq_along(droplevels(nsc.adult.sugar.mass$Mortality14)), 25, b$n, cex=.7)
lm.mort.fruc <- lm(fructose.mg.g ~ Mortality14, data = nsc.adult.sugar.mass)
anova(lm.mort.fruc)
summary(lm.mort.fruc)
#plot(lm.mort.fruc)

####Glucose####
b <- boxplot(glucose.mg.g ~ droplevels(Mortality14), data = nsc.adult.sugar.mass)
boxplot(glucose.mg.g ~ droplevels(Mortality14), data = nsc.adult.sugar.mass,
        ylab = expression(paste("Glucose mg g"^"-1")),ylim=c(0,55))
#text(seq_along(droplevels(nsc.adult.sugar.mass$Mortality14)), 25, b$n, cex=.7)
lm.mort.gluc <- lm(glucose.mg.g ~ Mortality14, data = nsc.adult.sugar.mass)
anova(lm.mort.gluc)
summary(lm.mort.gluc)
#plot(lm.mort.fruc)

#sample sizes
table(nsc.adult.sugar.mass$Mortality14)

####Growth####
#NSCs and average growth
nsc.adult.sugar.mass$DBH09 <- as.numeric(nsc.adult.sugar.mass$DBH09)
nsc.adult.sugar.mass$DBH11 <- as.numeric(nsc.adult.sugar.mass$DBH11)
nsc.adult.sugar.mass$DBH12 <- as.numeric(nsc.adult.sugar.mass$DBH12)
nsc.adult.sugar.mass$DBH13 <- as.numeric(nsc.adult.sugar.mass$DBH13)
nsc.adult.sugar.mass$DBH14 <- as.numeric(nsc.adult.sugar.mass$DBH14)

#avg growth #2 or 3 years is basically the same and sample size is very small starting in 09. n = 1353
nsc.adult.sugar.mass$avg.growth <- (nsc.adult.sugar.mass$DBH14 - nsc.adult.sugar.mass$DBH12)/2


####Sucrose####
#Growth
plot(nsc.adult.sugar.mass$avg.growth,nsc.adult.sugar.mass$sucrose.mg.g,
     xlab = "Avg. Radial growth 2011-2014", ylab = expression(paste("Sucrose mg g"^"-1")))
lm.growth.suc <- lm(nsc.adult.sugar.mass$sucrose.mg.g ~ nsc.adult.sugar.mass$avg.growth)
abline(lm.growth.suc)
summary(lm.growth.suc)
plot(lm.growth.suc)

####Fructose####
plot(nsc.adult.sugar.mass$avg.growth,nsc.adult.sugar.mass$fructose.mg.g,
     xlab = "Avg. Radial growth 2011-2014", ylab = expression(paste("Fructose mg g"^"-1")))
lm.growth.fruc <- lm(nsc.adult.sugar.mass$fructose.mg.g ~ nsc.adult.sugar.mass$avg.growth)
abline(lm.growth.fruc)
summary(lm.growth.fruc)
plot(lm.growth.fruc)

####Glucose####
plot(nsc.adult.sugar.mass$avg.growth,nsc.adult.sugar.mass$glucose.mg.g,
     xlab = "Avg. Radial growth 2011-2014", ylab = expression(paste("Glucose mg g"^"-1")))
lm.growth.gluc <- lm(nsc.adult.sugar.mass$glucose.mg.g ~ nsc.adult.sugar.mass$avg.growth)
abline(lm.growth.gluc)
summary(lm.growth.gluc)
plot(lm.growth.gluc)


####Sucrose####
#NSCs and canopy status n = 800
nsc.adult.sugar.mass$Canopy15 <- as.factor(nsc.adult.sugar.mass$Canopy15)
b <- boxplot(sucrose.mg.g ~ Canopy15, data = nsc.adult.sugar.mass, 
             ylab = expression(paste("Sucrose mg g"^"-1")), xlab= "canopy status")
boxplot(sucrose.mg.g ~ Canopy15, data = nsc.adult.sugar.mass, 
        ylab = expression(paste("Sucrose mg g"^"-1")), xlab= "canopy status",ylim=c(0,55))
text(seq_along(droplevels(nsc.adult.sugar.mass$Mortality14)), 53, b$n, cex=.7)
lm.canopy.suc <- lm(sucrose.mg.g ~ Canopy15, data = nsc.adult.sugar.mass)
anova(lm.canopy.suc)
summary(lm.canopy.suc)
#plot(lm.canopy.suc)

####Fructose####
b <- boxplot(fructose.mg.g ~ Canopy15, data = nsc.adult.sugar.mass, 
             ylab = expression(paste("Fructose mg g"^"-1")), xlab= "canopy status")
boxplot(fructose.mg.g ~ Canopy15, data = nsc.adult.sugar.mass, 
        ylab = expression(paste("Fructose mg g"^"-1")), xlab= "canopy status",ylim=c(0,28))
text(seq_along(droplevels(nsc.adult.sugar.mass$Mortality14)), 28, b$n, cex=.7)
lm.canopy.fruc <- lm(fructose.mg.g ~ Canopy15, data = nsc.adult.sugar.mass)
anova(lm.canopy.fruc)
summary(lm.canopy.fruc)

####Glucose####
b <- boxplot(glucose.mg.g ~ Canopy15, data = nsc.adult.sugar.mass, 
             ylab = expression(paste("Glucose mg g"^"-1")), xlab= "canopy status")
boxplot(glucose.mg.g ~ Canopy15, data = nsc.adult.sugar.mass, 
        ylab = expression(paste("Glucose mg g"^"-1")), xlab= "canopy status",ylim=c(0,55))
text(seq_along(droplevels(nsc.adult.sugar.mass$Mortality14)), 53, b$n, cex=.7)
lm.canopy.gluc <- lm(glucose.mg.g ~ as.factor(Canopy15), data = nsc.adult.sugar.mass)
anova(lm.canopy.gluc)
summary(lm.canopy.gluc)

#canopy status sample size
table(nsc.adult.sugar.mass$Canopy15)


############################Relationships btwn Sugar, Starch, Total############################

