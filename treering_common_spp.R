#identify common species for treering analysis
#identify species with most individuals and how many sites
#they are found in. Looking for enough individuals to xdate
#and ability to look across scales - region, site, plot, sub,
#individual

#read adult data
adult <- read.csv("/Users/Josh/Dropbox/Dietze_Lab_Undergrads/JAM - Xsite/Field Data Entry/Data Sheets/Adult_Field_Data_JAM_MCD.csv")

#subset for trees > 15cm DBH
big <- which(adult$DBH15 >= 15)
adult.big <- adult[big,]

#subset for stem 1 #may miss a couple of trees that were cored but were stem 2
#just check those by hand
leader <- which(adult.big$Num == 1)
adult.cored <- adult.big[leader,]

#number of cores per site
site.num <- count(adult.cored, vars = "Site")

#identify spp numbers
spp.num <- count(adult.cored, vars = "Spp")

#identify species by site
spp.site.num <- count(adult.cored, vars = c("Spp","Site"))
print(spp.site.num)
write.csv(spp.site.num, file = "/Users/josh/Dropbox/Dissertation/CH1_Treerings/Data/species_by_site.csv")

#high priority species:
#ACSA3, BEPA, CAGL8, CATO6, FAGR, FRAM2, LITU, QUAL, QUMO4, QURU, QUVE, TSCA

#TODO identify priority species using a rank-abundance curve or something
#similar. Make graph too.

