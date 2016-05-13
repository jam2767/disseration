#subsetting data by species and reordering folders for xdating
#crossdate species at each site individually using cofecha
#create species folders within sites
#identify missing trees from census data and tucson file matches

#match tucson files and adult data for each site

#Read all inventory data
trees.xdate <- read.csv("/Users/josh/Dropbox/Dissertation/CH1_Treerings/Data/Adult_Field_Data_JAM_MCD.csv")

#Harvard Forest
rings.xdate.h <- Read_Tuscon("/Users/josh/Dropbox/Dissertation/CH1_Treerings/Data/Windendro_text/Harvard (H)/Tucson/")
combined.xdate.h <- matchInventoryRings(trees=trees.xdate,rings=rings.xdate.h)
write.csv(x=combined.xdate.h,file="/Users/josh/Dropbox/Dissertation/CH1_Treerings/Data/H_Treerings.csv")

#Baskett 
rings.xdate.m <- Read_Tuscon("/Users/josh/Dropbox/Dissertation/CH1_Treerings/Data/Windendro_text/Baskett (M)/Tucson/")
combined.xdate.m <- matchInventoryRings(trees=trees.xdate,rings=rings.xdate.m)
write.csv(x=combined.xdate.m,file="/Users/josh/Dropbox/Dissertation/CH1_Treerings/Data/M_Treerings.csv")

#Smithsonian
rings.xdate.s <- Read_Tuscon("/Users/josh/Dropbox/Dissertation/CH1_Treerings/Data/Windendro_text/Smithsonian (S)/Tucson/")
combined.xdate.s <- matchInventoryRings(trees=trees.xdate,rings=rings.xdate.s)
write.csv(x=combined.xdate.s,file="/Users/josh/Dropbox/Dissertation/CH1_Treerings/Data/S_Treerings.csv")

#Bartlett
rings.xdate.b <- Read_Tuscon("/Users/josh/Dropbox/Dissertation/CH1_Treerings/Data/Windendro_text/Bartlett (B)/Tucson/")
combined.xdate.b <- matchInventoryRings(trees=trees.xdate,rings=rings.xdate.b)
write.csv(x=combined.xdate.b,file="/Users/josh/Dropbox/Dissertation/CH1_Treerings/Data/B_Treerings.csv")

#UNDERC
rings.xdate.w <- Read_Tuscon("/Users/josh/Dropbox/Dissertation/CH1_Treerings/Data/Windendro_text/UNDERC (W)/Tucson/")
combined.xdate.w <- matchInventoryRings(trees=trees.xdate,rings=rings.xdate.w)
write.csv(x=combined.xdate.w,file="/Users/josh/Dropbox/Dissertation/CH1_Treerings/Data/W_Treerings.csv")

#VRO
rings.xdate.v <- Read_Tuscon("/Users/josh/Dropbox/Dissertation/CH1_Treerings/Data/Windendro_text/Vermillion (V)/Tucson/")
combined.xdate.v <- matchInventoryRings(trees=trees.xdate,rings=rings.xdate.v)
write.csv(x=combined.xdate.v,file="/Users/josh/Dropbox/Dissertation/CH1_Treerings/Data/V_Treerings.csv")

#Oak Ridge
rings.xdate.t <- Read_Tuscon("/Users/josh/Dropbox/Dissertation/CH1_Treerings/Data/Windendro_text/ORNL (T)/Tucson/")
combined.xdate.t <- matchInventoryRings(trees=trees.xdate,rings=rings.xdate.t)
write.csv(x=combined.xdate.t,file="/Users/josh/Dropbox/Dissertation/CH1_Treerings/Data/T_Treerings.csv")

##need to convert files to Tucscon before I can xdate
#Duke Forest
rings.xdate.d <- Read_Tuscon("/Users/josh/Dropbox/Dissertation/CH1_Treerings/Data/Windendro_text/Duke (D) 2013/Tucson/")
combined.xdate.d <- matchInventoryRings(trees=trees.xdate,rings=rings.xdate.d)
write.csv(x=combined.xdate.d,file="/Users/josh/Dropbox/Dissertation/CH1_Treerings/Data/D_Treerings.csv")

##need to convert files to Tucscon before I can xdate
#Coweeta
rings.xdate.c <- Read_Tuscon("/Users/josh/Dropbox/Dissertation/CH1_Treerings/Data/Windendro_text/Coweeta (C)/Tucson/")
combined.xdate.c <- matchInventoryRings(trees=trees.xdate,rings=rings.xdate.c)
write.csv(x=combined.xdate.c,file="/Users/josh/Dropbox/Dissertation/CH1_Treerings/Data/C_Treerings.csv")

#Ordway Swisher
rings.xdate.f <- Read_Tuscon("/Users/josh/Dropbox/Dissertation/CH1_Treerings/Data/Windendro_text/Ordway Swisher (F)/Tucson/")
combined.xdate.f <- matchInventoryRings(trees=trees.xdate,rings=rings.xdate.f)
write.csv(x=combined.xdate.f,file="/Users/josh/Dropbox/Dissertation/CH1_Treerings/Data/F_Treerings.csv")


