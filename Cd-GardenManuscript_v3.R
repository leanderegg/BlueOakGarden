######################################################
##             Code for:                            ##
##      "Plasticity, geographic variation and trait 
##      coordination in blue oak drought physiology"
######################################################

# Analyzes hydraulic and morphological trait data from a common garden and the 7 wild source populations
# of Quercus douglasii

#load packages
library(RColorBrewer)
library(dplyr)
library(reshape)
library(nlme)
library(lme4)
library(lmerTest)
library(MuMIn)
library(ggplot2)
library(stringr)
require(lmodel2)
#require(rgdal) # getting retired 2023 and can no longer download
require(car)
require(cowplot)
require(quantreg)
require(tidyr)
require(ggcorrplot)
library(sjPlot)
library(png)
library(mgcv)

mypal <- c(brewer.pal(n=9, "Set1"), brewer.pal(n=8, "Dark2"))
palette(mypal)
pal2 <- brewer.pal(n=8, "Set2")
pal2light <- paste0(pal2,"33")

# also load to colorblind friendly palettes
# The palette with grey:
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

se <- function(x){
  se <- sd(x, na.rm=T)/sqrt(length(x[!is.na(x)]))
  return(se)
}
# function for plotting upper and lower error bars

error_bars<- function(xvar, yvar, upper, lower=upper, errordata, color="black", length=0.1, lwd=1, plotlower=T, plotupper=T,angle=90,...){
  if(plotupper==T){
    # upper error bar
    arrows(x0 = errordata[,xvar], y0=errordata[,yvar], 
           x1=errordata[,xvar], y1=errordata[,yvar]+errordata[,upper], #add the sterr to the mean value
           angle = angle, col=color, length = length, lwd=lwd,...)
  }
  if(plotlower==T){
    # lower error bar
    arrows(x0 = errordata[,xvar], y0=errordata[,yvar], 
           x1=errordata[,xvar], y1=errordata[,yvar]-errordata[,lower], #subtract the sterr from the mean value
           angle=angle, col=color, length = length, lwd=lwd,...)
  }
}

# set your working directory if desired
# setwd("")


save.figures <- T # whether to save figure pdfs

#results.version <- "v240105" # updated for revision
#results.version <- "v240210" # renumbered figures
#results.version <- "v240628" # corrected kstem/Ks units
#results.version <- "v240821" # updated for resubmission
results.version <- "v250505" # for NP submission
results.dir <- paste0("Results_",results.version)
if(save.figures == T) { dir.create(results.dir)}


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#######   BEGIN: LOAD DATA ###################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

popclim <- read.csv("Data/PopulationClimate_240830.csv")


#______________________________________________________________________
#######   * Blue Oak herbarium records ###################################################
#______________________________________________________________________
# extracted from Baldwin et al. 2017 
# Baldwin BG, Thornhill AH, Freyman WA, Ackerly DD, Kling MM, Morueta-Holme N, Mishler BD. 2017. Species richness and endemism in the native flora of California. American Journal of Botany 104: 487–501.

qudo <- read.csv("Data/California_QuercusDouglasii_specimens_20180318.csv")

# clean up the collection date
first <- str_extract(qudo$early_julian_day, "\\b[:digit:]+") # grab the first part of the date string
# this is either month (mm/dd/yyyy) or year (yyyy-mm-dd)
last <- str_extract(qudo$early_julian_day, "[:digit:]+$\\b")
year <- rep(NA, times=nrow(qudo))
year[which(nchar(first)==4)] <- first[which(nchar(first)==4)]
year[which(nchar(last)==4)] <- last[which(nchar(last)==4)]
qudo$year <- as.numeric(year)






#______________________________________________________________________
#######   * Garden growth data (from survey winter 2017) ###################################################
#______________________________________________________________________


## import hopland tree sizes
# from pre-harvest survey, (2017)
hoptrees <- read.csv("Data/Hop_trees_20240830.csv")
hoptrees$bioest <- (hoptrees$diam_50cm/2)^2 * pi * hoptrees$height  # currently as a cylinder 
hoptrees$ind <- factor(hoptrees$sub)
levels(hoptrees$ind) <- list(a = "1", b="2", c="3")
hoptrees$cwd <- popclim$cwd1951_1980[match(hoptrees$pop, popclim$Garden_pop)]
#hoptrees$soil_depth <- popinfo2$soil_depth[match(hoptrees$pop, popinfo2$Garden_pop)]
hoptrees$aet <- popclim$aet1951_1980[match(hoptrees$pop, popclim$Garden_pop)]
hoptrees$pet <- popclim$pet1951_1980[match(hoptrees$pop, popclim$Garden_pop)]
hoptrees$ppt <- popclim$ppt1951_1980[match(hoptrees$pop, popclim$Garden_pop)]
hoptrees$tmn <- popclim$tmn1951_1980[match(hoptrees$pop, popclim$Garden_pop)]
hoptrees$tmx <- popclim$tmx1951_1980[match(hoptrees$pop, popclim$Garden_pop)]

hoptrees$log.bio <- log(hoptrees$bioest, base=10)

## calcualte Transfer distance between garden and source locations
hoptrees$cwd.td <- hoptrees$cwd - popclim$cwd1951_1980[which(popclim$Site=="HREC Garden")]
hoptrees$aet.td <- hoptrees$aet - popclim$aet1951_1980[which(popclim$Site=="HREC Garden")]
hoptrees$pet.td <- hoptrees$pet - popclim$pet1951_1980[which(popclim$Site=="HREC Garden")]
hoptrees$ppt.td <- hoptrees$ppt - popclim$ppt1951_1980[which(popclim$Site=="HREC Garden")]
hoptrees$tmn.td <- hoptrees$tmn - popclim$tmn1951_1980[which(popclim$Site=="HREC Garden")]
hoptrees$tmx.td <- hoptrees$tmx - popclim$tmx1951_1980[which(popclim$Site=="HREC Garden")]


hoppop <- hoptrees %>% group_by(pop) %>% summarise(Bio = mean(bioest, na.rm=T), sdBio=sd(bioest,na.rm=T), seBio=se(bioest),
                                                   log.Bio =mean(log.bio, na.rm=T), sdlog.Bio=sd(log.bio,na.rm=T), selog.Bio=se(log.bio),
                                                       mort = length(which(is.na(height))), mort2 = length(which(stems==0)), mort3 = 18 - n() + length(which(is.na(height))),
                                                       Height = mean(height, na.rm=T), sdHeight=sd(height,na.rm=T), seHeight=se(height),
                                                       Diam = mean(diam_50cm, na.rm=T), sdDiam=sd(diam_50cm, na.rm=T),seDaim=se(diam_50cm),
                                                       cwd = unique(cwd), ppt = unique(ppt), aet=unique(aet), pet=unique(pet), tmn=unique(tmn), tmx=unique(tmx),
                                                   cwd.td = unique(cwd.td), ppt.td = unique(ppt.td), aet.td=unique(aet.td), pet.td=unique(pet.td), tmn.td=unique(tmn.td), tmx.td=unique(tmx.td))  



#______________________________________________________________________
#######   * Wild growth data (from tree cores Oct 2018) ###################################################
#______________________________________________________________________

wildgrowth <- read.csv("Data/WildGrowth_Oct2018_20230428.csv", row.names=1)

# calculate Basal Area Increment from 5yr radial growth and tree DBH
wildgrowth$BAI_cm2 <- (wildgrowth$DBH_cm/2)^2 * pi - (wildgrowth$DBH_cm/2 - wildgrowth$growth5yr_mm/10)^2 * pi
  # area of whole tree - area of tree 5 years ago

# -- NOT RUN -- 
# # visualize BAI and quantile regressions to select the right max growth quantile
# ggplot(wildgrowth, aes(x=DBH_cm, y=growth5yr_mm, col=Site)) + geom_point() # expected decreasing radial growth with DBH
# ggplot(wildgrowth, aes(x=DBH_cm, y=BAI_cm2, col=Site)) + geom_point() # expected increase in BAI with DBH, looks like linear upper boundary
# 
# # visualize quantile regressions of BAI~DBH to verify which quantile to use for calculating max size-specific BAI
# plot(BAI_cm2~DBH_cm, wildgrowth,type="n",xlab="DBH (cm)", ylab="Basal Area Increment (cm2)")
# points(BAI_cm2~DBH_cm, wildgrowth, col=factor(Site), pch=16)
# points(BAI_cm2~DBH_cm, wildgrowth[which(wildgrowth$Site %in% c("SJR", "LYN","SMR","SMT","SON", "SRD")),], pch=1, cex=1.1) # highlight the trees actually from this study
# # median
# abline(rq(BAI_cm2~DBH_cm, wildgrowth,tau=.5),col="blue")
# # other range of quantiles
# taus <- c(.05,.1,.25,.75,.90,.95)
# for( i in 1:length(taus)){
#   abline(rq(BAI_cm2~DBH_cm, wildgrowth,tau=taus[i]),col="gray")
# }
# # looks like 90th and 95th percentile are almost identical, so will use 90th percentile



# calculate raw BAI suppression from 90th percentile BAI growth, and %max BAI
maxgr.90 <- rq(BAI_cm2~DBH_cm, wildgrowth, tau=0.9) # fit the quantile regression


# write a function to calculate how much growth is supressed from maximum, both in cm2 BAI and % max BAI for a given DBH
BAI_suppression <- function(BAI, DBH, maxgrowth){
  # DBH column, BAI column, and a quantile regression object for max growth
  predgrowth <- maxgrowth$coefficients[1] + maxgrowth$coefficients[2] * DBH 
  growthsup <- predgrowth - BAI
  perc_BAI <- BAI/predgrowth
  return(cbind(growthsup, perc_BAI))
}

wildgrowth$growthsup_cm2 <- BAI_suppression(BAI=wildgrowth$BAI_cm2, DBH=wildgrowth$DBH_cm, maxgrowth=maxgr.90)[,1]
wildgrowth$perc_maxBAI <- BAI_suppression(BAI=wildgrowth$BAI_cm2, DBH=wildgrowth$DBH_cm, maxgrowth=maxgr.90)[,2]

# select down to only trees from the source populations
wg <- wildgrowth %>% filter(Site %in% c("SJR", "LYN","SMR","SMT","SON", "SRD"))

# Adding site Climate normals
wg$cwd <- popclim$cwd1951_1980[match(wg$Site, popclim$Code)]
wg$aet <- popclim$aet1951_1980[match(wg$Site, popclim$Code)]
wg$pet <- popclim$pet1951_1980[match(wg$Site, popclim$Code)]
wg$ppt <- popclim$ppt1951_1980[match(wg$Site, popclim$Code)]
wg$tmn <- popclim$tmn1951_1980[match(wg$Site, popclim$Code)]
wg$tmx <- popclim$tmx1951_1980[match(wg$Site, popclim$Code)]
# 2018 Wy climate
wg$cwd.2018wy <- popclim$cwd.2018wy[match(wg$Site, popclim$Code)]
wg$aet.2018wy <- popclim$aet.2018wy[match(wg$Site, popclim$Code)]
wg$pet.2018wy <- popclim$pet.2018wy[match(wg$Site, popclim$Code)]
wg$ppt.2018wy <- popclim$ppt.2018wy[match(wg$Site, popclim$Code)]
wg$tmn.2018wy <- popclim$tmn.2018wy[match(wg$Site, popclim$Code)]
wg$tmx.2018wy <- popclim$tmx.2018wy[match(wg$Site, popclim$Code)]
# # throwing in dormant season values 
# wg$cwd.2018ds <- popclim$cwd.2018ds[match(wg$Site, popclim$Code)]
# wg$aet.2018ds <- popclim$aet.2018ds[match(wg$Site, popclim$Code)]
# wg$pet.2018ds <- popclim$pet.2018ds[match(wg$Site, popclim$Code)]
# wg$ppt.2018ds <- popclim$ppt.2018ds[match(wg$Site, popclim$Code)]
# wg$tmn.2018ds <- popclim$tmn.2018ds[match(wg$Site, popclim$Code)]
# wg$tmx.2018ds <- popclim$tmx.2018ds[match(wg$Site, popclim$Code)]

#_________________________________________________________________________
####### * Load all trait measurements for variance decomp  ##################
#_________________________________________________________________________

# cleaned kstem and Ks values
stemk.cl <- read.csv("Data/StemK_cleaned_20240627.csv", row.names=1)
# rename Ks
names(stemk.cl)[which(names(stemk.cl)=="kstem.sa.l")] <- "Ks"
# cleaned kleaf values
leafk.cl <- read.csv("Data/LeafK_cleaned_20230428.csv",row.names = 1)
# cleaned WD values
WD <- read.csv("Data/WD_cleaned_20240428.csv", row.names=1)
# leaf and stem traits
hub <- read.csv("Data/LeafTraits_cleaned_20230428.csv", row.names = 1)


#_________________________________________________________________________
####### * Load individual-averaged traits  ##################
#_________________________________________________________________________

# aggregated from above data to the individual (with some name cleaning) and with P50 data from Skelton et al. 2019 New Phyt added
all.ind <- read.csv("Data/AllIndividuals_H-W_traits_clim_20240830.csv")

# calculate 2018-Normals climate differentials
all.ind$aet.diff <- all.ind$aet.2018wy - all.ind$aet
all.ind$pet.diff <- all.ind$pet.2018wy - all.ind$pet
all.ind$cwd.diff <- all.ind$cwd.2018wy - all.ind$cwd
all.ind$ppt.diff <- all.ind$ppt.2018wy - all.ind$ppt

# rename Ks to be easily located
names(all.ind)[which(names(all.ind)=="mkstem.sa.l")] <- "mKs"



#_________________________________________________________________________
####### * Merge data and Average traits to Population  ##################
#_________________________________________________________________________




## newer merging (post all.ind creation)
all.pop0 <- data.frame(all.ind %>% group_by(site.pop, site, pop, pop.name) %>% summarise(
  WD = mean(mWD, na.rm=T), medWD = median(mWD, na.rm=T), seWD = se(mWD), lqWD=quantile(mWD, probs = .25, na.rm=T), uqWD = quantile(mWD, probs=.75, na.rm=T),
  SLA = mean(mSLA, na.rm=T), medSLA = median(mSLA, na.rm=T), seSLA = se(mSLA), lqSLA=quantile(mSLA, probs = .25, na.rm=T), uqSLA = quantile(mSLA, probs=.75, na.rm=T),
  LDMC = mean(mLDMC, na.rm=T), medLDMC = median(mLDMC, na.rm=T), seLDMC = se(mLDMC), lqLDMC=quantile(mLDMC, probs = .25, na.rm=T), uqLDMC = quantile(mLDMC, probs=.75, na.rm=T),
  leafsize = mean(mleafsize, na.rm=T), medleafsize = median(mleafsize, na.rm=T), seleafsize = se(mleafsize), lqleafsize=quantile(mleafsize, probs = .25, na.rm=T), uqleafsize = quantile(mleafsize, probs=.75, na.rm=T),
  ml_ms = mean(mml_ms, na.rm=T), medml_ms = median(mml_ms, na.rm=T), seml_ms = se(mml_ms), lqml_ms=quantile(mml_ms, probs = .25, na.rm=T), uqml_ms = quantile(mml_ms, probs=.75, na.rm=T),
  Al_As = mean(mAl_As, na.rm=T), medAl_As = median(mAl_As, na.rm=T), seAl_As = se(mAl_As), lqAl_As=quantile(mAl_As, probs = .25, na.rm=T), uqAl_As = quantile(mAl_As, probs=.75, na.rm=T),
  kstem = mean(mkstem, na.rm=T), medkstem = median(mkstem, na.rm=T), sekstem = se(mkstem), lqkstem=quantile(mkstem, probs = .25, na.rm=T), uqkstem = quantile(mkstem, probs=.75, na.rm=T),
  kleaf = mean(mkleaf, na.rm=T), medkleaf = median(mkleaf, na.rm=T), sekleaf = se(mkleaf), lqkleaf=quantile(mkleaf, probs = .25, na.rm=T), uqkleaf = quantile(mkleaf, probs=.75, na.rm=T),
  Ks = mean(mKs, na.rm=T), medKs = median(mKs, na.rm=T), seKs = se(mKs), lqKs=quantile(mKs, probs = .25, na.rm=T), uqKs = quantile(mKs, probs=.75, na.rm=T),
  popP50leaf = mean(P50leaf, na.rm=T), medP50leaf = median(P50leaf, na.rm=T), seP50leaf = se(P50leaf), lqP50leaf = quantile(P50leaf, probs=.25, na.rm=T), upP50leaf = quantile(P50leaf, .75, na.rm=T),
  popP50stem = mean(P50stem, na.rm=T), medP50stem = median(P50stem, na.rm=T), seP50stem = se(P50stem), lqP50stem = quantile(P50stem, probs=.25, na.rm=T), upP50stem = quantile(P50stem, .75, na.rm=T)

))

# add in all the source pop climate info just because it's easiest not to carry it through
#all.pop$cwd2km <- popinfo$cwd_2km_mean[match(all.pop$pop, popinfo$pop)]
all.pop0$cwd <- popclim$cwd1951_1980[match(all.pop0$pop, popclim$Garden_pop)]
#all.pop0$soil_depth <- popclim$soil_depth[match(all.pop0$pop, popclim$Garden_pop)]
all.pop0$aet <- popclim$aet1951_1980[match(all.pop0$pop, popclim$Garden_pop)]
all.pop0$pet <- popclim$pet1951_1980[match(all.pop0$pop, popclim$Garden_pop)]
all.pop0$ppt <- popclim$ppt1951_1980[match(all.pop0$pop, popclim$Garden_pop)]
all.pop0$tmn <- popclim$tmn1951_1980[match(all.pop0$pop, popclim$Garden_pop)]
all.pop0$tmx <- popclim$tmx1951_1980[match(all.pop0$pop, popclim$Garden_pop)]
#all.pop0$avail_soil_water <- popclim$avail_soil_water[match(all.pop0$pop, popclim$Garden_pop)]
#all.pop0$dieback_10km <- popclim$proportion_10km[match(all.pop0$pop, popclim$Garden_pop)]
## add in a reverse of cwd for plotting purposes
all.pop0$cw.surp <- all.pop0$cwd *-1

# summarize wild growth to population
wg.pop <- wg %>% group_by(Site, pet) %>% summarise(DBH = mean(DBH_cm), medDBH = median(DBH_cm, na.rm=T), seDBH = se(DBH_cm), lqDBH=quantile(DBH_cm, probs = .25, na.rm=T), uqDBH = quantile(DBH_cm, probs=.75, na.rm=T),
                                              BAI = mean(BAI_cm2), , medBAI = median(BAI_cm2, na.rm=T), seBAI = se(BAI_cm2), lqBAI=quantile(BAI_cm2, probs = .25, na.rm=T), uqBAI = quantile(BAI_cm2, probs=.75, na.rm=T),
                                              growthsup = mean(growthsup_cm2), , medgrowthsup = median(growthsup_cm2, na.rm=T), segrowthsup = se(growthsup_cm2), lqgrowthsup=quantile(growthsup_cm2, probs = .25, na.rm=T), uqgrowthsup = quantile(growthsup_cm2, probs=.75, na.rm=T),
                                              mperc_maxBAI=mean(perc_maxBAI), medperc_maxBAI = median(perc_maxBAI, na.rm=T), seperc_maxBAI = se(perc_maxBAI), lqperc_maxBAI=quantile(perc_maxBAI, probs = .25, na.rm=T), uqperc_maxBAI = quantile(perc_maxBAI, probs=.75, na.rm=T))
names(wg.pop)[grep("mperc", names(wg.pop))] <- "perc_maxBAI" # get rid of m added for summarising



all.pop.m <- melt(all.pop0,id.vars =c(1:4,60:66) ) 
all.pop.wide <- cast(all.pop.m, pop + pop.name + cwd + aet + pet + ppt + tmn + tmx + cw.surp ~ site + variable, value.var="value")



### add in Hopland garden growth data for FINAL population averages
popw0 <- left_join(all.pop.wide, hoppop %>% select(1:15, "pet.td"))
### add in wild growth data
popw <- left_join(popw0, wg.pop, by=c("pop.name"="Site"))




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
################## END: LOAD DATA ######################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++










#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
################ BEGIN: ANALYSIS #######################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#______________________________________________________________________
####### * Garden growth transfer distance #############################
#______________________________________________________________________


# test for hump-shaped relationship with climate transfer distance in height
gardengrowth1 <- gamm(height~s(pet.td),random = list(pop=~1), data = hoptrees)
gardengrowth2 <- gamm(height~s(tmn.td),random = list(pop=~1), data = hoptrees)
gardengrowth3 <- gamm(height~s(ppt.td),random = list(pop=~1), data = hoptrees)
gardengrowth4 <- gamm(height~s(aet.td),random = list(pop=~1), data = hoptrees)
gardengrowth5 <- gamm(height~s(cwd.td),random = list(pop=~1), data = hoptrees)
AIC(gardengrowth1,gardengrowth2,gardengrowth3,gardengrowth4,gardengrowth5)
gam.check(gardengrowth1$gam) # things look decent from model criticism
# PET is the best climate predictor, and peaks at a transfer distance quite a ways from 0


# same but with diamater @ 50cm instead of height
# gardengrowth1.diam <- gamm(diam_50cm~s(pet.td),random = list(pop=~1), data = tmp)
# gardengrowth2.diam <- gamm(diam_50cm~s(tmn.td),random = list(pop=~1), data = tmp)
# gardengrowth3.diam <- gamm(diam_50cm~s(ppt.td),random = list(pop=~1), data = tmp)
# gardengrowth4.diam <- gamm(diam_50cm~s(aet.td),random = list(pop=~1), data = tmp)
# gardengrowth5.diam <- gamm(diam_50cm~s(cwd.td),random = list(pop=~1), data = tmp)
# AIC(gardengrowth1.diam,gardengrowth2.diam,gardengrowth3.diam,gardengrowth4.diam,gardengrowth5.diam)
# gam.check(gardengrowth1.diam$gam)
# CWD is best climate variable for 50cm diam, but shows very similar patterns
# both PET and CWD show hump shaped relationships with growth peaking at substantially drier sites
# however, we have a bigger problem with large outliers in diam





#______________________________________________________________________
######## * Population ANOVAs #################
#______________________________________________________________________

# function to calculate eta2 (among population variance)
eta2 <- function(aovres, group = "pop"){
  eta2 <- aovres[[1]]$'Sum Sq'[grep(group, rownames(aovres[[1]]))]/sum(aovres[[1]]$`Sum Sq`)
  return(round(eta2,2))
}

# function to calculate omega2 (unbiased among pop variance)
omega2 <- function(aovres, group = "pop"){
  omega2 <- (aovres[[1]]$'Sum Sq'[grep(group, rownames(aovres[[1]]))] - (aovres[[1]]$Df[grep(group, rownames(aovres[[1]]))]*aovres[[1]]$`Mean Sq`[grep("Residuals", rownames(aovres[[1]]))]))/(sum(aovres[[1]]$`Sum Sq`)+aovres[[1]]$`Mean Sq`[grep("Residuals", rownames(aovres[[1]]))])
  return(round(omega2,2))
}

# extract significance of pop factor in anova
aovsig <- function(aovres, group = "pop"){
  sig <- aovres[[1]]$'Pr(>F)'[grep(group, rownames(aovres[[1]]))]
  return(round(sig,3))
}

# just to assess the role of genetic variation vs genetic + plastic variation for different traits
# use ANOVAs seperately for W and H, with Pop as factor

###### . Hopland ANOVAS #########
all.ind$pop <- factor(all.ind$pop)
all.ind$rep <- factor(all.ind$rep)
Hop.ind <- all.ind[which(all.ind$site=="H"),]

traits.to.test <- c("SLA","LDMC","WD","ml_ms","Al_As","leafsize","kleaf","kstem","Ks","P50stem","P50leaf","Growth")
HopAOV <- data.frame(Trait=traits.to.test
                     , n = rep(NA, length(traits.to.test))
                     , eta2 = rep(NA, length(traits.to.test))
                     , omega2 = rep(NA, length(traits.to.test))
                     , sig = rep(NA, length(traits.to.test))
                     , pop.rand = rep(NA, length(traits.to.test))
                     , bestclim = rep(NA, length(traits.to.test))
                     , climDeltaAIC = rep(NA, length(traits.to.test))
                     , climsig = rep(NA, length(traits.to.test))
                     , climvar=rep(NA, length(traits.to.test))
                     , climname = rep(NA, length(traits.to.test))
                     , CV=rep(NA, length(traits.to.test)))

## SLA
HopAOV$n[which(HopAOV$Trait=="SLA")] <- length(which(Hop.ind$mSLA>0))
HopAOV$eta2[which(HopAOV$Trait=="SLA")] <- eta2(summary(aov(mSLA~pop , Hop.ind)))
HopAOV$omega2[which(HopAOV$Trait=="SLA")] <- omega2(summary(aov(mSLA~pop , Hop.ind)))
HopAOV$CV[which(HopAOV$Trait=="SLA")] <- sd(Hop.ind$mSLA, na.rm=T)/mean(Hop.ind$mSLA, na.rm=T)
HopAOV$sig[which(HopAOV$Trait=="SLA")] <- aovsig(summary(aov(mSLA~pop , Hop.ind)))
## LDMC
HopAOV$n[which(HopAOV$Trait=="LDMC")] <- length(which(Hop.ind$mLDMC>0))
HopAOV$eta2[which(HopAOV$Trait=="LDMC")] <- eta2(summary(aov(mLDMC~pop , Hop.ind)))
HopAOV$omega2[which(HopAOV$Trait=="LDMC")] <- omega2(summary(aov(mLDMC~pop , Hop.ind)))
HopAOV$CV[which(HopAOV$Trait=="LDMC")] <- sd(Hop.ind$mLDMC, na.rm=T)/mean(Hop.ind$mLDMC, na.rm=T)
HopAOV$sig[which(HopAOV$Trait=="LDMC")] <- aovsig(summary(aov(mLDMC~pop , Hop.ind)))
## WD
HopAOV$n[which(HopAOV$Trait=="WD")] <- length(which(Hop.ind$mWD>0))
HopAOV$eta2[which(HopAOV$Trait=="WD")] <- eta2(summary(aov(mWD~pop , Hop.ind)))
HopAOV$omega2[which(HopAOV$Trait=="WD")] <- omega2(summary(aov(mWD~pop , Hop.ind)))
HopAOV$CV[which(HopAOV$Trait=="WD")] <- sd(Hop.ind$mWD, na.rm=T)/mean(Hop.ind$mWD, na.rm=T)
HopAOV$sig[which(HopAOV$Trait=="WD")] <- aovsig(summary(aov(mWD~pop , Hop.ind)))
## ml_ms
HopAOV$n[which(HopAOV$Trait=="ml_ms")] <- length(which(Hop.ind$mml_ms>0))
HopAOV$eta2[which(HopAOV$Trait=="ml_ms")] <- eta2(summary(aov(mml_ms~pop , Hop.ind)))
HopAOV$omega2[which(HopAOV$Trait=="ml_ms")] <- omega2(summary(aov(mml_ms~pop , Hop.ind)))
HopAOV$CV[which(HopAOV$Trait=="ml_ms")] <- sd(Hop.ind$mml_ms, na.rm=T)/mean(Hop.ind$mml_ms, na.rm=T)
HopAOV$sig[which(HopAOV$Trait=="ml_ms")] <- aovsig(summary(aov(mml_ms~pop , Hop.ind)))

## mAl_As
HopAOV$n[which(HopAOV$Trait=="Al_As")] <- length(which(Hop.ind$mAl_As>0))
HopAOV$eta2[which(HopAOV$Trait=="Al_As")] <- eta2(summary(aov(mAl_As~pop , Hop.ind)))
HopAOV$omega2[which(HopAOV$Trait=="Al_As")] <- omega2(summary(aov(mAl_As~pop , Hop.ind)))
HopAOV$CV[which(HopAOV$Trait=="Al_As")] <- sd(Hop.ind$mAl_As, na.rm=T)/mean(Hop.ind$mAl_As, na.rm=T)
HopAOV$sig[which(HopAOV$Trait=="Al_As")] <- aovsig(summary(aov(mAl_As~pop , Hop.ind)))
## leafsize - one  decently  large outlier remains
HopAOV$n[which(HopAOV$Trait=="leafsize")] <- length(which(Hop.ind$mleafsize>0))
HopAOV$eta2[which(HopAOV$Trait=="leafsize")] <- eta2(summary(aov(mleafsize~pop , Hop.ind))) # looks somewhat non-normal (huge high outliers). but eta only goes from .14 to .18 if logged and W doesn't really look logged. leaves 7, 28, 33 as massive oultiers
HopAOV$omega2[which(HopAOV$Trait=="leafsize")] <- omega2(summary(aov(mleafsize~pop , Hop.ind)))
HopAOV$CV[which(HopAOV$Trait=="leafsize")] <- sd(Hop.ind$mleafsize, na.rm=T)/mean(Hop.ind$mleafsize, na.rm=T)
HopAOV$sig[which(HopAOV$Trait=="leafsize")] <- aovsig(summary(aov(mleafsize~pop , Hop.ind)))
## kleaf
HopAOV$n[which(HopAOV$Trait=="kleaf")] <- length(which(Hop.ind$mkleaf>0))
HopAOV$eta2[which(HopAOV$Trait=="kleaf")] <- eta2(summary(aov(mkleaf~pop , Hop.ind)))
HopAOV$omega2[which(HopAOV$Trait=="kleaf")] <- omega2(summary(aov(mkleaf~pop , Hop.ind)))
HopAOV$CV[which(HopAOV$Trait=="kleaf")] <- sd(Hop.ind$mkleaf, na.rm=T)/mean(Hop.ind$mkleaf, na.rm=T)
HopAOV$sig[which(HopAOV$Trait=="kleaf")] <- aovsig(summary(aov(mkleaf~pop , Hop.ind)))
## kstem
HopAOV$n[which(HopAOV$Trait=="kstem")] <- length(which(Hop.ind$mkstem>0))
HopAOV$eta2[which(HopAOV$Trait=="kstem")] <- eta2(summary(aov(mkstem~pop , Hop.ind)))
HopAOV$omega2[which(HopAOV$Trait=="kstem")] <- omega2(summary(aov(mkstem~pop , Hop.ind)))
HopAOV$CV[which(HopAOV$Trait=="kstem")] <- sd(Hop.ind$mkstem, na.rm=T)/mean(Hop.ind$mkstem, na.rm=T)
HopAOV$sig[which(HopAOV$Trait=="kstem")] <- aovsig(summary(aov(mkstem~pop , Hop.ind)))
## Ks
HopAOV$n[which(HopAOV$Trait=="Ks")] <- length(which(Hop.ind$mKs>0))
HopAOV$eta2[which(HopAOV$Trait=="Ks")] <- eta2(summary(aov(mKs~pop , Hop.ind)))
HopAOV$omega2[which(HopAOV$Trait=="Ks")] <- omega2(summary(aov(mKs~pop , Hop.ind)))
HopAOV$CV[which(HopAOV$Trait=="Ks")] <- sd(Hop.ind$mKs, na.rm=T)/mean(Hop.ind$mKs, na.rm=T)
HopAOV$sig[which(HopAOV$Trait=="Ks")] <- aovsig(summary(aov(mKs~pop , Hop.ind)))
## P50stem
HopAOV$n[which(HopAOV$Trait=="P50stem")] <- length(which(Hop.ind$P50stem<0))
HopAOV$eta2[which(HopAOV$Trait=="P50stem")] <- eta2(summary(aov(P50stem~pop , Hop.ind)))
HopAOV$omega2[which(HopAOV$Trait=="P50stem")] <- omega2(summary(aov(P50stem~pop , Hop.ind)))
HopAOV$CV[which(HopAOV$Trait=="P50stem")] <- sd(Hop.ind$P50stem, na.rm=T)/mean(Hop.ind$P50stem, na.rm=T)*-1
HopAOV$sig[which(HopAOV$Trait=="P50stem")] <- aovsig(summary(aov(P50stem~pop , Hop.ind)))
## P50leaf
HopAOV$n[which(HopAOV$Trait=="P50leaf")] <- length(which(Hop.ind$P50leaf<0))
HopAOV$eta2[which(HopAOV$Trait=="P50leaf")] <- eta2(summary(aov(P50leaf~pop , Hop.ind)))
HopAOV$omega2[which(HopAOV$Trait=="P50leaf")] <- omega2(summary(aov(P50leaf~pop , Hop.ind)))
HopAOV$CV[which(HopAOV$Trait=="P50leaf")] <- sd(Hop.ind$P50leaf, na.rm=T)/mean(Hop.ind$P50leaf, na.rm=T)*-1
HopAOV$sig[which(HopAOV$Trait=="P50leaf")] <- aovsig(summary(aov(P50leaf~pop , Hop.ind)))

## Growth
  # leaf rep in here, because it's statistically significant and there's replication within reps (3 trees)
HopAOV$n[which(HopAOV$Trait=="Growth")] <- length(which(hoptrees$bioest>0))
HopAOV$eta2[which(HopAOV$Trait=="Growth")] <- eta2(summary(aov(log.bio~factor(pop)+ factor(rep), hoptrees)))
HopAOV$omega2[which(HopAOV$Trait=="Growth")] <- omega2(summary(aov(log.bio~factor(pop) + factor(rep), hoptrees)))
HopAOV$CV[which(HopAOV$Trait=="Growth")] <- sd(hoptrees$log.bio, na.rm=T)/mean(hoptrees$log.bio, na.rm=T)
HopAOV$sig[which(HopAOV$Trait=="Growth")] <- aovsig(summary(aov(log.bio~factor(pop) + factor(rep), hoptrees)))

HopAOV$CV <- round(HopAOV$CV, 3)

# --NOT RUN -- double checking whether the garden blocking factor 'rep' needs to be included
# HopAOVnorep <- HopAOV
# HopAOVrep <- HopAOV
#   # rerun above ANOVA code including +rep
#   # Note on rep effect: probably reasonable to include in garden
#   # - but doesn't change the qualitative answers
#   # - not statistically significant in any trait

###### . Wild ANOVAS #########
all.ind$pop <- factor(all.ind$pop)
all.ind$rep <- factor(all.ind$rep)
Wild.ind <- all.ind[which(all.ind$site=="W"),]

traits.to.test <- c("SLA","LDMC","WD","ml_ms","Al_As","leafsize","kleaf","kstem","Ks","P50stem","P50leaf","Growth")
WildAOV <- data.frame(Trait=traits.to.test
                     , n = rep(NA, length(traits.to.test))
                     , eta2 = rep(NA, length(traits.to.test))
                     , omega2 = rep(NA, length(traits.to.test))
                     , sig = rep(NA, length(traits.to.test))
                     , pop.rand = rep(NA, length(traits.to.test))
                     , bestclim = rep(NA, length(traits.to.test))
                     , climDeltaAIC = rep(NA, length(traits.to.test))
                     , climsig = rep(NA, length(traits.to.test))
                     , climvar=rep(NA, length(traits.to.test))
                     , climname = rep(NA, length(traits.to.test))
                     , CV=rep(NA, length(traits.to.test)))



## SLA
WildAOV$n[which(WildAOV$Trait=="SLA")] <- length(which(Wild.ind$mSLA>0))
WildAOV$eta2[which(WildAOV$Trait=="SLA")] <- eta2(summary(aov(mSLA~pop , Wild.ind)))
WildAOV$omega2[which(WildAOV$Trait=="SLA")] <- omega2(summary(aov(mSLA~pop , Wild.ind)))
WildAOV$CV[which(WildAOV$Trait=="SLA")] <- sd(Wild.ind$mSLA, na.rm=T)/mean(Wild.ind$mSLA, na.rm=T)
WildAOV$sig[which(WildAOV$Trait=="SLA")] <- aovsig(summary(aov(mSLA~pop , Wild.ind)))
## LDMC
WildAOV$n[which(WildAOV$Trait=="LDMC")] <- length(which(Wild.ind$mLDMC>0))
WildAOV$eta2[which(WildAOV$Trait=="LDMC")] <- eta2(summary(aov(mLDMC~pop , Wild.ind)))
WildAOV$omega2[which(WildAOV$Trait=="LDMC")] <- omega2(summary(aov(mLDMC~pop , Wild.ind)))
WildAOV$CV[which(WildAOV$Trait=="LDMC")] <- sd(Wild.ind$mLDMC, na.rm=T)/mean(Wild.ind$mLDMC, na.rm=T)
WildAOV$sig[which(WildAOV$Trait=="LDMC")] <- aovsig(summary(aov(mLDMC~pop , Wild.ind)))
## WD
WildAOV$n[which(WildAOV$Trait=="WD")] <- length(which(Wild.ind$mWD>0))
WildAOV$eta2[which(WildAOV$Trait=="WD")] <- eta2(summary(aov(mWD~pop , Wild.ind)))
WildAOV$omega2[which(WildAOV$Trait=="WD")] <- omega2(summary(aov(mWD~pop , Wild.ind)))
WildAOV$CV[which(WildAOV$Trait=="WD")] <- sd(Wild.ind$mWD, na.rm=T)/mean(Wild.ind$mWD, na.rm=T)
WildAOV$sig[which(WildAOV$Trait=="WD")] <- aovsig(summary(aov(mWD~pop , Wild.ind)))
## mml_ms
  #make sure it doesn't need log-transformed
# qqp(resid(aov(log(mml_ms)~pop , Wild.ind))) # residuals are a bit weird, but log doesn't help
WildAOV$n[which(WildAOV$Trait=="ml_ms")] <- length(which(Wild.ind$mml_ms>0))
WildAOV$eta2[which(WildAOV$Trait=="ml_ms")] <- eta2(summary(aov(mml_ms~pop , Wild.ind)))
WildAOV$omega2[which(WildAOV$Trait=="ml_ms")] <- omega2(summary(aov(mml_ms~pop , Wild.ind)))
WildAOV$CV[which(WildAOV$Trait=="ml_ms")] <- sd(Wild.ind$mml_ms, na.rm=T)/mean(Wild.ind$mml_ms, na.rm=T)
WildAOV$sig[which(WildAOV$Trait=="ml_ms")] <- aovsig(summary(aov(mml_ms~pop , Wild.ind)))

## mAl_As
# qqp(resid(aov(mAl_As~pop , Wild.ind))) # logging kinda helps, but resids aren't horrible without, so I'll stick with raw for interpretability
WildAOV$n[which(WildAOV$Trait=="Al_As")] <- length(which(Wild.ind$mAl_As>0))
WildAOV$eta2[which(WildAOV$Trait=="Al_As")] <- eta2(summary(aov(mAl_As~pop , Wild.ind)))
WildAOV$omega2[which(WildAOV$Trait=="Al_As")] <- omega2(summary(aov(mAl_As~pop , Wild.ind)))
WildAOV$CV[which(WildAOV$Trait=="Al_As")] <- sd(Wild.ind$mAl_As, na.rm=T)/mean(Wild.ind$mAl_As, na.rm=T)
WildAOV$sig[which(WildAOV$Trait=="Al_As")] <- aovsig(summary(aov(mAl_As~pop , Wild.ind)))
## leafsize
WildAOV$n[which(WildAOV$Trait=="leafsize")] <- length(which(Wild.ind$mleafsize>0))
WildAOV$eta2[which(WildAOV$Trait=="leafsize")] <- eta2(summary(aov(mleafsize~pop , Wild.ind))) # looks somewhat non-normal (huge high outliers). but eta only goes from .14 to .18 if logged and W doesn't really look logged
WildAOV$omega2[which(WildAOV$Trait=="leafsize")] <- omega2(summary(aov(mleafsize~pop , Wild.ind)))
WildAOV$CV[which(WildAOV$Trait=="leafsize")] <- sd(Wild.ind$mleafsize, na.rm=T)/mean(Wild.ind$mleafsize, na.rm=T)
WildAOV$sig[which(WildAOV$Trait=="leafsize")] <- aovsig(summary(aov(mleafsize~pop , Wild.ind)))
## kleaf
WildAOV$n[which(WildAOV$Trait=="kleaf")] <- length(which(Wild.ind$mkleaf>0))
WildAOV$eta2[which(WildAOV$Trait=="kleaf")] <- eta2(summary(aov(mkleaf~pop , Wild.ind)))
WildAOV$omega2[which(WildAOV$Trait=="kleaf")] <- omega2(summary(aov(mkleaf~pop , Wild.ind)))
WildAOV$CV[which(WildAOV$Trait=="kleaf")] <- sd(Wild.ind$mkleaf, na.rm=T)/mean(Wild.ind$mkleaf, na.rm=T)
WildAOV$sig[which(WildAOV$Trait=="kleaf")] <- aovsig(summary(aov(mkleaf~pop , Wild.ind)))
## kstem
WildAOV$n[which(WildAOV$Trait=="kstem")] <- length(which(Wild.ind$mkstem>0))
WildAOV$eta2[which(WildAOV$Trait=="kstem")] <- eta2(summary(aov(mkstem~pop , Wild.ind)))
WildAOV$omega2[which(WildAOV$Trait=="kstem")] <- omega2(summary(aov(mkstem~pop , Wild.ind)))
WildAOV$CV[which(WildAOV$Trait=="kstem")] <- sd(Wild.ind$mkstem, na.rm=T)/mean(Wild.ind$mkstem, na.rm=T)
WildAOV$sig[which(WildAOV$Trait=="kstem")] <- aovsig(summary(aov(mkstem~pop , Wild.ind)))
## Ks
WildAOV$n[which(WildAOV$Trait=="Ks")] <- length(which(Wild.ind$mKs>0))
WildAOV$eta2[which(WildAOV$Trait=="Ks")] <- eta2(summary(aov(mKs~pop , Wild.ind)))
WildAOV$omega2[which(WildAOV$Trait=="Ks")] <- omega2(summary(aov(mKs~pop , Wild.ind)))
WildAOV$CV[which(WildAOV$Trait=="Ks")] <- sd(Wild.ind$mKs, na.rm=T)/mean(Wild.ind$mKs, na.rm=T)
WildAOV$sig[which(WildAOV$Trait=="Ks")] <- aovsig(summary(aov(mKs~pop , Wild.ind)))
## P50stem
WildAOV$n[which(WildAOV$Trait=="P50stem")] <- length(which(Wild.ind$P50stem<0))
WildAOV$eta2[which(WildAOV$Trait=="P50stem")] <- eta2(summary(aov(P50stem~pop , Wild.ind)))
WildAOV$omega2[which(WildAOV$Trait=="P50stem")] <- omega2(summary(aov(P50stem~pop , Wild.ind)))
WildAOV$CV[which(WildAOV$Trait=="P50stem")] <- sd(Wild.ind$P50stem, na.rm=T)/mean(Wild.ind$P50stem, na.rm=T)*-1
WildAOV$sig[which(WildAOV$Trait=="P50stem")] <- aovsig(summary(aov(P50stem~pop , Wild.ind)))
## P50leaf
WildAOV$n[which(WildAOV$Trait=="P50leaf")] <- length(which(Wild.ind$P50leaf<0))
WildAOV$eta2[which(WildAOV$Trait=="P50leaf")] <- eta2(summary(aov(P50leaf~pop , Wild.ind)))
WildAOV$omega2[which(WildAOV$Trait=="P50leaf")] <- omega2(summary(aov(P50leaf~pop , Wild.ind)))
WildAOV$CV[which(WildAOV$Trait=="P50leaf")] <- sd(Wild.ind$P50leaf, na.rm=T)/mean(Wild.ind$P50leaf, na.rm=T)*-1
WildAOV$sig[which(WildAOV$Trait=="P50leaf")] <- aovsig(summary(aov(P50leaf~pop , Wild.ind)))

## Growth
WildAOV$n[which(WildAOV$Trait=="Growth")] <- length(which(wg$perc_maxBAI>0 ))
WildAOV$eta2[which(WildAOV$Trait=="Growth")] <- eta2(summary(aov(perc_maxBAI~pop, wg %>% select(perc_maxBAI,pop=Site)))) # have to rename Site for eta2 function
WildAOV$omega2[which(WildAOV$Trait=="Growth")] <- omega2(summary(aov(perc_maxBAI~pop, wg %>% select(perc_maxBAI,pop=Site))))
WildAOV$CV[which(WildAOV$Trait=="Growth")] <- sd(wg$perc_maxBAI)/mean(wg$perc_maxBAI)
WildAOV$sig[which(WildAOV$Trait=="Growth")] <- aovsig(summary(aov(perc_maxBAI~pop, wg %>% select(perc_maxBAI,pop=Site))))

WildAOV$CV <- round(WildAOV$CV, 3)


#______________________________________________________________________
### * Trait-Climate Analysis ###################
#______________________________________________________________________


### create a function that tests for univariate predictors of trait values in the wild and hopland populations
# using climate 

# ## function for visualization of linear models to check reasonableness
# test.trait <- function(trait, dataz) {
#   quartz(width=11, height=4.5)
#   par(mfcol=c(2, 7), mar=c(4,0,0,0), oma=c(0,4,2,1), mgp=c(2.2,1,0))
#   for(i in c("cwd","aet","pet","ppt","aet.diff","pet.diff","cwd.diff")){
#     modW <- lm(get(trait)~get(i), dataz[which(dataz$site=="W"),])
#     modH <- lm(get(trait)~get(i), dataz[which(dataz$site=="H"),])
#     plot(get(trait)~get(i), dataz[which(dataz$site=="H"),], pch=1,xlab="", yaxt="n", ylab=trait, col=ifelse(summary(modH)$coefficients["get(i)",4]<0.05,yes = "red","black") )
#     if(i=="cwd") {mtext(trait, side=2, line=2)
#       axis(2)
#       mtext("Hopland", side=3, adj = 0)}
#     #    mtext(text = paste("site=", round(test[[1]]$`Pr(>F)`[1],2), " pop=", round(test[[1]]$`Pr(>F)`[[2]])),side=3)
#     plot(get(trait)~get(i), dataz[which(dataz$site=="W"),], pch=3,xlab=i, yaxt="n", ylab=trait, col=ifelse(summary(modW)$coefficients["get(i)",4]<0.05,yes = "red","black")  )
#     if(i=="cwd") {mtext(trait, side=2, line=2)
#       axis(2)
#       mtext("Wild", side=3, adj = 0)}
#     
#   }
#   print("WILD:")
#   print(summary(aov(get(trait)~ pop, all.ind[which(all.ind$site=="W"),])) )
#   print("HOPLAND:")
#   print(summary(aov(get(trait)~ pop, all.ind[which(all.ind$site=="H"),])) )
# }
# 
# ## function for visualization of linear mixed modesl to check reasonableness
# test.mertrait <- function(trait, dataz) {
#   quartz(width=11, height=4.5)
#   par(mfcol=c(2, 7), mar=c(4,0,0,0), oma=c(0,4,2,1), mgp=c(2.2,1,0))
#   for(i in c("cwd","aet","pet","ppt","sitePD","siteMD","siteE.drop")){
#     modW <- lmer(get(trait)~get(i) + (1|pop.name), dataz[which(dataz$site=="W"),])
#     modH <- lmer(get(trait)~get(i)+ (1|pop.name), dataz[which(dataz$site=="H"),])
#     plot(get(trait)~get(i), dataz[which(dataz$site=="H"),], pch=1,xlab="", yaxt="n", ylab=trait, col=ifelse(summary(modH)$coefficients["get(i)",5]<0.05,yes = "red","black") )
#     if(i=="cwd") {mtext(trait, side=2, line=2)
#       axis(2)
#       mtext("Hopland", side=3, adj = 0)}
#     #    mtext(text = paste("site=", round(test[[1]]$`Pr(>F)`[1],2), " pop=", round(test[[1]]$`Pr(>F)`[[2]])),side=3)
#     plot(get(trait)~get(i), dataz[which(dataz$site=="W"),], pch=3,xlab=i, yaxt="n", ylab=trait, col=ifelse(summary(modW)$coefficients["get(i)",5]<0.05,yes = "red","black")  )
#     if(i=="cwd") {mtext(trait, side=2, line=2)
#       axis(2)
#       mtext("Wild", side=3, adj = 0)}
#     
#   }
#   print(test <- summary(aov(get(trait)~ pop, all.ind[which(all.ind$site=="W"),])) )
#   print(test <- summary(aov(get(trait)~ pop, all.ind[which(all.ind$site=="H"),])) )
# }

## function to test: Do we need pop random effects?
  # to double check whether we need to fit linear or linear mixed models
test.rand <- function(dataz, trait){
  require(nlme)
  datazz <- dataz[which(abs(dataz[,trait])>0),]
  form <- paste(as.character(trait),"~ppt + pet + tmn")
  null <- gls(model = as.formula(form), datazz)
  pop <- lme(fixed = as.formula(form), random= ~1|pop.name, data=datazz)
  res <-anova(null, pop)$`p-value`[2]
  if(res<0.05){return("yes")}
  else{return("no")}
}

# ## select best climate predictor for Wild data using linear models
#___________ Newer Function _____________

best.mermod <- function(trait, dataz, clim.vars) {
  results <- data.frame("mod"=rep(NA, times=length(clim.vars)+1)
                        , "n"=rep(NA, times=length(clim.vars)+1)
                        , "AICc"=rep(NA, times=length(clim.vars)+1)
                        , "AIC"=rep(NA, times=length(clim.vars)+1)
                        ,"BIC"=rep(NA, times=length(clim.vars)+1)
                        ,"df"=rep(NA, times=length(clim.vars)+1)
  )
  null <- lmer(get(trait)~ 1 + (1|pop), dataz, REML=F)
  results$mod[1] <- "null"
  results$n[1] <- nrow(null@frame)
  results$AICc[1] <- AICc(null)
  results$AIC[1] <- summary(null)$AICtab[1]
  results$BIC[1] <- summary(null)$AICtab[2]
  results$df[1] <- summary(null)$AICtab[5]
  for (i in 1:length(clim.vars)){
    mod <- lmer(get(trait)~ get(clim.vars[i]) + (1|pop), dataz, REML=F)
    results$mod[i+1] <- clim.vars[i]
    results$n[i+1] <- nrow(mod@frame)
    results$AICc[i+1] <- AICc(mod)
    results$AIC[i+1] <- summary(mod)$AICtab[1]
    results$BIC[i+1] <- summary(mod)$AICtab[2]
    results$df[i+1] <- summary(mod)$AICtab[5]
  }
  res <- list("table" = results %>% arrange(AICc), "deltaAIC" =  NA)
  res$deltaAIC <-round(res$table$AICc[1] - res$table$AICc[which(res$table$mod=="null")], 2)
  return(res)
}


best.mod <- function(trait, dataz, clim.vars) {
  results <- data.frame("mod"=rep(NA, times=length(clim.vars)+1)
                        , "n"=rep(NA, times=length(clim.vars)+1)
                        , "AICc"=rep(NA, times=length(clim.vars)+1)
                        , "AIC"=rep(NA, times=length(clim.vars)+1)
                        ,"BIC"=rep(NA, times=length(clim.vars)+1)
                        ,"df"=rep(NA, times=length(clim.vars)+1)
  )
  null <- lm(get(trait)~ 1 , dataz)
  results$mod[1] <- "null"
  results$n[1] <- nrow(null$model)
  results$AICc[1] <- AICc(null)
  results$AIC[1] <- AIC(null)
  results$BIC[1] <-BIC(null)
  results$df[1] <- null$df.residual
  for (i in 1:length(clim.vars)){
    mod <- lm(get(trait)~ get(clim.vars[i]), dataz)
    results$mod[i+1] <- clim.vars[i]
    results$n[i+1] <- nrow(mod$model)
    results$AICc[i+1] <- AICc(mod)
    results$AIC[i+1] <- AIC(mod)
    results$BIC[i+1] <- BIC(mod)
    results$df[i+1] <- mod$df.residual
  }
  res <- list("table" = results %>% arrange(AICc), "deltaAIC" =  NA)
  res$deltaAIC <-round(res$table$AICc[1] - res$table$AICc[which(res$table$mod=="null")], 2)
  return(res)
}


# climate variables for wild populations (climate normals, 2018 dormant season, 2018-normals anamoly)
clim.variablesW <- c("cwd","aet","pet","ppt","tmn","cwd.2018ds","aet.2018ds","pet.2018ds","ppt.2018ds","tmn.2018ds","tmx.2018ds", "cwd.diff","aet.diff","pet.diff","ppt.diff")
# climate variables for the garden (just climate normals)
clim.variablesH <- c("cwd","aet","pet","ppt", "tmn")


#_______________________________________________________________________________
############ . identify best Trait-climate relationships ############################
#_______________________________________________________________________________


##### determining if population random effects are needed
HopAOV$pop.rand[which(HopAOV$Trait=="SLA")] <- test.rand(dataz=Hop.ind, trait="mSLA")
HopAOV$pop.rand[which(HopAOV$Trait=="LDMC")] <- test.rand(dataz=Hop.ind, trait="mLDMC")
HopAOV$pop.rand[which(HopAOV$Trait=="WD")] <- test.rand(dataz=Hop.ind, trait="mWD")
HopAOV$pop.rand[which(HopAOV$Trait=="ml_ms")] <- test.rand(dataz=Hop.ind, trait="mml_ms")
HopAOV$pop.rand[which(HopAOV$Trait=="Al_As")] <- test.rand(dataz=Hop.ind, trait="mAl_As")
HopAOV$pop.rand[which(HopAOV$Trait=="leafsize")] <- test.rand(dataz=Hop.ind, trait="mleafsize")
HopAOV$pop.rand[which(HopAOV$Trait=="kleaf")] <- test.rand(dataz=Hop.ind, trait="mkleaf")
HopAOV$pop.rand[which(HopAOV$Trait=="kstem")] <- test.rand(dataz=Hop.ind, trait="mkstem")
HopAOV$pop.rand[which(HopAOV$Trait=="Ks")] <- test.rand(dataz=Hop.ind, trait="mKs")
HopAOV$pop.rand[which(HopAOV$Trait=="P50stem")] <- test.rand(dataz=Hop.ind, trait="P50stem")
HopAOV$pop.rand[which(HopAOV$Trait=="P50leaf")] <- test.rand(dataz=Hop.ind, trait="P50leaf")
HopAOV$pop.rand[which(HopAOV$Trait=="Growth")] <- "yes" # can't use function, have to do manually
# null <- gls(model = log.bio~ppt + pet + tmn, hoptrees[which(hoptrees$bioest>0),])
# pop <- lme(fixed = log.bio~ppt + pet + tmn, random= ~1|factor(pop), data=hoptrees[which(hoptrees$bioest>0),])
# res <-anova(null, pop)$`p-value`[2] # DEFO need pop random intercepts for log.bio

  # so no hopland traits require random intercepts for pop

WildAOV$pop.rand[which(WildAOV$Trait=="SLA")] <- test.rand(dataz=Wild.ind, trait="mSLA")
WildAOV$pop.rand[which(WildAOV$Trait=="LDMC")] <- test.rand(dataz=Wild.ind, trait="mLDMC")
WildAOV$pop.rand[which(WildAOV$Trait=="WD")] <- test.rand(dataz=Wild.ind, trait="mWD")
WildAOV$pop.rand[which(WildAOV$Trait=="ml_ms")] <- test.rand(dataz=Wild.ind, trait="mml_ms")
WildAOV$pop.rand[which(WildAOV$Trait=="Al_As")] <- test.rand(dataz=Wild.ind, trait="mAl_As")
WildAOV$pop.rand[which(WildAOV$Trait=="leafsize")] <- test.rand(dataz=Wild.ind, trait="mleafsize")
WildAOV$pop.rand[which(WildAOV$Trait=="kleaf")] <- test.rand(dataz=Wild.ind, trait="mkleaf")
WildAOV$pop.rand[which(WildAOV$Trait=="kstem")] <- test.rand(dataz=Wild.ind, trait="mkstem")
WildAOV$pop.rand[which(WildAOV$Trait=="Ks")] <- test.rand(dataz=Wild.ind, trait="mKs")
WildAOV$pop.rand[which(WildAOV$Trait=="P50stem")] <- test.rand(dataz=Wild.ind, trait="P50stem")
WildAOV$pop.rand[which(WildAOV$Trait=="P50leaf")] <- test.rand(dataz=Wild.ind, trait="P50leaf")
WildAOV$pop.rand[which(WildAOV$Trait=="Growth")] <- "yes" # definitely required

  ## test for Growth
# null <- gls(model = perc_maxBAI~ppt + pet + tmn, wg)
# pop <- lme(fixed = perc_maxBAI~ppt + pet + tmn, random= ~1|Site, data=wg)
# res <-anova(null, pop)$`p-value`[2] # need pop random intercepts for perc_maxBAI

  # all wild traits except ml_ms, Ks and P50stem/P50stem require random intercepts for pop

## Filling ANOVA tables:

##### . SLA #######
# Hopland
(HopSLA <- best.mermod("mSLA",dataz=Hop.ind,clim.vars = clim.variablesH))
HopAOV$climname[which(HopAOV$Trait=="SLA")] <- HopSLA$table$mod[1] # grab name of best climate variable from AIC table
HopAOV$climDeltaAIC[which(HopAOV$Trait=="SLA")] <- HopSLA$deltaAIC # grab name of best climate variable from AIC table
HopmodSLA <- lm(mSLA~pet, Hop.ind)
# plot(HopmodSLA)
# PET
HopAOV$bestclim[which(HopAOV$Trait=="SLA")] <- "PET[30yr]"
HopAOV$climvar[which(HopAOV$Trait=="SLA")] <- round(r.squaredLR(HopmodSLA),2)
HopAOV$climsig[which(HopAOV$Trait=="SLA")] <- round(summary(HopmodSLA)$coefficients[2,4],4)
# Wild
#best.mod("mSLA",dataz=Wild.ind)
#WildmodSLA <- lm(mSLA~tmx.2018ds, Wild.ind)
# plot(WildmodSLA)
(WildSLA<- best.mermod("mSLA",dataz=Wild.ind, clim.vars = clim.variablesW))
  # with linear model,tmn.2018ds is vaguely important but prob ns. with lmer, null is best
# tmx/tmin/ppt all marginally significant
WildAOV$climname[which(WildAOV$Trait=="SLA")] <- WildSLA$table$mod[1] # grab name of best climate variable from AIC table
WildAOV$climDeltaAIC[which(WildAOV$Trait=="SLA")] <- WildSLA$deltaAIC # grab delta from AIC table
WildmodSLA <- lmer(mSLA~pet + (1|pop), Wild.ind, REML=T)
WildmodSLA <- lmer(mSLA~pet*site + (1|pop), all.ind, REML=T)

# plot(WildmodSLA)
# PET
WildAOV$bestclim[which(WildAOV$Trait=="SLA")] <- "AET[gy]"
WildAOV$climvar[which(WildAOV$Trait=="SLA")] <- round(r.squaredGLMM(WildmodSLA)[1],2)
WildAOV$climsig[which(WildAOV$Trait=="SLA")] <- round(summary(WildmodSLA)$coefficients[2,5],4)


##### . LDMC #######
# Hopland
(HopLDMC <- best.mod("mLDMC",dataz=Hop.ind,clim.vars = clim.variablesH))
HopAOV$climname[which(HopAOV$Trait=="LDMC")] <- HopLDMC$table$mod[1] # grab name of best climate variable from AIC table
HopAOV$climDeltaAIC[which(HopAOV$Trait=="LDMC")] <- HopLDMC$deltaAIC # grab name of best climate variable from AIC table
#HopmodLDMC <- lm(mLDMC~pet, Hop.ind)
# plot(HopmodLDMC)
HopAOV$bestclim[which(HopAOV$Trait=="LDMC")] <- "none"
HopAOV$climvar[which(HopAOV$Trait=="LDMC")] <- 0
HopAOV$climsig[which(HopAOV$Trait=="LDMC")] <- 1

# Wild
(WildLDMC<- best.mermod("mLDMC",dataz=Wild.ind, clim.vars = clim.variablesW))
# with linear model,tmn.2018ds is vaguely important but prob ns. with lmer, null is best
# tmx/tmin/ppt all marginally significant
WildAOV$climname[which(WildAOV$Trait=="LDMC")] <- WildLDMC$table$mod[1] # grab name of best climate variable from AIC table
WildAOV$climDeltaAIC[which(WildAOV$Trait=="LDMC")] <- WildLDMC$deltaAIC # grab delta from AIC table
#WildmodLDMC <- lmer(mLDMC~pet + (1|pop), Wild.ind, REML=T)
# plot(WildmodLDMC)
# PET
WildAOV$bestclim[which(WildAOV$Trait=="LDMC")] <- "none"
WildAOV$climvar[which(WildAOV$Trait=="LDMC")] <- 0 #round(r.squaredGLMM(WildmodLDMC)[1],2)
WildAOV$climsig[which(WildAOV$Trait=="LDMC")] <- 1 #round(summary(WildmodLDMC)$coefficients[2,5],4)




##### . WD #######
# Hopland
(HopWD <- best.mermod("mWD",dataz=Hop.ind,clim.vars = clim.variablesH))
HopAOV$climname[which(HopAOV$Trait=="WD")] <- HopWD$table$mod[1] # grab name of best climate variable from AIC table
HopAOV$climDeltaAIC[which(HopAOV$Trait=="WD")] <- HopWD$deltaAIC # grab name of best climate variable from AIC table
HopmodWD <- lm(mWD~ppt, Hop.ind)
# plot(HopmodWD)
HopAOV$bestclim[which(HopAOV$Trait=="WD")] <- "PPT[30yr]"
HopAOV$climvar[which(HopAOV$Trait=="WD")] <- round(r.squaredLR(HopmodWD),2)
HopAOV$climsig[which(HopAOV$Trait=="WD")] <- round(summary(HopmodWD)$coefficients[2,4],4)
# Wild
#best.mod("mWD",dataz=Wild.ind)
#WildmodWD <- lm(mWD~tmx.2018ds, Wild.ind)
# plot(WildmodWD)
(WildWD<- best.mermod("mWD",dataz=Wild.ind, clim.vars = clim.variablesW))
# anomolies are best, but don't really make sense because the branches are multiple years of wood...
WildAOV$climname[which(WildAOV$Trait=="WD")] <- WildWD$table$mod[1] # grab name of best climate variable from AIC table
WildAOV$climDeltaAIC[which(WildAOV$Trait=="WD")] <- WildWD$deltaAIC # grab delta from AIC table
WildmodWD <- lmer(mWD~aet.diff + (1|pop), Wild.ind, REML=T)
# plot(WildmodWD)
# qqp(resid(WildmodWD))
#qqp(ranef(WildmodWD)$pop[[1]])
WildAOV$bestclim[which(WildAOV$Trait=="WD")] <- "AET[anom]"
WildAOV$climvar[which(WildAOV$Trait=="WD")] <- round(r.squaredGLMM(WildmodWD)[1],2)
WildAOV$climsig[which(WildAOV$Trait=="WD")] <- round(summary(WildmodWD)$coefficients[2,5],4)




##### . ml_ms #######
# Hopland
(Hopml_ms <- best.mod("mml_ms",dataz=Hop.ind,clim.vars = clim.variablesH))
HopAOV$climname[which(HopAOV$Trait=="ml_ms")] <- Hopml_ms$table$mod[1] # grab name of best climate variable from AIC table
HopAOV$climDeltaAIC[which(HopAOV$Trait=="ml_ms")] <- Hopml_ms$deltaAIC # grab name of best climate variable from AIC table
#Hopmodml_ms <- lm(mml_ms~pet, Hop.ind)
# plot(Hopmodml_ms)
HopAOV$bestclim[which(HopAOV$Trait=="ml_ms")] <- "none"
HopAOV$climvar[which(HopAOV$Trait=="ml_ms")] <- 0 #round(r.squaredLR(Hopmodml_ms),2)
HopAOV$climsig[which(HopAOV$Trait=="ml_ms")] <- 1 #round(summary(Hopmodml_ms)$coefficients[2,4],4)
# Wild
(Wildml_ms<- best.mermod("mml_ms",dataz=Wild.ind, clim.vars = clim.variablesW))
WildAOV$climname[which(WildAOV$Trait=="ml_ms")] <- Wildml_ms$table$mod[1] # grab name of best climate variable from AIC table
WildAOV$climDeltaAIC[which(WildAOV$Trait=="ml_ms")] <- Wildml_ms$deltaAIC # grab delta from AIC table
Wildmodml_ms <- lm(mml_ms~tmn.2018ds, Wild.ind)
# plot(Wildmodml_ms)
# PET
WildAOV$bestclim[which(WildAOV$Trait=="ml_ms")] <- "Tmin[gy]"
WildAOV$climvar[which(WildAOV$Trait=="ml_ms")] <- round(r.squaredGLMM(Wildmodml_ms)[1],2)
WildAOV$climsig[which(WildAOV$Trait=="ml_ms")] <- round(summary(Wildmodml_ms)$coefficients[2,4],4)



##### . Al_As #######
# Hopland
(HopAl_As <- best.mod("mAl_As",dataz=Hop.ind,clim.vars = clim.variablesH))
HopAOV$climname[which(HopAOV$Trait=="Al_As")] <- HopAl_As$table$mod[1] # grab name of best climate variable from AIC table
HopAOV$climDeltaAIC[which(HopAOV$Trait=="Al_As")] <- HopAl_As$deltaAIC # grab name of best climate variable from AIC table
#HopmodAl_As <- lm(mAl_As~pet, Hop.ind)
# plot(HopmodAl_As)

HopAOV$bestclim[which(HopAOV$Trait=="Al_As")] <- "none"
HopAOV$climvar[which(HopAOV$Trait=="Al_As")] <- 0 #round(r.squaredLR(HopmodAl_As),2)
HopAOV$climsig[which(HopAOV$Trait=="Al_As")] <- 1 #round(summary(HopmodAl_As)$coefficients[2,4],4)
# Wild
(WildAl_As<- best.mermod("mAl_As",dataz=Wild.ind, clim.vars = clim.variablesW))
WildAOV$climname[which(WildAOV$Trait=="Al_As")] <- WildAl_As$table$mod[1] # grab name of best climate variable from AIC table
WildAOV$climDeltaAIC[which(WildAOV$Trait=="Al_As")] <- WildAl_As$deltaAIC # grab delta from AIC table
WildmodAl_As <- lmer(mAl_As~cwd.diff + (1|pop), Wild.ind, REML=T)
# plot(WildmodAl_As)
WildAOV$bestclim[which(WildAOV$Trait=="Al_As")] <- "CWD[anom]"
WildAOV$climvar[which(WildAOV$Trait=="Al_As")] <- round(r.squaredGLMM(WildmodAl_As)[1],2)
WildAOV$climsig[which(WildAOV$Trait=="Al_As")] <- round(summary(WildmodAl_As)$coefficients[2,5],4)



##### . leafsize #######
# Hopland
(Hopleafsize <- best.mermod("mleafsize",dataz=Hop.ind,clim.vars = clim.variablesH))
HopAOV$climname[which(HopAOV$Trait=="leafsize")] <- Hopleafsize$table$mod[1] # grab name of best climate variable from AIC table
HopAOV$climDeltaAIC[which(HopAOV$Trait=="leafsize")] <- Hopleafsize$deltaAIC # grab name of best climate variable from AIC table
Hopmodleafsize <- lm(mleafsize~ppt, Hop.ind)
# plot(Hopmodleafsize)

HopAOV$bestclim[which(HopAOV$Trait=="leafsize")] <- "PPT[30yr]"
HopAOV$climvar[which(HopAOV$Trait=="leafsize")] <- round(r.squaredLR(Hopmodleafsize),2)
HopAOV$climsig[which(HopAOV$Trait=="leafsize")] <- round(summary(Hopmodleafsize)$coefficients[2,4],4)
# Wild
(Wildleafsize<- best.mermod("mleafsize",dataz=Wild.ind, clim.vars = clim.variablesW))
WildAOV$climname[which(WildAOV$Trait=="leafsize")] <- Wildleafsize$table$mod[1] # grab name of best climate variable from AIC table
WildAOV$climDeltaAIC[which(WildAOV$Trait=="leafsize")] <- Wildleafsize$deltaAIC # grab delta from AIC table
Wildmodleafsize <- lmer(mleafsize~tmn.2018ds + (1|pop), Wild.ind, REML=T)
# plot(Wildmodleafsize)
WildAOV$bestclim[which(WildAOV$Trait=="leafsize")] <- "Tmin[gy]"
WildAOV$climvar[which(WildAOV$Trait=="leafsize")] <- round(r.squaredGLMM(Wildmodleafsize)[1],2)
WildAOV$climsig[which(WildAOV$Trait=="leafsize")] <- round(summary(Wildmodleafsize)$coefficients[2,5],4)




##### . kleaf #######
# Hopland
(Hopkleaf <- best.mod("mkleaf",dataz=Hop.ind,clim.vars = clim.variablesH))
HopAOV$climname[which(HopAOV$Trait=="kleaf")] <- Hopkleaf$table$mod[1] # grab name of best climate variable from AIC table
HopAOV$climDeltaAIC[which(HopAOV$Trait=="kleaf")] <- Hopkleaf$deltaAIC # grab name of best climate variable from AIC table
Hopmodkleaf <- lm(mkleaf~tmn, Hop.ind)
# plot(Hopmodkleaf)

HopAOV$bestclim[which(HopAOV$Trait=="kleaf")] <- "Tmin[30yr]"
HopAOV$climvar[which(HopAOV$Trait=="kleaf")] <- round(r.squaredLR(Hopmodkleaf),2)
HopAOV$climsig[which(HopAOV$Trait=="kleaf")] <- round(summary(Hopmodkleaf)$coefficients[2,4],4)
# Wild
(Wildkleaf<- best.mermod("mkleaf",dataz=Wild.ind, clim.vars = clim.variablesW))
WildAOV$climname[which(WildAOV$Trait=="kleaf")] <- Wildkleaf$table$mod[1] # grab name of best climate variable from AIC table
WildAOV$climDeltaAIC[which(WildAOV$Trait=="kleaf")] <- Wildkleaf$deltaAIC # grab delta from AIC table
Wildmodkleaf <- lmer(mkleaf~cwd.2018ds + (1|pop), Wild.ind, REML=T)
# plot(Wildmodkleaf)
WildAOV$bestclim[which(WildAOV$Trait=="kleaf")] <- "CWD[gy]"
WildAOV$climvar[which(WildAOV$Trait=="kleaf")] <- round(r.squaredGLMM(Wildmodkleaf)[1],2)
WildAOV$climsig[which(WildAOV$Trait=="kleaf")] <- round(summary(Wildmodkleaf)$coefficients[2,5],4)



# ##### . kstem #######
# Hopland
(Hopkstem <- best.mod("mkstem",dataz=Hop.ind,clim.vars = clim.variablesH))
HopAOV$climname[which(HopAOV$Trait=="kstem")] <- Hopkstem$table$mod[1] # grab name of best climate variable from AIC table
HopAOV$climDeltaAIC[which(HopAOV$Trait=="kstem")] <- Hopkstem$deltaAIC # grab name of best climate variable from AIC table
Hopmodkstem <- lm(mkstem~pet, Hop.ind)
# plot(Hopmodkstem)

HopAOV$bestclim[which(HopAOV$Trait=="kstem")] <- "PET[30yr]"
HopAOV$climvar[which(HopAOV$Trait=="kstem")] <- round(r.squaredLR(Hopmodkstem),2)
HopAOV$climsig[which(HopAOV$Trait=="kstem")] <- round(summary(Hopmodkstem)$coefficients[2,4],4)
# Wild
(Wildkstem<- best.mermod("mkstem",dataz=Wild.ind, clim.vars = clim.variablesW))
WildAOV$climname[which(WildAOV$Trait=="kstem")] <- Wildkstem$table$mod[1] # grab name of best climate variable from AIC table
WildAOV$climDeltaAIC[which(WildAOV$Trait=="kstem")] <- Wildkstem$deltaAIC # grab delta from AIC table
Wildmodkstem <- lmer(mkstem~tmn.2018ds + (1|pop), Wild.ind, REML=T)
# plot(Wildmodkstem)
WildAOV$bestclim[which(WildAOV$Trait=="kstem")] <- "Tmin[gy]"
WildAOV$climvar[which(WildAOV$Trait=="kstem")] <- round(r.squaredGLMM(Wildmodkstem)[1],2)
WildAOV$climsig[which(WildAOV$Trait=="kstem")] <- round(summary(Wildmodkstem)$coefficients[2,5],4)




##### . Ks #######
# Hopland
(HopKs <- best.mod("mKs",dataz=Hop.ind,clim.vars = clim.variablesH))
HopAOV$climname[which(HopAOV$Trait=="Ks")] <- HopKs$table$mod[1] # grab name of best climate variable from AIC table
HopAOV$climDeltaAIC[which(HopAOV$Trait=="Ks")] <- HopKs$deltaAIC # grab name of best climate variable from AIC table
#HopmodKs <- lm(mKs~pet, Hop.ind)
# plot(HopmodKs)
HopAOV$bestclim[which(HopAOV$Trait=="Ks")] <- "none"
HopAOV$climvar[which(HopAOV$Trait=="Ks")] <- 0 #round(r.squaredLR(HopmodKs),2)
HopAOV$climsig[which(HopAOV$Trait=="Ks")] <- 1 #round(summary(HopmodKs)$coefficients[2,4],4)
# Wild
(WildKs<- best.mod("mKs",dataz=Wild.ind, clim.vars = clim.variablesW))
WildAOV$climname[which(WildAOV$Trait=="Ks")] <- WildKs$table$mod[1] # grab name of best climate variable from AIC table
WildAOV$climDeltaAIC[which(WildAOV$Trait=="Ks")] <- WildKs$deltaAIC # grab delta from AIC table
#WildmodKs <- lm(mKs~tmn.2018ds, Wild.ind)
# plot(WildmodKs)

WildAOV$bestclim[which(WildAOV$Trait=="Ks")] <- "none"
WildAOV$climvar[which(WildAOV$Trait=="Ks")] <- 0 #round(r.squaredGLMM(WildmodKs)[1],2)
WildAOV$climsig[which(WildAOV$Trait=="Ks")] <- 1 #round(summary(WildmodKs)$coefficients[2,4],4)



##### . P50stem #######
##### . P50stem #######
# Hopland
(HopP50stem <- best.mod("P50stem",dataz=Hop.ind,clim.vars = clim.variablesH))
HopAOV$climname[which(HopAOV$Trait=="P50stem")] <- HopP50stem$table$mod[1] # grab name of best climate variable from AIC table
HopAOV$climDeltaAIC[which(HopAOV$Trait=="P50stem")] <- HopP50stem$deltaAIC # grab name of best climate variable from AIC table
#HopmodP50stem <- lm(P50stem~tmn, Hop.ind)
# plot(HopmodP50stem)
HopAOV$bestclim[which(HopAOV$Trait=="P50stem")] <- "none"
HopAOV$climvar[which(HopAOV$Trait=="P50stem")] <- 0 #round(r.squaredLR(HopmodP50stem),2)
HopAOV$climsig[which(HopAOV$Trait=="P50stem")] <- 1 #round(summary(HopmodP50stem)$coefficients[2,4],4)
# Wild
(WildP50stem<- best.mod("P50stem",dataz=Wild.ind, clim.vars = clim.variablesW))
WildAOV$climname[which(WildAOV$Trait=="P50stem")] <- WildP50stem$table$mod[1] # grab name of best climate variable from AIC table
WildAOV$climDeltaAIC[which(WildAOV$Trait=="P50stem")] <- WildP50stem$deltaAIC # grab delta from AIC table
WildmodP50stem <- lm(P50stem~aet, Wild.ind)
# plot(WildmodP50stem) # 16 is a bit of an outlier, but still ms without it
WildAOV$bestclim[which(WildAOV$Trait=="P50stem")] <- "AET[30yr]"
WildAOV$climvar[which(WildAOV$Trait=="P50stem")] <- round(r.squaredGLMM(WildmodP50stem)[1],2)
WildAOV$climsig[which(WildAOV$Trait=="P50stem")] <- round(summary(WildmodP50stem)$coefficients[2,4],4)



##### . P50leaf #######
# Hopland
(HopP50leaf <- best.mod("P50leaf",dataz=Hop.ind,clim.vars = clim.variablesH))
HopAOV$climname[which(HopAOV$Trait=="P50leaf")] <- HopP50leaf$table$mod[1] # grab name of best climate variable from AIC table
HopAOV$climDeltaAIC[which(HopAOV$Trait=="P50leaf")] <- HopP50leaf$deltaAIC # grab name of best climate variable from AIC table
#HopmodP50leaf <- lm(mP50leaf~pet, Hop.ind)
# plot(HopmodP50leaf)
HopAOV$bestclim[which(HopAOV$Trait=="P50leaf")] <- "none"
HopAOV$climvar[which(HopAOV$Trait=="P50leaf")] <- 0 #round(r.squaredLR(HopmodP50leaf),2)
HopAOV$climsig[which(HopAOV$Trait=="P50leaf")] <- 1 #round(summary(HopmodP50leaf)$coefficients[2,4],4)
# Wild
(WildP50leaf<- best.mod("P50leaf",dataz=Wild.ind, clim.vars = clim.variablesW))
WildAOV$climname[which(WildAOV$Trait=="P50leaf")] <- WildP50leaf$table$mod[1] # grab name of best climate variable from AIC table
WildAOV$climDeltaAIC[which(WildAOV$Trait=="P50leaf")] <- WildP50leaf$deltaAIC # grab delta from AIC table
# WildmodP50leaf <- lm(mP50leaf~tmn.2018ds, Wild.ind)
# plot(WildmodP50leaf)

WildAOV$bestclim[which(WildAOV$Trait=="P50leaf")] <- "none"
WildAOV$climvar[which(WildAOV$Trait=="P50leaf")] <- 0 #round(r.squaredGLMM(WildmodP50leaf)[1],2)
WildAOV$climsig[which(WildAOV$Trait=="P50leaf")] <- 1 #round(summary(WildmodP50leaf)$coefficients[2,4],4)


##### . Growth #######
# Hopland
hoptrees$pop.name <- factor(hoptrees$pop)

modnull <- lmer(log.bio~1 + (1|pop.name),hoptrees, REML=F)
modcwd <- lmer(log.bio~cwd+ (1|pop.name), hoptrees, REML=F)
modaet <- lmer(log.bio~aet+ (1|pop.name), hoptrees, REML=F)
modpet <- lmer(log.bio~pet+ (1|pop.name), hoptrees, REML=F)
modppt <- lmer(log.bio~ppt+ (1|pop.name), hoptrees, REML=F)
modtmn <- lmer(log.bio~tmn+ (1|pop.name), hoptrees, REML=F)
HopGrowth <- data.frame(AICc(modnull, modcwd, modaet, modpet, modppt, modtmn)
           , n=c(length(resid(modnull)), length(resid(modcwd)), length(resid(modaet)), length(resid(modpet)), length(resid(modppt)), length(resid(modtmn))))[order(AICc(modnull, modcwd, modaet, modpet, modppt, modtmn)[,2]),]


HopmodGrowth<- lmer(log.bio~pet + (1|pop.name), hoptrees) 
#qqp(resid(HopmodGrowth)) #some major low outliers. but generally not horrible
#qqp(random.effects(HopmodGrowth)$pop.name[,1]) # random effects are decent
# pet
HopAOV$climDeltaAIC[which(HopAOV$Trait=="Growth")] <- HopGrowth$AICc[1] - HopGrowth$AICc[which(rownames(HopGrowth)=="modnull")]
HopAOV$climname[which(HopAOV$Trait=="Growth")] <- "pet"
HopAOV$bestclim[which(HopAOV$Trait=="Growth")] <- "PET[30yr]"
HopAOV$climvar[which(HopAOV$Trait=="Growth")] <- round(r.squaredGLMM(HopmodGrowth)[1],2)
HopAOV$climsig[which(HopAOV$Trait=="Growth")] <- round(summary(HopmodGrowth)$coefficients[2,5],4)


# Wild
  # for growth, current year values less relevant because we're using 5yrs worth of growth data
modnull <- lmer(perc_maxBAI~1 + (1|Site),wg, REML=F)
modcwd <- lmer(perc_maxBAI~cwd+ (1|Site), wg, REML=F)
modaet <- lmer(perc_maxBAI~aet+ (1|Site), wg, REML=F)
modpet <- lmer(perc_maxBAI~pet+ (1|Site), wg, REML=F)
modppt <- lmer(perc_maxBAI~ppt+ (1|Site), wg, REML=F)
modtmn <- lmer(perc_maxBAI~tmn+ (1|Site), wg, REML=F)
WildGrowth <- data.frame(AICc(modnull, modcwd, modaet, modpet, modppt, modtmn)
           , n=c(length(resid(modnull)), length(resid(modcwd)), length(resid(modaet)), length(resid(modpet)), length(resid(modppt)), length(resid(modtmn))))[order(AICc(modnull, modcwd, modaet, modpet, modppt, modtmn)[,2]),]

WildAOV$climDeltaAIC[which(WildAOV$Trait=="Growth")] <- 0
WildAOV$climname[which(WildAOV$Trait=="Growth")] <- "null"
WildAOV$bestclim[which(WildAOV$Trait=="Growth")] <- "none"
WildAOV$climvar[which(WildAOV$Trait=="Growth")] <- 0 #round(r.squaredGLMM(WildmodGrowth)[1],2)
WildAOV$climsig[which(WildAOV$Trait=="Growth")] <- 1 #round(summary(WildmodGrowth)$coefficients[2,4],4)



### write csvs with results
write.csv(HopAOV, paste0(results.dir,"/HopAOVresults.csv"))
write.csv(WildAOV, paste0(results.dir,"/WildAOVresults.csv"))



#______________________________________________________________________
######## * Variance Decomp #################
#______________________________________________________________________
# variance decomp of full measurement dataset (to capture within-tree variance as well)
# must use the individual trait datasets before they were averaged to .ind versions

# SLA, LDMC, leafsize: hub
# WD: WD
# Al_As: need to combine hub and stemk somehow
# kleaf: leafk.cl 
# kstem/Ks: stemk
# P50: p50wide -  no w/in tree
# growth: hoptrees or wg - no w/in tree



# note: this ends up being somewhat conservative for pop effects. 
#     - also, you get a lot of singular fits (boundary (singular) fit errors), particularly for the Garden. This is because the pop effect is often estimated to be 0

#++++ Garden ++++++++++

#**** Note: variance estimates are essentially identical without (1|rep) random effect. So I left it out for better fitting

# HopSLAvd <- lmer(SLA~ 1 + (1|site.pop) + (1|ind.id), data = hub, subset=site=="H")
HopSLAvd <- lmer(SLA~ 1  + (1|site.pop/ind.id), data = hub, subset=site=="H")
# HopSLAvd <- aov(SLA~site.pop/ind.id + rep, hub[which(hub$site=="H"),]) #huh, this overestimates the pop variance compared to lmer
# HopSLAvd <- aov(SLA~ rep+site.pop/ind.id , hub[which(hub$site=="H"),]) #whereas this is more similar to lmer
HopLDMCvd <- lmer(LDMC~ 1  + (1|site.pop/ind.id), data = hub, subset=site=="H")
Hopleafsizevd <- lmer(leafsize~ 1  + (1|site.pop/ind.id), data = hub, subset=site=="H")
Hopml_msvd <- lmer(ml_ms~ 1  + (1|site.pop/ind.id), data = hub, subset=site=="H")
Hoplogml_msvd <- lmer(log.ml_ms~ 1  + (1|site.pop/ind.id), data = hub, subset=site=="H")
  # resids look way better logged.
HopAl_Asvd <- lmer(Al_As~1 + (1|site.pop/ind.id), data = hub, subset=site=="H")
HoplogAl_Asvd <- lmer(log.Al_As~1 + (1|site.pop/ind.id), data = hub, subset=site=="H")
  # high resids with raw, low resids with logged.
HopWDvd <- lmer(WD~ 1  + (1|site.pop/ind.id), data = WD, subset=site=="H")
  #note: WD only typically has 1-2 reps w/in tree in H, but 1-3 in W
Hopkleafvd <- lmer(kleaf~ 1  + (1|site.pop/ind.id), data = leafk.cl, subset=site=="H")
Hopkstemvd <- lmer(kstem~ 1  + (1|site.pop/ind.id), data = stemk.cl, subset=site=="H")
HopKsvd <- lmer(Ks~ 1  + (1|site.pop/ind.id), data = stemk.cl, subset=site=="H")

## now use VarCorr extractor to pull out the variance and sd (cols 4:5 of this df)
HopSLAvar <- data.frame(VarCorr(HopSLAvd))
HopLDMCvar <- data.frame(VarCorr(HopLDMCvd))
Hopleafsizevar <- data.frame(VarCorr(Hopleafsizevd))
Hopml_msvar <- data.frame(VarCorr(Hopml_msvd))
Hoplogml_msvar <- data.frame(VarCorr(Hoplogml_msvd))
  # raw, w.in pop = .64 resid, log -> w/in pop is .75 resid
HopAl_Asvar <- data.frame(VarCorr(HopAl_Asvd))
HoplogAl_Asvar <- data.frame(VarCorr(HoplogAl_Asvd))
  # raw w.in pop = .65 resid, log -> w/in pop is .52 resid
Hopkleafvar <- data.frame(VarCorr(Hopkleafvd))
HopWDvar <- data.frame(VarCorr(HopWDvd))
Hopkstemvar <- data.frame(VarCorr(Hopkstemvd))
HopKsvar <- data.frame(VarCorr(HopKsvd))
  # I'm going to stick with raw for ease of interpretation.
  # Using the raw Al_As and ml_ms values slightly decreases the proportion of w/in pop (i.e. btw tree) variance
  # but doesn't meaningfully change the among population variance (which is ~0)



# make a data frame, rows are initially: between tree (within ppop, between pop, rep (remove/ignore), and within tree
variancesHop <- data.frame(HopSLAvar[,4],HopLDMCvar[,4],HopWDvar[,4],Hopml_msvar[,4],HopAl_Asvar[,4],Hopleafsizevar[,4],Hopkleafvar[,4],Hopkstemvar[,4],HopKsvar[,4], row.names = c("wiPop", "btwPop","wiTree"))

colnames(variancesHop) <- c("SLA","LDMC","WD","ml_ms","Al_As","leafsize","kleaf","kstem","Ks")
# and reorder variance componentsto make them nice plots
variancesHop <- variancesHop[c(2,1,3),] # reorder to go btw pop, w/in pop, w/in tree

scaledvariancesHop <- apply(variancesHop, MARGIN=2, FUN= function(x){x/sum(x)})


## without w.in tree for P50s and growth
HopP50stemvd <- lmer(P50stem~ 1  + (1|site.pop), data = Hop.ind)
HopP50leafvd <- lmer(P50leaf~ 1  + (1|site.pop), data = Hop.ind)
HopGrowthvd <- lmer(log.bio~ 1  + (1|pop), data = hoptrees)

HopP50stemvar <- data.frame(VarCorr(HopP50stemvd))
HopP50leafvar <- data.frame(VarCorr(HopP50leafvd))
HopGrowthvar <- data.frame(VarCorr(HopGrowthvd))

# make a data frame, rows are initially: between tree (within ppop, between pop, rep (remove/ignore), and within tree
variancesHop2 <- data.frame(c(HopP50stemvar[,4],0),c(HopP50leafvar[,4],0),c(HopGrowthvar[,4],0), row.names = c("btwPop","wiPop","wiTree"))

colnames(variancesHop2) <- c("P50stem","P50leaf","Growth")

scaledvariancesHop2 <- apply(variancesHop2, MARGIN=2, FUN= function(x){x/sum(x)})


#++++ Wild ++++++++++


WildSLAvd <- lmer(SLA~ 1  + (1|site.pop/ind.id), data = hub, subset=site=="W")
WildLDMCvd <- lmer(LDMC~ 1  + (1|site.pop/ind.id), data = hub, subset=site=="W")
Wildleafsizevd <- lmer(leafsize~ 1  + (1|site.pop/ind.id), data = hub, subset=site=="W") 
Wildml_msvd <- lmer(ml_ms~1 + (1|site.pop/ind.id), data = hub, subset=site=="W")
WildAl_Asvd <- lmer(Al_As~1 + (1|site.pop/ind.id), data = hub, subset=site=="W")
WildWDvd <- lmer(WD~ 1  + (1|site.pop/ind.id), data = WD, subset=site=="W")
#note: WD only typically has 1-2 reps w/in tree in H, but 1-3 in W
Wildkleafvd <- lmer(kleaf~ 1  + (1|site.pop/ind.id), data = leafk.cl, subset=site=="W")
Wildkstemvd <- lmer(kstem~ 1  + (1|site.pop/ind.id), data = stemk.cl, subset=site=="W")
WildKsvd <- lmer(Ks~ 1  + (1|site.pop/ind.id), data = stemk.cl, subset=site=="W")

## now use VarCorr extractor to pull out the variance and sd (cols 4:5 of this df)
WildSLAvar <- data.frame(VarCorr(WildSLAvd))
WildLDMCvar <- data.frame(VarCorr(WildLDMCvd))
Wildleafsizevar <- data.frame(VarCorr(Wildleafsizevd))
Wildml_msvar <- data.frame(VarCorr(Wildml_msvd))
WildAl_Asvar <- data.frame(VarCorr(WildAl_Asvd))
Wildkleafvar <- data.frame(VarCorr(Wildkleafvd))
WildWDvar <- data.frame(VarCorr(WildWDvd))
Wildkstemvar <- data.frame(VarCorr(Wildkstemvd))
WildKsvar <- data.frame(VarCorr(WildKsvd))




# make a data frame, rows are initially: between tree (within ppop, between pop, rep (remove/ignore), and within tree
variancesWild <- data.frame(WildSLAvar[,4],WildLDMCvar[,4],WildWDvar[,4],Wildml_msvar[,4],WildAl_Asvar[,4],Wildleafsizevar[,4],Wildkleafvar[,4],Wildkstemvar[,4],WildKsvar[,4], row.names = c("wiPop", "btwPop","wiTree"))

colnames(variancesWild) <- c("SLA","LDMC","WD","ml_ms","Al_As","leafsize","kleaf","kstem","Ks")
# and reorder variance componentsto make them nice plots
variancesWild <- variancesWild[c(2,1,3),]

scaledvariancesWild <- apply(variancesWild, MARGIN=2, FUN= function(x){x/sum(x)})



## without w.in tree for P50s and growth
WildP50stemvd <- lmer(P50stem~ 1  + (1|site.pop), data = Wild.ind)
WildP50leafvd <- lmer(P50leaf~ 1  + (1|site.pop), data = Wild.ind)
WildGrowthvd <- lmer(perc_maxBAI~1 + (1|Site), wg)

WildP50stemvar <- data.frame(VarCorr(WildP50stemvd))
WildP50leafvar <- data.frame(VarCorr(WildP50leafvd))
WildGrowthvar <- data.frame(VarCorr(WildGrowthvd))

# make a data frame, rows are initially: between tree (within ppop, between pop, rep (remove/ignore), and within tree
variancesWild2 <- data.frame(c(WildP50stemvar[,4],0),c(WildP50leafvar[,4],0),c(WildGrowthvar[,4],0), row.names = c("btwPop","wiPop","wiTree"))

colnames(variancesWild2) <- c("P50stem","P50leaf","Growth")

scaledvariancesWild2 <- apply(variancesWild2, MARGIN=2, FUN= function(x){x/sum(x)})


## combined variances of traits with within/tree variation and those without 
variancesHop.comb <- cbind(variancesHop, variancesHop2)
variancesWild.comb <- cbind(variancesWild, variancesWild2)

scaledvariancesHop.comb <- cbind(scaledvariancesHop, scaledvariancesHop2)
scaledvariancesWild.comb <- cbind(scaledvariancesWild, scaledvariancesWild2)

### make dataframe of among-pop variance component only for Fig 3 (raw)
comparison <- data.frame( "Garden" = t(variancesHop.comb)[,1], "Wild" = t(variancesWild.comb)[,1])
comparison$perc_wild <- comparison$Garden/comparison$Wild
comparison$perc_wild[which(comparison$perc_wild == Inf)] <- NA








#______________________________________________________________________
######## ** FIG 1: Map + Growth #################
#______________________________________________________________________
mypallight <- paste0(brewer.pal(n=8, "Dark2"), "33")


# turn sampling locations into spatial data for plotting

latlon <- sp::CRS("+proj=longlat +datum=WGS84")
sampledlocs <- sp::SpatialPoints(coords = data.frame(popclim$Lon.dd, popclim$Lat.dd), proj4string = latlon)

plotlocs <- popclim[which(popclim$Garden_pop %in% c(1,3,15,16,18,22,26)),] %>% select(Garden_pop, Code, Lon.dd, Lat.dd)
plotlocs$label.Lon <- plotlocs$Lon.dd
plotlocs$label.Lon[5] <- -120.9
plotlocs$label.Lon[6] <- -119.9

plotlocs$label.Lat <- plotlocs$Lat.dd+.24
plotlocs$label.Lat[4] <- plotlocs$Lat.dd[4]-0.24
plotlocs$label.Lat[5] <- plotlocs$Lat.dd[5]+.15
oakPicture<-readPNG(source="Data/BlueOakPicture.png")

#jpeg(filename="figures/FIG1_SamplingMap_v1.jpg",width =5, height=5, units = "in", res = 600)
quartz(width=3.8,height=4.5)
#par(mar=c(0,0,0,0))
CA <- maps::map('state', region = c('california'), col="black") 
t <- maps::map.scale(x=-123.5, y=33.9,ratio=F,srt=45)
maps::map.axes()
mtext("Latitude (°)",side=2, line=2.3)
mtext("Longitude (°)",side=1, line=2.3)
points(latitude~longitude, qudo[which(qudo$latitude>33.5),], pch=16, col="lightgrey")#mypallight[1])
legend('topright',legend = c("herbarium sample", 'common garden','source pops'), pch=c(16,17,16), col = c("lightgrey","black","black"),pt.bg=mypal[1], pt.cex = c(1,1.5,1), bty="n", cex=.8)
#points(gardenlocs.latlon)
#points(sampledlocs[which(popclim$Garden_pop %in% c(1,3,15,16,18,22,26))], pch=16, cex=.8, col="black")#mypal[4])
points(Lat.dd~Lon.dd, plotlocs, pch=16, cex=1, col=factor(Code))
text(x=plotlocs$label.Lon
     ,y=plotlocs$label.Lat
     ,labels = popclim$Code[which(popclim$Garden_pop %in% c(1,3,15,16,18,22,26))]
     ,cex=.7)
# text(garden.sampledlocs.latlon, labels=gardens.sampled$Garden_pop)
points(sampledlocs[which(popclim$Code =="HREC")], pch=17,col="black", bg=mypal[1], cex=1.6)
rasterImage(oakPicture,xleft = -117.5, ybottom = 37.5, xright=-118+3.75, ytop=37+3.75*(540/720), interpolate = T)
mtext("a)", side=3, line=0.3, adj=0)
# add a north arrow
#narrow <- sp::layout.north.arrow(type=1)
# shift coordinates to be bottom left corner
#Narrow <- maptools::elide(Narrow, shift=c(extent()))
# #text(gardenlocs.latlon, labels=gardens$Code)
# points(tosample.locs.latlon, pch=3, col="black")
## -- add in water potential locations
#points(sampledlocs[which(popclim$Code %in% wp.site$Site)], pch=15, cex=1, col=4)
#legend('topright',legend = c("herbarium sample", 'common garden','source pops', 'water potentials'), pch=c(16,24,15, 3), col = c(mypallight[1],"black",mypal[4], "black"),pt.bg=mypal[1], pt.cex = c(1,1.5,1), bty="n")
#legend('topright',legend = c("herbarium sample", 'sampled sites'), pch=c(16,15), col = c(mypallight[1],4),pt.bg=mypal[1], pt.cex = c(1,1), bty="n")
#dev.off()

if(save.figures==T){
  quartz.save(file=paste0(results.dir,"/Fig1_Map_v3.pdf"), type="pdf")
}


### Panel b) Garden Height:
quartz(width=3.2,height=4.5)
par(mfrow=c(2,1), mgp=c(2,.6,0), mar=c(3,5,1,1))

# fit a gam and make predictions
tmp <- hoptrees %>% filter(height>0) %>% arrange(pet.td)
gardengrowth.best <- gamm(height~s(pet.td),random = list(pop=~1), data = tmp)
test <- predict(gardengrowth.best, newdata = data.frame("pet.td"=seq(from=min(tmp$pet.td),to=max(tmp$pet.td), by=1)),se.fit = T)

# plot predictions ober observations
plot(height~pet.td, hoptrees, pch=16, col="#00000022", ylab="Garden Growth\n(height, m)", xlab="PET transfer distance (mm)")
# all individual heights
points(Height~pet.td, popw, col=factor(pop.name), pch=16, cex=1.2)
# pop average heights for 7 focal pops
lines(test[[1]]~seq(from=min(tmp$pet.td),to=max(tmp$pet.td), by=1), col=mypal[1], lwd=2)
# gam prediction
abline(v=0, lty=2)
# location of the common garden
mtext("b)", side=3, line=0.2, adj=0)
mtext(paste0("R^2=",round(summary(gardengrowth.best$gam)$r.sq,2)), side=3, line=-1.2)

plot(perc_maxBAI~pet, wg, pch=16, col="#00000022", ylab="Wild Growth\n(% max BAI)", xlab="PET (mm)")
points(perc_maxBAI~pet, wg.pop, pch=16, col=factor(Site), cex=1.2)
mtext("c)", side=3, line=0.2, adj=0)

if(save.figures==T){
  quartz.save(file=paste0(results.dir,"/Fig1_Map_Growthbc_v1.pdf"), type="pdf")
}






#______________________________________________________________________
######## ** FIG 2: Var Decomp #################
#______________________________________________________________________

## set color choices
#pal.vardecomp <- brewer.pal(n=9, "Set1")
pal.vardecomp <- c(paste0(brewer.pal(n=3, "Set1")[1],"66"), brewer.pal(n=9, "Blues")[c(8,5,2)])
palette(pal.vardecomp)
# set the colors that wil denote within-species, within-genus, within family and across CWMs
colchoices <- c(1,2,4,3,6)


quartz(width=5.3,height=5.4)
#jpeg(file=paste0("./","/Fig_VarainceDecomp_scaled.jpg"), width=4.3, height=5.4, units="in", res=600)
par(mfrow=c(2,1), mar=c(1,3.6,1,3.6), mgp=c(2.3,1,0), oma=c(3.6,0,3,0), cex.lab=1.4, cex.axis=1.1)
p<-barplot(height=as.matrix(scaledvariancesHop)
           , beside=F, names.arg=rep(NA, times=9)
           , col = pal.vardecomp
           , legend.text= c("Btw Pops", "Btw Trees", "Within Tree")
           , args.legend=list(bty="n", x=11.5, y=1.5, xpd=NA, cex=1,xjust=0.5, ncol=4)
           , ylab="% Var in Garden", las=3
           , xlim=c(1.1,18.2), width=1, space=1)#,  yaxt="n" # ylim=c(0,0.008),
axis(4,at=c(0,.4,.8), labels = c(0,.4,.6)*round(max(HopAOV$CV, na.rm=T),1), xpd=NA,col="#333333", col.axis="#333333")
mtext("Trait CV", side = 4, line =2.3)
barplot(height=as.matrix(t(HopAOV$CV[-c(10,11,12)]/max(HopAOV$CV[-c(10,11,12)], na.rm=T))), add=T
        , space=c(4,3,3,3,3,3,3,3,3), width=.5, names.arg=rep(NA, times=9))
#text(x=p,y=colSums(scaledvariancesWD)+.001, labels=climsens, srt=90, xpd=T, font=c(2,2,1,1,2,2,2,1))
mtext(text = "a)", side=3, adj=-0.12, line=.2)

p<-barplot(height=as.matrix(scaledvariancesWild)
           , beside=F, names.arg=rep(NA, times=9)
           , col = pal.vardecomp
           , legend.text= NULL
           , args.legend=list(bty="n", x=mean(c(1.1,16)), y=1.5, xpd=NA, cex=1,xjust=0.5, ncol=4)
           , ylab="% Var in Wild", las=3
           , xlim=c(1.1,18.2), width=1, space=1)#,  yaxt="n" # ylim=c(0,0.008),
axis(4,at=c(0,.4,.8), labels = c(0,.4,.6)*round(max(WildAOV$CV, na.rm=T),1), xpd=NA,col="#333333", col.axis="#333333")
mtext("Trait CV", side = 4, line =2.3)
barplot(height=as.matrix(t(WildAOV$CV[-c(10,11,12)]/max(WildAOV$CV[-c(10,11,12)], na.rm=T))), add=T
        , space=c(4,3,3,3,3,3,3,3,3), width=.5, names.arg=rep(NA, times=9))
#text(x=p,y=colSums(scaledvariancesWD)+.001, labels=climsens, srt=90, xpd=T, font=c(2,2,1,1,2,2,2,1))
mtext(text = "b)", side=3, adj=-0.12, line=.2)

axis(1,at=p,labels= c("SLA","LDMC","WD","Ml:Ms","Al:As","Leaf Size","k[leaf]","k[stem]","K[s]"),font=1, cex.axis = 1.1, tick=F,las=2, mgp=c(2,.2,0))
#dev.off()
palette(mypal)
if(save.figures==T){
  quartz.save(file=paste0(results.dir,"/Fig1_VarainceDecomp_scaled_v5.pdf"),type = "pdf" )
}

#_________________
####### **... Fig 1 continued: second bar plot for 3 traits without within-canopy ######
#_________________

colchoices <- c(1,2,4,3,6)
quartz(width=3,height=5.4)
#jpeg(file=paste0("./","/Fig_VarainceDecomp_scaled.jpg"), width=4.3, height=5.4, units="in", res=600)
par(mfrow=c(2,1), mar=c(1,3.6,1,1), mgp=c(2.3,1,0), oma=c(3.6,0,3,0), cex.lab=1.4, cex.axis=1.1)
p<-barplot(height=as.matrix(scaledvariancesHop2)
           , beside=F, names.arg=rep(NA, times=3)
           , col = pal.vardecomp
           #, legend.text= c("Btw Pops", "Btw Trees", "Within Tree")
           #, args.legend=list(bty="n", x=11.5, y=1.5, xpd=NA, cex=1,xjust=0.5, ncol=4)
           , ylab="% Var in Garden", las=3
           , xlim=c(1.1,6.6), width=1, space=1)#,  yaxt="n" # ylim=c(0,0.008),
#axis(4,at=c(0,.4,.8), labels = c(0,.4,.6)*round(max(HopAOV$CV, na.rm=T),1), xpd=NA,col="#333333", col.axis="#333333")
#mtext("Trait CV", side = 4, line =2.3)
#barplot(height=as.matrix(t(HopAOV$CV[-c(10,11)]/max(HopAOV$CV, na.rm=T))), add=T
#        , space=c(4,3,3,3,3,3,3,3,3), width=.5, names.arg=rep(NA, times=9))
#text(x=p,y=colSums(scaledvariancesWD)+.001, labels=climsens, srt=90, xpd=T, font=c(2,2,1,1,2,2,2,1))
#mtext(text = "a)", side=3, adj=-0.12, line=.2)

p<-barplot(height=as.matrix(scaledvariancesWild2)
           , beside=F, names.arg=rep(NA, times=3)
           , col = pal.vardecomp
           , legend.text= NULL
           , args.legend=list(bty="n", x=mean(c(1.1,16)), y=1.5, xpd=NA, cex=1,xjust=0.5, ncol=4)
           , ylab="% Var in Wild", las=3
           , xlim=c(1.1,6.6), width=1, space=1)#,  yaxt="n" # ylim=c(0,0.008),
#axis(4,at=c(0,.4,.8), labels = c(0,.4,.6)*round(max(WildAOV$CV, na.rm=T),1), xpd=NA,col="#333333", col.axis="#333333")
#mtext("Trait CV", side = 4, line =2.3)
#barplot(height=as.matrix(t(WildAOV$CV[-c(10,11)]/max(WildAOV$CV, na.rm=T))), add=T
#        , space=c(4,3,3,3,3,3,3,3,3), width=.5, names.arg=rep(NA, times=9))
#text(x=p,y=colSums(scaledvariancesWD)+.001, labels=climsens, srt=90, xpd=T, font=c(2,2,1,1,2,2,2,1))
#mtext(text = "b)", side=3, adj=-0.12, line=.2)

axis(1,at=p,labels= c("P50stem","P50leaf","Growth"),font=1, cex.axis = 1.1, tick=F,las=2, mgp=c(2,.2,0))
#dev.off()
palette(mypal)
if(save.figures==T){
  quartz.save(file=paste0(results.dir,"/Fig1_VarainceDecompp50Growth_scaled_v1.pdf"),type = "pdf" )
}

# 
# #_________________
# ###### Fig 1 V2: all bars together ##
# #_________________
# 
# quartz(width=6.3,height=5.4) 
# #jpeg(file=paste0("./","/Fig_VarainceDecomp_scaled.jpg"), width=4.3, height=5.4, units="in", res=600)
# par(mfrow=c(2,1), mar=c(1,3.6,1,3.6), mgp=c(2.3,1,0), oma=c(3.6,0,3,0), cex.lab=1.4, cex.axis=1.1)
# p<-barplot(height=as.matrix(scaledvariancesHop.comb)
#            , beside=F, names.arg=rep(NA, times=12)
#            , col = pal.vardecomp
#            , legend.text= c("Btw Pops", "Btw Trees", "Within Tree")
#            , args.legend=list(bty="n", x=11.5, y=1.5, xpd=NA, cex=1,xjust=0.5, ncol=4)
#            , ylab="% Var in Garden", las=3
#            , xlim=c(1.1,24.2), width=1, space=1)#,  yaxt="n" # ylim=c(0,0.008),
# axis(4,at=c(0,.4,.8), labels = c(0,.4,.6)*round(max(HopAOV$CV, na.rm=T),1), xpd=NA,col="#333333", col.axis="#333333")
# mtext("Trait CV", side = 4, line =2.3)
# barplot(height=as.matrix(t(HopAOV$CV/max(HopAOV$CV, na.rm=T))), add=T
#         , space=c(4,3,3,3,3,3,3,3,3,3,3,3), width=.5, names.arg=rep(NA, times=12))
# #text(x=p,y=colSums(scaledvariancesWD)+.001, labels=climsens, srt=90, xpd=T, font=c(2,2,1,1,2,2,2,1))
# mtext(text = "a)", side=3, adj=-0.12, line=.2)
# 
# p<-barplot(height=as.matrix(scaledvariancesWild.comb)
#            , beside=F, names.arg=rep(NA, times=12)
#            , col = pal.vardecomp
#            , legend.text= NULL
#            , args.legend=list(bty="n", x=mean(c(1.1,16)), y=1.5, xpd=NA, cex=1,xjust=0.5, ncol=4)
#            , ylab="% Var in Wild", las=3
#            , xlim=c(1.1,24.2), width=1, space=1)#,  yaxt="n" # ylim=c(0,0.008),
# axis(4,at=c(0,.4,.8), labels = c(0,.4,.6)*round(max(WildAOV$CV, na.rm=T),1), xpd=NA,col="#333333", col.axis="#333333")
# mtext("Trait CV", side = 4, line =2.3)
# barplot(height=as.matrix(t(WildAOV$CV/max(WildAOV$CV, na.rm=T))), add=T
#         , space=c(4,3,3,3,3,3,3,3,3,3,3,3), width=.5, names.arg=rep(NA, times=12))
# #text(x=p,y=colSums(scaledvariancesWD)+.001, labels=climsens, srt=90, xpd=T, font=c(2,2,1,1,2,2,2,1))
# mtext(text = "b)", side=3, adj=-0.12, line=.2)
# 
# axis(1,at=p,labels= c("SLA","LDMC","WD","Ml:Ms","Al:As","Leaf Size","k[leaf]","k[stem]","K[s]","P50[stem]","P50[leaf]","Growth")
#      ,font=1, cex.axis = 1.1, tick=F,las=2, mgp=c(2,.2,0))
# #dev.off()
# palette(mypal)
# if(save.figures==T){
#   quartz.save(file=paste0(results.dir,"/Fig2_VarainceDecomp_scaled_combined_v1.pdf"),type = "pdf" )
# }




#________________________________________________________________________
############## ** FIG 3: G vs G+E trait variation ###############
#_______________________________________________________________________


# Make a column for plotting significant traits as filled and ns as open
HopAOV$Plotting.Sig <- 1
HopAOV$Plotting.Sig[which(HopAOV$sig<= 0.05)]<-16

WildAOV$Plotting.Sig <- 1
WildAOV$Plotting.Sig[which(WildAOV$sig<= 0.05)]<-16
WildAOV$Plotting.Clim.Sig <- 1
WildAOV$Plotting.Clim.Sig[which(WildAOV$climsig<= 0.05)]<-16
palette(mypal)

quartz(width=7, height=3.5)
layout(matrix(c(1,2), nrow=1),widths = c(1.3,1) )

par(mar=c(3.5,3.5,2,6), mgp=c(2.5,1,0))
plot(WildAOV$eta2*100~I(HopAOV$eta2*100), col=factor(WildAOV$Trait), pch=WildAOV$Plotting.Sig, cex=1, lwd=1.8,
     xlim=c(0,.7)*100, ylim=c(0,.7)*100,xaxs="i", yaxs="i",
     xlab="Garden: Btw-Pop Variation (%)", ylab="Wild: Btw-Site Variation (%)")
abline(a=0,b=1)
# plot the significant garden pop differences with a black circle
points(WildAOV$eta2[which(HopAOV$sig<0.05)]*100~I(HopAOV$eta2[which(HopAOV$sig<0.05)]*100), pch=3, cex=1.5, lwd=1.8, col=factor(HopAOV$Trait)[which(HopAOV$sig<0.1)])
legend.names <- c("SLA"
                  ,"LDMC"
                  ,"WD"
                  , "Ml:Ms"
                  , "Al:As"
                  , "Leaf Size"
                  , "k[leaf]"
                  , "k[stem]"
                  , "Ks"
                  , "P50[stem]"
                  , "P50[leaf]"
                  , "Growth")

legend("bottomright", legend=c("Garden sig", "Wild sig", "ns","","",""), pch=c(3,16,1,1,1),col=c("black","black","black",rep("white",3)), bty="n", xpd=F, cex=.8)

legend(x=75,y=70, legend=legend.names, col=mypal[c(11,6,12,8,1,7,3,5,4,10,9,2)]
       , pch=15, xpd=NA,horiz = F, bty='n')
mtext(text = "a)", side=3, adj=-0.12, line=.4)

textlocs.x <- HopAOV$eta2*100
textlocs.y <- WildAOV$eta2*100
textlocs.y[which(textlocs.y<51.5)] <- textlocs.y[which(textlocs.y<51.5)]- 2.7
textlocs.y[which(textlocs.y>51.5)] <- textlocs.y[which(textlocs.y>51.5)] + 2.7
textlocs.x[2] <- 17 # shift LDMC left
textlocs.x[5] <- 28 # shift Al:As right
textlocs.x[7] <- 12 # shift kleaf right
text(y=textlocs.y, x=textlocs.x, legend.names, cex=.6)

par(mar=c(3.5,3.5,2,1), mgp=c(2.5,1,0))
plot(I(climvar*100)~I(eta2*100), data=WildAOV[which(WildAOV$sig<= 0.05 | WildAOV$climsig<=0.05),], pch=Plotting.Clim.Sig,cex=1,col=factor(WildAOV$Trait)[which(WildAOV$sig<= 0.05 | WildAOV$climsig<=0.05)],
     xlim=c(0,.7)*100, ylim=c(0,.7)*100,xaxs="i",yaxs="i",
     lwd=1.8,
     xlab="Btw-Site or Btw-Pop Variation (%)", ylab="Variation Explained by Climate")
points(I(climvar*100)~I(eta2*100), data=HopAOV[which(HopAOV$climsig<=0.05),], pch=3,cex=1.3,col=factor(HopAOV$Trait)[which(HopAOV$climsig<=0.05)], lwd=1.8)
abline(a=0,b=1)
mtext(text = "b)", side=3, adj=-0.12, line=.4)
# trait labels for Wild
textlocs.x <- WildAOV$eta2[which(WildAOV$sig<= 0.05 | WildAOV$climsig<=0.05)]*100
textlocs.y <- WildAOV$climvar[which(WildAOV$sig<= 0.05 | WildAOV$climsig<=0.05)]*100 + 2.7
textlocs.x[8] <- 51+4  # shift kstem right
textlocs.x[5] <- 55+5 # shift Al:As right
text(y=textlocs.y, x=textlocs.x, legend.names[which(WildAOV$sig<= 0.05 | WildAOV$climsig<=0.05)], cex=.6)

# trait labels for garden
with(HopAOV[which(HopAOV$climsig<=0.05),], text(y=climvar*100-2.5, x= eta2*100-5, cex=.6, font=2, c("SLA","WD","Leaf Size","Growth") ))


#legend('topleft', legend=c("Garden","Wild"), pch=c(3,17), bty="n")
if(save.figures==T){
  quartz.save(file=paste0(results.dir,"/Fig3_G-vs-GE_variation_v6.pdf"), type="pdf")
}









#________________________________________________________________________
########## ** FIG 4: Bootstrapping variance comparison #######################
#________________________________________________________________________

var.boot <- function(trait, boots=1000){
  set.seed(42)
  variances <- data.frame(Hop = rep(NA, boots), Wild = rep(NA, boots))
  for(i in 1:boots){
    # only bootstrap to the smallest dataset size to control for different sample sizes
    tosamp <- min(length(which(!is.na(Hop.ind[,trait]))), length(which(!is.na(Wild.ind[,trait]))))
    variances$Hop[i] <- var(Hop.ind[which(!is.na(Hop.ind[,trait])),trait][sample(1:length(which(!is.na(Hop.ind[,trait]))), tosamp, TRUE)])
    variances$Wild[i] <- var(Wild.ind[which(!is.na(Wild.ind[,trait])),trait][sample(1:length(which(!is.na(Wild.ind[,trait]))), tosamp, TRUE)])
  }
  varsumsHop <- quantile(variances$Hop, c(.05,.25,.5,.75,.95))
  names(varsumsHop) <- c("Hop.05", "Hop.25", "Hop.50", "Hop.75", "Hop.95")
  varsumsWild <- quantile(variances$Wild, c(.05,.25,.5,.75,.95))
  names(varsumsWild) <- c("Wild.05", "Wild.25", "Wild.50", "Wild.75", "Wild.95")
  varratio <- quantile(variances$Wild/variances$Hop, c(.05, .1,.5,.9,.95))
  names(varratio) <- c("ratio.05", "ratio.10", "ratio.50", "ratio.90", "ratio.95") # wild/garden
  varperc <- quantile(variances$Hop/variances$Wild, c(.05, .1,.5,.9,.95)) # garden/wild
  names(varperc) <- c("perc.05", "perc.10", "perc.50", "perc.90", "perc.95")
  varsums <- c(varsumsHop, varsumsWild,varratio, varperc)
  return(varsums)
}

vars <- c("mSLA","mLDMC","mWD","mml_ms","mAl_As","mleafsize","mkleaf", "mkstem", "mKs", "P50stem","P50leaf")



varcomps <- var.boot(vars[1])
for(i in 2:length(vars)){
  varcomps <- rbind(varcomps, var.boot(vars[i]))
}
varcomps <- data.frame(varcomps)
rownames(varcomps) <- vars
#varcomps$sig



## Showing the % of Wild variance in the Garden
quartz(width=6, height=3)
par(mar=c(5,4,1,1))
tmp <- varcomps$Hop.50/varcomps$Wild.50
statest <- varcomps$Wild.05/varcomps$Hop.50
bp <- barplot(varcomps$perc.50, las=2
              , names.arg = c("SLA","LDMC","WD","Ml:Ms","Al:As","Leaf size","k[leaf]","k[stem]", "Ks","P50[stem]","P50[leaf]")
              , ylab="Garden Var : Wild Var", ylim=c(0,1.9)
              , border = c(rep("white",3), rep("black", 8)), # "white",rep("black",3))
              , col = c(rep("lightgray",3), rep("white", 8) )) #"lightgray",rep("white",3))) #with marginal significance
#arrows(x0=bp, y0=varcomps$ratio.10, y1=varcomps$ratio.90, lwd=4, length = 0)
arrows(x0=bp, y0=varcomps$perc.05, y1=varcomps$perc.95, lwd=4, length = 0, col="darkgrey")
#arrows(x0=bp, y0=varcomps$ratio.10, y1=varcomps$ratio.90, lwd=4, length = 0, col="darkgray")
legend("topleft",legend=c("Total Var","btw Pop Var"), pch=c(NA, 16), fill = c("white",NA), merge=T, border = c("black","white"), bty="n")

abline(h=1, lty=2, lwd=2)

points(comparison$perc_wild[1:11]~bp, cex=2, pch=16)

# text(c("btw Pop\nVar", "total Var"), x=c(1.9, 7.9), y=c(.6, 2.2))
# arrows(x0=c(2.2,7.9), x1=c(2.95,9), y0=c(0.6,2.08), y1=c(0.22,1.7 ), length=.1)

if(save.figures==T){
  quartz.save(file=paste0(results.dir,"/Fig4_GvWvariance_comparison_v1.pdf"),type = "pdf")
}




#_________________________________________________________________________________
################# ** FIG 5: Trait Integration correlation plots #################
#_________________________________________________________________________________

#### . Calculate Correlations########

# define the variables we want to include

# V8 moved Ml:Ms back down
# V9, removed growth
varsH <- c("mleafsize","mAl_As","mkleaf", "mkstem", "mKs","mml_ms","mSLA","mLDMC","mWD", "P50stem","P50leaf")#, "Height")
varsW <- c("mleafsize","mAl_As","mkleaf", "mkstem", "mKs","mml_ms","mSLA","mLDMC","mWD", "P50stem","P50leaf")#, "perc_maxBAI")



# calculate trait-trait correlations for H
tmpH <- data.frame(all.ind[which(all.ind$site=="H"),]) %>% select(all_of(varsH))
# remove outlier stem P50 in H
tmpH$P50stem[which(tmpH$P50stem> -2)] <- NA
corH <- round(cor(tmpH, use="pairwise.complete.obs"),2)
pmatH <- cor_pmat(tmpH) #matrix of p-values

# calculate trait-trait cors in W
tmpW <- data.frame(all.ind[which(all.ind$site=="W"),]) %>% select(all_of(varsW))
corW <- round(cor(tmpW, use="pairwise.complete.obs"),2)
pmatW <- cor_pmat(tmpW)

# make a master correlation table for plotting heatmap
corAll <- corH 
pmatAll <- pmatH
# make upper tri for corW
corAll[upper.tri(corAll)] <- corW[upper.tri(corW)]
pmatAll[upper.tri(pmatAll)] <- pmatW[upper.tri(pmatW)]


corAll.m <-  melt(corAll)




########## Make two panesl of heat map (edited outside of R)
quartz(width=4, height=4)
#corrplot::corrplot(corAll.named, p.mat=1/pmatAll*.01, sig.level=0.1, insig="pch",pch=16, pch.cex=.8
#                   , type="lower", diag=T, tl.pos="lt", tl.col="black", cl.pos="r", addgrid.col = NA)#, title="Wild = Upper, Garden = Lower")
corAll.named <- corAll
colnames(corAll.named) <- c("Leaf size","Al:As","k[leaf]","k[stem]", "Ks","Ml:Ms","SLA","LDMC","WD", "P50[stem]","P50[leaf]")
rownames(corAll.named) <- c("Leaf size","Al:As","k[leaf]","k[stem]", "Ks","Ml:Ms","SLA","LDMC","WD", "P50[stem]","P50[leaf]")


corrplot::corrplot(corAll.named, p.mat=1/pmatAll*.01, sig.level=0.1, insig="pch",pch=16, pch.cex=.8
                   , type="lower", diag=F, add=F, cl.pos="n", method="ellipse", outline=F, addgrid.col = NA, tl.col="black"
)
if(save.figures==T){
  quartz.save(file=paste0(results.dir,"/Fig5_Correlation_Heatmap_Hopland_v9.", type="pdf"))
}



quartz(width=4, height=4)
#corrplot::corrplot(corAll.named, p.mat=1/pmatAll*.01, sig.level=0.1, insig="pch",pch=16, pch.cex=.8
#                   , type="lower", diag=T, tl.pos="lt", tl.col="black", cl.pos="r", addgrid.col = NA)#, title="Wild = Upper, Garden = Lower")
corW.named <- corW
colnames(corW.named) <- c("Leaf size","Al:As","k[leaf]","k[stem]", "Ks","Ml:Ms","SLA","LDMC","WD", "P50[stem]","P50[leaf]")
rownames(corW.named) <- c("Leaf size","Al:As","k[leaf]","k[stem]", "Ks","Ml:Ms","SLA","LDMC","WD", "P50[stem]","P50[leaf]")

corrplot::corrplot(corW.named, p.mat=1/pmatW*.01, sig.level=0.1, insig="pch",pch=16, pch.cex=.8
                   , type="lower", add=F, diag=F, cl.pos="n", tl.col="black"
                   , method="ellipse", outline = F, addgrid.col = NA)
if(save.figures==T){
  quartz.save(file=paste0(results.dir,"/Fig5_Correlation_Heatmap_Wild_v9.", type="pdf"))
}





#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
################ BEGIN: SUPPLEMENTAL FIGURES #######################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#______________________________________________________________________
######## ** Table S1: Site Info #################
#______________________________________________________________________

siteinfo <- popclim[which(popclim$Garden_pop %in% c(1,3,15,16,18,22,26,109)),] %>% select(Code,elevation,tmn1951_1980, tmx1951_1980, ppt1951_1980,pet1951_1980,cwd1951_1980)
siteinfo.round <- apply(siteinfo[,-1], MARGIN=2, FUN=round, digits=1)
siteinfo.clean <- data.frame("Site"=siteinfo$Code, siteinfo.round)
siteinfo.clean$elevation[which(siteinfo.clean$Site=="HREC")] <- 270 # pulled elevation for the common garden from Google Earth
write.csv(siteinfo.clean, paste0(results.dir,"/TableS1_SiteCharacteristics.csv"))



#_________________________________________________________________________________
################# ** FIG S1: BAI standardization to size #################
#_________________________________________________________________________________

quartz(width=4, height=4)
plot(BAI_cm2~DBH_cm, wildgrowth,type="n",xlab="DBH (cm)", ylab="Basal Area Increment (cm2)")
points(BAI_cm2~DBH_cm, wildgrowth, col=factor(Site), pch=16)
points(BAI_cm2~DBH_cm, wildgrowth[which(wildgrowth$Site %in% c("SJR", "LYN","SMR","SMT","SON", "SRD")),], pch=1, cex=1.1) # highlight the trees actually from this study
# median
abline(rq(BAI_cm2~DBH_cm, wildgrowth,tau=.9),col="blue")
# # other range of quantiles
# taus <- c(.05,.1,.25,0.5,.75,.90,.95)
# 
# for( i in 1:length(taus)){
#   abline(rq(BAI_cm2~DBH_cm, wildgrowth,tau=taus[i]),col="gray")
# }
legend('topleft', legend=c("6 study sites", "5 additional sites"), pch=c(21,16), col=c("black","grey"), pt.bg = "gray")
if(save.figures==T){
  quartz.save(file=paste0(results.dir,"/FigS1_BAI-v-DBH_v1.pdf"), type="pdf")
}




#_________________________________________________________________________________
################ ** FIG S2: Wild vs Garden Growth ################################
#_________________________________________________________________________________

# demonstrates that performance (and presumably any trait signals that are linked to performance) in the garden 
# are not likely to be driven by maternal effects. 
# plot(Height~perc_maxBAI, all.pop.wide, pch=16, cex=1.2, ylab='Height (m)', xlab="% max BAI\n(size-standardized)")
# 
# plot(Diam~DBH, all.pop.wide, pch=16, cex=1.2, ylab='Height (m)', xlab="% max BAI\n(size-standardized)")


growthcomparison <- popw %>% select(pop, pop.name, Bio, log.Bio, Height, Diam, DBH, BAI, perc_maxBAI)
quartz(width=3.5, height=3.5)
par(mar=c(4,5,1,1))
plot(Height~perc_maxBAI, growthcomparison,cex=2, col=factor(pop.name), pch=16
     , ylab="Garden Height (m)"
     , xlab="Wild % max BAI\n(size standardized growth)")
legend("bottomright", legend=levels(factor(growthcomparison$pop.name))[-7], col=mypal[1:6], pch=16, cex=.8)
if(save.figures==T){
  quartz.save(file=paste0(results.dir,"/FigS2_Garden-v-WildGrowth_v1.pdf"), type="pdf")
}






#_________________________________________________________________________________
################# ** FIG S3&S4: Wild vs Garden pop trait means #################
#_________________________________________________________________________________

WildAOV$climsig.plotting <- 0 # don't plot ns results
WildAOV$climsig.plotting[which(WildAOV$climsig<0.1 & WildAOV$climsig>0.05)] <- 2 # plot marg as dotted lines
WildAOV$climsig.plotting[which(WildAOV$climsig<0.05)] <- 1 # plot sig as solid lines

HopAOV$climsig.plotting <- 0 # don't plot ns results
HopAOV$climsig.plotting[which(HopAOV$climsig<0.1 & HopAOV$climsig>0.05)] <- 2 # plot marg as dotted lines
HopAOV$climsig.plotting[which(HopAOV$climsig<0.05)] <- 1 # plot sig as solid lines


quartz(width=8, height=4)
par(mfcol=c(2,6), mar=c(3,3,1.5,1), mgp=c(2,.8,0))

#SLA
# vars <- c("mSLA","mLDMC","mWD","mml_ms","mAl_As","mleafsize","mkleaf", "mkstem", "mKs", "P50stem","P50leaf")
# for (i in vars){}
plot(mSLA~get(HopAOV$climname[which(HopAOV$Trait=="SLA")]), Hop.ind, ylab="SLA (cm2g-1)", xlab=HopAOV$bestclim[which(HopAOV$Trait=="SLA")], col=factor(pop.name))
abline(HopmodSLA, lty=HopAOV$climsig.plotting[which(HopAOV$Trait=="SLA")])
mtext("a)", adj=0)

plot(mSLA~get(WildAOV$climname[which(WildAOV$Trait=="SLA")]), Wild.ind, ylab="SLA (cm2g-1)", xlab=WildAOV$bestclim[which(WildAOV$Trait=="SLA")], col=factor(pop.name))
abline(lm(mSLA~aet.2018ds, Wild.ind), lty=WildAOV$climsig.plotting[which(WildAOV$Trait=="SLA")])
mtext("g)", adj=0)


#LDMC
#plot(mLDMC~get(HopAOV$climname[which(HopAOV$Trait=="LDMC")]), Hop.ind, ylab="LDMC (g g-1)", xlab=HopAOV$bestclim[which(HopAOV$Trait=="LDMC")])
plot(mLDMC~ppt, Hop.ind,ylab="LDMC (g g-1)", xlab="PPT[30yr]", col="grey" )
#abline(HopmodLDMC, lty=HopAOV$climsig.plotting[which(HopAOV$Trait=="LDMC")])
mtext("b)", adj=0)

#plot(mLDMC~get(WildAOV$climname[which(WildAOV$Trait=="LDMC")]), Wild.ind, ylab="LDMC (g g-1)", xlab=WildAOV$bestclim[which(WildAOV$Trait=="LDMC")])
plot(mLDMC~ppt, Wild.ind,ylab="LDMC (g g-1)", xlab="PPT[30yr]", col="grey" )
abline(lm(mLDMC~aet.2018ds, Wild.ind), lty=WildAOV$climsig.plotting[which(WildAOV$Trait=="LDMC")])
mtext("h)", adj=0)


#WD
plot(mWD~get(HopAOV$climname[which(HopAOV$Trait=="WD")]), Hop.ind, ylab="WD (g cm-2)", xlab=HopAOV$bestclim[which(HopAOV$Trait=="WD")], col=factor(pop.name))
abline(HopmodWD, lty=HopAOV$climsig.plotting[which(HopAOV$Trait=="WD")])
mtext("c)", adj=0)

plot(mWD~get(WildAOV$climname[which(WildAOV$Trait=="WD")]), Wild.ind, ylab="WD (g cm-2)", xlab=WildAOV$bestclim[which(WildAOV$Trait=="WD")], col=factor(pop.name))
abline(lm(mWD~aet.diff, Wild.ind), lty=WildAOV$climsig.plotting[which(WildAOV$Trait=="WD")])
mtext("i)", adj=0)


#ml_ms
#plot(mml_ms~get(HopAOV$climname[which(HopAOV$Trait=="ml_ms")]), Hop.ind, ylab="Ml:Ms (g g-1)", xlab=HopAOV$bestclim[which(HopAOV$Trait=="ml_ms")])
plot(mml_ms~ppt, Hop.ind,ylab="Ml:Ms (g g-1)", xlab="PPT[30yr]", col="grey" )
#abline(Hopmodml_ms, lty=HopAOV$climsig.plotting[which(HopAOV$Trait=="ml_ms")])
mtext("d)", adj=0)

plot(mml_ms~get(WildAOV$climname[which(WildAOV$Trait=="ml_ms")]), Wild.ind, ylab="Ml:Ms (g g-1)", xlab=WildAOV$bestclim[which(WildAOV$Trait=="ml_ms")], col=factor(pop.name))
abline(Wildmodml_ms, lty=WildAOV$climsig.plotting[which(WildAOV$Trait=="ml_ms")])
mtext("j)", adj=0)


#Al_As
#plot(mAl_As~get(HopAOV$climname[which(HopAOV$Trait=="Al_As")]), Hop.ind, ylab="Ml:Ms (g g-1)", xlab=HopAOV$bestclim[which(HopAOV$Trait=="Al_As")])
plot(mAl_As~ppt, Hop.ind,ylab="Al:As (cm2 mm-2)", xlab="PPT[30yr]", col="grey" )
#abline(HopmodAl_As, lty=HopAOV$climsig.plotting[which(HopAOV$Trait=="Al_As")])
mtext("e)", adj=0)

plot(mAl_As~get(WildAOV$climname[which(WildAOV$Trait=="Al_As")]), Wild.ind, ylab="Al:As (cm2 mm-2)", xlab=WildAOV$bestclim[which(WildAOV$Trait=="Al_As")], col=factor(pop.name))
abline(lm(mAl_As~get(WildAOV$climname[which(WildAOV$Trait=="Al_As")]), Wild.ind), lty=WildAOV$climsig.plotting[which(WildAOV$Trait=="Al_As")])
mtext("k)", adj=0)


# leafsize
plot(mleafsize~get(HopAOV$climname[which(HopAOV$Trait=="leafsize")]), Hop.ind, ylab="leafsize (cm2)", xlab=HopAOV$bestclim[which(HopAOV$Trait=="leafsize")], col=factor(pop.name))
abline(Hopmodleafsize, lty=HopAOV$climsig.plotting[which(HopAOV$Trait=="leafsize")])
mtext("f)", adj=0)

plot(mleafsize~get(WildAOV$climname[which(WildAOV$Trait=="leafsize")]), Wild.ind, ylab="leafsize (cm2)", xlab=WildAOV$bestclim[which(WildAOV$Trait=="leafsize")], col=factor(pop.name))
abline(lm(mleafsize~aet.2018ds, Wild.ind), lty=WildAOV$climsig.plotting[which(WildAOV$Trait=="leafsize")])
mtext("l)", adj=0)


if(save.figures==T){
  quartz.save(file=paste0(results.dir,"/FigS3_TraitClim1_v2.pdf"),type = "pdf" )
}





quartz(width=8, height=4)
par(mfcol=c(2,6), mar=c(3,3,1.5,1), mgp=c(2,.8,0))

#kleaf
plot(mkleaf~get(HopAOV$climname[which(HopAOV$Trait=="kleaf")]), Hop.ind, ylab="kleaf", xlab=HopAOV$bestclim[which(HopAOV$Trait=="kleaf")], col=factor(pop.name))
abline(Hopmodkleaf, lty=HopAOV$climsig.plotting[which(HopAOV$Trait=="kleaf")])
mtext("a)", adj=0)

plot(mkleaf~get(WildAOV$climname[which(WildAOV$Trait=="kleaf")]), Wild.ind, ylab="kleaf", xlab=WildAOV$bestclim[which(WildAOV$Trait=="kleaf")], col=factor(pop.name))
abline(lm(mkleaf~get(WildAOV$climname[which(WildAOV$Trait=="kleaf")]), Wild.ind), lty=WildAOV$climsig.plotting[which(WildAOV$Trait=="kleaf")])
mtext("g)", adj=0)


#kstem
plot(mkstem~get(HopAOV$climname[which(HopAOV$Trait=="kstem")]), Hop.ind, ylab="kstem", xlab=HopAOV$bestclim[which(HopAOV$Trait=="kstem")], col=factor(pop.name))
abline(Hopmodkstem, lty=HopAOV$climsig.plotting[which(HopAOV$Trait=="kstem")])
mtext("b)", adj=0)

plot(mkstem~get(WildAOV$climname[which(WildAOV$Trait=="kstem")]), Wild.ind, ylab="kstem", xlab=WildAOV$bestclim[which(WildAOV$Trait=="kstem")], col=factor(pop.name))
abline(lm(mkstem~get(WildAOV$climname[which(WildAOV$Trait=="kstem")]), Wild.ind), lty=WildAOV$climsig.plotting[which(WildAOV$Trait=="kstem")])
mtext("h)", adj=0)


#Ks
# plot(mKs~get(HopAOV$climname[which(HopAOV$Trait=="Ks")]), Hop.ind, ylab="Ks", xlab=HopAOV$bestclim[which(HopAOV$Trait=="Ks")])
plot(mKs~ppt, Hop.ind,ylab="Ks", xlab="PPT[30yr]", col="grey" )
#abline(HopmodKs, lty=HopAOV$climsig.plotting[which(HopAOV$Trait=="Ks")])
mtext("c)", adj=0)

# plot(mKs~get(WildAOV$climname[which(WildAOV$Trait=="Ks")]), Wild.ind, ylab="Ks (g cm-2)", xlab=WildAOV$bestclim[which(WildAOV$Trait=="Ks")])
plot(mKs~ppt, Wild.ind,ylab="Ks", xlab="PPT[30yr]", col="grey" )
# abline(lm(mKs~aet.diff, Wild.ind), lty=WildAOV$climsig.plotting[which(WildAOV$Trait=="Ks")])
mtext("i)", adj=0)


#P50stem
#plot(P50stem~get(HopAOV$climname[which(HopAOV$Trait=="P50stem")]), Hop.ind, ylab="Ml:Ms (g g-1)", xlab=HopAOV$bestclim[which(HopAOV$Trait=="P50stem")])
plot(P50stem~ppt, Hop.ind,ylab="P50stem", xlab="PPT[30yr]", col="grey" )
#abline(HopmodP50stem, lty=HopAOV$climsig.plotting[which(HopAOV$Trait=="P50stem")])
mtext("d)", adj=0)

plot(P50stem~get(WildAOV$climname[which(WildAOV$Trait=="P50stem")]), Wild.ind, ylab="P50stem", xlab=WildAOV$bestclim[which(WildAOV$Trait=="P50stem")], col=factor(pop.name))
abline(WildmodP50stem, lty=WildAOV$climsig.plotting[which(WildAOV$Trait=="P50stem")])
mtext("j)", adj=0)


#P50leaf
#plot(P50leaf~get(HopAOV$climname[which(HopAOV$Trait=="P50leaf")]), Hop.ind, ylab="Ml:Ms (g g-1)", xlab=HopAOV$bestclim[which(HopAOV$Trait=="P50leaf")])
plot(P50leaf~ppt, Hop.ind,ylab="P50leaf", xlab="PPT[30yr]", col="grey" )
#abline(HopmodP50leaf, lty=HopAOV$climsig.plotting[which(HopAOV$Trait=="P50leaf")])
mtext("e)", adj=0)

#plot(P50leaf~get(WildAOV$climname[which(WildAOV$Trait=="P50leaf")]), Wild.ind, ylab="Al:As (cm2 mm-2)", xlab=WildAOV$bestclim[which(WildAOV$Trait=="P50leaf")])
plot(P50leaf~ppt, Wild.ind,ylab="P50leaf", xlab="PPT[30yr]", col="grey" )
# abline(lm(P50leaf~get(WildAOV$climname[which(WildAOV$Trait=="P50leaf")]), Wild.ind), lty=WildAOV$climsig.plotting[which(WildAOV$Trait=="P50leaf")])
mtext("k)", adj=0)


# Growth
plot(log.bio~get(HopAOV$climname[which(HopAOV$Trait=="Growth")]), hoptrees, ylab="Growth (log10(stem biomass))", xlab=HopAOV$bestclim[which(HopAOV$Trait=="Growth")], col=factor(pop.name))
abline(lm(log.bio~get(HopAOV$climname[which(HopAOV$Trait=="Growth")]), hoptrees), lty=HopAOV$climsig.plotting[which(HopAOV$Trait=="Growth")])
mtext("f)", adj=0)

plot(perc_maxBAI~ppt, wg, ylab="Growth (% max BAI)", xlab="PPT[30yr", col="grey")
#abline(lm(mleafsize~aet.2018ds, Wild.ind), lty=WildAOV$climsig.plotting[which(WildAOV$Trait=="leafsize")])
mtext("l)", adj=0)


if(save.figures==T){
  quartz.save(file=paste0(results.dir,"/FigS4_TraitClim2_v2.pdf"),type = "pdf" )
}





#_________________________________________________________________________________
################# ** FIG S5: Norm of Reaction Plots #################
#_________________________________________________________________________________


all.pop.m$site.numeric <- as.numeric(as.factor(all.pop.m$site))
# add in a column for PPT30yr with garden's climate for plotting norms of reaction
all.pop.m$ppt.actual <- all.pop.m$ppt
all.pop.m$ppt.actual[which(all.pop.m$site=="H")] <- popclim$ppt1951_1980[which(popclim$Code=="HOP")]
all.pop.m$pet.actual <- all.pop.m$pet
all.pop.m$pet.actual[which(all.pop.m$site=="H")] <- popclim$pet1951_1980[which(popclim$Code=="HOP")]


## Plotting just 'Garden' vs 'Wild' (not true norm of reaction)

traits <- c("SLA","LDMC","WD","ml_ms","Al_As","leafsize","kleaf","kstem","Ks")#,"P50stem")
labs <- c("SLA","LDMC","WD","Ml:Ms","Al:As","Leaf size","k[leaf]","k[stem]","Ks")#,"P50")


quartz(width=6.5, height=6.5)
par(mfrow=c(3,3), mar=c(3.2,4,0,0), oma=c(0,2,1,1), mgp=c(2,1,0), cex.lab=1.5 )

for (j in 1:length(traits)){
  tr <- traits[j]
  dataz <- all.pop.m[which(all.pop.m$variable==tr),]
  dataz$pop.name <- factor(dataz$pop.name)
  plot(value~site.numeric, data=dataz, type="p"
       , xlim=c(0.75,2.25)
       , ylim=c(min(dataz$value,na.rm=T), max(dataz$value, na.rm=T) + (max(dataz$value, na.rm=T)-min(dataz$value, na.rm=T))*0.1)
       # make space for significance stars in y-axis
       , pch=16, col=factor(pop.name), cex=1.5, xaxt="n", xlab="", ylab=labs[j])
  axis(1, at=c(1,2), labels=c("Garden","Wild"))
  #mtext(labs[j], side=3,line = 2)
  for(i in levels(dataz$pop.name)){
    lines(value~site.numeric, data=dataz[which(dataz$pop.name==i),], col=pop.name, lwd=1.3)
  }
  # add in Garden Significance
  if(tr %in% c("SLA","leafsize")){
    text("*", x=1,y=max(dataz$value, na.rm=T) + (max(dataz$value, na.rm=T)-min(dataz$value, na.rm=T))*0.075, cex=3)
  }
  # add in wild significance (all traits except Ks)
  if(tr != "Ks"){
    text("*", x=2,y=max(dataz$value, na.rm=T) + (max(dataz$value, na.rm=T)-min(dataz$value, na.rm=T))*0.075, cex=3)
  }
  
  
}

if(save.figures==T){
  quartz.save(file=paste0(results.dir,"/FigS5_Norm_of_Reaction_v1.pdf"), type="pdf")
}


## plotting true norm of reaction based on pop PPT

quartz(width=6.5, height=6.5)
par(mfrow=c(3,3), mar=c(3.2,4,0,0), oma=c(2,2,1,1), mgp=c(2,1,0), cex.lab=1.5 )
traits <- c("SLA","LDMC","WD","ml_ms","Al_As","leafsize","kleaf","kstem","Ks")#,"P50stem")
labs <- c("SLA","LDMC","WD","Ml:Ms","Al:As","Leaf size","k[leaf]","k[stem]","Ks")#,"P50")


for (j in 1:length(traits)){
  tr <- traits[j]
  dataz <- all.pop.m[which(all.pop.m$variable==tr),]
  dataz$pop.name <- factor(dataz$pop.name)
  plot(value~ppt.actual, data=dataz, type="p"
       , ylim=c(min(dataz$value,na.rm=T), max(dataz$value, na.rm=T) + (max(dataz$value, na.rm=T)-min(dataz$value, na.rm=T))*0.1)
       # make space for significance stars in y-axis
       , pch=(site.numeric-1)*15+1, col=factor(pop.name), cex=1.8, xaxt="n", xlab="", ylab=labs[j])
  axis(1)
  #mtext(labs[j], side=3,line = 2)
  for(i in levels(dataz$pop.name)){
    lines(value~ppt.actual, data=dataz[which(dataz$pop.name==i),], col=pop.name, lwd=1.3)
  }
  # add in Garden Significance
  if(tr %in% c("SLA","leafsize")){
    text("*", x=1,y=max(dataz$value, na.rm=T) + (max(dataz$value, na.rm=T)-min(dataz$value, na.rm=T))*0.075, cex=3)
  }
  # add in wild significance (all traits except Ks)
  if(tr != "Ks"){
    text("*", x=2,y=max(dataz$value, na.rm=T) + (max(dataz$value, na.rm=T)-min(dataz$value, na.rm=T))*0.075, cex=3)
  }
  if(tr =="ml_ms"){
    legend("topleft", legend=levels(as.factor(all.pop.m$pop.name)), pch=16, col=1:7    )
  }
  if(tr=="kleaf"){legend("topleft", legend = c("W","G"), pch=c(16,1))}
  if(tr=="kstem") {mtext("PPT[30yr] (mm)", side=1, line=3, cex=1)}
  
}


if(save.figures==T){
  quartz.save(file=paste0(results.dir,"/FigS5_TrueNorm_of_Reaction_ppt_v1.pdf"), type="pdf")
}


## plotting true norm of reaction based on pop PPT

quartz(width=6.5, height=6.5)
par(mfrow=c(3,3), mar=c(3.2,4,0,0), oma=c(2,2,1,1), mgp=c(2,1,0), cex.lab=1.5 )
traits <- c("SLA","LDMC","WD","ml_ms","Al_As","leafsize","kleaf","kstem","Ks")#,"P50stem")
labs <- c("SLA","LDMC","WD","Ml:Ms","Al:As","Leaf size","k[leaf]","k[stem]","Ks")#,"P50")


for (j in 1:length(traits)){
  tr <- traits[j]
  dataz <- all.pop.m[which(all.pop.m$variable==tr),]
  dataz$pop.name <- factor(dataz$pop.name)
  plot(value~pet.actual, data=dataz, type="p"
       , ylim=c(min(dataz$value,na.rm=T), max(dataz$value, na.rm=T) + (max(dataz$value, na.rm=T)-min(dataz$value, na.rm=T))*0.1)
       # make space for significance stars in y-axis
       , pch=(site.numeric-1)*15+1, col=factor(pop.name), cex=1.8, xaxt="n", xlab="", ylab=labs[j])
  axis(1)
  #mtext(labs[j], side=3,line = 2)
  for(i in levels(dataz$pop.name)){
    lines(value~pet.actual, data=dataz[which(dataz$pop.name==i),], col=pop.name, lwd=1.3)
  }
  # add in Garden Significance
  if(tr %in% c("SLA","leafsize")){
    text("*", x=1,y=max(dataz$value, na.rm=T) + (max(dataz$value, na.rm=T)-min(dataz$value, na.rm=T))*0.075, cex=3)
  }
  # add in wild significance (all traits except Ks)
  if(tr != "Ks"){
    text("*", x=2,y=max(dataz$value, na.rm=T) + (max(dataz$value, na.rm=T)-min(dataz$value, na.rm=T))*0.075, cex=3)
  }
  if(tr =="ml_ms"){
    legend("topright", legend=levels(as.factor(all.pop.m$pop.name)), pch=16, col=1:7    )
  }
  if(tr=="kleaf"){legend("topright", legend = c("W","G"), pch=c(16,1))}
  if(tr=="kstem") {mtext("PET[30yr] (mm)", side=1, line=3, cex=1)}
}



if(save.figures==T){
  quartz.save(file=paste0(results.dir,"/FigS5_TrueNorm_of_Reaction_pet_v1.pdf"), type="pdf")
}


#_________________________________________________________________________________
################# ** FIG S6: Wild vs Garden pop trait means #################
#_________________________________________________________________________________


quartz(width=6, height=6)
par(mfrow=c(3,3), mar=c(3.5,3.5,0,0), oma=c(0,0,1,1), mgp=c(2.1,1,0), cex.lab=1.2)
traits <- c("SLA","LDMC","WD","ml_ms","Al_As","leafsize","kleaf","kstem","Ks")#,"P50stem")
labs <- c("SLA","LDMC","WD","Ml:Ms","Al:As","Leaf size","k[leaf]","k[stem]","Ks")#,"P50")


for(i in 1:9){
  # subset data down to make it easier to grep colnames
  tmp <- popw %>% select(grep("Height", colnames(popw)), grep(traits[i], colnames(popw)))
  # make plots appropriately large for error bars
  xlimit <- c(min(tmp[,paste0("H_",traits[i])],na.rm=T) - 1.2* max(tmp[,paste0("H_se",traits[i])], na.rm=T), 
              max(tmp[,paste0("H_",traits[i])],na.rm=T) + 1.2* max(tmp[,paste0("H_se",traits[i])], na.rm=T) )
  ylimit <- c(min(tmp[,paste0("W_",traits[i])], na.rm=T) - 1.2* max(tmp[,paste0("W_se",traits[i])], na.rm=T), 
              max(tmp[,paste0("W_",traits[i])], na.rm=T) + 1.2* max(tmp[,paste0("W_se",traits[i])], na.rm=T) )
  # plot with xlimit & ylimit to focus data in center
  sqlimit <- c(min(xlimit[1], ylimit[1]), max(xlimit[2], ylimit[2]))
  # plot with sqlimit to have 1:1 line center
  
  plot(get(paste0("W_", traits[i]))~get(paste0("H_", traits[i])), popw
       , xlim=sqlimit, ylim=sqlimit, pch=16, cex=2
       , xlab=paste0("Garden ",labs[i]), ylab=paste0("Wild ", labs[i]))  
  # add in y error bars (height)
  arrows(x0=tmp[,paste0("H_", traits[i])], x1=tmp[,paste0("H_", traits[i])],
         y0=tmp[,paste0("W_", traits[i])], y1=tmp[,paste0("W_", traits[i])]+tmp[,paste0("W_se", traits[i])], length =0, lwd=2, col="grey")
  arrows(x0=tmp[,paste0("H_", traits[i])], x1=tmp[,paste0("H_", traits[i])],
         y0=tmp[,paste0("W_", traits[i])], y1=tmp[,paste0("W_", traits[i])]-tmp[,paste0("W_se", traits[i])], length =0, lwd=2, col="grey")
  # add in x error bars (trait)
  arrows(x=tmp[,paste0("H_", traits[i])], x1=tmp[,paste0("H_", traits[i])] + tmp[,paste0("H_se", traits[i])],
         y0=tmp[,paste0("W_", traits[i])], lwd=2, length=0, col="grey")
  arrows(x=tmp[,paste0("H_", traits[i])], x1=tmp[,paste0("H_", traits[i])] - tmp[,paste0("H_se", traits[i])],
         y0=tmp[,paste0("W_", traits[i])], lwd=2, length=0, col="grey")
  abline(a=0,b=1, lwd=2, col="grey")
  #mod <- lm(get(paste0("W_", traits[i]))~get(paste0("H_", traits[i])), tmp, weights=1/sqrt(tmp[,paste0("H_se", traits[i])]))
  points(get(paste0("W_", traits[i]))~get(paste0("H_", traits[i])), popw
         , xlim=sqlimit, ylim=sqlimit, pch=16, cex=2, col=factor(pop.name))
  mod <- lmodel2(get(paste0("W_", traits[i]))~get(paste0("H_", traits[i])), tmp)
  if(mod$P.param<0.1){
    abline(a=mod$regression.results$Intercept[3], b=mod$regression.results$Slope[3])
    #abline(a=mod$regression.results$Intercept[2], b=mod$regression.results$Slope[2])
    mtext(paste("p=", round(mod$P.param,3)),side = 1, line=-1, adj=.1)
  }
  
}

if(save.figures==T){
  quartz.save(file=paste0(results.dir,"/FigS6_Garden-v-Wild_TraitValues_v5.pdf"), type="pdf")
}











#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
########### ** FIG S7: Supplemental Leaf size vs number on Al_As ** #####################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

quartz(width=6, height=3.2)
#psize <- ggplot(hub[-which(hub$site=="SF"),], aes(x=leafsize, y=Al_As, col=site)) + geom_point() + geom_smooth(method='lm')
#pnum <- ggplot(hub[-which(hub$site=="SF"),], aes(x=log(leaf_number,10), y=log(Al_As,10), col=site)) + geom_point() + geom_smooth(method='lm')
#psize
#pnum
palette(cbp2)
par(mfrow=c(1,2), mar=c(4,4,2,1), mgp=c(2.5,1,0))

# leaf size vs Al:As
plot(Al_As~leafsize, hub[-which(hub$site=="SF"),], col=factor(site), pch=16
     , ylab=expression(paste(A[L]:A[S]," (",cm^2*mm^-2,")"))
     , xlab=expression(paste("average Leaf Size (",cm^2,")")))
mW <- lmodel2(Al_As~leafsize, hub[which(hub$site=="W"),], nperm = 1000)
mH <- lmodel2(Al_As~leafsize, hub[which(hub$site=="H"),], nperm = 1000)
# add SMA regression
abline(a = mW$regression.results$Intercept[3], b=mW$regression.results$Slope[3], lwd=2, col=2)
abline(a = mH$regression.results$Intercept[3], b=mH$regression.results$Slope[3], lwd=2, col=1)
# add R2 
mtext(paste("R2=",round(mH$rsquare,2)),side=3,adj=0, col=1)
mtext(paste("R2=",round(mW$rsquare,2)),side=3,adj=1, col=2)

# leaf number vs Al:As
plot(Al_As~leaf_number, hub[-which(hub$site=="SF"),], col=factor(site), pch=16, xlim=c(0,30)
     , ylab=expression(paste(A[L]:A[S]," (",cm^2*mm^-2,")"))
     , xlab="Leaves per branch")
mW <- lmodel2(Al_As~leaf_number, hub[which(hub$site=="W"),], nperm = 1000)
mH <- lmodel2(Al_As~leaf_number, hub[which(hub$site=="H"),], nperm = 1000)
# add SMA regression
abline(a = mW$regression.results$Intercept[3], b=mW$regression.results$Slope[3], lwd=2, col=2)
abline(a = mH$regression.results$Intercept[3], b=mH$regression.results$Slope[3], lwd=2, col=1, lty=2)
legend("topright",legend = c("Garden","Wild"), col=c(1,2), pch=16)
mtext(paste("R2=",round(mH$rsquare,2)),side=3,adj=0, col=1)
mtext(paste("R2=",round(mW$rsquare,2)),side=3,adj=1, col=2)


if(save.figures==T){
  quartz.save(file=paste0(results.dir,"/FigS7_LeafSizevsLeafNum-Al_As_v2.pdf"), type="pdf")
}




#______________________________________________________________________________
########### ** FIG S8: Bootstrapping Trait-Trait correlations ################
#______________________________________________________________________________



# Doing 10000 bootstrap samples from Wild vs Garden trait-trait correlations
W.H.test <- function(trait1, trait2, boots=1000, sig.type="alpha"){
  corrW <- corrH <- rep(0, boots)
  set.seed(42) # set seed to make answer repeatable
  
  for(i in 1:boots){
    # only bootstrap to the smallest dataset size to control for different sample sizes
    tosamp <- min(length(which(complete.cases(Hop.ind[,c(trait1,trait2)])==T)), length(which(complete.cases(Wild.ind[,c(trait1,trait2)])==T)))
    corrW[i] <- cor(Wild.ind[which(complete.cases(Wild.ind[,c(trait1,trait2)])==T),c(trait1,trait2)][sample(1:length(which(complete.cases(Wild.ind[,c(trait1,trait2)])==T)), tosamp, TRUE),])[1,2]
    corrH[i] <- cor(Hop.ind[which(complete.cases(Hop.ind[,c(trait1,trait2)])==T),c(trait1,trait2)][sample(1:length(which(complete.cases(Hop.ind[,c(trait1,trait2)])==T)), tosamp, TRUE),])[1,2]
  }
  diffs <- corrH - corrW
  # store the significance, ns = alpha >0.2, with alpha < 0.2, alpha <0.1, and alpha <0.05 as possibilities
  sig <-0.3
  if(quantile(diffs, 0.1)>0 | quantile(diffs,0.9)<0) sig <- 0.2
  if(quantile(diffs, 0.05)>0 | quantile(diffs,0.95)<0) sig <- 0.1
  if(quantile(diffs, 0.025)>0 | quantile(diffs,0.975)<0) sig <- 0.05
  # store a reverse scaled significance [0,1] for corrplotting with alpha <0.5 =>1
  revsig <-0
  if(quantile(diffs, 0.1)>0 | quantile(diffs,0.9)<0) revsig <- 0.33
  if(quantile(diffs, 0.05)>0 | quantile(diffs,0.95)<0) revsig <- 0.66
  if(quantile(diffs, 0.025)>0 | quantile(diffs,0.975)<0) revsig <- 1
  if(sig.type=="alpha") return(sig)
  if(sig.type=="rev.scaled") return(revsig)
}

traitz <- varsW
hwcomp <- matrix(NA, nrow=length(traitz), ncol=length(traitz))
colnames(hwcomp) <- traitz
rownames(hwcomp) <- traitz
hwcomprevscale <- hwcomp
for(i in 1:length(traitz)){
  if(i < length(traitz)){
    for(j in (i+1):length(traitz)){
      hwcomp[j,i] <- W.H.test(traitz[i], traitz[j], boots=5000)
      # create a version scaled so that alpha <0.05 = 1 and >0.2 = 0
      hwcomprevscale[j,i] <- W.H.test(traitz[i], traitz[j], boots=5000, sig.type = "rev.scaled")
    }
  }  
}



#colnames(hwcomprevscale) <- c("SLA","LDMC","WD","Ml:Ms","Al:As","Leaf size","k[leaf]","k[stem]", "Ks","P50[stem]","P50[leaf]")
#rownames(hwcomprevscale)  <-c("SLA","LDMC","WD","Ml_Ms","Al:As","Leaf size","k[leaf]","k[stem]", "Ks","P50[stem]","P50[leaf]")
# with reordered rows (v3)
# colnames(hwcomprevscale) <- c("Leaf size","Al:As","Ml:Ms", "k[leaf]","k[stem]", "Ks","SLA","LDMC","WD", "P50[stem]","P50[leaf]")
# rownames(hwcomprevscale) <- c("Leaf size","Al:As","Ml:Ms","k[leaf]","k[stem]", "Ks","SLA","LDMC","WD", "P50[stem]","P50[leaf]")
# moving Ml:As down (v4)
colnames(hwcomprevscale) <- c("Leaf size","Al:As","k[leaf]","k[stem]", "Ks","Ml:Ms","SLA","LDMC","WD", "P50[stem]","P50[leaf]")
rownames(hwcomprevscale) <- c("Leaf size","Al:As","k[leaf]","k[stem]", "Ks","Ml:Ms","SLA","LDMC","WD", "P50[stem]","P50[leaf]")

corcolors <- corrplot::COL2(n = 7) # pull the colors for making a legend

# # get rid of the alpha<0.2 correlation changes (v5)
# hwcomprevscale_trm <- hwcomprevscale
# hwcomprevscale_trm[which(hwcomprevscale_trm<0.6)] <- 0

quartz(width=4, height=4)
corrplot::corrplot(hwcomprevscale, p.mat=NULL, sig.level=0.04, insig="pch",pch=16, pch.cex=.8
                   , type="lower", add=F, diag=F, cl.pos="n", tl.col="black"
                   , method="color", outline = T, addgrid.col = F)
# v4 w 0.2 alpha cutoff
legend("topright", bty="n", legend=c("alpha <0.2","alpha <0.1", "alpha <0.05"), pch=15, pt.cex = 2, col=corcolors[5:7], xpd=NA)
# v5 w/out .2, hwcomprevscaled_trm
# legend("topright", bty="n", legend=c("alpha <0.1", "alpha <0.05"), pch=15, pt.cex = 2, col=corcolors[6:7], xpd=NA)


if(save.figures==T){
  quartz.save(file=paste0(results.dir,"/FigS8_Correlation_SigFig_v5.pdf"), type="pdf")
}







#_________________________________________________________________________________
################# ** FIG S9&S10: Population Average Trait-Growth relationships #################
#_________________________________________________________________________________

### Loop through each trait and plot with error bars

quartz(width=6.5, height=6.5)
par(mfrow=c(3,3), mar=c(3.2,2,0,0), oma=c(0,2,1,1), mgp=c(2,1,0), cex.lab=1.5 )
traits <- c("SLA","LDMC","WD","ml_ms","Al_As","leafsize","kleaf","kstem","Ks","P50stem")
labs <- c("SLA","LDMC","WD","Ml:Ms","Al:As","Leaf size","k[leaf]","k[stem]","Ks","P50")

#traits <- c("P50leaf","P50stem","kstem","Ks","kleaf","kleaf.sa.temp","Al_As","WD","SLA")
for(i in 1:9){
  # subset data down to make it easier to grep colnames
  tmp <- popw %>% select(grep("Height", colnames(popw)), grep(traits[i], colnames(popw)))
  # make plots appropriately large for error bars
  xlimit <- c(min(tmp[,paste0("H_",traits[i])]) - 1.2* max(tmp[,paste0("H_se",traits[i])], na.rm=T), 
              max(tmp[,paste0("H_",traits[i])]) + 1.2* max(tmp[,paste0("H_se",traits[i])], na.rm=T) )
  ylimit <- c(min(tmp$Height)-1.2*max(tmp$seHeight), max(tmp$Height) + 1.2 * max(tmp$seHeight))
  plot(Height~get(paste0("H_", traits[i])), popw
       , xlim=xlimit, ylim=ylimit, pch=16, cex=2
       , xlab=labs[i], ylab="")
  if(i==4){mtext("30 yr height (m)",side=2, line=2)}
  # add in y error bars (height)
  arrows(x0=tmp[,paste0("H_", traits[i])], x1=tmp[,paste0("H_", traits[i])],
         y0=tmp$Height, y1=tmp$Height+tmp$seHeight, length =0, lwd=2, col="grey")
  arrows(x0=tmp[,paste0("H_", traits[i])], x1=tmp[,paste0("H_", traits[i])],
         y0=tmp$Height, y1=tmp$Height-tmp$seHeight, length =0, lwd=2, col="grey")
  # add in x error bars (trait)
  arrows(x=tmp[,paste0("H_", traits[i])], x1=tmp[,paste0("H_", traits[i])] + tmp[,paste0("H_se", traits[i])],
         y0=tmp$Height, lwd=2, length=0, col="grey")
  arrows(x=tmp[,paste0("H_", traits[i])], x1=tmp[,paste0("H_", traits[i])] - tmp[,paste0("H_se", traits[i])],
         y0=tmp$Height, lwd=2, length=0, col="grey")
  points(Height~get(paste0("H_", traits[i])), popw, pch=16, cex=2, col=factor(pop.name))
  mod <- lm(Height~get(paste0("H_", traits[i])), tmp, weights=1/sqrt(tmp[,paste0("H_se", traits[i])]))
  if(summary(mod)$coefficients[2,4]<0.1){
    abline(mod, lty = ifelse(summary(mod)$coefficients[2,4]<0.05,1,2))
    mtext(paste("p=", round(summary(mod)$coefficients[2,4],3)),side = 1, line=-1, adj=.9)
  }
  
}

if(save.figures==T){
  quartz.save(file=paste0(results.dir,"/FigS9_Trait-Growth_GardenHeight_v1.pdf"),type = "pdf")
}


### same figure with volume rather than height

quartz(width=6.5, height=6.5)
par(mfrow=c(3,3), mar=c(3.2,2,0,0), oma=c(0,2,1,1), mgp=c(2,1,0), cex.lab=1.5 )
traits <- c("SLA","LDMC","WD","ml_ms","Al_As","leafsize","kleaf","kstem","Ks")#,"P50stem")
labs <- c("SLA","LDMC","WD","Ml:Ms","Al:As","Leaf size","k[leaf]","k[stem]","Ks")#,"P50")

#traits <- c("P50leaf","P50stem","kstem","Ks","kleaf","kleaf.sa.temp","Al_As","WD","SLA")
for(i in 1:9){
  # subset data down to make it easier to grep colnames
  tmp <- popw %>% select(grep("log.Bio", colnames(popw)), grep(traits[i], colnames(popw)))
  # make plots appropriately large for error bars
  xlimit <- c(min(tmp[,paste0("H_",traits[i])]) - 1.2* max(tmp[,paste0("H_se",traits[i])], na.rm=T), 
              max(tmp[,paste0("H_",traits[i])]) + 1.2* max(tmp[,paste0("H_se",traits[i])], na.rm=T) )
  ylimit <- c(min(tmp$log.Bio)-1.2*max(tmp$selog.Bio), max(tmp$log.Bio) + 1.2 * max(tmp$selog.Bio))
  plot(log.Bio~get(paste0("H_", traits[i])), popw
       , xlim=xlimit, ylim=ylimit, pch=16, cex=2
       , xlab=labs[i], ylab="")
  if(i==4){mtext(expression(paste(log[10],"(30 yr stem volume)")),side=2, line=2)}
  # add in y error bars (log.Bio)
  arrows(x0=tmp[,paste0("H_", traits[i])], x1=tmp[,paste0("H_", traits[i])],
         y0=tmp$log.Bio, y1=tmp$log.Bio+tmp$selog.Bio, length =0, lwd=2, col="grey")
  arrows(x0=tmp[,paste0("H_", traits[i])], x1=tmp[,paste0("H_", traits[i])],
         y0=tmp$log.Bio, y1=tmp$log.Bio-tmp$selog.Bio, length =0, lwd=2, col="grey")
  # add in x error bars (trait)
  arrows(x=tmp[,paste0("H_", traits[i])], x1=tmp[,paste0("H_", traits[i])] + tmp[,paste0("H_se", traits[i])],
         y0=tmp$log.Bio, lwd=2, length=0, col="grey")
  arrows(x=tmp[,paste0("H_", traits[i])], x1=tmp[,paste0("H_", traits[i])] - tmp[,paste0("H_se", traits[i])],
         y0=tmp$log.Bio, lwd=2, length=0, col="grey")
  points(log.Bio~get(paste0("H_", traits[i])), popw, pch=16, cex=2, col=factor(pop.name))
  mod <- lm(log.Bio~get(paste0("H_", traits[i])), tmp, weights=1/sqrt(tmp[,paste0("H_se", traits[i])]))
  if(summary(mod)$coefficients[2,4]<0.1){
    abline(mod, lty = ifelse(summary(mod)$coefficients[2,4]<0.05,1,2))
    mtext(paste("p=", round(summary(mod)$coefficients[2,4],3)),side = 1, line=-1, adj=.9)
  }
  
}

if(save.figures==T){
  quartz.save(file=paste0(results.dir,"/FigS9_Trait-Growth_GardenVolume_v1.pdf"),type = "pdf")
}




### Wild version



#_______________________________________________________________________
#################### ** FIG S11: Supplemental Safety-Efficiency Tradeoff ##############
#_______________________________________________________________________

# visualize whether there is any tradeoff between 
# P50stem and Ks or kstem
# P50leaf and kleaf


palette(cbp2)

quartz(width=7, height=3)
par(mar=c(3.5,4,1.5,1), mfrow=c(1,3), mgp=c(2.5,1,0),cex=1)
plot(mkstem~P50stem, all.ind, pch=16, col=factor(site)
     , ylab=expression(paste(k[stem], " (",mmol*m^-2*s^-1*MPa^-1, ")"))
     , xlab=expression(paste(P50[stem]," (MPa)")))

mtext("a) branch conductance", side=3, line=0, adj=0)



plot(mKs~P50stem, all.ind, pch=16, col=factor(site)
     , ylab=expression(paste(K[s], " (",mmol*m^-1*s^-1*kPa^-1, ")"))
     , xlab=expression(paste(P50[stem]," (MPa)")))
mtext("b) branch conductivity", side=3, line=0, adj=0)


plot(mkleaf~P50leaf, all.ind, pch=16, col=factor(site)
     , ylab=expression(paste(k[leaf], " (",mmol*m^-2*s^-1*MPa^-1, ")"))
     , xlab=expression(paste(P50[leaf]," (MPa)")))
mtext("c) leaf conductance", side=3, line=0, adj=0)
legend("topleft",legend = c("Garden","Wild"), col=c(1,2), pch=16, cex=.8, bty='n')


if(save.figures==T){
  quartz.save(file=paste0(results.dir,"/Fig11_SafetyEfficiencyTradeoff_v1.pdf"), type="pdf")
}




###########################################
#-------- Height tests ---------------------




