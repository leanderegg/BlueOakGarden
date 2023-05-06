######################################################
##             Code for Anderegg et al.              ##
##      "Plasticity drives geographic variation and trait 
##      coordination in blue oak drought physiology"
######################################################

# contact: Leander DL Anderegg
# landeregg@ucsb.edu

# Analyzes hydraulic and morphological trait data from a common garden and the 7 wild source populations
# of Quercus douglasii

##### NOTE: STILL NEED Ml_Ms data from W22!!!



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
#results.version <- "v190422" 
#results.version <- "v220817" # directory version for saving figures/results
#results.version <- "v221218"
#results.version <- "v230113" # updated the SLA from W22
#results.version <- "v230418" # updated leaf size and whatnot from W22
results.version <- "v230505" # full new version after code update
results.dir <- paste0("Results_",results.version)
if(save.figures == T) { dir.create(results.dir)}


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#######   BEGIN: LOAD DATA ###################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

popclim <- read.csv("Data/PopulationClimate_230418.csv")[,-1]
  # NOTE: THIS CURRENTLY HAS Canyon Ranch climate data for LPD Los Padres new site


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
# from pre-harvest survey, (2017?)
hoptrees <- read.csv("Data/Hop_trees_20180730.csv")
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
hoptrees$jjaT <- popclim$jja1951_1980[match(hoptrees$pop, popclim$Garden_pop)]
hoptrees$djfT <- popclim$djf1951_1980[match(hoptrees$pop, popclim$Garden_pop)]
hoptrees$log.bio <- log(hoptrees$bioest, base=10)


hoppop <- hoptrees %>% group_by(pop) %>% summarise(Bio = mean(bioest, na.rm=T), sdBio=sd(bioest,na.rm=T), seBio=se(bioest),
                                                   log.Bio =mean(log.bio, na.rm=T), sdlog.Bio=sd(log.bio,na.rm=T), selog.Bio=se(log.bio),
                                                       mort = length(which(is.na(height))), mort2 = length(which(stems==0)), mort3 = 18 - n() + length(which(is.na(height))),
                                                       Height = mean(height, na.rm=T), sdHeight=sd(height,na.rm=T), seHeight=se(height),
                                                       Diam = mean(diam_50cm, na.rm=T), sdDiam=sd(diam_50cm, na.rm=T),seDaim=se(diam_50cm),
                                                       cwd = unique(cwd), ppt = unique(ppt), aet=unique(aet), pet=unique(pet), tmn=unique(tmn), tmx=unique(tmx)) 



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
# throwing in dormant season values 
wg$cwd.2018ds <- popclim$cwd.2018ds[match(wg$Site, popclim$Code)]
wg$aet.2018ds <- popclim$aet.2018ds[match(wg$Site, popclim$Code)]
wg$pet.2018ds <- popclim$pet.2018ds[match(wg$Site, popclim$Code)]
wg$ppt.2018ds <- popclim$ppt.2018ds[match(wg$Site, popclim$Code)]
wg$tmn.2018ds <- popclim$tmn.2018ds[match(wg$Site, popclim$Code)]
wg$tmx.2018ds <- popclim$tmx.2018ds[match(wg$Site, popclim$Code)]

#_________________________________________________________________________
####### * Load all trait measurements for variance decomp  ##################
#_________________________________________________________________________

# cleaned kstem and Ks values
stemk.cl <- read.csv("Data/StemK_cleaned_20230428.csv", row.names=1)
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
all.ind <- read.csv("Data/AllIndividuals_H-W_traits_clim_20230505.csv")[,-1]


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
  kstem.sa.l = mean(mkstem.sa.l, na.rm=T), medkstem.sa.l = median(mkstem.sa.l, na.rm=T), sekstem.sa.l = se(mkstem.sa.l), lqkstem.sa.l=quantile(mkstem.sa.l, probs = .25, na.rm=T), uqkstem.sa.l = quantile(mkstem.sa.l, probs=.75, na.rm=T),
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
all.pop0$avail_soil_water <- popclim$avail_soil_water[match(all.pop0$pop, popclim$Garden_pop)]
all.pop0$dieback_10km <- popclim$proportion_10km[match(all.pop0$pop, popclim$Garden_pop)]
## add in a reverse of cwd for plotting purposes
all.pop0$cw.surp <- all.pop0$cwd *-1

# summarize wild growth to population
wg.pop <- wg %>% group_by(Site) %>% summarise(DBH = mean(DBH_cm), medDBH = median(DBH_cm, na.rm=T), seDBH = se(DBH_cm), lqDBH=quantile(DBH_cm, probs = .25, na.rm=T), uqDBH = quantile(DBH_cm, probs=.75, na.rm=T),
                                              BAI = mean(BAI_cm2), , medBAI = median(BAI_cm2, na.rm=T), seBAI = se(BAI_cm2), lqBAI=quantile(BAI_cm2, probs = .25, na.rm=T), uqBAI = quantile(BAI_cm2, probs=.75, na.rm=T),
                                              growthsup = mean(growthsup_cm2), , medgrowthsup = median(growthsup_cm2, na.rm=T), segrowthsup = se(growthsup_cm2), lqgrowthsup=quantile(growthsup_cm2, probs = .25, na.rm=T), uqgrowthsup = quantile(growthsup_cm2, probs=.75, na.rm=T),
                                              mperc_maxBAI=mean(perc_maxBAI), medperc_maxBAI = median(perc_maxBAI, na.rm=T), seperc_maxBAI = se(perc_maxBAI), lqperc_maxBAI=quantile(perc_maxBAI, probs = .25, na.rm=T), uqperc_maxBAI = quantile(perc_maxBAI, probs=.75, na.rm=T))
names(wg.pop)[grep("mperc", names(wg.pop))] <- "perc_maxBAI" # get rid of m added for summarising



all.pop.m <- melt(all.pop0,id.vars =c(1:4,60:68) ) 
all.pop.wide <- cast(all.pop.m, pop + pop.name + cwd + aet + pet + ppt + tmn + tmx + avail_soil_water + dieback_10km + cw.surp ~ site + variable, value.var="value")



### add in Hopland garden growth data for FINAL population averages
popw0 <- left_join(all.pop.wide, hoppop[,c(1:15)])
### add in wild growth data
popw <- left_join(popw0, wg.pop, by=c("pop.name"="Site"))


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
################## END: LOAD DATA ######################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++










#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
################ BEGIN: ANALYSIS #######################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



#______________________________________________________________________
######## * Population ANOVAs #################
#______________________________________________________________________

# function to calculate eta2 (among population variance)
eta2 <- function(aovres, group = "pop"){
  eta2 <- aovres[[1]]$'Sum Sq'[grep(group, rownames(aovres[[1]]))]/sum(aovres[[1]]$`Sum Sq`)
  return(eta2)
}

# function to calculate omega2 (unbiased among pop variance)
omega2 <- function(aovres, group = "pop"){
  omega2 <- (aovres[[1]]$'Sum Sq'[grep(group, rownames(aovres[[1]]))] - (aovres[[1]]$Df[grep(group, rownames(aovres[[1]]))]*aovres[[1]]$`Mean Sq`[grep("Residuals", rownames(aovres[[1]]))]))/(sum(aovres[[1]]$`Sum Sq`)+aovres[[1]]$`Mean Sq`[grep("Residuals", rownames(aovres[[1]]))])
  return(omega2)
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
HopAOV <- data.frame(Trait=traits.to.test, n = rep(NA, length(traits.to.test)), eta2 = rep(NA, length(traits.to.test)), omega2 = rep(NA, length(traits.to.test)),sig = rep(NA, length(traits.to.test)), pop.rand = rep(NA, length(traits.to.test)), bestclim = rep(NA, length(traits.to.test)), climvar=rep(NA, length(traits.to.test)), CV=rep(NA, length(traits.to.test)))

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
HopAOV$n[which(HopAOV$Trait=="Ks")] <- length(which(Hop.ind$mkstem.sa.l>0))
HopAOV$eta2[which(HopAOV$Trait=="Ks")] <- eta2(summary(aov(mkstem.sa.l~pop , Hop.ind)))
HopAOV$omega2[which(HopAOV$Trait=="Ks")] <- omega2(summary(aov(mkstem.sa.l~pop , Hop.ind)))
HopAOV$CV[which(HopAOV$Trait=="Ks")] <- sd(Hop.ind$mkstem.sa.l, na.rm=T)/mean(Hop.ind$mkstem.sa.l, na.rm=T)
HopAOV$sig[which(HopAOV$Trait=="Ks")] <- aovsig(summary(aov(mkstem.sa.l~pop , Hop.ind)))
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
WildAOV <- data.frame(Trait=traits.to.test, n = rep(NA, length(traits.to.test)), eta2 = rep(NA, length(traits.to.test)), omega2 = rep(NA, length(traits.to.test)),sig = rep(NA, length(traits.to.test)),pop.rand=rep(NA, length(traits.to.test)), bestclim = rep(NA, length(traits.to.test)), climvar=rep(NA, length(traits.to.test)), CV=rep(NA, length(traits.to.test)))

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
WildAOV$n[which(WildAOV$Trait=="Ks")] <- length(which(Wild.ind$mkstem.sa.l>0))
WildAOV$eta2[which(WildAOV$Trait=="Ks")] <- eta2(summary(aov(mkstem.sa.l~pop , Wild.ind)))
WildAOV$omega2[which(WildAOV$Trait=="Ks")] <- omega2(summary(aov(mkstem.sa.l~pop , Wild.ind)))
WildAOV$CV[which(WildAOV$Trait=="Ks")] <- sd(Wild.ind$mkstem.sa.l, na.rm=T)/mean(Wild.ind$mkstem.sa.l, na.rm=T)
WildAOV$sig[which(WildAOV$Trait=="Ks")] <- aovsig(summary(aov(mkstem.sa.l~pop , Wild.ind)))
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




#______________________________________________________________________
############### * Trait-Climate Analysis ###################
#______________________________________________________________________


### create a function that tests for univariate predictors of trait values in the wild and hopland populations
# using climate 

## function for visualization of linear models to check reasonableness
test.trait <- function(trait, dataz) {
  quartz(width=11, height=4.5)
  par(mfcol=c(2, 7), mar=c(4,0,0,0), oma=c(0,4,2,1), mgp=c(2.2,1,0))
  for(i in c("cwd","aet","pet","ppt","sitePD","siteMD","siteE.drop")){
    modW <- lm(get(trait)~get(i), dataz[which(dataz$site=="W"),])
    modH <- lm(get(trait)~get(i), dataz[which(dataz$site=="H"),])
    plot(get(trait)~get(i), dataz[which(dataz$site=="H"),], pch=1,xlab="", yaxt="n", ylab=trait, col=ifelse(summary(modH)$coefficients["get(i)",4]<0.05,yes = "red","black") )
    if(i=="cwd") {mtext(trait, side=2, line=2)
      axis(2)
      mtext("Hopland", side=3, adj = 0)}
    #    mtext(text = paste("site=", round(test[[1]]$`Pr(>F)`[1],2), " pop=", round(test[[1]]$`Pr(>F)`[[2]])),side=3)
    plot(get(trait)~get(i), dataz[which(dataz$site=="W"),], pch=3,xlab=i, yaxt="n", ylab=trait, col=ifelse(summary(modW)$coefficients["get(i)",4]<0.05,yes = "red","black")  )
    if(i=="cwd") {mtext(trait, side=2, line=2)
      axis(2)
      mtext("Wild", side=3, adj = 0)}
    
  }
  print("WILD:")
  print(summary(aov(get(trait)~ pop, all.ind[which(all.ind$site=="W"),])) )
  print("HOPLAND:")
  print(summary(aov(get(trait)~ pop, all.ind[which(all.ind$site=="H"),])) )
}

## function for visualization of linear mixed modesl to check reasonableness
test.mertrait <- function(trait, dataz) {
  quartz(width=11, height=4.5)
  par(mfcol=c(2, 7), mar=c(4,0,0,0), oma=c(0,4,2,1), mgp=c(2.2,1,0))
  for(i in c("cwd","aet","pet","ppt","sitePD","siteMD","siteE.drop")){
    modW <- lmer(get(trait)~get(i) + (1|pop.name), dataz[which(dataz$site=="W"),])
    modH <- lmer(get(trait)~get(i)+ (1|pop.name), dataz[which(dataz$site=="H"),])
    plot(get(trait)~get(i), dataz[which(dataz$site=="H"),], pch=1,xlab="", yaxt="n", ylab=trait, col=ifelse(summary(modH)$coefficients["get(i)",5]<0.05,yes = "red","black") )
    if(i=="cwd") {mtext(trait, side=2, line=2)
      axis(2)
      mtext("Hopland", side=3, adj = 0)}
    #    mtext(text = paste("site=", round(test[[1]]$`Pr(>F)`[1],2), " pop=", round(test[[1]]$`Pr(>F)`[[2]])),side=3)
    plot(get(trait)~get(i), dataz[which(dataz$site=="W"),], pch=3,xlab=i, yaxt="n", ylab=trait, col=ifelse(summary(modW)$coefficients["get(i)",5]<0.05,yes = "red","black")  )
    if(i=="cwd") {mtext(trait, side=2, line=2)
      axis(2)
      mtext("Wild", side=3, adj = 0)}
    
  }
  print(test <- summary(aov(get(trait)~ pop, all.ind[which(all.ind$site=="W"),])) )
  print(test <- summary(aov(get(trait)~ pop, all.ind[which(all.ind$site=="H"),])) )
}

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

## select best climate predictor for Wild data using linear models
# if no pop random effect needed
best.mod <- function(trait, dataz) {
  modnull <- lm(get(trait)~1 ,dataz)
  modcwd <- lm(get(trait)~cwd, dataz)
  modaet <- lm(get(trait)~aet, dataz)
  modpet <- lm(get(trait)~pet, dataz)
  modppt <- lm(get(trait)~ppt, dataz)
  modtmn <- lm(get(trait)~tmn, dataz)
  modtmn.2018ds <- lm(get(trait)~tmn.2018ds, dataz)
  modtmx.2018ds <- lm(get(trait)~tmx.2018ds, dataz)
  modppt.2018ds <- lm(get(trait)~ppt.2018ds, dataz)
  modcwd.2018ds <- lm(get(trait)~cwd.2018ds, dataz)
  #modPD <- lm(get(trait)~sitePD, dataz) # field water potentials, but only for 6 sites
  #modMD <- lm(get(trait)~siteMD, dataz)
  #modE.drop <- lm(get(trait)~siteE.drop, dataz)
  data.frame(AICc(modnull, modcwd, modaet, modpet, modppt, modtmn, modtmn.2018ds, modtmx.2018ds, modppt.2018ds, modcwd.2018ds)
             , n=c(length(resid(modnull)), length(resid(modcwd)), length(resid(modaet)), length(resid(modpet)), length(resid(modppt)), length(resid(modtmn)), length(resid(modtmn.2018ds)), length(resid(modtmx.2018ds)), length(resid(modppt.2018ds)), length(resid(modcwd.2018ds))))[order(AIC(modnull, modcwd, modaet, modpet, modppt, modtmn, modtmn.2018ds,modtmx.2018ds, modppt.2018ds,modcwd.2018ds)[,2]),]
  
}

## Seperate climate selection function for the garden, with only climate normals
best.modH <- function(trait, dataz) {
  modnull <- lm(get(trait)~1 ,dataz)
  modcwd <- lm(get(trait)~cwd, dataz)
  modaet <- lm(get(trait)~aet, dataz)
  modpet <- lm(get(trait)~pet, dataz)
  modppt <- lm(get(trait)~ppt, dataz)
  modtmn <- lm(get(trait)~tmn, dataz)
  #modtmn.2018ds <- lm(get(trait)~tmn.2018ds, dataz)
  #modtmx.2018ds <- lm(get(trait)~tmx.2018ds, dataz)
  #modppt.2018ds <- lm(get(trait)~ppt.2018ds, dataz)
  #modcwd.2018ds <- lm(get(trait)~cwd.2018ds, dataz)
  #modPD <- lm(get(trait)~sitePD, dataz)
  #modMD <- lm(get(trait)~siteMD, dataz)
  #modE.drop <- lm(get(trait)~siteE.drop, dataz)
  data.frame(AICc(modnull, modcwd, modaet, modpet, modppt, modtmn)
             , n=c(length(resid(modnull)), length(resid(modcwd)), length(resid(modaet)), length(resid(modpet)), length(resid(modppt)), length(resid(modtmn))))[order(AICc(modnull, modcwd, modaet, modpet, modppt, modtmn)[,2]),]
  
}

## climate variable selection function with mixed models
# use if pop random effect needed
best.mermod <- function(trait, dataz) {
  modnull <- lmer(get(trait)~1 + (1|pop.name),dataz, REML=F)
  modcwd <- lmer(get(trait)~cwd+ (1|pop.name), dataz, REML=F)
  modaet <- lmer(get(trait)~aet+ (1|pop.name), dataz, REML=F)
  modpet <- lmer(get(trait)~pet+ (1|pop.name), dataz, REML=F)
  modppt <- lmer(get(trait)~ppt+ (1|pop.name), dataz, REML=F)
  modtmn <- lmer(get(trait)~tmn+ (1|pop.name), dataz, REML=F)
  modtmn.2018ds <- lmer(get(trait)~tmn.2018ds+ (1|pop.name), dataz, REML=F)
  modtmx.2018ds <- lmer(get(trait)~tmx.2018ds+ (1|pop.name), dataz, REML=F)
  modppt.2018ds <- lmer(get(trait)~ppt.2018ds+ (1|pop.name), dataz, REML=F)
  modcwd.2018ds <- lmer(get(trait)~cwd.2018ds+ (1|pop.name), dataz, REML=F)
  #modPD <- lmer(get(trait)~sitePD+ (1|pop.name), dataz, REML=F)
  #modMD <- lmer(get(trait)~siteMD+ (1|pop.name), dataz, REML=F)
  #modE.drop <- lmer(get(trait)~siteE.drop+ (1|pop), dataz, REML=F)
  data.frame(AICc(modnull, modcwd, modaet, modpet, modppt, modtmn, modtmn.2018ds, modtmx.2018ds, modppt.2018ds, modcwd.2018ds)
             , n=c(length(resid(modnull)), length(resid(modcwd)), length(resid(modaet)), length(resid(modpet)), length(resid(modppt)), length(resid(modtmn)), length(resid(modtmn.2018ds)), length(resid(modtmx.2018ds)), length(resid(modppt.2018ds)), length(resid(modcwd.2018ds))))[order(AICc(modnull, modcwd, modaet, modpet, modppt, modtmn, modtmn.2018ds,modtmx.2018ds, modppt.2018ds,modcwd.2018ds)[,2]),]
  
}
## climate variable selection for garden with pop random effect.
best.mermodH <- function(trait, dataz) {
  modnull <- lmer(get(trait)~1 + (1|pop.name),dataz, REML=F)
  modcwd <- lmer(get(trait)~cwd+ (1|pop.name), dataz, REML=F)
  modaet <- lmer(get(trait)~aet+ (1|pop.name), dataz, REML=F)
  modpet <- lmer(get(trait)~pet+ (1|pop.name), dataz, REML=F)
  modppt <- lmer(get(trait)~ppt+ (1|pop.name), dataz, REML=F)
  modtmn <- lmer(get(trait)~tmn+ (1|pop.name), dataz, REML=F)
  #modtmn.2018ds <- lmer(get(trait)~tmn.2018ds+ (1|pop.name), dataz, REML=F)
  #modtmx.2018ds <- lmer(get(trait)~tmx.2018ds+ (1|pop.name), dataz, REML=F)
  #modppt.2018ds <- lmer(get(trait)~ppt.2018ds+ (1|pop.name), dataz, REML=F)
  #modcwd.2018ds <- lmer(get(trait)~cwd.2018ds+ (1|pop.name), dataz, REML=F)
  #modPD <- lmer(get(trait)~sitePD+ (1|pop.name), dataz, REML=F)
  #modMD <- lmer(get(trait)~siteMD+ (1|pop.name), dataz, REML=F)
  #modE.drop <- lmer(get(trait)~siteE.drop+ (1|pop), dataz, REML=F)
  data.frame(AICc(modnull, modcwd, modaet, modpet, modppt, modtmn)
             , n=c(length(resid(modnull)), length(resid(modcwd)), length(resid(modaet)), length(resid(modpet)), length(resid(modppt)), length(resid(modtmn))))[order(AICc(modnull, modcwd, modaet, modpet, modppt, modtmn)[,2]),]
  
}

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
HopAOV$pop.rand[which(HopAOV$Trait=="Ks")] <- test.rand(dataz=Hop.ind, trait="mkstem.sa.l")
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
WildAOV$pop.rand[which(WildAOV$Trait=="Ks")] <- test.rand(dataz=Wild.ind, trait="mkstem.sa.l")
WildAOV$pop.rand[which(WildAOV$Trait=="P50stem")] <- test.rand(dataz=Wild.ind, trait="P50stem")
WildAOV$pop.rand[which(WildAOV$Trait=="P50leaf")] <- test.rand(dataz=Wild.ind, trait="P50leaf")
WildAOV$pop.rand[which(WildAOV$Trait=="Growth")] <- "yes"

  ## test for Growth
# null <- gls(model = perc_maxBAI~ppt + pet + tmn, wg)
# pop <- lme(fixed = perc_maxBAI~ppt + pet + tmn, random= ~1|Site, data=wg)
# res <-anova(null, pop)$`p-value`[2] # need pop random intercepts for perc_maxBAI

  # all wild traits except ml_ms, Ks and P50stem/P50stem require random intercepts for pop

## Filling ANOVA tables:

##### . SLA #######
# Hopland
 best.modH("mSLA",dataz=Hop.ind)
HopmodSLA <- lm(mSLA~pet, Hop.ind)
# plot(HopmodSLA)
# PET
HopAOV$bestclim[which(HopAOV$Trait=="SLA")] <- "PET"
HopAOV$climvar[which(HopAOV$Trait=="SLA")] <- r.squaredLR(HopmodSLA)

# Wild
#best.mod("mSLA",dataz=Wild.ind)
#WildmodSLA <- lm(mSLA~tmx.2018ds, Wild.ind)
# plot(WildmodSLA)
best.mermod("mSLA",dataz=Wild.ind)
  # with linear model,tmn.2018ds is vaguely important but prob ns. with lmer, null is best
# tmx/tmin/ppt all marginally significant
WildAOV$bestclim[which(WildAOV$Trait=="SLA")] <- "none"
WildAOV$climvar[which(WildAOV$Trait=="SLA")] <- 0 #r.squaredLR(WildmodSLA)


##### . LDMC #######
# Hopland
best.modH("mLDMC",dataz=Hop.ind)
# HopmodLDMC <- lm(mLDMC~pet, Hop.ind)
## plot(HopmodLDMC)
# null
HopAOV$bestclim[which(HopAOV$Trait=="LDMC")] <- "none"
HopAOV$climvar[which(HopAOV$Trait=="LDMC")] <- 0

# Wild
# best.mod("mLDMC",dataz=Wild.ind)
#  WildmodLDMC <- lm(mLDMC~aet, Wild.ind)
# plot(WildmodLDMC)
best.mermod("mLDMC",dataz=Wild.ind)

# tmx/tmin/ppt all marginally significant with linear, but ns with mixed effects
WildAOV$bestclim[which(WildAOV$Trait=="LDMC")] <- "none" #"AET"
WildAOV$climvar[which(WildAOV$Trait=="LDMC")] <- 0 #r.squaredLR(WildmodLDMC)




##### . WD #######
# Hopland
best.modH("mWD",dataz=Hop.ind)
HopmodWD <- lm(mWD~ppt, Hop.ind)
# plot(HopmodWD)
# ppt
HopAOV$bestclim[which(HopAOV$Trait=="WD")] <- "PPT"
HopAOV$climvar[which(HopAOV$Trait=="WD")] <- r.squaredLR(HopmodWD)

# Wild
# best.mod("mWD",dataz=Wild.ind)
# WildmodWD <- lm(mWD~aet, Wild.ind)
# plot(WildmodWD)
best.mermod("mWD",dataz=Wild.ind)
WildmodWD <- lmer(mWD~aet + (1|pop.name), Wild.ind) #  totally ns
# aet was sig with linear, but not with mixed effects
WildAOV$bestclim[which(WildAOV$Trait=="WD")] <- "none" #"AET"
WildAOV$climvar[which(WildAOV$Trait=="WD")] <- 0 #r.squaredLR(WildmodWD)


##### . ml_ms #######
# Hopland
best.modH("mml_ms",dataz=Hop.ind)
Hopmodml_ms <- lm(mml_ms~ppt, Hop.ind)
# plot(Hopmodml_ms)
# null
HopAOV$bestclim[which(HopAOV$Trait=="ml_ms")] <- "none"
HopAOV$climvar[which(HopAOV$Trait=="ml_ms")] <- 0

# Wild
 best.mod("mml_ms",dataz=Wild.ind)
Wildmodml_ms <- lm(mml_ms~tmn.2018ds, Wild.ind) # ml_ms goes down with higher Tmin 2018, goes up with PPT
#plot(Wildmodml_ms)
#qqp(resid(Wildmodml_ms)) # problematic high risids, log-transform?

Wild.ind$log.ml_ms <- log(Wild.ind$mml_ms, base=10)
best.mod("log.ml_ms",dataz=Wild.ind)
  #but it doesn't fundamentally change any of the inferences (same order for predictors, still highly sig)

# Tmin and aet pretty similar (as is Tmin2018) were significant with linear model, but not mixed
WildAOV$bestclim[which(WildAOV$Trait=="ml_ms")] <- "winter Tmin"
WildAOV$climvar[which(WildAOV$Trait=="ml_ms")] <-  r.squaredLR(Wildmodml_ms)




##### . Al_As #######
# Hopland
best.modH("mAl_As",dataz=Hop.ind)
HopmodAl_As <- lm(mAl_As~tmn, Hop.ind)
# plot(HopmodAl_As)
# null
HopAOV$bestclim[which(HopAOV$Trait=="Al_As")] <- "none"
HopAOV$climvar[which(HopAOV$Trait=="Al_As")] <- 0

# Wild
# best.mod("mAl_As",dataz=Wild.ind)
# WildmodAl_As <- lm(mAl_As~tmn, Wild.ind) # Al_As goes down with higher Tmin, goes up with larger AET
# plot(WildmodAl_As)
best.mermod("mAl_As",dataz=Wild.ind)

# Tmin and aet pretty similar (as is Tmin2018) were significant with linear model, but not mixed
WildAOV$bestclim[which(WildAOV$Trait=="Al_As")] <- "none" #"Tmin"
WildAOV$climvar[which(WildAOV$Trait=="Al_As")] <- 0 # r.squaredLR(WildmodAl_As)



##### . leafsize #######
# Hopland
best.modH("mleafsize",dataz=Hop.ind)
Hopmodleafsize <- lm(mleafsize~ppt, Hop.ind) #increases with ppt
# plot(Hopmodleafsize)
# ppt, with 2 decent large outliers. more significant without outliers
HopAOV$bestclim[which(HopAOV$Trait=="leafsize")] <- "PPT"
HopAOV$climvar[which(HopAOV$Trait=="leafsize")] <- r.squaredLR(Hopmodleafsize)

# Wild
best.mermod("mleafsize",dataz=Wild.ind[-c(2,12),]) # marginal at best with or without two large outliers
#Wildmodleafsize <- lmer(mleafsize~tmn.2018ds + (1|pop.name), Wild.ind[-c(2,12),])
#plot(Wildmodleafsize) #ns, with two large outliers
#qqp(resid(Wildmodleafsize))
# Null best, with both mermod and lm, and with and without two large outliers
  # visual patterns very unconvincing as well  

WildAOV$bestclim[which(WildAOV$Trait=="leafsize")] <- "none"
WildAOV$climvar[which(WildAOV$Trait=="leafsize")] <- 0






##### . kleaf #######
# Hopland
best.modH("mkleaf",dataz=Hop.ind) 
Hopmodkleaf <- lm(mkleaf~tmn, Hop.ind) # delta AIC only 0.3 over null
#plot(Hopmodkleaf)
# tmn ns (p=0.11)
HopAOV$bestclim[which(HopAOV$Trait=="kleaf")] <- "none" #"nsTmin"
HopAOV$climvar[which(HopAOV$Trait=="kleaf")] <- 0 #r.squaredLR(Hopmodkleaf)

# Wild
# best.mod("mkleaf",dataz=Wild.ind)
# Wildmodkleaf <- lm(mkleaf~cwd.2018ds, Wild.ind) #cwd 2018ds and normal cwd best
#  plot(Wildmodkleaf)
best.mermod("mkleaf",dataz=Wild.ind)
Wildmodkleaf <- lmer(mkleaf~cwd.2018ds + (1|pop.name), Wild.ind) #cwd 2018ds and normal cwd best
# plot(Wildmodkleaf)
# qqp(resid(Wildmodkleaf))
# summary(Wildmodkleaf)
#cwd 2018ds and normal cwd best, but only ms
WildAOV$bestclim[which(WildAOV$Trait=="kleaf")] <- "winter CWD"
WildAOV$climvar[which(WildAOV$Trait=="kleaf")] <- r.squaredGLMM(Wildmodkleaf)[1]



##### . kstem #######
# Hopland
best.modH("mkstem",dataz=Hop.ind) 
Hopmodkstem <- lm(mkstem~pet, Hop.ind) # deltaAIC ~0 over null
# plot(Hopmodkstem)
# pet is nonsignificant, null is next best
HopAOV$bestclim[which(HopAOV$Trait=="kstem")] <- "none"
HopAOV$climvar[which(HopAOV$Trait=="kstem")] <-0

# Wild
# best.mod("mkstem",dataz=Wild.ind)
# Wildmodkstem <- lm(mkstem~tmn.2018ds, Wild.ind) 
#  plot(Wildmodkstem)
best.mermod("mkstem",dataz=Wild.ind)
Wildmodkstem <- lmer(mkstem~tmn.2018ds + (1|pop.name), Wild.ind) 
#plot(Wildmodkstem)
#qqp(resid(Wildmodkstem))
#summary(Wildmodkstem)
#tmn.2018ds and cwd 2018ds both pretty close, higher tmn, cwd, higher Kstem, similar with lm and lmer
WildAOV$bestclim[which(WildAOV$Trait=="kstem")] <- "winter Tmin"
WildAOV$climvar[which(WildAOV$Trait=="kstem")] <- r.squaredGLMM(Wildmodkstem)[1]



##### . Ks #######
# Hopland
best.modH("mkstem.sa.l",dataz=Hop.ind) 
#HopmodKs <- lm(mkstem.sa.l~tmx, Hop.ind) #null best
# plot(HopmodKs)
# null  best
HopAOV$bestclim[which(HopAOV$Trait=="Ks")] <- "none"
HopAOV$climvar[which(HopAOV$Trait=="Ks")] <-0

# Wild - doesnt need pop random intercept
best.mod("mkstem.sa.l",dataz=Wild.ind)
# WildmodKs <- lm(mkstem.sa.l~cwd.2018ds, Wild.ind) #null is best, slightly more significant without large outlier, but still ns
# plot(WildmodKs)
# qqp(resid(WildmodKs))

#null
WildAOV$bestclim[which(WildAOV$Trait=="Ks")] <- "none"
WildAOV$climvar[which(WildAOV$Trait=="Ks")] <- 0



##### . P50stem #######
# Hopland
best.modH("P50stem",dataz=Hop.ind) 
HopmodP50stem <- lm(P50stem~tmn, Hop.ind) #all equally bad
# plot(HopmodP50stem)
# tmn p=0.1468
HopAOV$bestclim[which(HopAOV$Trait=="P50stem")] <- "none"
HopAOV$climvar[which(HopAOV$Trait=="P50stem")] <- 0

# Wild
best.mod("P50stem",dataz=Wild.ind)
WildmodP50stem <- lm(P50stem~aet, Wild.ind) #aet best 
# plot(WildmodP50stem)# 16 is a bit of a high outlier with some leverage, but still ms without outlier

#aet and cwd best
WildAOV$bestclim[which(WildAOV$Trait=="P50stem")] <- "AET"
WildAOV$climvar[which(WildAOV$Trait=="P50stem")] <- r.squaredLR(WildmodP50stem) # very in line with results of Skelton et al. 2019



##### . P50leaf #######
# Hopland
best.modH("P50leaf",dataz=Hop.ind) 
# HopmodP50leaf <- lm(P50leaf~tmn, Hop.ind) #all equally bad
# plot(HopmodP50stem)
# tmn p=0.173
HopAOV$bestclim[which(HopAOV$Trait=="P50leaf")] <- "none"
HopAOV$climvar[which(HopAOV$Trait=="P50leaf")] <- 0

# Wild
best.mod("P50leaf",dataz=Wild.ind)
#WildmodP50stem <- lm(P50leaf~aet, Wild.ind) #null clearly the best 
# plot(WildmodP50leaf)

#aet and cwd best
WildAOV$bestclim[which(WildAOV$Trait=="P50leaf")] <- "none"
WildAOV$climvar[which(WildAOV$Trait=="P50leaf")] <- 0



##### . Growth #######
# Hopland
hoptrees$pop.name <- factor(hoptrees$pop)

modnull <- lmer(log.bio~1 + (1|pop.name),hoptrees, REML=F)
modcwd <- lmer(log.bio~cwd+ (1|pop.name), hoptrees, REML=F)
modaet <- lmer(log.bio~aet+ (1|pop.name), hoptrees, REML=F)
modpet <- lmer(log.bio~pet+ (1|pop.name), hoptrees, REML=F)
modppt <- lmer(log.bio~ppt+ (1|pop.name), hoptrees, REML=F)
modtmn <- lmer(log.bio~tmn+ (1|pop.name), hoptrees, REML=F)
data.frame(AICc(modnull, modcwd, modaet, modpet, modppt, modtmn)
           , n=c(length(resid(modnull)), length(resid(modcwd)), length(resid(modaet)), length(resid(modpet)), length(resid(modppt)), length(resid(modtmn))))[order(AICc(modnull, modcwd, modaet, modpet, modppt, modtmn)[,2]),]


HopmodGrowth<- lmer(log.bio~pet + (1|pop.name), hoptrees) 
#qqp(resid(HopmodGrowth)) #some major low outliers. but generally not horrible
#qqp(random.effects(HopmodGrowth)$pop.name[,1]) # random effects are decent
# pet
HopAOV$bestclim[which(HopAOV$Trait=="Growth")] <- "PET"
HopAOV$climvar[which(HopAOV$Trait=="Growth")] <- r.squaredGLMM(HopmodGrowth)[,1]

# Wild
  # for growth, current year values less relevant because we're using 5yrs worth of growth data
modnull <- lmer(perc_maxBAI~1 + (1|Site),wg, REML=F)
modcwd <- lmer(perc_maxBAI~cwd+ (1|Site), wg, REML=F)
modaet <- lmer(perc_maxBAI~aet+ (1|Site), wg, REML=F)
modpet <- lmer(perc_maxBAI~pet+ (1|Site), wg, REML=F)
modppt <- lmer(perc_maxBAI~ppt+ (1|Site), wg, REML=F)
modtmn <- lmer(perc_maxBAI~tmn+ (1|Site), wg, REML=F)
data.frame(AICc(modnull, modcwd, modaet, modpet, modppt, modtmn)
           , n=c(length(resid(modnull)), length(resid(modcwd)), length(resid(modaet)), length(resid(modpet)), length(resid(modppt)), length(resid(modtmn))))[order(AICc(modnull, modcwd, modaet, modpet, modppt, modtmn)[,2]),]


#null is the best
WildAOV$bestclim[which(WildAOV$Trait=="Growth")] <- "none"
WildAOV$climvar[which(WildAOV$Trait=="Growth")] <- 0

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
HopKsvd <- lmer(kstem.sa.l~ 1  + (1|site.pop/ind.id), data = stemk.cl, subset=site=="H")

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
variancesHop2 <- data.frame(HopP50stemvar[,4],HopP50leafvar[,4],HopGrowthvar[,4], row.names = c("btwPop","wiPop"))

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
WildKsvd <- lmer(kstem.sa.l~ 1  + (1|site.pop/ind.id), data = stemk.cl, subset=site=="W")

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
variancesWild2 <- data.frame(WildP50stemvar[,4],WildP50leafvar[,4],WildGrowthvar[,4], row.names = c("btwPop","wiPop"))

colnames(variancesWild2) <- c("P50stem","P50leaf","Growth")

scaledvariancesWild2 <- apply(variancesWild2, MARGIN=2, FUN= function(x){x/sum(x)})






#______________________________________________________________________
######## ** FIG1: Var Decomp #################
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
  quartz.save(file=paste0(results.dir,"/Fig_VarainceDecomp_scaled_v5.pdf"),type = "pdf" )
}

#_________________
####### Fig 1 continued: second bar plot for 3 traits without within-canopy ######
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
  quartz.save(file=paste0(results.dir,"/Fig_VarainceDecompp50Growth_scaled_v1.pdf"),type = "pdf" )
}




#________________________________________________________________________
############## ** FIG2: G vs G+E trait variation ###############
#_______________________________________________________________________


# Make a column for plotting significant traits as filled and ns as open
HopAOV$Plotting.Sig <- 1
HopAOV$Plotting.Sig[which(HopAOV$sig<= 0.05)]<-16

WildAOV$Plotting.Sig <- 1
WildAOV$Plotting.Sig[which(WildAOV$sig<= 0.05)]<-16

palette(mypal)

quartz(width=7, height=3.5)
layout(matrix(c(1,2), nrow=1),widths = c(1.3,1) )

par(mar=c(3.5,3.5,2,6), mgp=c(2.5,1,0))
plot(WildAOV$eta2*100~I(HopAOV$eta2*100), col=factor(WildAOV$Trait), pch=WildAOV$Plotting.Sig, cex=2, lwd=1.8,
     xlim=c(0,.7)*100, ylim=c(0,.7)*100,xaxs="i", yaxs="i",
     xlab="Garden: Btw-Pop Variation (%)", ylab="Wild: Btw-Site Variation (%)")
abline(a=0,b=1)
# plot the significant garden pop differences with a black circle
points(WildAOV$eta2[which(HopAOV$sig<0.05)]*100~I(HopAOV$eta2[which(HopAOV$sig<0.05)]*100), pch=3, cex=2, lwd=1.8, col=factor(HopAOV$Trait)[which(HopAOV$sig<0.1)])
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
legend("bottomright", legend=c("Garden sig", "Wild sig", "ns",""), pch=c(3,16,1,1),col=c("black","black","black","white"), bty="n", xpd=F)

legend(x=75,y=70, legend=legend.names, col=mypal[c(11,6,12,8,1,7,3,5,4,10,9,2)]
       , pch=15, xpd=NA,horiz = F, bty='n')
mtext(text = "a)", side=3, adj=-0.12, line=.4)


par(mar=c(3.5,3.5,2,1), mgp=c(2.5,1,0))
plot(I(climvar*100)~I(eta2*100), data=WildAOV[which(WildAOV$sig<= 0.05),], pch=16,cex=1.5,col=factor(WildAOV$Trait)[which(WildAOV$sig<= 0.05)],
     xlim=c(0,.7)*100, ylim=c(0,.7)*100,xaxs="i",
     xlab="Btw-Site or Btw-Pop Variation (%)", ylab="Variation Explained by Climate")
points(I(climvar*100)~I(eta2*100), data=HopAOV[which(HopAOV$sig<=0.05),], pch=3,cex=1.5,col=factor(HopAOV$Trait)[which(HopAOV$sig<=0.05)])
abline(a=0,b=1)
mtext(text = "b)", side=3, adj=-0.12, line=.4)

#legend('topleft', legend=c("Garden","Wild"), pch=c(3,17), bty="n")
if(save.figures==T){
  quartz.save(file=paste0(results.dir,"/Fig_G-vs-GE_variation_v5.pdf"), type="pdf")
}







#________________________________________________________________________
########## ** FIG3: Bootstrapping variance comparison #######################
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
  names(varratio) <- c("ratio.05", "ratio.10", "ratio.50", "ratio.90", "ratio.95")
  varsums <- c(varsumsHop, varsumsWild,varratio)
  return(varsums)
}

vars <- c("mSLA","mLDMC","mWD","mml_ms","mAl_As","mleafsize","mkleaf", "mkstem", "mkstem.sa.l", "P50stem","P50leaf")



varcomps <- var.boot(vars[1])
for(i in 2:length(vars)){
  varcomps <- rbind(varcomps, var.boot(vars[i]))
}
varcomps <- data.frame(varcomps)
rownames(varcomps) <- vars
#varcomps$sig

quartz(width=6, height=3)
par(mar=c(5,4,1,1))
tmp <- varcomps$Wild.50/varcomps$Hop.50
statest <- varcomps$Wild.05/varcomps$Hop.50
bp <- barplot(varcomps$ratio.50, las=2
              , names.arg = c("SLA","LDMC","WD","Ml:Ms","Al:As","Leaf size","k[leaf]","k[stem]", "Ks","P50[stem]","P50[leaf]")
              , ylab="Wild Var : Garden Var", ylim=c(0,2.9)
              , border = c(rep("white",3), rep("black", 4), "white",rep("black",3))
              , col = c(rep("lightgray",3), rep("white", 4), "lightgray",rep("white",3)))
#arrows(x0=bp, y0=varcomps$ratio.10, y1=varcomps$ratio.90, lwd=4, length = 0)
arrows(x0=bp, y0=varcomps$ratio.05, y1=varcomps$ratio.95, lwd=2, length = 0, col="grey")
arrows(x0=bp, y0=varcomps$ratio.10, y1=varcomps$ratio.90, lwd=4, length = 0, col="darkgray")


abline(h=1, lty=2)

if(save.figures==T){
  quartz.save(file=paste0(results.dir,"/Fig_WvGvariance_comparison_v2.pdf"),type = "pdf")
}





#_________________________________________________________________________________
################# ** FIG4: Trait Integration correlation plots #################
#_________________________________________________________________________________

#### . Calculate Correlations########

# define the variables we want to include

# V8 moved Ml:Ms back down
# V9, removed growth
varsH <- c("mleafsize","mAl_As","mkleaf", "mkstem", "mkstem.sa.l","mml_ms","mSLA","mLDMC","mWD", "P50stem","P50leaf")#, "Height")
varsW <- c("mleafsize","mAl_As","mkleaf", "mkstem", "mkstem.sa.l","mml_ms","mSLA","mLDMC","mWD", "P50stem","P50leaf")#, "perc_maxBAI")



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
#colnames(corAll.named) <- c("SLA","LDMC","WD","Ml:Ms","Al:As","Leaf size","k[leaf]","k[stem]", "Ks","P50[stem]","P50[leaf]")
#rownames(corAll.named)  <-c("SLA","LDMC","WD","Ml_Ms","Al:As","Leaf size","k[leaf]","k[stem]", "Ks","P50[stem]","P50[leaf]")
#colnames(corAll.named) <- c("SLA","LDMC","WD","Leaf size","Ml:Ms","Al:As","k[leaf]","k[stem]", "Ks","P50[stem]","P50[leaf]")
#rownames(corAll.named)  <-c("SLA","LDMC","WD","Leaf size","Ml_Ms","Al:As","k[leaf]","k[stem]", "Ks","P50[stem]","P50[leaf]")
# with reordered columns and growth
# colnames(corAll.named) <- c("Leaf size","Al:As","Ml:Ms","k[leaf]","k[stem]", "Ks","SLA","LDMC","WD", "P50[stem]","P50[leaf]", "30yr Height")
# rownames(corAll.named) <- c("Leaf size","Al:As","Ml:Ms","k[leaf]","k[stem]", "Ks","SLA","LDMC","WD", "P50[stem]","P50[leaf]", "30yr Height")

## reordered to highlight leaf size and hydraulics
# colnames(corAll.named) <- c("Leaf size","Al:As","k[leaf]","k[stem]", "Ks","Ml:Ms","SLA","LDMC","WD", "P50[stem]","P50[leaf]", "30yr Height")
# rownames(corAll.named) <- c("Leaf size","Al:As","k[leaf]","k[stem]", "Ks","Ml:Ms","SLA","LDMC","WD", "P50[stem]","P50[leaf]", "30yr Height")
## removed growth
colnames(corAll.named) <- c("Leaf size","Al:As","k[leaf]","k[stem]", "Ks","Ml:Ms","SLA","LDMC","WD", "P50[stem]","P50[leaf]")
rownames(corAll.named) <- c("Leaf size","Al:As","k[leaf]","k[stem]", "Ks","Ml:Ms","SLA","LDMC","WD", "P50[stem]","P50[leaf]")


corrplot::corrplot(corAll.named, p.mat=1/pmatAll*.01, sig.level=0.1, insig="pch",pch=16, pch.cex=.8
                   , type="lower", diag=F, add=F, cl.pos="n", method="ellipse", outline=F, addgrid.col = NA, tl.col="black"
)
if(save.figures==T){
  quartz.save(file=paste0(results.dir,"/Fig_Correlation_Heatmap_Hopland_v9.", type="pdf"))
}



quartz(width=4, height=4)
#corrplot::corrplot(corAll.named, p.mat=1/pmatAll*.01, sig.level=0.1, insig="pch",pch=16, pch.cex=.8
#                   , type="lower", diag=T, tl.pos="lt", tl.col="black", cl.pos="r", addgrid.col = NA)#, title="Wild = Upper, Garden = Lower")
corW.named <- corW
# colnames(corW.named) <- c("SLA","LDMC","WD","Ml:Ms","Al:As","Leaf size","k[leaf]","k[stem]", "Ks","P50[stem]","P50[leaf]")
# rownames(corW.named)  <-c("SLA","LDMC","WD","Ml:Ms","Al:As","Leaf size","k[leaf]","k[stem]", "Ks","P50[stem]","P50[leaf]")
#colnames(corW.named) <- c("SLA","LDMC","WD","Leaf size","Ml:Ms","Al:As","k[leaf]","k[stem]", "Ks","P50[stem]","P50[leaf]")
#rownames(corW.named)  <-c("SLA","LDMC","WD","Leaf size","Ml_Ms","Al:As","k[leaf]","k[stem]", "Ks","P50[stem]","P50[leaf]")

# with reordered columns and growth
# colnames(corW.named) <- c("Leaf size","Al:As","Ml:Ms","k[leaf]","k[stem]", "Ks","SLA","LDMC","WD", "P50[stem]","P50[leaf]", "% max BAI")
# rownames(corW.named) <- c("Leaf size","Al:As","Ml:Ms","k[leaf]","k[stem]", "Ks","SLA","LDMC","WD", "P50[stem]","P50[leaf]", "% max BAI")

# # move Ml:Ms down with SLA
# colnames(corW.named) <- c("Leaf size","Al:As","k[leaf]","k[stem]", "Ks","Ml:Ms","SLA","LDMC","WD", "P50[stem]","P50[leaf]", "% max BAI")
# rownames(corW.named) <- c("Leaf size","Al:As","k[leaf]","k[stem]", "Ks","Ml:Ms","SLA","LDMC","WD", "P50[stem]","P50[leaf]", "% max BAI")

# move Ml:Ms down with SLA
# removed growth
colnames(corW.named) <- c("Leaf size","Al:As","k[leaf]","k[stem]", "Ks","Ml:Ms","SLA","LDMC","WD", "P50[stem]","P50[leaf]")
rownames(corW.named) <- c("Leaf size","Al:As","k[leaf]","k[stem]", "Ks","Ml:Ms","SLA","LDMC","WD", "P50[stem]","P50[leaf]")

corrplot::corrplot(corW.named, p.mat=1/pmatW*.01, sig.level=0.1, insig="pch",pch=16, pch.cex=.8
                   , type="lower", add=F, diag=F, cl.pos="n", tl.col="black"
                   , method="ellipse", outline = F, addgrid.col = NA)
if(save.figures==T){
  quartz.save(file=paste0(results.dir,"/Fig_Correlation_Heatmap_Wild_v9.", type="pdf"))
}






#______________________________________________________________________________
########### ** FIG5: Bootstrapping Trait-Trait correlations ################
#______________________________________________________________________________



# Doing 10000 bootstrap samples from Wild vs Garden trait-trait correlations
W.H.test <- function(trait1, trait2, boots=1000, sig.type="alpha"){
  corrW <- corrH <- rep(0, boots)
  set.seed(42)
  
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
      hwcomp[j,i] <- W.H.test(traitz[i], traitz[j], boots=1000)
      # create a version scaled so that alpha <0.05 = 1 and >0.2 = 0
      hwcomprevscale[j,i] <- W.H.test(traitz[i], traitz[j], boots=1000, sig.type = "rev.scaled")
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

# get rid of the alpha<0.2 correlation changes (v5)
hwcomprevscale_trm <- hwcomprevscale
hwcomprevscale_trm[which(hwcomprevscale_trm<0.6)] <- 0

quartz(width=4, height=4)
corrplot::corrplot(hwcomprevscale_trm, p.mat=NULL, sig.level=0.04, insig="pch",pch=16, pch.cex=.8
                   , type="lower", add=F, diag=F, cl.pos="n", tl.col="black"
                   , method="color", outline = T, addgrid.col = F)
# v4 w 0.2 alpha cutoff
# legend("topright", bty="n", legend=c("alpha <0.2","alpha <0.1", "alpha <0.05"), pch=15, pt.cex = 2, col=corcolors[5:7], xpd=NA)
# v5 w/out .2, hwcomprevscaled_trm
legend("topright", bty="n", legend=c("alpha <0.1", "alpha <0.05"), pch=15, pt.cex = 2, col=corcolors[6:7], xpd=NA)


if(save.figures==T){
  quartz.save(file=paste0(results.dir,"/Fig_Correlation_SigFig_v5.pdf"), type="pdf")
}








#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
################ BEGIN: SUPPLEMENTAL FIGURES #######################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++






#______________________________________________________________________
######## ** FIG S1: Map #################
#______________________________________________________________________
mypallight <- paste0(brewer.pal(n=8, "Dark2"), "55")


# turn sampling locations into spatial data for plotting

latlon <- sp::CRS("+proj=longlat +datum=WGS84")
sampledlocs <- sp::SpatialPoints(coords = data.frame(popclim$Lon.dd, popclim$Lat.dd), proj4string = latlon)

#jpeg(filename="figures/FIG1_SamplingMap_v1.jpg",width =5, height=5, units = "in", res = 600)
quartz(width=3.8,height=5)
#par(mar=c(0,0,0,0))
CA <- maps::map('state', region = c('california'), col="black") 

points(latitude~longitude, qudo[which(qudo$latitude>33.5),], pch=16, col=mypallight[1])
legend('topright',legend = c("herbarium sample", 'common garden','source pops'), pch=c(16,17,15), col = c(mypallight[1],"black",mypal[4]),pt.bg=mypal[1], pt.cex = c(1,1.5,1), bty="n")
#points(gardenlocs.latlon)
points(sampledlocs[which(popclim$Garden_pop %in% c(1,3,15,16,18,22,26))], pch=15, cex=1, col=mypal[4])
# text(garden.sampledlocs.latlon, labels=gardens.sampled$Garden_pop)
points(sampledlocs[which(popclim$Code =="HREC")], pch=17,col="black", bg=mypal[1], cex=1.6)
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
  quartz.save(file=paste0(results.dir,"/Fig_Map_v1.pdf"), type="pdf")
}





#_________________________________________________________________________________
################# ** FIG S2: Wild vs Garden pop trait means #################
#_________________________________________________________________________________


quartz(width=6, height=6)
par(mfrow=c(3,3), mar=c(3.5,3.5,0,0), oma=c(0,0,1,1), mgp=c(2.1,1,0), cex.lab=1.2)
traits <- c("SLA","LDMC","WD","ml_ms","Al_As","leafsize","kleaf","kstem","kstem.sa.l")#,"P50stem")
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
         , xlim=sqlimit, ylim=sqlimit, pch=16, cex=2)
  mod <- lmodel2(get(paste0("W_", traits[i]))~get(paste0("H_", traits[i])), tmp)
  if(mod$P.param<0.1){
    abline(a=mod$regression.results$Intercept[3], b=mod$regression.results$Slope[3])
    #abline(a=mod$regression.results$Intercept[2], b=mod$regression.results$Slope[2])
    mtext(paste("p=", round(mod$P.param,3)),side = 1, line=-1, adj=.1)
  }
  
}

if(save.figures==T){
  quartz.save(file=paste0(results.dir,"/Fig_Garden-v-Wild_TraitValues_v4.pdf"), type="pdf")
}











#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
########### ** FIG S3: Supplemental Leaf size vs number on Al_As ** #####################
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
  quartz.save(file=paste0(results.dir,"/FigSX_LeafSizevsLeafNum-Al_As_v2.pdf"), type="pdf")
}





#_________________________________________________________________________________
################# ** FIG S4&S5: Population Average Trait-Growth relationships #################
#_________________________________________________________________________________

### Loop through each trait and plot with error bars

quartz(width=6.5, height=6.5)
par(mfrow=c(3,3), mar=c(3.2,2,0,0), oma=c(0,2,1,1), mgp=c(2,1,0), cex.lab=1.5 )
traits <- c("SLA","LDMC","WD","ml_ms","Al_As","leafsize","kleaf","kstem","kstem.sa.l","P50stem")
labs <- c("SLA","LDMC","WD","Ml:Ms","Al:As","Leaf size","k[leaf]","k[stem]","Ks","P50")

#traits <- c("P50leaf","P50stem","kstem","kstem.sa.l","kleaf","kleaf.sa.temp","Al_As","WD","SLA")
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
  points(Height~get(paste0("H_", traits[i])), popw, pch=16, cex=2)
  mod <- lm(Height~get(paste0("H_", traits[i])), tmp, weights=1/sqrt(tmp[,paste0("H_se", traits[i])]))
  if(summary(mod)$coefficients[2,4]<0.1){
    abline(mod, lty = ifelse(summary(mod)$coefficients[2,4]<0.05,1,2))
    mtext(paste("p=", round(summary(mod)$coefficients[2,4],3)),side = 1, line=-1, adj=.9)
  }
  
}

if(save.figures==T){
  quartz.save(file=paste0(results.dir,"/Fig_Trait-Growth_GardenHeight_v1.pdf"),type = "pdf")
}


### same figure with volume rather than height

quartz(width=6.5, height=6.5)
par(mfrow=c(3,3), mar=c(3.2,2,0,0), oma=c(0,2,1,1), mgp=c(2,1,0), cex.lab=1.5 )
traits <- c("SLA","LDMC","WD","ml_ms","Al_As","leafsize","kleaf","kstem","kstem.sa.l")#,"P50stem")
labs <- c("SLA","LDMC","WD","Ml:Ms","Al:As","Leaf size","k[leaf]","k[stem]","Ks")#,"P50")

#traits <- c("P50leaf","P50stem","kstem","kstem.sa.l","kleaf","kleaf.sa.temp","Al_As","WD","SLA")
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
  points(log.Bio~get(paste0("H_", traits[i])), popw, pch=16, cex=2)
  mod <- lm(log.Bio~get(paste0("H_", traits[i])), tmp, weights=1/sqrt(tmp[,paste0("H_se", traits[i])]))
  if(summary(mod)$coefficients[2,4]<0.1){
    abline(mod, lty = ifelse(summary(mod)$coefficients[2,4]<0.05,1,2))
    mtext(paste("p=", round(summary(mod)$coefficients[2,4],3)),side = 1, line=-1, adj=.9)
  }
  
}

if(save.figures==T){
  quartz.save(file=paste0(results.dir,"/Fig_Trait-Growth_GardenVolume_v1.pdf"),type = "pdf")
}




### Wild version
quartz(width=6.5, height=6.5)
par(mfrow=c(3,3), mar=c(3.2,2,0,0), oma=c(0,2,1,1), mgp=c(2,1,0), cex.lab=1.5 )
traits <- c("SLA","LDMC","WD","ml_ms","Al_As","leafsize","kleaf","kstem","kstem.sa.l","P50stem")
labs <- c("SLA","LDMC","WD","Ml:Ms","Al:As","Leaf size","k[leaf]","k[stem]","Ks","P50")

#traits <- c("P50leaf","P50stem","kstem","kstem.sa.l","kleaf","kleaf.sa.temp","Al_As","WD","SLA")
for(i in 1:9){
  # subset data down to make it easier to grep colnames
  tmp <- popw %>% select(grep("perc_maxBAI", colnames(popw)), grep(traits[i], colnames(popw)))
  # make plots appropriately large for error bars
  xlimit <- c(min(tmp[,paste0("W_",traits[i])], na.rm=T) - 1.2* max(tmp[,paste0("W_se",traits[i])], na.rm=T), 
              max(tmp[,paste0("W_",traits[i])], na.rm=T) + 1.2* max(tmp[,paste0("W_se",traits[i])], na.rm=T) )
  ylimit <- c(min(tmp$perc_maxBAI, na.rm=T, na.rm=T)-1.2*max(tmp$seperc_maxBAI,na.rm=T),
              max(tmp$perc_maxBAI,na.rm=T, na.rm=T) + 1.2 * max(tmp$seperc_maxBAI, na.rm=T))
  #  ylimit <- c(min(tmp$perc_maxBAI, na.rm=T), max(tmp$perc_maxBAI, na.rm=T))
  plot(perc_maxBAI~get(paste0("W_", traits[i])), popw
       , xlim=xlimit, ylim=ylimit, pch=16, cex=2
       , xlab=labs[i], ylab="")
  if(i==4){mtext("% of maximum size-specific BAI",side=2, line=2)}
  # add in y error bars (perc_maxBAI)
  arrows(x0=tmp[,paste0("W_", traits[i])], x1=tmp[,paste0("W_", traits[i])],
         y0=tmp$perc_maxBAI, y1=tmp$perc_maxBAI+tmp$seperc_maxBAI, length =0, lwd=2, col="grey")
  arrows(x0=tmp[,paste0("W_", traits[i])], x1=tmp[,paste0("W_", traits[i])],
         y0=tmp$perc_maxBAI, y1=tmp$perc_maxBAI-tmp$seperc_maxBAI, length =0, lwd=2, col="grey")
  # add in x error bars (trait)
  arrows(x=tmp[,paste0("W_", traits[i])], x1=tmp[,paste0("W_", traits[i])] + tmp[,paste0("W_se", traits[i])],
         y0=tmp$perc_maxBAI, lwd=2, length=0, col="grey")
  arrows(x=tmp[,paste0("W_", traits[i])], x1=tmp[,paste0("W_", traits[i])] - tmp[,paste0("W_se", traits[i])],
         y0=tmp$perc_maxBAI, lwd=2, length=0, col="grey")
  points(perc_maxBAI~get(paste0("W_", traits[i])), popw, pch=16, cex=2)
  mod <- lm(perc_maxBAI~get(paste0("W_", traits[i])), tmp, weights=1/sqrt(tmp[,paste0("W_se", traits[i])]))
  if(summary(mod)$coefficients[2,4]<0.1){
    abline(mod)
    mtext(paste("p=", round(summary(mod)$coefficients[2,4],3)),side = 1, line=-1, adj=.9)
  }
  
}

if(save.figures==T){
  quartz.save(file=paste0(results.dir,"/Fig_Trait-Growth_Wild_v1.pdf"),type = "pdf")
}





#_______________________________________________________________________
#################### ** FIG S6: Supplemental Safety-Efficiency Tradeoff ##############
#_______________________________________________________________________

# visualize whether there is any tradeoff between 
# P50stem and Ks or kstem
# P50leaf and kleaf


palette(cbp2)

quartz(width=7, height=3)
par(mar=c(3.5,4,1.5,1), mfrow=c(1,3), mgp=c(2.5,1,0),cex=1)
plot(mkstem~P50stem, all.ind, pch=16, col=factor(site)
     , ylab=expression(paste(k[stem], " (",mmol*cm^-2*s^-1*kPa^-1, ")"))
     , xlab=expression(paste(P50[stem]," (MPa)")))

mtext("a) branch conductance", side=3, line=0, adj=0)



plot(mkstem.sa.l~P50stem, all.ind, pch=16, col=factor(site)
     , ylab=expression(paste(K[s], " (",mmol*cm*s^-1*kPa^-1, ")"))
     , xlab=expression(paste(P50[stem]," (MPa)")))
mtext("b) branch conductivity", side=3, line=0, adj=0)


plot(mkleaf~P50leaf, all.ind, pch=16, col=factor(site)
     , ylab=expression(paste(k[leaf], " (",mmol*m^-2*s^-1*MPa^-1, ")"))
     , xlab=expression(paste(P50[leaf]," (MPa)")))
mtext("c) leaf conductance", side=3, line=0, adj=0)
legend("topleft",legend = c("Garden","Wild"), col=c(1,2), pch=16, cex=.8, bty='n')


if(save.figures==T){
  quartz.save(file=paste0(results.dir,"/FigSXY_SafetyEfficiencyTradeoff_v1.pdf"), type="pdf")
}










