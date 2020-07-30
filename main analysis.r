#---------------------------------------------------------------------------------#
#                                                                                 #
# Script to analyse how many and which above- and belowground traits are needed   #
# to accurately predict levels of ecosystem functioning in the Jena main          #
# experiment.                                                                     #
#                                                                                 #
# This is the main script, where:                                                 #
# 1. FI and FD metrics are calculated                                             #
# 2. all 41 traits are analysed as predictors                                     #
# 3. the above analysis is re-done with false discovery rates                     #
# 4. overlap in significant trait predictors is calculated                        #
# 5. step 4 is redone, but based on PCA variables instead of raw traits           #
# 6. species richness instead of traits is analyzed                               #
#                                                                                 #
# Authors: Thomas Schroeder-Georgi and Fons van der Plas                          #
#                                                                                 #
#---------------------------------------------------------------------------------#





### Switches ----------------------------------------------------------------------

analyse.richness <- 0 # if zero: don't include species richness in analysis, if 1: do include
FDR <- 0 # if zero: no false discover rate performed. If 1: false discovery rate
abundances.included <- 1 # if zero: FD and CWMs based on presence/absence data. If 1: based on abundances



#---------------------------------------------------------------------------------#
#                               Open R packages                                   #
#---------------------------------------------------------------------------------#

# set working directory
setwd("C:/2019.01.21.full.BU/64 Schroeder FD CWM/R files/")


# load required packages
require(lm.beta)
library(MuMIn)
library(lme4)
library(lmerTest)
library(sjstats)
library(car)
library(FD)
require(missMDA)





#---------------------------------------------------------------------------------#
#                               Some basic functions                              #
#---------------------------------------------------------------------------------#

# Function to calculate Rao's quadratic entropy (Botta-Dukat, 2005)
# a = community matrix: rows=sides, columns=species, rownames=sideslables
# d = any kind of distance matrix
# species names in "a" have to match species names in "d" being in the same order.

Raosq <- function(a, d){ 
  d <- as.matrix(d)
  fdq <- rep(NA, nrow(a))
  names(fdq) <- rownames(a)
  for(i in 1:nrow(a)){
    com <- as.matrix(a[i, is.na(a[i,])==F])
    dis <- as.matrix(d)[rownames(com), rownames(com)]
    p <- as.vector(com/sum(com))
    pp <- p %o% p
    if(sum(dis) > 0){
      fdq[i] <- sum(dis * pp) / 2
    }else{
      fdq[i] <- 0
    }
  }
  return(fdq)
}

# Function to calculate community-weighted means of traits
CWM <- function(x, a){ 
  cwm <- matrix(NA, nrow(a), ncol(x))
  rownames(cwm) <- rownames(a)
  colnames(cwm) <- colnames(x)
  for(i in 1:nrow(cwm)){
    com <- as.matrix(a[i,is.na(a[i,])==F])
    t <- as.matrix(x[rownames(com),])
    p <- as.vector(com/sum(com))
    pt <- t * p
    cwm[i,] <- apply(pt,2,sum)
  }
  return(cwm)
}

### scale reflection
scale.reflection <- function(x){
  x2 <- x * (-1)
  x2 <- x2 + max(x, na.rm=T)
  return(x2)
}

### function to standardize data by their maximum
standardize <- function(x){
  minx <- min(x, na.rm=T)
  if(minx<0){
    x <- x + abs(minx)
  }
  s <- x / max(x, na.rm=T)
  return(s)
}

# MFw = FT/F x ST/S x mean(fT) = FT/F x ST/S x sum(fT)/FT
# -> MFw = ST/S * sum(fT)/F
MFw <- function(x,S,TH){ 
  dat <- na.omit(data.frame(S,x))
  S <- length(unique(dat$S))
  Fun <- nrow(dat)
  dat <- dat[dat$x >= TH,]
  ST <- length(unique(dat$S))
  FT <- nrow(dat)
  MF <- 0
  if(FT>0){
    MF <- ST/S*sum(dat$x)/Fun
  }
  MF <- round(MF, 4)
  return(MF)
}





#---------------------------------------------------------------------------------#
#                         Load trait and function data                            #
#---------------------------------------------------------------------------------#

traits.st <- read.table("traits.st.csv", header=T, sep=",")
row.names(traits.st) <- traits.st$sp
functions4.st <- read.table("functions4.st.csv", header=T, sep=",")





#---------------------------------------------------------------------------------#
# Load and prepare community data                                                 #
#---------------------------------------------------------------------------------#

# load species cover data
cover <- read.csv("CoverMain-JenaExp2003-2015all.csv", header=T, sep=",") 
head(cover)

# replace missing values by 0
cover[is.na(cover)] <- 0
cover[cover=="NULL"] <- 0 # change covers of NULL to 0
cover[,11]
for(i in 4:ncol(cover)){
  if(is.character(cover[,i])){
    droplevels(cover[,i])
  }
  cover[,i] <- as.numeric(as.character(cover[,i]))
  print(i)
}
cover[is.na(cover)] <- 0
cover[cover==-9999] <- 0 # change covers of -9999 to 0
cover[cover=="NULL"] <- 0 # change covers of -9999 to 0

# correct some values
cover[cover==1] <- 0.5 # change covers of 1 to 0.5
cover[cover==2] <- 3 # change covers of 2 to 3

# add columns
cover2 <- data.frame("plotcode"=cover$Plot, "year"=cover$Year,
                     "season"=cover$Month, "sowndiv"=rep(NA, nrow(cover)),
                     "target.biomass"=rep(NA, nrow(cover)))
cover2 <- cbind(cover2, cover[,c(4:63)])

# change month data into season data
levels(cover2$season) <- c("august", "may", "sept", "spring", "summer")
cover2$season[which(cover2$season=="may")] <- "spring"
cover2$season[which(cover2$season=="august")] <- "summer"
cover2$season[which(cover2$season=="sept")] <- "summer"
cover2 <- cover2[which(cover2$year %in% c(2002:2012)),]

# load target biomass data
community <- read.csv("biomass02_12.csv", header=T, sep=",") 

# replace -9999 by NA
community2 <- replace(community, community == -9999, NA)
summary(community2)

# monoculture plots:
mono <- paste(community2[community2$sowndiv==1, 1],
              community2[community2$sowndiv==1, 5], sep="_")
mono <- mono[!duplicated(mono)]

# remove monocultures:
community3 <- community2[community2$sowndiv>1, ] # without monocultures

# show plots, without any target species and target.biomass = -9999:
plots_target_0 <- community3[community3$target.biomass==0,][,c(1,5)]
plots_target_0 <- paste(plots_target_0[,1],plots_target_0[,2],sep="_")

# remove plots without any target species:
community3 <- community2[community2$target.biomass>0,]
community3 <- community3[,c(1,4,5,6,12,84:144)] # the communities of interest
summary(community3)

# mean of subsamples
cnames <- colnames(community3)
community4 <- aggregate(as.matrix(community3[,6:ncol(community3)]) ~
                        community3$plotcode + community3$year +
                        community3$month + community3$sowndiv,
                        FUN="median", na.action=na.pass, na.rm=T) 
colnames(community4) <- cnames[-2]

# change levels of community3$month ( "Aug"  "Jun"  "May"  "Sept")
# to "spring" and "summer":
levels(community4$month) <- c("summer","spring","spring","summer")
com <- community4
colnames(com)[3] <- "season"
rownames(com) <- paste(com$plotcode, com$year, com$season, sep=".")
dim(com)
names(cover2)[c(6:65)] <- names(com)[c(6:65)]

# real species richness:
lenx <- function(x){length(x[!is.na(x) & x>1])}
SR <- apply(com[,6:65], 1, lenx)

# remove plots with less than 1 target species present
com <- com[SR>=1,]

# real species richness:
SR2 <- apply(cover2[,6:65], 1, lenx)
cover2 <- cover2[SR2>=1,]

# harmonize column names
names(cover2)[6:65] <- names(com)[6:65]

# add correct sowndiversity and target biomass data
cover3 <- merge(cover2, com[,c(1,2,3,4,5)], by=c("plotcode","year","season"),
                all=T)
cover4 <- cover3[,c(-66,-67)]
cover4[,c(4,5)] <- cover3[,c(66,67)]
names(cover4)[c(4,5)] <- c("sowndiv","target.biomass")

# fill out missing sowndiv levels
sowndiv.levels <- aggregate(sowndiv ~ plotcode, data=cover4, mean, na.rm=T)

cover5 <- cover4
for(i in 1:length(cover4$sowndiv[is.na(cover4$sowndiv)])){
  cover5$sowndiv[is.na(cover4$sowndiv)][i] <- sowndiv.levels[which(sowndiv.levels[,1]==cover4$plotcode[is.na(cover4$sowndiv)][i]),2] 
}

# remove sowndiv levels of 0 and 60
cover5 <- cover5[-which(cover5$sowndiv==60),]
cover5 <- cover5[-which(cover5$sowndiv==0),]

# drop redundant season and plotcode levels
cover5$season <- droplevels(cover5$season)
cover5$plotcode <- droplevels(cover5$plotcode)

# add rownames
rownames(cover5) <- paste(cover5$plotcode, cover5$year, cover5$season, sep=".")

# idenity monoculture plots:
mono <- paste(cover5[cover5$sowndiv==1, 1], cover5[cover5$sowndiv==1, 2], sep="_")
mono <- mono[!duplicated(mono)]

# calculate realized species richness:
lenx <- function(x){length(x[!is.na(x) & x>1])}
SR <- apply(cover5[,6:65], 1, lenx)

# remove plots with less than 1 target species present
cover5 <- cover5[SR>=1,]

# rename cover dataframe
com <- cover5 





#---------------------------------------------------------------------------------#
# Calculate FD and CWM data based on raw trait values                                                 #
#---------------------------------------------------------------------------------#

x <- traits.st[,-1]

# Josse, J & Husson, F. (2013). Handling missing values in exploratory
# multivariate data analysis methods. Journal de la SFdS. 153 (2), pp. 79-99.
nb <- estim_ncpPCA(traits.st[,-1], ncp.min=0, ncp.max=10, method.cv="Kfold",
                   nbsim=50) # estimate the number of dimensions needed for NA imputation in impute.PCA
imputed <- imputePCA(traits.st[,-1], ncp=nb$ncp, method="Regularized") # NA imputation, creating a complete dataset using PCA predictions
x <- imputed$completeObs

# community data:
a <- com[,6:ncol(com)]

if(abundances.included == 1){
  a <- a
}else{
  for(i in 1:ncol(a)){
    a[which(a[,i]>0),i] <- 1
  }
}

# in summer 2004 there are all target.biomass = 0
# order species in x and a alphabetically:
x <- x[order(rownames(x)),]
a <- as.matrix(a[,rownames(x)])
species.names <- colnames(a)
a <- as.numeric(as.character(a))
a <- matrix(a, ncol=59)
colnames(a) <- species.names

# Calculate RaosQ as a measure of FD:
FD_raw <- com[,1:3]

for(i in colnames(x)){
  x.i <- as.data.frame(x[,i])
  rownames(x.i) <- rownames(x)
  assign(paste("d", i, sep="."), gowdis(x.i))
  assign(paste("FDq", i, sep="."), Raosq(a,get(paste("d", i, sep="."))))
  FD_raw <- cbind(FD_raw, get(paste("FDq", i, sep=".")))
  colnames(FD_raw)[ncol(FD_raw)] <- paste("FD", i, sep="_")
}
FD_raw

# calculate CWMs
cwm_raw <- CWM(x, a)
colnames(cwm_raw) <- paste("FI", colnames(cwm_raw), sep="_")
cwm_raw <- data.frame(com[,1:3], cwm_raw)

# merge FD and CWM data
predictors_raw <- merge(FD_raw, cwm_raw)
predictors_raw <- predictors_raw[,-c(4,6,48,50)] # take away traits leading to circularity
head(predictors_raw)
dim(predictors_raw)

# change colnames of predictors ("plotcode" to "plot")
colnames(predictors_raw)[1] <- "plot"
write.table(predictors_raw, "predictors_raw.csv", row.names=F, sep=",")
predictors_raw <- read.table("predictors_raw.csv", header=T, sep=",")

# calculate trait correlations
correlations.traits <- round(cor(x), 2)
correlations.traits <- round(cor(x[,-c(1,3,29)]), 2)
row.names(correlations.traits)
correlations.traits2 <- correlations.traits[c(16,21,10,2,3,15,23,14,25,12,
                                              13,26,20,19,22,24,17,18,11,30,
                                              31,36,32,35,33,34,29,1,9,8,
                                              4,5,27,6,28,40,38,41,37,39,
                                              7),
                                            c(16,21,10,2,3,15,23,14,25,12,
                                               13,26,20,19,22,24,17,18,11,30,
                                               31,36,32,35,33,34,29,1,9,8,
                                               4,5,27,6,28,40,38,41,37,39,
                                               7)]
write.table(correlations.traits2, "correlations.traits.csv", row.names=T, sep=",")





#---------------------------------------------------------------------------------#
# test for best model structure: year as fixed or random factor?                  #
#---------------------------------------------------------------------------------#

# include block in predictors
block <- as.numeric(as.character(substr(predictors_raw$plot,2,2)))
predictors_raw <- data.frame(predictors_raw[,c(1,2,3)], block,
                             predictors_raw[,-c(1,2,3)])

# exclude redundant traits
names(predictors_raw)
predictors_raw <- predictors_raw[,-c(31, 73)] # leaf C:N

# count number of functions
no.functions <- ncol(functions4.st) - 4

if(analyse.richness==0){
  significant.predictors <- as.data.frame(matrix(rep(0,(no.functions*
                            (ncol(predictors_raw[,-c(1,2,3,4)])))),
                             nrow=no.functions, ncol=(ncol(predictors_raw[
                             ,-c(1,2,3,4)]))))
  names(significant.predictors) <- 
  c("FD_srr", "FD_Nacr", "FD_h1", "FD_lbpr", "FD_la", "FD_lSize", "FD_lt",
    "FD_sla", "FD_lsd", "FD_lar", "FD_fc", "FD_ldmc", "FD_leafC", "FD_leafN",
    "FD_cond", "FD_pp", "FD_sd", "FD_ssd", "FD_erect",
    "FD_BM_Target_Inflorescence", "FD_NrInflorescence_Target_Shoot",
    "FD_DurationFlowreing_Target_Phenology", "FD_SeedsProjArea",
    "FD_NrSeedling_Target", "FD_SeedsTCW", "FD_SeedsWLRatio", "FD_sdmc", "FD_ra",
    "FD_RootingDepth_Target", "FD_k.rad", "FD_sra", "FD_srl", "FD_rtd", "FD_RNU",
    "FD_RCN", "FD_rootP", "FD_rootK", "FD_rootS", "FD_rootCa", "FD_rootNa",
    "FD_nute", "FI_srr", "FI_Nacr", "FI_h1", "FI_lbpr", "FI_la", "FI_lSize",
    "FI_lt", "FI_sla", "FI_lsd", "FI_lar", "FI_fc", "FI_ldmc", "FI_leafC",
    "FI_leafN", "FI_cond", "FI_pp", "FI_sd", "FI_ssd", "FI_erect",
    "FI_BM_Target_Inflorescence", "FI_NrInflorescence_Target_Shoot", 
    "FI_DurationFlowreing_Target_Phenology", "FI_SeedsProjArea",
    "FI_NrSeedling_Target", "FI_SeedsTCW", "FI_SeedsWLRatio", "FI_sdmc", "FI_ra",
    "FI_RootingDepth_Target", "FI_k.rad", "FI_sra", "FI_srl", "FI_rtd", "FI_RNU",
    "FI_RCN", "FI_rootP", "FI_rootK", "FI_rootS", "FI_rootCa", "FI_rootNa",
    "FI_nute")
}else if(analyse.richness==1){
  significant.predictors <- as.data.frame(matrix(rep(0,(no.functions*
                            (ncol(predictors_raw[,-c(1,2,3,4)])+1))),
                            nrow = no.functions, ncol=(1+ncol(predictors_raw[,-c(1,2,3,4)]))))
  names(significant.predictors) <- 
    c("sowndiv","FD_srr", "FD_Nacr", "FD_h1", "FD_lbpr", "FD_la", "FD_lSize",
      "FD_lt", "FD_sla", "FD_lsd", "FD_lar", "FD_fc", "FD_ldmc", "FD_leafC",
      "FD_leafN", "FD_cond", "FD_pp", "FD_sd", "FD_ssd", "FD_erect",
      "FD_BM_Target_Inflorescence", "FD_NrInflorescence_Target_Shoot",
      "FD_DurationFlowreing_Target_Phenology", "FD_SeedsProjArea",
      "FD_NrSeedling_Target", "FD_SeedsTCW", "FD_SeedsWLRatio", "FD_sdmc", "FD_ra",
      "FD_RootingDepth_Target", "FD_k.rad", "FD_sra", "FD_srl", "FD_rtd", "FD_RNU",
      "FD_RCN", "FD_rootP", "FD_rootK", "FD_rootS", "FD_rootCa", "FD_rootNa",
      "FD_nute", "FI_srr", "FI_Nacr", "FI_h1", "FI_lbpr", "FI_la", "FI_lSize",
      "FI_lt", "FI_sla", "FI_lsd", "FI_lar", "FI_fc", "FI_ldmc", "FI_leafC",
      "FI_leafN", "FI_cond", "FI_pp", "FI_sd", "FI_ssd", "FI_erect",
      "FI_BM_Target_Inflorescence", "FI_NrInflorescence_Target_Shoot", 
      "FI_DurationFlowreing_Target_Phenology", "FI_SeedsProjArea",
      "FI_NrSeedling_Target", "FI_SeedsTCW", "FI_SeedsWLRatio", "FI_sdmc", "FI_ra",
      "FI_RootingDepth_Target", "FI_k.rad", "FI_sra", "FI_srl", "FI_rtd", "FI_RNU",
      "FI_RCN", "FI_rootP", "FI_rootK", "FI_rootS", "FI_rootCa", "FI_rootNa",
      "FI_nute")
}

effect.sizes <- significant.predictors
effect.sizes[effect.sizes==0] <- NA

nr_predictors <- c()
R2.model <- c()
R2.model.cond <- c()
function.names <- c()
year.analyis <- c()
number.analysis <- 0
predictors_vector <- ncol(predictors_raw[,-c(1,2,3,4)])

# test how to best treat year
AIC.values <- c()
for(i in 1:(ncol(functions4.st)-4)){
  j <- i + 4
  response <- functions4.st[,c(j,1:4)]
  nrow(response)
  
  predictors_raw.2 <-
    aggregate(cbind(FD_srr, FD_Nacr, FD_h1, FD_lbpr, FD_la, FD_lSize, FD_lt,
                    FD_sla, FD_lsd, FD_lar, FD_fc, FD_ldmc, FD_leafC, FD_leafN,
                    FD_rootP, FD_rootK, FD_rootS, FD_rootCa, FD_rootNa, FD_cond,
                    FD_pp, FD_sd, FD_ssd, FD_erect, FD_BM_Target_Inflorescence,
                    FD_NrInflorescence_Target_Shoot, 
                    FD_DurationFlowreing_Target_Phenology, FD_SeedsProjArea,
                    FD_NrSeedling_Target, FD_SeedsTCW, FD_SeedsWLRatio, FD_sdmc,
                    FD_ra, FD_RootingDepth_Target, FD_k.rad, FD_sra, FD_srl,
                    FD_rtd, FD_RNU, FD_RCN, FD_nute, FI_srr, FI_Nacr, FI_h1,
                    FI_lbpr, FI_la, FI_lSize, FI_lt, FI_sla, FI_lsd, FI_lar,
                    FI_fc, FI_ldmc, FI_leafC, FI_leafN, FI_rootP, FI_rootK,
                    FI_rootS, FI_rootCa, FI_rootNa, FI_cond, FI_pp, FI_sd, FI_ssd,
                    FI_erect, FI_BM_Target_Inflorescence,
                    FI_NrInflorescence_Target_Shoot, 
                    FI_DurationFlowreing_Target_Phenology, FI_SeedsProjArea,
                    FI_NrSeedling_Target, FI_SeedsTCW, FI_SeedsWLRatio, FI_sdmc,
                    FI_ra, FI_RootingDepth_Target, FI_k.rad, FI_sra, FI_srl,
                    FI_rtd, FI_RNU, FI_RCN, FI_nute)~
                    year+plot, predictors_raw, mean)
  
  nrow(predictors_raw.2)
  all.data <- merge(response, predictors_raw.2, by=c("plot","year"), all.y=T, all.x=T)
  names(all.data)[3] <- "y"
  predictors.test <- all.data[complete.cases(all.data$y),c(-4)][c(1,2,3,4)]
  model.formula.1 <- as.formula(
    paste("y", "~", paste(names(predictors.test)[c(-1,-2,-3,-4)], collapse=" + "), " + ", "(1|year) + (1|block)", sep="")
  )
  m1 <- lmer(model.formula.1, data=predictors.test)
  model.formula.2 <- as.formula(
    paste("y", "~", paste(names(predictors.test)[c(-1,-2,-3,-4)], collapse=" + "), " + year + ", "(1|block)", sep="")
  )
  m2 <- lmer(model.formula.2, data=predictors.test)
  AIC.values[i] <- which(c(AIC(m1),AIC(m2))==min(c(AIC(m1),AIC(m2))))
}
AIC.values # generally best to treat year as a random factor (not fixed)





#---------------------------------------------------------------------------------#
# Main analysis: which trait relate to which functions? Depending on switch       #
# setting (line 26) either conventional P value based, or based on False          #
# Discovery Rates                                                                 #
#---------------------------------------------------------------------------------#

# loop for single functions to select the most parsimonious model
q.star <- 0.05 # for false discovery rates: initial q.star = 0.05
for(i in 1:(ncol(functions4.st)-4)){
  j <- i + 4
  response <- functions4.st[,c(j,1:4)]
  nrow(response)
  predictors_raw.2 <-
    aggregate(cbind(FD_srr, FD_Nacr, FD_h1, FD_lbpr, FD_la, FD_lSize, FD_lt,
              FD_sla, FD_lsd, FD_lar, FD_fc, FD_ldmc, FD_leafC, FD_leafN,
              FD_rootP, FD_rootK, FD_rootS, FD_rootCa, FD_rootNa, FD_cond, FD_pp,
              FD_sd, FD_ssd, FD_erect, FD_BM_Target_Inflorescence,
              FD_NrInflorescence_Target_Shoot,
              FD_DurationFlowreing_Target_Phenology, FD_SeedsProjArea,
              FD_NrSeedling_Target, FD_SeedsTCW, FD_SeedsWLRatio, FD_sdmc, FD_ra,
              FD_RootingDepth_Target, FD_k.rad, FD_sra, FD_srl,FD_rtd, FD_RNU,
              FD_RCN, FD_nute, FI_srr, FI_Nacr, FI_h1, FI_lbpr, FI_la, FI_lSize,
              FI_lt, FI_sla, FI_lsd, FI_lar, FI_fc, FI_ldmc, FI_leafC, FI_leafN,
              FI_rootP, FI_rootK, FI_rootS, FI_rootCa, FI_rootNa, FI_cond, FI_pp,
              FI_sd, FI_ssd, FI_erect, FI_BM_Target_Inflorescence,
              FI_NrInflorescence_Target_Shoot, 
              FI_DurationFlowreing_Target_Phenology, FI_SeedsProjArea,
              FI_NrSeedling_Target, FI_SeedsTCW, FI_SeedsWLRatio, FI_sdmc, FI_ra,
              FI_RootingDepth_Target, FI_k.rad, FI_sra, FI_srl, FI_rtd, FI_RNU,
              FI_RCN, FI_nute)~year+plot, predictors_raw, mean)
 
  nrow(predictors_raw.2)
  all.data <- merge(response, predictors_raw.2, by=c("plot","year"),
                    all.y=T, all.x=T)
  nrow(all.data)
  if(analyse.richness==0){
    all.data <- all.data[,c(1:3,5:87)]
  }else if(analyse.richness==1){
    all.data <- all.data[,c(1:3,5,4,6:87)]
  }
  
  P.value <- 0
  n.predictors.left <- ncol(all.data[,c(5:ncol(all.data))])
  predictors <- c()
  predictors.include2 <- predictors.include3 <- c()
  names(all.data)[3] <- "y"
  count <- 0
  alpha <- q.star

  while(P.value < alpha && n.predictors.left > 0){
    count <- count + 1
    P.values <- c()
    VIF.values <- c()
    predictors.test.range <- c()
    for(k in 5:(n.predictors.left+4)){
      if(length(predictors)==0){
        predictors.test <- all.data[complete.cases(all.data$y),][c(1,2,3,4,k)]
        model.formula <- as.formula(
          paste("y", "~", paste(names(predictors.test)[c(-1,-2,-3,-4)], collapse=" + "), " + ", "(1|year) + (1|plot)", sep="")
        )
        m1 <- lmer(model.formula, data=predictors.test)
        P.values[k] <- summary(m1)$coefficients[-1,5]
        VIF.values[k] <- 0
        predictors.test.range[k] <- paste(names(predictors.test)[c(-1,-2,-3,-4)], sep=",")
      }else{
        all.data.3 <- cbind(all.data[,c("y", "year", "plot", "block", predictors.include3)], all.data.2[,k])
        names(all.data.3)[ncol(all.data.3)] <- names(all.data.2)[k]
        predictors.test <- all.data.3[complete.cases(all.data$y),]
        model.formula <- as.formula(
          paste("y", "~", paste(names(predictors.test)[c(-1,-2,-3,-4)], collapse=" + "), " + ", "(1|year) + (1|plot)", sep="")
        )
        m1 <- lmer(model.formula, data=predictors.test)
        P.values[k] <- max(summary(m1)$coefficients[-1,5])
        VIF.values[k] <- max(vif(m1))
        predictors.test.range[k] <- paste(names(predictors.test)[c(-1,-2,-3,-4)], collapse=",")
      }
    }
    P.values[which(VIF.values>5)] <- 1 # to avoid that collinear factors are added
    P.value <- min(P.values, na.rm=T)
    variable.include <- which(P.values == min(P.values, na.rm=T))
    predictors.include2 <- predictors.test.range[variable.include]
    predictors.include3 <- predictors.include2
    
    if(length(predictors)==0){
      predictors <- all.data[,predictors.include2]
      n.predictors.left <- n.predictors.left - 1
      if(P.value < 0.05){
        predictors.include.numbers <- which(names(all.data) %in% predictors.include2)
        all.data.2 <- all.data[, -predictors.include.numbers]
        predictors.include.numbers.all <- c(3,2,1,4,predictors.include.numbers)
        predictors.test <- all.data[complete.cases(all.data$y), predictors.include.numbers.all]
        model.formula <- as.formula(
          paste("y", "~", paste(names(predictors.test)[c(-1,-2,-3,-4)], collapse=" + "), " + ", "(1|year) + (1|plot)", sep="")
        )
        best.model <- lmer(model.formula, data=predictors.test)        
      }else{
        model.formula <- as.formula(
          paste("y", "~", 1, " + ", "(1|year) + (1|plot)", sep="")
        )
        best.model <- lmer(model.formula, data=predictors.test)
      }
    }else{
      predictors.include3 <- as.character(strsplit(predictors.include2, split=",", fixed=TRUE)[[1]])
      predictors <- all.data[,predictors.include3] 
      n.predictors.left <- n.predictors.left - 1
      if(FDR==0){
        alpha <- q.star
      }else{
        alpha <- (ncol(all.data[,c(5:ncol(all.data))]) - n.predictors.left) /
          ncol(all.data[,c(5:ncol(all.data))]) * q.star
      }
      if(P.value < alpha){
        predictors.include.numbers <- which(names(all.data) %in% predictors.include3)
        all.data.2 <- all.data[, -predictors.include.numbers]
        predictors.test <- all.data[complete.cases(all.data$y),c("y", "year", "plot", "block", predictors.include3)]
        model.formula <- as.formula(
          paste("y", "~", paste(names(predictors.test)[c(-1,-2,-3,-4)], collapse=" + "), " + ", "(1|year) + (1|plot)", sep="")
        )
        best.model <- lmer(model.formula, data=predictors.test)
        if(max(summary(best.model)$coefficients[-1,5])>0.05){stop("error: accepted model contains non-significant paths")}
      }
    }
    print(paste("i = ", i, "k = ", k, "count = ", count))
  }
  best.model2 <- lmer(model.formula, data=predictors.test)
  if(nrow(summary(best.model2)$coefficients)>1){
    if(length(which(summary(best.model2)$coefficients[-1,5] < 0.05))>0){
      model.formula <- as.formula(
        paste("y", "~", paste(names(which(summary(best.model2)$coefficients[-1,5]
                                          < 0.05)), collapse=" + "), " + ", "(1|year) + (1|plot)", sep="")
      )
      predictors.test <- predictors.test[,c(1:4,
                                            which(names(predictors.test) %in%
                                                    names(which(summary(best.model2)$coefficients[-1,5] < 0.05))))]
      best.model <- lmer(model.formula, data=predictors.test)
      significant.predictors[i, which(names(significant.predictors) %in% row.names(summary(best.model)$coefficients)[-1])] <- 1
      effect.sizes[i, which(names(significant.predictors) %in% row.names(summary(best.model)$coefficients)[-1])] <- std_beta(best.model)[-1,2]
    }else{
      model.formula <- as.formula(
        paste("y", "~", paste(names(which(summary(best.model2)$coefficients[-1,5] !=
                                            max(summary(best.model2)$coefficients[-1,5]))), collapse=" + "), " + ", "(1|year) + (1|plot)", sep="")
      )
      predictors.test <- predictors.test[,c(1:4,
                                            which(names(predictors.test) %in%
                                                    names(which(summary(best.model2)$coefficients[-1,5]!=
                                                                  max(summary(best.model2)$coefficients[-1,5])))))]
      best.model <- lmer(model.formula, data=predictors.test)
      significant.predictors[i, which(names(significant.predictors) %in% row.names(summary(best.model)$coefficients)[-1])] <- 1
      effect.sizes[i, which(names(significant.predictors) %in% row.names(summary(best.model)$coefficients)[-1])] <- std_beta(best.model)[-1,2]
    }
  }
  nr_predictors[i] <- length(row.names(summary(best.model)$coefficients)[-1]) 
  R2.model[i] <- r.squaredGLMM(best.model)[1]
  R2.model.cond[i] <- r.squaredGLMM(best.model)[2]
  function.names[i] <- names(functions4.st)[j]
  print(i)
}

write.table(significant.predictors, "significant.predictors.new.txt", sep="\t")

# identify significant predictors
if(analyse.richness==1){
  significant.predictors.2 <- cbind(significant.predictors[,1],
            significant.predictors[,c(2:42)] + significant.predictors[,c(43:83)])
  
  for(i in 1:nrow(significant.predictors.2)){
    significant.predictors.2[i,which(significant.predictors.2[i,]==2)] <- 1
  }
  
  names(significant.predictors.2)
  row.names(significant.predictors.2)
  order.functions <- c(23, 29, 4, 1, 8, 9, 2, 5, 6, 7, 32, 24, 22, 40, 30, 31,  
                       41, 39, 19, 3, 18, 42, 43, 13, 14, 25, 26, 35, 36, 27, 28,
                       20, 37, 38, 10, 11, 15, 12, 16, 17, 34, 33, 21)
  significant.predictors.2 <- significant.predictors.2[order.functions,]
  
  names(significant.predictors.2) <-
    c("sowndiv","shoot root ratio", "shoot:root N ratio", "plant height",
      "leaf biomass production rate", "total leaf area", "leaf area",
      "leaf thickness", "specific leaf area", "leaf specific density",
      "leaf area ratio", "leaf form coefficient", "leaf dry matter content",
      "leaf C content", "leaf N content", "conductance", "toughness",
      "stem diameter", "stem specific density", "erectedness",
      "biomass inflorescence", "inflorences per shoot", "duration flowering",
      "seeds projected area", "nr seedlings", "seed weight",
      "seed width length ratio", "seed dry matter content", "root area",
      "rooting depth", "root area distribution", "specific root area",
      "specific root length", "root tissue density", "root nitrogen uptake",
      "root CN ratio", "root P content", "root K content", "root S content",
      "root Ca content", "root Na content", "nutrient uptake efficiency")
  row.names(significant.predictors.2) <- function.names[order.functions]
  
  significant.predictors.2[significant.predictors.2==2] <- 1
}

# if a traits contributes both with FI and FD to a function, it's still only
# one trait that counts
significant.predictors.2 <- significant.predictors[,c(1:41)] + significant.predictors[,c(42:82)]
for(i in 1:nrow(significant.predictors.2)){
  significant.predictors.2[i,which(significant.predictors.2[i,]==2)] <- 1
}

names(significant.predictors.2)
row.names(significant.predictors.2)
order.functions <- c(22, 28, 4, 1, 8, 9, 2, 5, 6, 7, 31, 23, 21, 39, 29, 30,  
                     40, 38, 19, 3, 18, 41, 42, 13, 14, 24, 25, 34, 35, 26, 27, 
                     36, 37, 10, 11, 15, 12, 16, 17, 33, 32, 20)
significant.predictors.2 <- significant.predictors.2[order.functions,]
   
names(significant.predictors.2) <-
  c("shoot root ratio", "shoot:root N ratio", "plant height",
    "leaf biomass production rate", "total leaf area", "leaf area",
    "leaf thickness", "specific leaf area", "leaf specific density",
    "leaf area ratio", "leaf form coefficient", "leaf dry matter content",
    "leaf C content", "leaf N content", "conductance", "toughness",
    "stem diameter", "stem specific density", "erectedness",
    "biomass inflorescence", "inflorences per shoot", "duration flowering",
    "seeds projected area", "nr seedlings", "seed weight",
    "seed width length ratio", "seed dry matter content", "root area",
    "rooting depth", "root area distribution", "specific root area",
    "specific root length", "root tissue density", "root nitrogen uptake",
    "root CN ratio", "root P content", "root K content", "root S content",
    "root Ca content", "root Na content", "nutrient uptake efficiency")
row.names(significant.predictors.2) <- function.names[order.functions]

significant.predictors.2[significant.predictors.2==2] <- 1

write.table(significant.predictors.2, "significant.predictors.3.txt", sep="\t")
write.table(effect.sizes, "effect.sizes.original.txt", sep="\t")

effect.sizes.2 <- effect.sizes[,c(1:41)] 
for(i in 1:nrow(significant.predictors.2)){
  for(j in 1:ncol(significant.predictors.2)){
    candidates <- effect.sizes[i,c(j,j+41)]
    effect.sizes.2[i,j] <- sum(as.numeric(as.character(candidates)), na.rm=T)
    if(effect.sizes.2[i,j]==0){effect.sizes.2[i,j] <- NA}
  }
}

# process output on effect sizes, and identify whether FD or FI has strongest effect
significant.predictors.3 <- significant.predictors[,c(1:41)] +
  significant.predictors[,c(42:82)]
which.strongest <- significant.predictors[,c(1:41)]
for(i in 1:nrow(significant.predictors.2)){
  for(j in 1:ncol(significant.predictors.2)){
    if(significant.predictors.3[i,j]==2){
      candidates <- effect.sizes[i,c(j,j+41)]
      effect.sizes.2[i,j] <- candidates[which(abs(candidates)==
                             max(abs(candidates)))]
      candidates[which(is.na(candidates))] <- 0
      which.strongest[i,j] <- which(abs(candidates)==max(abs(candidates)))
    }else if(significant.predictors.3[i,j]==0){
      which.strongest[i,j] <- NA
    }else{
      candidates <- effect.sizes[i,c(j,j+41)]
      candidates[which(is.na(candidates))] <- 0
      which.strongest[i,j] <- which(abs(candidates)==max(abs(candidates)))
    }
  }
}

names(effect.sizes.2)
row.names(effect.sizes.2)
order.functions <- c(order.functions)

which.strongest <- which.strongest[order.functions,]

effect.sizes.2 <- effect.sizes.2[order.functions,]

names(effect.sizes.2) <- names(which.strongest) <-
  c("shoot root ratio", "shoot:root N ratio", "plant height",
    "leaf biomass production rate", "total leaf area", "leaf area",
    "leaf thickness", "specific leaf area", "leaf specific density",
    "leaf area ratio", "leaf form coefficient", "leaf dry matter content",
    "leaf C content", "leaf N content", "conductance", "toughness",
    "stem diameter", "stem specific density", "erectedness",
    "biomass inflorescence", "inflorences per shoot", "duration flowering",
    "seeds projected area", "nr seedlings", "seed weight",
    "seed width length ratio", "seed dry matter content", "root area",
    "rooting depth", "root area distribution", "specific root area",
    "specific root length", "root tissue density", "root nitrogen uptake",
    "root CN ratio", "root P content", "root K content", "root S content",
    "root Ca content", "root Na content", "nutrient uptake efficiency")
row.names(effect.sizes.2) <- function.names[c(order.functions)]
row.names(which.strongest) <- function.names[c(order.functions)]

write.table(effect.sizes.2, "effect.sizes.2.txt", sep="\t")
write.table(which.strongest, "which.strongest.txt", sep="\t")


names(significant.predictors.2) <-
  c("shoot root ratio", "shoot:root N ratio", "plant height",
    "leaf biomass production rate", "total leaf area", "leaf area",
    "leaf thickness", "specific leaf area", "leaf specific density",
    "leaf area ratio", "leaf form coefficient", "leaf dry matter content",
    "leaf C content", "leaf N content", "conductance", "toughness",
    "stem diameter", "stem specific density", "erectedness",
    "biomass inflorescence", "inflorences per shoot", "duration flowering",
    "seeds projected area", "nr seedlings", "seed weight",
    "seed width length ratio", "seed dry matter content", "root area",
    "rooting depth", "root area distribution", "specific root area",
    "specific root length", "root tissue density", "root nitrogen uptake",
    "root CN ratio", "root P content", "root K content", "root S content",
    "root Ca content", "root Na content", "nutrient uptake efficiency")

names(sort(colMeans(significant.predictors.2)))
write.table(names(sort(colMeans(significant.predictors.2))), "traits.ordered.txt", sep="\t")

library("colourlovers")
palette1 <- clpalette('cadetblue3')

# make barplot indicating how often traits affected a function
colour.codes <- rep("lightsteelblue",ncol(significant.predictors.2))
colour.codes[which(names(sort(100*(colSums(significant.predictors.2)))) %in%
            c("specific leaf area", "plant height", "leaf N content",
              "leaf dry matter content", "leaf area", "seed weight",
              "stem specific density", "specific root length"))] <- 'royalblue3'
barplot(sort(colSums(significant.predictors.2)), col=colour.codes, ylim=c(0,14))
box()

write.table(names(sort(100*(colSums(significant.predictors.2)))),
            "functions.sorted.csv")   

# make histogram with R2 values of models
a <- hist(R2.model, col="grey")
box()
mean(R2.model)
quantile(R2.model,0.5)

write.table(R2.model[order.functions], "R2.model.txt", sep="\t")
write.table(R2.model.cond[order.functions], "R2.model.cond.txt", sep="\t")

# make histogram with number of significant predictors per function
a <- hist(nr_predictors)
rowSums(significant.predictors.2)
a <- hist(rowSums(significant.predictors.2), col="grey", breaks=seq(-0.5,15.5,1))
box()

mean(rowSums(significant.predictors.2))
mean(colSums(significant.predictors.2))





#---------------------------------------------------------------------------------#
# Overlap analysis: which proportion of significant predictors do different       #
# functions have in common?                                                       #
#---------------------------------------------------------------------------------#

significant.predictors.2
overlap <- matrix(nrow=nrow(significant.predictors.2),
                  ncol=nrow(significant.predictors.2))
for(i in 1:nrow(significant.predictors.2)){
  for(j in 1:nrow(significant.predictors.2)){
    row.A <- significant.predictors.2[i,]
    row.B <- significant.predictors.2[j,]
    both <- length(which((row.A + row.B)==2))
    combined <- 0.5 * sum(row.A,row.B)
    overlap[i,j] <- both / combined
  }
  print(i)
}

overlap
overlap2 <- overlap

# replace values on diagonal by NA
for(i in 1:nrow(significant.predictors.2)){
  for(j in 1:nrow(significant.predictors.2)){
    function.names.general.A <- function.names[i]
    function.names.general.B <- function.names[j]
    if(function.names.general.A==function.names.general.B){
      overlap2[i,j] <- NA
    }
  }
}
overlap2
mean(overlap2, na.rm=T) # mean overlap value





#---------------------------------------------------------------------------------#
# Repeat Main analysis and overlap analysis, but this time with PCA variables     #
# instead of raw trait data.                                                      #
#---------------------------------------------------------------------------------#

predictors_raw
pca.analysis <- prcomp(predictors_raw[,-c(1,2,3,4)])
which(pca.analysis$sdev>1)
pca.analysis[1]
which(summary(pca.analysis)$importance[2,] > 1/82) # first 15 pc's explain more
# than 1/82 of variance
PCAs <- predict(pca.analysis)[,c(1:15)]
predictors_PCA <- data.frame(predictors_raw[,c(1,2,3,4)],PCAs)

significant.predictors <- as.data.frame(matrix(rep(0,(no.functions*
                          ncol(predictors_PCA[,-c(1,2,3,4)]))),
                          nrow=no.functions, ncol=
                          ncol(predictors_PCA[,-c(1,2,3,4)])))
names(significant.predictors) <- c("PC1","PC2","PC3","PC4","PC5","PC6",
                                   "PC7","PC8","PC9","PC10","PC11","PC12",
                                   "PC13","PC14","PC15")
effect.sizes <- as.data.frame(matrix(rep(NA,(no.functions*ncol(predictors_PCA[,
                -c(1,2,3,4)]))), nrow = no.functions, ncol=ncol(predictors_PCA[,
                -c(1,2,3,4)])))
names(effect.sizes) <- names(significant.predictors)

nr_predictors <- c()
R2.model <- c()
R2.model.cond <- c()
function.names <- c()
year.analyis <- c()
number.analysis <- 0

for(i in 1:(ncol(functions4.st)-4)){
  j <- i + 4
  response <- functions4.st[,c(j,1:4)]
  nrow(response)
  predictors_raw.2 <- aggregate(cbind(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,
                                      PC11,PC12,PC13,PC14,PC15)~
                                  year+plot, predictors_PCA, mean)
  nrow(predictors_raw.2)
  all.data <- merge(response, predictors_raw.2, by=c("plot","year"), all.y=T, all.x=T)
  nrow(all.data)
  
  P.value <- 0
  n.predictors.left <- ncol(all.data[,c(6:ncol(all.data))])
  predictors <- c()
  predictors.include2 <- predictors.include3 <- c()
  names(all.data)[3] <- "y"
  count <- 0
  
  while(P.value < 0.05 && n.predictors.left > 0){
    count <- count + 1
    P.values <- c()
    VIF.values <- c()
    predictors.test.range <- c()
    for(k in 5:(n.predictors.left+4)){
      if(length(predictors)==0){
        predictors.test <- all.data[complete.cases(all.data$y),c(-4)][c(1,2,3,4,k)]
        model.formula <- as.formula(
          paste("y", "~", paste(names(predictors.test)[c(-1,-2,-3,-4)], collapse=" + "), " + ", "(1|year) + (1|plot)", sep="")
        )
        m1 <- lmer(model.formula, data=predictors.test)
        P.values[k] <- summary(m1)$coefficients[-1,5]
        VIF.values[k] <- 0
        predictors.test.range[k] <- paste(names(predictors.test)[c(-1,-2,-3,-4)], sep=",")
      }else{
        all.data.3 <- cbind(all.data[,c("y", "year", "plot", "block", predictors.include3)], all.data.2[,k+1])
        names(all.data.3)[ncol(all.data.3)] <- names(all.data.2)[k+1]
        predictors.test <- all.data.3[complete.cases(all.data$y),]
        model.formula <- as.formula(
          paste("y", "~", paste(names(predictors.test)[c(-1,-2,-3,-4)], collapse=" + "), " + ", "(1|year) + (1|plot)", sep="")
        )
        m1 <- lmer(model.formula, data=predictors.test)
        P.values[k] <- max(summary(m1)$coefficients[-1,5])
        VIF.values[k] <- max(vif(m1))
        predictors.test.range[k] <- paste(names(predictors.test)[c(-1,-2,-3,-4)], collapse=",")
      }
    }
    P.values[which(VIF.values>5)] <- 1 # to avoid that collinear factors are added
    P.value <- min(P.values, na.rm=T)
    variable.include <- which(P.values == min(P.values, na.rm=T))
    predictors.include2 <- predictors.test.range[variable.include]
    predictors.include3 <- predictors.include2
    
    if(length(predictors)==0){
      predictors <- all.data[,predictors.include2]
      n.predictors.left <- n.predictors.left - 1
      if(P.value < 0.05){
        predictors.include.numbers <- which(names(all.data) %in% predictors.include2)
        all.data.2 <- all.data[, -predictors.include.numbers]
      }else{
        model.formula <- as.formula(
          paste("y", "~", paste(names(predictors.test)[c(-1,-2,-3,-4)], collapse=" + "), " + ", "(1|year) + (1|plot)", sep="")
        )
        best.model <- lmer(model.formula, data=predictors.test)
      }
    }else{
      predictors.include3 <- as.character(strsplit(predictors.include2, split=",", fixed=TRUE)[[1]])
      predictors <- all.data[,predictors.include3] 
      n.predictors.left <- n.predictors.left - 1
      if(P.value < 0.05){
        predictors.include.numbers <- which(names(all.data) %in% predictors.include3)
        all.data.2 <- all.data[, -predictors.include.numbers]
        predictors.test <- all.data[complete.cases(all.data$y),c("y", "year", "plot", "block", predictors.include3)]
        model.formula <- as.formula(
          paste("y", "~", paste(names(predictors.test)[c(-1,-2,-3,-4)], collapse=" + "), " + ", "(1|year) + (1|plot)", sep="")
        )
        best.model <- lmer(model.formula, data=predictors.test)
        if(max(summary(best.model)$coefficients[-1,5])>0.05){stop("errorr: accepted model contains non-significant paths")}
      }
    }
    print(paste("i = ", i, "k = ", k, "count = ", count))
  }
  best.model2 <- lmer(model.formula, data=predictors.test)
  if(nrow(summary(best.model2)$coefficients)>1 && nchar(model.formula)[3] > 24){
    if(length(which(summary(best.model2)$coefficients[-1,5] < 0.05))>0){
      model.formula <- as.formula(
        paste("y", "~", paste(names(which(summary(best.model2)$coefficients[-1,5]
                                          < 0.05)), collapse=" + "), " + ", "(1|year) + (1|plot)", sep="")
      )
      predictors.test <- predictors.test[,c(1:4,
                                            which(names(predictors.test) %in%
                                                    names(which(summary(best.model2)$coefficients[-1,5] < 0.05))))]
      best.model <- lmer(model.formula, data=predictors.test)
      significant.predictors[i, which(names(significant.predictors) %in% row.names(summary(best.model)$coefficients)[-1])] <- 1
      effect.sizes[i, which(names(significant.predictors) %in% row.names(summary(best.model)$coefficients)[-1])] <- std_beta(best.model)[-1,2]
    }else{
      model.formula <- as.formula(
        paste("y", "~", paste(names(which(summary(best.model2)$coefficients[-1,5] !=
                                            max(summary(best.model2)$coefficients[-1,5]))), collapse=" + "), " + ", "(1|year) + (1|plot)", sep="")
      )
      predictors.test <- predictors.test[,c(1:4,
                                            which(names(predictors.test) %in%
                                            names(which(summary(best.model2)$coefficients[-1,5]!=
                                            max(summary(best.model2)$coefficients[-1,5])))))]
      best.model <- lmer(model.formula, data=predictors.test)
      significant.predictors[i, which(names(significant.predictors) %in% row.names(summary(best.model)$coefficients)[-1])] <- 1
      if(nrow(summary(best.model)$coefficients)>1){
        effect.sizes[i, which(names(significant.predictors) %in% row.names(summary(best.model)$coefficients)[-1])] <- std_beta(best.model)[-1,2]
      }else{
        effect.sizes[i, ] <- 0
      }
    }
  }
  nr_predictors[i] <- length(row.names(summary(best.model)$coefficients)[-1]) 
  R2.model[i] <- r.squaredGLMM(best.model)[1]
  R2.model.cond[i] <- r.squaredGLMM(best.model)[2]
  function.names[i] <- names(functions4.st)[j]
  print(i)
}

overlap <- matrix(nrow=nrow(significant.predictors), ncol=nrow(significant.predictors))
for(i in 1:nrow(significant.predictors)){
  for(j in 1:nrow(significant.predictors)){
    row.A <- significant.predictors[i,]
    row.B <- significant.predictors[j,]
    both <- length(which((row.A + row.B)==2))
    combined <- 0.5 * sum(row.A,row.B)
    if(combined>0){
      overlap[i,j] <- both / combined
    }else{
      overlap[i,j] <- NA
    }
  }
  print(i)
}

overlap
overlap2 <- overlap

for(i in 1:nrow(significant.predictors)){
  for(j in 1:nrow(significant.predictors)){
    function.names.general.A <- function.names[i]
    function.names.general.B <- function.names[j]
    if(function.names.general.A==function.names.general.B){
      overlap2[i,j] <- NA
    }
  }
}
overlap2
mean(overlap2, na.rm=T)





#---------------------------------------------------------------------------------#
# How much does the R2 of final models decrease if only the 6 most frequently     #
# assessed traits are analysed?                                                   #
#---------------------------------------------------------------------------------#

traits.included <- 6
no.functions <- ncol(functions4.st) - 4
significant.predictors <- as.data.frame(matrix(rep(0,(no.functions*traits.included)), nrow = no.functions, ncol=traits.included))
names(significant.predictors) <- paste("var", c(1:traits.included), sep="")
nr_predictors <- c()
R2.model <- c()
function.names <- c()
year.analyis <- c()
number.analysis <- 0

predictors_vector <- ncol(predictors_raw.2[,-c(1,2,3,4)])

# Big loop to run all models
# loop for single functions to select the most parsimonious model
for(i in 1:(ncol(functions4.st)-4)){
  j <- i + 4
  response <- functions4.st[,c(j,1:4)]
  nrow(response)
  predictors_raw.2 <- aggregate(cbind(
    FD_srr, FD_Nacr, FD_h1, FD_lbpr, FD_la, FD_lSize, FD_lt, FD_sla, FD_lsd,
    FD_lar, FD_fc, FD_ldmc, FD_leafC, FD_leafN, FD_cond, FD_pp, FD_sd, FD_ssd,
    FD_erect, FD_BM_Target_Inflorescence, FD_NrInflorescence_Target_Shoot, 
    FD_DurationFlowreing_Target_Phenology, FD_SeedsProjArea, FD_NrSeedling_Target,
    FD_SeedsTCW, FD_SeedsWLRatio, FD_sdmc, FD_ra, FD_RootingDepth_Target, 
    FD_k.rad, FD_sra, FD_srl, FD_rtd, FD_RNU, FD_RCN, FD_rootP, FD_rootK,
    FD_rootS, FD_rootCa, FD_rootNa, FD_nute, FI_srr, FI_Nacr, FI_h1, FI_lbpr,
    FI_la, FI_lSize, FI_lt, FI_sla, FI_lsd, FI_lar, FI_fc, FI_ldmc, FI_leafC,
    FI_leafN, FI_cond, FI_pp, FI_sd, FI_ssd, FI_erect, FI_BM_Target_Inflorescence,
    FI_NrInflorescence_Target_Shoot, FI_DurationFlowreing_Target_Phenology,
    FI_SeedsProjArea, FI_NrSeedling_Target, FI_SeedsTCW, FI_SeedsWLRatio,
    FI_sdmc, FI_ra, FI_RootingDepth_Target, FI_k.rad, FI_sra, FI_srl, FI_rtd,
    FI_RNU, FI_RCN, FI_rootP, FI_rootK, FI_rootS, FI_rootCa, FI_rootNa, FI_nute)~
    year+plot, predictors_raw, mean)
  most.common.traits <- c(5, 10, 20, 14, 16, 8, 41+c(5, 10, 20, 14, 16, 8))
  predictors_raw.2 <- predictors_raw.2[,c(1,2,most.common.traits)]
  names(predictors_raw.2)[c(3:ncol(predictors_raw.2))] <- paste("var", c(1:(2*traits.included)), sep="")
  nrow(predictors_raw.2)
  all.data <- merge(response, predictors_raw.2, by=c("plot","year"), all.y=T, all.x=T)
  nrow(all.data)
  
  P.value <- 0
  n.predictors.left <- ncol(all.data[,c(6:ncol(all.data))])
  predictors <- c()
  predictors.include2 <- predictors.include3 <- c()
  names(all.data)[3] <- "y"
  count <- 0
  
  while(P.value < 0.05 && n.predictors.left > 0){
    count <- count + 1
    P.values <- c()
    VIF.values <- c()
    predictors.test.range <- c()
    for(k in 5:(n.predictors.left+4)){
      if(length(predictors)==0){
        predictors.test <- all.data[complete.cases(all.data$y),c(-4)][c(1,2,3,4,k)]
        model.formula <- as.formula(
          paste("y", "~", paste(names(predictors.test)[c(-1,-2,-3,-4)],
          collapse=" + "), " + ", "(1|year) + (1|plot)", sep="")
        )
        m1 <- lmer(model.formula, data=predictors.test)
        P.values[k] <- summary(m1)$coefficients[-1,5]
        VIF.values[k] <- 0
        predictors.test.range[k] <- paste(names(predictors.test)[c(-1,-2,-3,-4)], sep=",")
      }else{
        all.data.3 <- cbind(all.data[,c("y", "year", "plot", "block",
                                        predictors.include3)], all.data.2[,k+1])
        names(all.data.3)[ncol(all.data.3)] <- names(all.data.2)[k+1]
        predictors.test <- all.data.3[complete.cases(all.data$y),]
        model.formula <- as.formula(
          paste("y", "~", paste(names(predictors.test)[c(-1,-2,-3,-4)],
          collapse=" + "), " + ", "(1|year) + (1|plot)", sep="")
        )
        m1 <- lmer(model.formula, data=predictors.test)
        P.values[k] <- max(summary(m1)$coefficients[-1,5])
        VIF.values[k] <- max(vif(m1))
        predictors.test.range[k] <- paste(names(predictors.test)[c(-1,-2,-3,-4)], collapse=",")
      }
    }
    P.values[which(VIF.values>5)] <- 1 # to avoid that collinear factors are added
    P.value <- min(P.values, na.rm=T)
    variable.include <- which(P.values == min(P.values, na.rm=T))
    predictors.include2 <- predictors.test.range[variable.include]
    predictors.include3 <- predictors.include2
    
    if(length(predictors)==0){
      predictors <- all.data[,predictors.include2]
      n.predictors.left <- n.predictors.left - 1
      if(P.value < 0.05){
        predictors.include.numbers <- which(names(all.data) %in% predictors.include2)
        all.data.2 <- all.data[, -predictors.include.numbers]
        predictors.include.numbers.all <- c(3,2,1,5,predictors.include.numbers)
        predictors.test <- all.data[complete.cases(all.data$y),predictors.include.numbers.all]
        model.formula <- as.formula(
          paste("y", "~", paste(names(predictors.test)[c(-1,-2,-3,-4)],
          collapse=" + "), " + ", "(1|year) + (1|plot)", sep="")
        )
        best.model <- lmer(model.formula, data=predictors.test)
      }else{
        model.formula <- as.formula(
          paste("y", "~", 1, " + ", "(1|year) + (1|plot)", sep="")
        )
        best.model <- lmer(model.formula, data=predictors.test)
      }
    }else{
      predictors.include3 <- as.character(strsplit(predictors.include2, split=",", fixed=TRUE)[[1]])
      predictors <- all.data[,predictors.include3] 
      n.predictors.left <- n.predictors.left - 1
      if(P.value < 0.05){
        predictors.include.numbers <- which(names(all.data) %in% predictors.include3)
        all.data.2 <- all.data[, -predictors.include.numbers]
        included.factors <- c(3,2,1,5,predictors.include.numbers)
        predictors.test <- subset(all.data, complete.cases(all.data$y), 
               select=included.factors)
        model.formula <- as.formula(
          paste("y", "~", paste(names(predictors.test)[c(-1,-2,-3,-4)],
          collapse=" + "), " + ", "(1|year) + (1|plot)", sep="")
        )
        best.model <- lmer(model.formula, data=predictors.test)
        if(max(summary(best.model)$coefficients[-1,5])>0.05){stop("errorr: accepted model contains non-significant paths")}
      }
    }
    print(paste("i = ", i, "k = ", k, "count = ", count))
  }
  significant.predictors[i, which(names(significant.predictors) %in% row.names(summary(best.model)$coefficients)[-1])] <- 1
  nr_predictors[i] <- length(row.names(summary(best.model)$coefficients)[-1]) 
  R2.model[i] <- r.squaredGLMM(best.model)[1]
  function.names[i] <- names(functions4.st)[j]
  print(i)
}

a <- hist(R2.model, col="grey")
mean(R2.model)
R2.model.common <- mean(R2.model)

# make histogram with R2 values
a <- hist(nr_predictors, col="grey", breaks=c(-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5))
box()
nr_predictors.common <- mean(nr_predictors)





#---------------------------------------------------------------------------------#
# Analyze species richness instead of traits                                      #
#---------------------------------------------------------------------------------#

R2.model.random <- c()
nr_predictors.random <- matrix(nrow=no.functions, ncol=n.randomizations)
names.predictors.random <- matrix(nrow=n.randomizations, ncol=traits.included)
R2.model.all <- matrix(nrow=no.functions, ncol=n.randomizations)
SR <- SR[SR>=1]
com2 <- cbind(com, SR)
names(com2)[1] <- "plot"
all.data4 <- merge(all.data, com2[,c(1,2,66)], by=c("plot", "year"), all.x=T)

R2.model <- c()
R2.model.random <- c()
for(i in 1:(ncol(functions4.st)-4)){
  j <- i + 4
  response <- functions4.st[,c(j,1:4)]
  nrow(response)
  all.data <- merge(response, all.data4[,c(1,2,18)],
                    by=c("plot","year"), all.y=T, all.x=T)
  nrow(all.data)
  
  names(all.data)[3] <- "y"
  m1 <- lmer(y ~ SR + (1|year) + (1|plot), data=all.data)
  summary(m1)
  P.value <- summary(m1)$coefficients[2,5]
  R2.model[i] <- r.squaredGLMM(m1)[1]
  R2.model.random[i] <- r.squaredGLMM(m1)[2]
}

mean(R2.model)

