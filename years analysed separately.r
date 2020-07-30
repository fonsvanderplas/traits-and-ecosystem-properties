#---------------------------------------------------------------------------------#
#                                                                                 # 
# Script to investigate how much variation in ecosystem functioning can be        #                                  
# explained by traits within years.                                               #
#                                                                                 #
# Authors: Thomas Schroeder-Georgi and Fons van der Plas                          #
#                                                                                 #
#---------------------------------------------------------------------------------#





### Switches ----------------------------------------------------------------------

analyse.richness <- 0 # if zero: don't include species richness in analysis, if 1: do include
FDR <- 0 # if zero: no false discover rate operformed. If 1: false discovery rate
abundances.included <- 1 # if zero: FD and CWMs based on presence/absence data. If 1: based on abundances



### Open data and R packages ------------------------------------------------------

# set working directory
setwd("C:/2019.01.21.full.BU/64 Schroeder FD CWM/R files/")


### load required packages
require(lm.beta)
library(MuMIn)
library(lme4)
library(lmerTest)
library(sjstats)
library(car)





#---------------------------------------------------------------------------------#
#                                     Load data                                   #
#---------------------------------------------------------------------------------#

traits.st <- read.table("traits.st.csv", header=T, sep=",")
row.names(traits.st) <- traits.st$sp
predictors_raw <- read.table("predictors_raw.csv", header=T, sep=",")
functions4.st <- read.table("functions4.st.csv", header=T, sep=",")

#---------------------------------------------------------------------------------#









#---------------------------------------------------------------------------------#
#                                Analyse all data                                 #
#---------------------------------------------------------------------------------#

# include block in predictors
block <- as.numeric(as.character(substr(predictors_raw$plot,2,2)))
predictors_raw <- data.frame(predictors_raw[,c(1,2,3)], block, predictors_raw[,-c(1,2,3)])

# exclude redundant traits
names(predictors_raw)
predictors_raw <- predictors_raw[,-c(31, 73)] # leaf C:N

no.functions <- ncol(functions4.st) - 4
significant.predictors <- as.data.frame(matrix(rep(0,(no.functions*ncol(predictors_raw[,-c(1,2,3,4)]))),
                                               nrow = no.functions, ncol=ncol(predictors_raw[,-c(1,2,3,4)])))
names(significant.predictors) <- c("FD_srr", "FD_Nacr", "FD_h1", "FD_lbpr",
                                   "FD_la", "FD_lSize", "FD_lt", "FD_sla",
                                   "FD_lsd", "FD_lar", "FD_fc", "FD_ldmc",
                                   "FD_leafC", "FD_leafN", "FD_cond", "FD_pp",
                                   "FD_sd", "FD_ssd", "FD_erect",
                                   "FD_BM_Target_Inflorescence",
                                   "FD_NrInflorescence_Target_Shoot", 
                                   "FD_DurationFlowreing_Target_Phenology",
                                   "FD_SeedsProjArea", "FD_NrSeedling_Target",
                                   "FD_SeedsTCW", "FD_SeedsWLRatio", "FD_sdmc",
                                   "FD_ra", "FD_RootingDepth_Target", "FD_k.rad",
                                   "FD_sra", "FD_srl", "FD_rtd", "FD_RNU",
                                   "FD_RCN", "FD_rootP", "FD_rootK", "FD_rootS",
                                   "FD_rootCa", "FD_rootNa", "FD_nute", "FI_srr",
                                   "FI_Nacr", "FI_h1", "FI_lbpr", "FI_la",
                                   "FI_lSize", "FI_lt", "FI_sla", "FI_lsd",
                                   "FI_lar", "FI_fc", "FI_ldmc", "FI_leafC",
                                   "FI_leafN", "FI_cond", "FI_pp", "FI_sd",
                                   "FI_ssd", "FI_erect",
                                   "FI_BM_Target_Inflorescence",
                                   "FI_NrInflorescence_Target_Shoot", 
                                   "FI_DurationFlowreing_Target_Phenology",
                                   "FI_SeedsProjArea", "FI_NrSeedling_Target",
                                   "FI_SeedsTCW", "FI_SeedsWLRatio", "FI_sdmc",
                                   "FI_ra", "FI_RootingDepth_Target", "FI_k.rad",
                                   "FI_sra", "FI_srl", "FI_rtd", "FI_RNU",
                                   "FI_RCN", "FI_rootP", "FI_rootK", "FI_rootS",
                                   "FDI_rootCa", "FI_rootNa", "FI_nute")
effect.sizes <- as.data.frame(matrix(rep(NA,(no.functions*ncol(predictors_raw[,-c(1,2,3,4)]))),
                                     nrow = no.functions, ncol=ncol(predictors_raw[,-c(1,2,3,4)])))
names(effect.sizes) <- names(significant.predictors)
nr_predictors <- c()
R2.model <- c()
R2.model.cond <- c()
function.names <- c()
year.analyis <- c()
number.analysis <- 0

predictors_vector <- ncol(predictors_raw[,-c(1,2,3,4)])

# Big loop to run all models
# loop for single functions to select the most parsimonious model
for(i in 1:(ncol(functions4.st)-4)){
  j <- i + 4
  response <- functions4.st[,c(j,1:4)]
  nrow(response)
  predictors_raw.2 <- aggregate(cbind(FD_srr, FD_Nacr, FD_h1, FD_lbpr, FD_la, FD_lSize,
                                      FD_lt, FD_sla, FD_lsd, FD_lar, FD_fc, FD_ldmc,
                                      FD_leafC, FD_leafN, FD_cond,
                                      FD_pp, FD_sd, FD_ssd, FD_erect, FD_BM_Target_Inflorescence,
                                      FD_NrInflorescence_Target_Shoot, 
                                      FD_DurationFlowreing_Target_Phenology, FD_SeedsProjArea,
                                      FD_NrSeedling_Target, FD_SeedsTCW, FD_SeedsWLRatio,
                                      FD_sdmc, FD_ra, FD_RootingDepth_Target, 
                                      FD_k.rad, FD_sra, FD_srl, FD_rtd, FD_RNU, FD_RCN, FD_rootP, FD_rootK, FD_rootS,
                                      FD_rootCa, FD_rootNa, FD_nute, 
                                      FI_srr, FI_Nacr, FI_h1, FI_lbpr, FI_la, FI_lSize,
                                      FI_lt, FI_sla, FI_lsd, FI_lar, FI_fc, FI_ldmc,
                                      FI_leafC, FI_leafN, FI_cond,
                                      FI_pp, FI_sd, FI_ssd, FI_erect, FI_BM_Target_Inflorescence,
                                      FI_NrInflorescence_Target_Shoot, 
                                      FI_DurationFlowreing_Target_Phenology, FI_SeedsProjArea,
                                      FI_NrSeedling_Target, FI_SeedsTCW, FI_SeedsWLRatio,
                                      FI_sdmc, FI_ra, FI_RootingDepth_Target, 
                                      FI_k.rad, FI_sra, FI_srl, FI_rtd, FI_RNU, FI_RCN, FI_rootP, FI_rootK, FI_rootS,
                                      FI_rootCa, FI_rootNa, FI_nute)~
                                      year+plot, predictors_raw, mean)
  nrow(predictors_raw.2)
  all.data <- merge(response, predictors_raw.2, by=c("plot","year"), all.y=T, all.x=T)
  nrow(all.data)
  
  all.data <- all.data[which(complete.cases(all.data[,3])),]
  nr.years <- length(unique(all.data$year))
  
  R2.model.year <- c()
  for(l in 1:nr.years){
    year.data <- all.data[all.data$year==unique(all.data$year)[l],]
    P.value <- 0
    n.predictors.left <- ncol(year.data[,c(6:ncol(year.data))])
    predictors <- c()
    predictors.include2 <- predictors.include3 <- c()
    names(year.data)[3] <- "y"
    count <- 0
    
    while(P.value < 0.05 && n.predictors.left > 0){
      count <- count + 1
      P.values <- c()
      VIF.values <- c()
      predictors.test.range <- c()
      for(k in 5:(n.predictors.left+4)){
        if(length(predictors)==0){
          predictors.test <- year.data[complete.cases(year.data$y),c(-4)][c(1,2,3,4,k)]
          model.formula <- as.formula(
            paste("y", "~", paste(names(predictors.test)[c(-1,-2,-3,-4)], collapse=" + "), sep="")
          )
          m1 <- lm(model.formula, data=predictors.test)
          P.values[k] <- summary(m1)$coefficients[-1,4]
          VIF.values[k] <- 0
          predictors.test.range[k] <- paste(names(predictors.test)[c(-1,-2,-3,-4)], sep=",")
        }else{
          all.data.3 <- cbind(year.data[,c("y", "year", "plot", "block", predictors.include3)], all.data.2[,k+1])
          names(all.data.3)[ncol(all.data.3)] <- names(all.data.2)[k+1]
          predictors.test <- all.data.3[complete.cases(year.data$y),]
          model.formula <- as.formula(
            paste("y", "~", paste(names(predictors.test)[c(-1,-2,-3,-4)], collapse=" + "), sep=""))
          m1 <- lm(model.formula, data=predictors.test)
          P.values[k] <- max(summary(m1)$coefficients[-1,4])
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
          predictors.include.numbers <- which(names(year.data) %in% predictors.include2)
          all.data.2 <- year.data[, -predictors.include.numbers]
          predictors.include.numbers.all <- c(3,2,1,5,predictors.include.numbers)
          predictors.test <- year.data[complete.cases(year.data$y), predictors.include.numbers.all]
          model.formula <- as.formula(
            paste("y", "~", paste(names(predictors.test)[c(-1,-2,-3,-4)], collapse=" + "), sep="")
          )
          best.model <- lm(model.formula, data=predictors.test)        
          
        }else{
          model.formula <- as.formula(
            paste("y", "~", 1, sep="")
          )
          best.model <- lm(model.formula, data=predictors.test)
        }
      }else{
        predictors.include3 <- as.character(strsplit(predictors.include2, split=",", fixed=TRUE)[[1]])
        predictors <- year.data[,predictors.include3] 
        n.predictors.left <- n.predictors.left - 1
        if(P.value < 0.05){
          predictors.include.numbers <- which(names(year.data) %in% predictors.include3)
          all.data.2 <- year.data[, -predictors.include.numbers]
          predictors.test <- year.data[complete.cases(year.data$y),c("y", "year", "plot", "block", predictors.include3)]
          model.formula <- as.formula(
            paste("y", "~", paste(names(predictors.test)[c(-1,-2,-3,-4)], collapse=" + "), sep="")
          )
          best.model <- lm(model.formula, data=predictors.test)
          if(max(summary(best.model)$coefficients[-1,4])>0.05){stop("errorr: accepted model contains non-significant paths")}
        }
      }
    }
  R2.model.year[l] <- summary(best.model)$adj.r.squared
  }
  R2.model[i] <- mean(R2.model.year)
  print(i)
}

mean(R2.model)


