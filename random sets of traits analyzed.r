#---------------------------------------------------------------------------------#
#                                                                                 #
# Script to analyse how many and which above- and belowground traits are needed   #
# to accurately predict levels of ecosystem functioning in the Jena main          #
# experiment.                                                                     #
#                                                                                 #
# Here, random sets of 2-39 traits are analyzed                                   #
#                                                                                 #
# Authors: Thomas Schroeder-Georgi and Fons van der Plas                          #
#                                                                                 #
#---------------------------------------------------------------------------------#

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

#---------------------------------------------------------------------------------#





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
block <- as.numeric(as.character(substr(predictors_raw$plot, 2, 2)))
predictors_raw <- data.frame(predictors_raw[,c(1,2,3)], block,
                             predictors_raw[,-c(1,2,3)])

# exclude redundant traits
names(predictors_raw)
predictors_raw <- predictors_raw[,-c(31, 73)] # leaf C:N

# count number of functions
no.functions <- ncol(functions4.st) - 4

# make empty dataframes to store output stats
significant.predictors <- as.data.frame(matrix(rep(0,(no.functions*
                          ncol(predictors_raw[,-c(1,2,3,4)]))), nrow=no.functions,
                          ncol=ncol(predictors_raw[,-c(1,2,3,4)])))
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
effect.sizes <- as.data.frame(matrix(rep(NA,(no.functions*ncol(
                              predictors_raw[,-c(1,2,3,4)]))),
                              nrow=no.functions,
                              ncol=ncol(predictors_raw[,-c(1,2,3,4)])))
names(effect.sizes) <- names(significant.predictors)

# analyse an x number of random traits
n.randomizations <- 100
traits.included <- 2
R2.model.random <- c()
nr_predictors.random <- matrix(nrow=no.functions, ncol=n.randomizations)
names.predictors.random <- matrix(nrow=n.randomizations, ncol=traits.included)
R2.model.all <- matrix(nrow=no.functions, ncol=n.randomizations)

predictors_raw.2 <- aggregate(cbind(FD_srr, FD_Nacr, FD_h1, FD_lbpr, FD_la,
                              FD_lSize, FD_lt, FD_sla, FD_lsd, FD_lar, FD_fc,
                              FD_ldmc, FD_leafC, FD_leafN, FD_rootP, FD_rootK,
                              FD_rootS, FD_rootCa, FD_rootNa, FD_cond, FD_pp,
                              FD_sd, FD_ssd, FD_erect, FD_BM_Target_Inflorescence,
                              FD_NrInflorescence_Target_Shoot, 
                              FD_DurationFlowreing_Target_Phenology,
                              FD_SeedsProjArea, FD_NrSeedling_Target, FD_SeedsTCW,
                              FD_SeedsWLRatio, FD_sdmc, FD_ra,
                              FD_RootingDepth_Target, FD_k.rad, FD_sra, FD_srl,
                              FD_rtd, FD_RNU, FD_RCN, FD_nute, FI_srr, FI_Nacr,
                              FI_h1, FI_lbpr, FI_la, FI_lSize, FI_lt, FI_sla,
                              FI_lsd, FI_lar, FI_fc, FI_ldmc, FI_leafC, FI_leafN,
                              FI_rootP, FI_rootK, FI_rootS, FI_rootCa, FD_rootNa,
                              FI_cond, FI_pp, FI_sd, FI_ssd, FI_erect,
                              FI_BM_Target_Inflorescence,
                              FI_NrInflorescence_Target_Shoot, 
                              FI_DurationFlowreing_Target_Phenology,
                              FI_SeedsProjArea, FI_NrSeedling_Target,
                              FI_SeedsTCW, FI_SeedsWLRatio, FI_sdmc, FI_ra,
                              FI_RootingDepth_Target, FI_k.rad, FI_sra, FI_srl,
                              FI_rtd, FI_RNU, FI_RCN, FI_nute)~
                              year+plot, predictors_raw, mean)

for(n in 1:n.randomizations){
  no.functions <- ncol(functions4.st) - 4
  significant.predictors <- as.data.frame(matrix(rep(0,(no.functions*traits.included)), nrow = no.functions, ncol=traits.included))
  names(significant.predictors) <- paste("var", c(1:traits.included), sep="")
  nr_predictors <- c()
  R2.model <- c()
  function.names <- c()
  year.analyis <- c()
  number.analysis <- 0
  
  predictors_vector <- ncol(predictors_raw.2[,-c(1,2,3,4)])
  random.traits <- sample(c(5:45), traits.included, replace=F) # check!
  
  names.predictors.random[n,] <- c("FD_srr", "FD_Nacr", "FD_h1", "FD_lbpr",
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
                                   "FD_rootCa", "FD_rootNa", "FD_nute")[
                                   random.traits-4]
  
  # Big loop to run all models
  # loop for single functions to select the most parsimonious model
  for(i in 1:(ncol(functions4.st)-4)){
    j <- i + 4
    response <- functions4.st[,c(j,1:4)]
    nrow(response)
    predictors_raw.2 <- aggregate(cbind(FD_srr, FD_Nacr, FD_h1, FD_lbpr, FD_la,
                                  FD_lSize, FD_lt, FD_sla, FD_lsd, FD_lar, FD_fc,
                                  FD_ldmc, FD_leafC, FD_leafN, FD_cond, FD_pp,
                                  FD_sd, FD_ssd, FD_erect,
                                  FD_BM_Target_Inflorescence,
                                  FD_NrInflorescence_Target_Shoot, 
                                  FD_DurationFlowreing_Target_Phenology,
                                  FD_SeedsProjArea, FD_NrSeedling_Target,
                                  FD_SeedsTCW, FD_SeedsWLRatio, FD_sdmc, FD_ra,
                                  FD_RootingDepth_Target, FD_k.rad, FD_sra, FD_srl,
                                  FD_rtd, FD_RNU, FD_RCN, FD_rootP, FD_rootK,
                                  FD_rootS, FD_rootCa, FD_rootNa, FD_nute, FI_srr,
                                  FI_Nacr, FI_h1, FI_lbpr, FI_la, FI_lSize, FI_lt,
                                  FI_sla, FI_lsd, FI_lar, FI_fc, FI_ldmc, FI_leafC,
                                  FI_leafN, FI_cond, FI_pp, FI_sd, FI_ssd,
                                  FI_erect, FI_BM_Target_Inflorescence,
                                  FI_NrInflorescence_Target_Shoot, 
                                  FI_DurationFlowreing_Target_Phenology,
                                  FI_SeedsProjArea, FI_NrSeedling_Target,
                                  FI_SeedsTCW, FI_SeedsWLRatio, FI_sdmc, FI_ra,
                                  FI_RootingDepth_Target, FI_k.rad, FI_sra,
                                  FI_srl, FI_rtd, FI_RNU, FI_RCN, FI_rootP,
                                  FI_rootK, FI_rootS, FI_rootCa, FI_rootNa,
                                  FI_nute)~year+plot, predictors_raw, mean)
    most.common.traits <- c(random.traits-2, 41+random.traits-2)
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
          predictors.include.numbers.all <- c(3,2,1,5,predictors.include.numbers)
          predictors.test <- all.data[complete.cases(all.data$y), predictors.include.numbers.all]
          model.formula <- as.formula(
            paste("y", "~", paste(names(predictors.test)[c(-1,-2,-3,-4)], collapse=" + "), " + ", "(1|year) + (1|plot)", sep="")
          )
          best.model <- lmer(model.formula, data=predictors.test)
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
      print(paste("i = ", i, "n = ", n))
    }
    significant.predictors[i, which(names(significant.predictors) %in% row.names(summary(best.model)$coefficients)[-1])] <- 1
    nr_predictors[i] <- length(row.names(summary(best.model)$coefficients)[-1]) 
    R2.model[i] <- r.squaredGLMM(best.model)[1]
    function.names[i] <- names(functions4.st)[j]
    print(i)
  }
  R2.model.random[n] <- mean(R2.model)
  nr_predictors.random[,n] <- mean(nr_predictors)
  R2.model.all[,n] <- R2.model
}


filename <- paste("R2.model.random.test", traits.included, ".txt", sep="")
write.table(R2.model.random, filename, sep="\t")

random.output <- data.frame(R2.model.random, names.predictors.random)
filename <- paste("random.output.test", traits.included, ".txt", sep="")
write.table(random.output, filename, sep="\t")

