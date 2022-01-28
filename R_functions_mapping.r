# ****************************************************************************************
# Trial name             : 
# Protocol No.           : 
# ----------------------------------------------------------------------------------------
# SAP No.                : 
# ----------------------------------------------------------------------------------------
# Title                  : R functions for Mapping EORTC QLQ-C30 onto EQ-5D-5L
# ----------------------------------------------------------------------------------------
# Program name           : R_functions_mapping
# Programmer / Location  : Yasuhiro Hagiwara 
# (Department of Biostatistics, University of Tokyo) 
# Date                   : 
# Pgm Modifier / Location: 
# Date / Comment         : 
# ----------------------------------------------------------------------------------------
# Attention:
# ****************************************************************************************


########################################################
# Packages
########################################################

# Before running this code, you should install the packages below.

library(tidyverse)
library(Matrix)
library(xgboost)


#######################################################
# Define functions for direct mapping
#######################################################

# Caution: Direct mapping is best suited to generating EQ-5D-5L index based on the Japanese value set


# Function for direct mapping based on an the gradient boosted tree approach
index.jp.gbt <- function(data) {
  # Data preparation
  data    <- data %>% select(ID, AGE, SEX, QLQPF, QLQRF, QLQEF, QLQCF, QLQSF, QLQGH,
                             QLQFA, QLQNV, QLQPA, QLQDY, QLQSL, QLQAP, QLQCO, QLQDI, QLQFI)
  modelmx <- sparse.model.matrix(ID ~ -1 + ., data)
  datamx  <- xgb.DMatrix(modelmx, label= data$ID)
  
  # Prediction of index
  load("GBTindex.rda")
  index.expected <- predict(gbt.train, datamx)
  
  return(index.expected)
}


# Function for direct mapping based on the two-part beta regression approach
index.jp.tpb <- function(data, seed) {
  # Calculation of parameters
  index.jp.max <- 1 - 0.060924 -0.043632
  index.jp.min <- -0.025449
  pr.full      <- plogis(-11.64060 + 0.03395*data$QLQPF + 0.03742*data$QLQRF + 0.02694*data$QLQEF + 0.03182*data$QLQGH - 0.06337*data$QLQPA)
  mu.notfull   <- plogis(-1.66963 + 0.15271*data$SEX + 0.02475*data$QLQPF + 0.00656*data$QLQRF + 0.00540*data$QLQEF + 0.01167*data$QLQGH)
  alpha        <- 6.05147 * mu.notfull
  beta         <- 6.05147 * (1 - mu.notfull)
  
  # Random sampling
  set.seed(seed)  
  full          <- rbinom(n = length(pr.full), size = 1, prob = pr.full) 
  index.notfull <- rbeta(n = length(mu.notfull), shape1 = alpha, shape2 = beta) 
  
  # Index (it is necessary to transform not-full-health part)
  index.expected <- 1 * pr.full + ((index.jp.max - index.jp.min) * mu.notfull     + index.jp.min) * (1 - pr.full)
  index.random   <- 1 * full    + ((index.jp.max - index.jp.min) * index.notfull  + index.jp.min) * (1 - full)
  
  return(list(index.expected = index.expected, index.random = index.random))
}


###########################################################
# Define functions for indirect mapping
###########################################################


# Function for indirect mapping based on an the gradient boosted tree approach.
index.response.gbt <- function(data,
                               tariff.mob,
                               tariff.sel,
                               tariff.usu,
                               tariff.pai,
                               tariff.anx,
                               tariff.notfull
                               ) {
  # Data preparation
  data    <- data %>% select(ID, AGE, SEX, QLQPF, QLQRF, QLQEF, QLQCF, QLQSF, QLQGH,
                             QLQFA, QLQNV, QLQPA, QLQDY, QLQSL, QLQAP, QLQCO, QLQDI, QLQFI)
  modelmx <- sparse.model.matrix(ID ~ -1 + ., data)
  datamx  <- xgb.DMatrix(modelmx, label= data$ID)
  
  # Calculation of Probabilities
  load("GBTresponse1.rda")
  load("GBTresponse2.rda")
  load("GBTresponse3.rda")
  load("GBTresponse4.rda")
  load("GBTresponse5.rda")
  load("GBTresponse6.rda")
  prob.mob <- matrix(predict(gbt.train1, datamx), ncol=5, byrow=T)
  prob.sel <- matrix(predict(gbt.train2, datamx), ncol=5, byrow=T)
  prob.usu <- matrix(predict(gbt.train3, datamx), ncol=5, byrow=T)
  prob.pai <- matrix(predict(gbt.train4, datamx), ncol=5, byrow=T)
  prob.anx <- matrix(predict(gbt.train5, datamx), ncol=5, byrow=T)
  prob.notfull <- predict(gbt.train6, datamx)
  
  # Expected Disutilities
  disutil.mob <- prob.mob %*% tariff.mob
  disutil.sel <- prob.sel %*% tariff.sel
  disutil.usu <- prob.usu %*% tariff.usu
  disutil.pai <- prob.pai %*% tariff.pai
  disutil.anx <- prob.anx %*% tariff.anx
  disutil.notfull <- prob.notfull * tariff.notfull

  # Index
  index.expected <- as.vector(1 + disutil.notfull + 
                                  disutil.mob + disutil.sel + disutil.usu + disutil.pai + disutil.anx)
  
  return(list(prob.mob = prob.mob, prob.sel = prob.sel, prob.usu = prob.usu, prob.pai = prob.pai, prob.anx = prob.anx,
              prob.notfull = prob.notfull, index.expected = index.expected))
}



# Function for indirect mapping based on the ordinal logistic regression approach
index.response.ol <- function(data, 
                              seed,
                              tariff.mob,
                              tariff.sel,
                              tariff.usu,
                              tariff.pai,
                              tariff.anx,
                              tariff.notfull
                              ) {
  # Mobility
  mob.lp    <- - 0.02469*data$AGE + 0.06934*data$QLQPF + 0.01498*data$QLQRF + 0.01767*data$QLQGH
  prob.mob1 <- plogis(-5.91277 + mob.lp)
  prob.mob2 <- plogis(-3.79148 + mob.lp) - prob.mob1
  prob.mob3 <- plogis(-1.85222 + mob.lp) - prob.mob1 - prob.mob2
  prob.mob4 <- plogis( 0.92106 + mob.lp) - prob.mob1 - prob.mob2 - prob.mob3
  prob.mob5 <- 1 - prob.mob1 - prob.mob2 - prob.mob3 - prob.mob4
  prob.mob  <- cbind(prob.mob1, prob.mob2, prob.mob3, prob.mob4, prob.mob5)
  
  # Self-care
  sel.lp    <- 0.60145*data$SEX + 0.06035*data$QLQPF + 0.01912*data$QLQRF
  prob.sel1 <- plogis(-4.24535 + sel.lp)
  prob.sel2 <- plogis(-2.56552 + sel.lp) - prob.sel1
  prob.sel3 <- plogis(-1.23110 + sel.lp) - prob.sel1 - prob.sel2
  prob.sel4 <- plogis( 0.17985 + sel.lp) - prob.sel1 - prob.sel2 - prob.sel3
  prob.sel5 <- 1 - prob.sel1 - prob.sel2 - prob.sel3 - prob.sel4
  prob.sel  <- cbind(prob.sel1, prob.sel2, prob.sel3, prob.sel4, prob.sel5)
  
  # Usual activities
  usu.lp    <- 0.33291*data$SEX + 0.04812*data$QLQPF + 0.03914*data$QLQRF + 0.00770*data$QLQSF + 0.02419*data$QLQGH
  prob.usu1 <- plogis(-9.28712 + usu.lp)
  prob.usu2 <- plogis(-6.21474 + usu.lp) - prob.usu1
  prob.usu3 <- plogis(-4.16801 + usu.lp) - prob.usu1 - prob.usu2
  prob.usu4 <- plogis(-1.46187 + usu.lp) - prob.usu1 - prob.usu2 - prob.usu3
  prob.usu5 <- 1 - prob.usu1 - prob.usu2 - prob.usu3 - prob.usu4
  prob.usu  <- cbind(prob.usu1, prob.usu2, prob.usu3, prob.usu4, prob.usu5)
  
  # Pain/discomfort
  pai.lp    <- 0.02343*data$QLQGH - 0.07683*data$QLQPA
  prob.pai1 <- plogis(-0.82400 + pai.lp)
  prob.pai2 <- plogis( 2.75434 + pai.lp) - prob.pai1
  prob.pai3 <- plogis( 5.35925 + pai.lp) - prob.pai1 - prob.pai2
  prob.pai4 <- plogis( 7.64850 + pai.lp) - prob.pai1 - prob.pai2 - prob.pai3
  prob.pai5 <- 1 - prob.pai1 - prob.pai2 - prob.pai3 - prob.pai4;
  prob.pai  <- cbind(prob.pai1, prob.pai2, prob.pai3, prob.pai4, prob.pai5)
  
  # Anxiety/Depression
  anx.lp    <- 0.02119*data$AGE + 0.06947*data$QLQEF
  prob.anx1 <- plogis(-6.65962 + anx.lp);
  prob.anx2 <- plogis(-4.31221 + anx.lp) - prob.anx1;
  prob.anx3 <- plogis(-2.81811 + anx.lp) - prob.anx1 - prob.anx2;
  prob.anx4 <- plogis(-0.86976 + anx.lp) - prob.anx1 - prob.anx2 - prob.anx3;
  prob.anx5 <- 1 - prob.anx1 - prob.anx2 - prob.anx3 - prob.anx4;
  prob.anx  <- cbind(prob.anx1, prob.anx2, prob.anx3, prob.anx4, prob.anx5)
  
  # Not full health
  prob.notfull <- 1 - (prob.mob[, 1] * prob.sel[, 1] * prob.usu[, 1] * prob.pai[, 1] * prob.anx[, 1])
  
  # Random sampling
  set.seed(seed)
  pred.mob <- matrix(nrow = nrow(prob.mob), ncol = 5)
  pred.sel <- matrix(nrow = nrow(prob.sel), ncol = 5)
  pred.usu <- matrix(nrow = nrow(prob.usu), ncol = 5)
  pred.pai <- matrix(nrow = nrow(prob.pai), ncol = 5)
  pred.anx <- matrix(nrow = nrow(prob.anx), ncol = 5)
  for (i in 1:nrow(prob.mob)) {
    pred.mob[i, ] <- t(rmultinom(n = 1, size = 1, prob = prob.mob[i, ]))
    pred.sel[i, ] <- t(rmultinom(n = 1, size = 1, prob = prob.sel[i, ]))
    pred.usu[i, ] <- t(rmultinom(n = 1, size = 1, prob = prob.usu[i, ]))
    pred.pai[i, ] <- t(rmultinom(n = 1, size = 1, prob = prob.pai[i, ]))
    pred.anx[i, ] <- t(rmultinom(n = 1, size = 1, prob = prob.anx[i, ]))
  }
  mult.pred     <- pred.mob * pred.sel * pred.usu * pred.pai * pred.anx
  pred.notfull <-  1 - pmin(1, mult.pred[, 1])
  
  # Expected disutilities
  exp.disutil.mob <- prob.mob %*% tariff.mob
  exp.disutil.sel <- prob.sel %*% tariff.sel
  exp.disutil.usu <- prob.usu %*% tariff.usu
  exp.disutil.pai <- prob.pai %*% tariff.pai
  exp.disutil.anx <- prob.anx %*% tariff.anx
  exp.disutil.notfull <- prob.notfull * tariff.notfull
  
  # Predicted disutilities
  pred.disutil.mob <- pred.mob %*% tariff.mob
  pred.disutil.sel <- pred.sel %*% tariff.sel
  pred.disutil.usu <- pred.usu %*% tariff.usu
  pred.disutil.pai <- pred.pai %*% tariff.pai
  pred.disutil.anx <- pred.anx %*% tariff.anx
  pred.disutil.notfull <- pred.notfull * tariff.notfull
  
  # Index
  index.expected <- as.vector(1 + exp.disutil.notfull + exp.disutil.mob + exp.disutil.sel +
                                  exp.disutil.usu     + exp.disutil.pai + exp.disutil.anx)
  index.random   <- as.vector(1 + pred.disutil.notfull + pred.disutil.mob + pred.disutil.sel +
                              pred.disutil.usu     + pred.disutil.pai + pred.disutil.anx)
  
  return(list(prob.mob = prob.mob, prob.sel = prob.sel, prob.usu = prob.usu, prob.pai = prob.pai, prob.anx = prob.anx,
              prob.notfull = prob.notfull, index.expected = index.expected, index.random = index.random))
}


#######################################################
# Instructions for preparing EORTC QLQ-C30 data
#######################################################

# 1: Please prepare your data in the same form as the sample data "EORTCQLQC30.csv".
#    Especially, you should use the same variable names (caution: R is case-sensitive).      
# 2: Our mapping algorithms did NOT learn missing values. 
#    You should use only complete records that have age, sex, and all 15 subscale scores in EORTC QLQ-C30.
# 3: Age should be on the year-scale and sex should be coded as 1 for male and 2 for female. 


############################################################
# Examples of mapping
############################################################

# Reading sample data
input <- read_csv("EORTCQLQC30.csv")

# Direct mapping
index.expected.gbt <- index.jp.gbt(input)
index.expected.tpb <- index.jp.tpb(input, 1234)$index.expected
index.random.tpb   <- index.jp.tpb(input, 1234)$index.random

# Example: Japanese value set
tariff.jp.mob <- c(0, -0.063865, -0.112618, -0.179043, -0.242916)
tariff.jp.sel <- c(0, -0.043632, -0.076660, -0.124265, -0.159659)
tariff.jp.usu <- c(0, -0.050407, -0.091131, -0.147929, -0.174786)
tariff.jp.pai <- c(0, -0.044545, -0.068178, -0.131436, -0.191203)
tariff.jp.anx <- c(0, -0.071779, -0.110496, -0.168171, -0.195961)
tariff.jp.notfull <- -0.060924

# Example: US value set
tariff.us.mob <- c(0, -0.096, -0.122, -0.237, -0.322)
tariff.us.sel <- c(0, -0.089, -0.107, -0.220, -0.261)
tariff.us.usu <- c(0, -0.068, -0.101, -0.255, -0.255)
tariff.us.pai <- c(0, -0.060, -0.098, -0.318, -0.414) 
tariff.us.anx <- c(0, -0.057, -0.123, -0.299, -0.321)
tariff.us.notfull <- 0

# Indirect mapping for the Japanese value set
index.expected.response.gbt.jp <- index.response.gbt(input, tariff.jp.mob, tariff.jp.sel, tariff.jp.usu, tariff.jp.pai, tariff.jp.anx, tariff.jp.notfull)$index.expected

index.expected.response.ol.jp  <- index.response.ol(input, 1234, tariff.jp.mob, tariff.jp.sel, tariff.jp.usu, tariff.jp.pai, tariff.jp.anx, tariff.jp.notfull)$index.expected
index.random.response.ol.jp    <- index.response.ol(input, 1234, tariff.jp.mob, tariff.jp.sel, tariff.jp.usu, tariff.jp.pai, tariff.jp.anx, tariff.jp.notfull)$index.random

# Indirect mapping for the US value set
index.expected.response.gbt.us <- index.response.gbt(input, tariff.us.mob, tariff.us.sel, tariff.us.usu, tariff.us.pai, tariff.us.anx, tariff.us.notfull)$index.expected

index.expected.response.ol.us  <- index.response.ol(input, 1234, tariff.us.mob, tariff.us.sel, tariff.us.usu, tariff.us.pai, tariff.us.anx, tariff.us.notfull)$index.expected
index.random.response.ol.us    <- index.response.ol(input, 1234, tariff.us.mob, tariff.us.sel, tariff.us.usu, tariff.us.pai, tariff.us.anx, tariff.us.notfull)$index.random

# EQ-5D-5L index
index.expected.gbt 
index.expected.tpb 
index.random.tpb   

index.expected.response.gbt.jp 
index.expected.response.ol.jp  
index.random.response.ol.jp    

index.expected.response.gbt.us 
index.expected.response.ol.us  
index.random.response.ol.us    

#Probabilities of 5 levels of, for example, mobility
index.response.gbt(input, tariff.jp.mob, tariff.jp.sel, tariff.jp.usu, tariff.jp.pai, tariff.jp.anx, tariff.jp.notfull)$prob.mob
index.response.ol(input, 1234, tariff.jp.mob, tariff.jp.sel, tariff.jp.usu, tariff.jp.pai, tariff.jp.anx, tariff.jp.notfull)$prob.mob

