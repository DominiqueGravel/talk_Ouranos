# ----------------------
# load calibration data 
# ----------------------

data = read.csv("fit_model/transitionsFourState.csv")
#data = read.csv("../STModel-Data/out_files/statesFourState.csv")
head(data)
dim(data)

# ----------------------
### choice of variables
# ----------------------

selectedVars = c("annual_mean_temp2","annual_pp2")

datSel = data[,c("state2",selectedVars)]

rm(data)

# ----------------------
# Clean data
# ---------------------

# Clean Undefined state
str(datSel)
datSel_wo_U <- subset(datSel, state2 != "U")
datSel_wo_U$state <- droplevels(datSel_wo_U$state2)

# ----------------------
# models 
# ----------------------
# evaluation statistics
# HK <- function (Pred, Obs) 
# {

# 	Misc = table(Pred, Obs)
	
#     if (nrow(Misc)!=ncol(Misc)) stop("wrong misclassification table")
#     Misc <- unclass(Misc)
#     k  <- ncol(Misc)
#     Nobs <- apply(Misc, 2, sum)
#     Npred <- apply(Misc, 1, sum)
#     N <- sum(Nobs)
  

#    HK <- (sum(diag(Misc))/N - sum(as.numeric(Nobs)*as.numeric(Npred))/N/N ) / ( 1 - sum(as.numeric(Nobs)*as.numeric(Nobs))/N/N )

#     return(HK)
# }

# cross validation - separation of the dataset
#sampl = sample(1:nrow(datSel_wo_U), 2*nrow(datSel_wo_U)/3)
#calib = datSel_wo_U[sampl,c("state2",selectedVars)] 
#valid = datSel_wo_U[-sampl,c("state2",selectedVars)]

# Run the models

# random Forest
# calib
#library(randomForest)
#rs = runif(1,0,1)
#set.seed(rs)
#SDM2 = randomForest(state2 ~ . , data = calib, ntree = 500)
#SDM2
#save(SDM2,rs,sampl,file= "RandomForest_temp.rObj")

# valid
#set.seed(rs)
#pred2 = predict(SDM2,new=datSel_wo_U,"response", OOB=TRUE)
#(HK2 = HK(pred2[-sampl], valid$state2)) 

# multimodal
#calib
library(nnet)
SDM1 = multinom(state2 ~ .^2 + I(annual_mean_temp2^2) + I(annual_mean_temp2^3) + I(annual_pp2^2) + I(annual_pp2^3), data = datSel_wo_U[,c("state2",selectedVars)], maxit =300)

#summary(SDM1)
#save(SDM1,file= "Multinom_temp.rObj")

#valid
#pred1 = predict(SDM1, new=valid,"class")
#(HK1 = HK(pred1, valid$state)) 

# ----------------------
# projection 
# ----------------------
## ----recap data
#load("../data/Multinom_temp.rObj")
#load("RandomForest_temp.rObj")
#selectedVars = c("annual_mean_temp2",  "annual_pp2")
#---------------

#dataProj = read.csv("transitionsFourState.csv")
#head(dataProj)

# means between climates of year of state 0 and year of state 1
dataProj$annual_mean_temp2 = apply(dataProj[,c("annual_mean_temp1", "annual_mean_temp2")], 1, mean)
dataProj$annual_pp2 = apply(dataProj[,c("annual_pp1", "annual_pp2")], 1, mean)

datProjSel = dataProj[,selectedVars]

# projection
set.seed(rs)
projProba = predict(SDM2,new=datProjSel,"prob", OOB=TRUE)

head(projProba)

# sauvegarde
write.table(projProba, file = "projection_neigbor_rf_temp.txt", quote=F, row.names=FALSE)






