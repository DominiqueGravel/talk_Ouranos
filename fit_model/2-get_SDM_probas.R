# ----------------------
# projection 
# ----------------------
## ----recap data
#load("../data/Multinom_temp.rObj")
load("RandomForest_temp.rObj")
selectedVars = c("annual_mean_temp2",  "annual_pp2")
#---------------

dataProj = read.csv("transitionsFourState.csv")
head(dataProj)

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




T = seq(-2,5,0.01)
P = numeric(length(T)) + 1000
proj = data.frame(annual_mean_temp2=T,annual_pp2=P)
projProba = predict(SDM2,new=proj,"prob", OOB=TRUE)