# Load libraries
library(randomForest)
library(nnet)

# Load data

# Load the model


# Draw states
draw = function(p) {
	states = c("B","T","M","R")
	draw = rmultinom(n=1,size=1,prob=p)
	states[which(draw==1)]
}

# Get predictions

#################
# Plot the TP map


# Compute the output of the SDM
set.seed(23)
T = seq(-4,10,0.01)
P = seq(-4,10,0.01)
newTP = expand.grid(T,P)
pred_TP = predict(SDM,new=data.frame(T=newTP$T, newTP$P),"prob")

# Get predicted states
p_states = 
pred_states = apply(p_states,2,draw) 

# Plot the results
Z = matrix(coexist,nr = length(X), nc = length(Y))
quartz(width = 6, height = 6)
layout(matrix(c(1,2),nr=2,nc=1,byrow=TRUE),heights = c(1,6))
par(mar=c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend("center",legend = c("Instable","Boréal","Tempéré","Mixte"),fill = c("black","darkcyan","orange","palegreen3"),bty = "n",horiz = TRUE,cex = 1)
par(mar=c(5,5,0,2))
image(T,P*1000,Z,xlab = "Température moyenne annuelle", ylab = "Précipitations annuelles (mm)", cex.lab = 1.5, cex.axis = 1.25, col = c("black","darkcyan","orange","palegreen3"))



#################
# Plot the Québec map




