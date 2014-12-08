rm(list = ls())

# Read parameters
setwd("/Users/DGravel/Documents/Projects_On_Going/Maple_migration/talk_Ouranos")
source("scripts/read_pars.R")
pars = as.list(read_pars("data/pars_v0.txt"))
#pars = as.list(read_pars("data/pars_v1.txt"))
#pars = as.list(read_pars("data/pars_all.txt"))
#pars = as.list(read_pars("data/GenSA_test_sub25000.txt"))

#################################
## SOLVE THE MODEL             ##
#################################

source("scripts/get_eq.R")

T = seq(-2,4,0.025)
P = seq(-2,2,0.025)
clim_space = expand.grid(T,P)

prob = matrix(nr = nrow(clim_space), nc = 4)
maxEig = numeric(nrow(clim_space))

# Equilibrium and eigen values
for(i in 1:nrow(clim_space)) {
  EQ = get_eq(ENV1 = clim_space[i,1],ENV2 = clim_space[i,2], pars)
  prob[i,] = EQ[1:4]
  maxEig[i] = max(EQ[5:7])
  cat(round(i/nrow(clim_space)*100,2),'\n')
}

write.table(STM2000,"data/SDM2080_clim_space_parsv0.txt")
write.table(maxEig,"data/maxEig_clim_space_parsv0.txt")

#write.table(STM2000,"data/SDM2080_clim_space_parsv1.txt")
#write.table(maxEig,"data/maxEig_clim_space_parsv1.txt")

#write.table(STM2000,"data/SDM2080_clim_space_parssub25000.txt")
#write.table(maxEig,"data/maxEig_clim_space_parssub25000.txt")

#################################
## PANEL A: CLIMATE SPACE      ##
#################################

# Draw states
draw = function(p) {
  states = c(1,2,3,0)
  draw = rmultinom(n=1,size=1,prob=p)
  states[which(draw==1)]
}

States = apply(prob,1,draw)

# Plot the results
Z = matrix(States,nr = length(T), nc = length(P))
#Z = matrix(coexist,nr = length(T), nc = length(P))
quartz(height = 6, width = 6)
layout(matrix(c(1,2),nr=2,nc=1,byrow=TRUE),heights = c(1,6))
par(mar=c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend("center",legend = c("Regeneration","Boreal","Temperate","Mixt"),fill = c("black","darkcyan","orange","palegreen3"),bty = "n",horiz = TRUE,cex = 1)
par(mar=c(5,5,0,2))
image(T*1.891697+1.529991 ,(P*132.858+1044.194),Z,xlab = "Annual mean temperature", ylab = "Precipitations (mm)", cex.lab = 1.5, cex.axis = 1.25, col = c("black","darkcyan","orange","palegreen3"))

dev.copy2pdf(file = "Figs/STM_clim_space.pdf")

#################################
## EIGEN VALUES
#################################

source("scripts/get_eq.R")

# Transform
logEig = log(-maxEig)

Z = matrix(logEig,nr = length(T), nc = length(P))
quartz(height = 6, width = 6)
layout(matrix(c(1,2),nr=2,nc=1,byrow=TRUE),heights = c(1,6))
par(mar=c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
par(mar=c(5,5,0,2))
image(T*1.891697+1.529991 ,(P*132.858+1044.194),Z,xlab = "Annual mean temperature", ylab = "Precipitations (mm)", cex.lab = 1.5, cex.axis = 1.25,col=rainbow(10000,end = 0.7))
#image(T*1.891697+1.529991 ,(P*132.858+1044.194),Z,xlab = "Annual mean temperature", ylab = "Precipitations (mm)", cex.lab = 1.5, cex.axis = 1.25,col=brewer.pal(11,"Spectral"))

dev.copy2pdf(file = "Figs/Eig_clim_space.pdf")




