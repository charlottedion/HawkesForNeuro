rm(list=ls())

#setwd("HawkesForNeuro")

# packages 
library(ggplot2)
library(gridExtra)
library(fda)
library(readr)
library(hawkes) 

source("Functions.R")


#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################

### Synthetic data (from Ptyhon)

#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################



# We imagine that the spiking of the neurons around the central neuron has been observed on [0,T] 
# with $T=100$. 
# We simulate now the trajectory of the membrane potential driven by the jumps on [0,6]

Delta = 0.0048 
grid <- seq(0, 6, by= Delta) 
M = 8 

# Spike trains of trial 1
trial1 <- read_csv("data/events1.csv", col_names = FALSE)

spikesneurons = as.list(1:M)
for (i in 1: M){
  vect <- trial1[i,which(trial1[i,]!=0)]
  ind = which((vect>0) & (vect <max(grid)))
  spikesneurons[[i]] = vect[ind]
}

jumptimes <- sort(unlist(spikesneurons))


# Simulation of the jump-diffusion 

isjump1 <- matrix(0,length(grid),length(jumptimes))
for (i in 1:length(grid)-1){
  for (j in 1:length(jumptimes)){
    isjump1[i,j] <- (grid[i] < jumptimes[j] & jumptimes[j] <= grid[i+1])
  }
}

isjumpN <- apply(isjump1,1,sum)

# simulation of the potential, definition of the coefficients 
bdrift = function(x){
  return(-20*x-1080)
}
sig = function(x){
  return(rep(11,length(x)))
}

ajump = function(x){
  return(-0.07*x-2)
}

set.seed(123)
X0 = runif(1,-55, -45)
X = simu_jumpdiff(X0=X0, grid= grid, bfunc= bdrift, sigfunc= sig, afunc = ajump, Jumps = jumptimes, isjumpN = isjumpN)

n = length(X)-2

npas = 200
q1 <- min(X)
q2 <- max(X) 
gridx = seq(q1, q2,length = npas) 


# Graph
op = par(mfrow = c(1, 2), mar = c(1.5, 1.5, 1, 1), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0),
         cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)
plot(grid, X,type='l', xlab="",ylab="" ,ylim = c(-70, -30))

plot(NA, NA, xlab='', ylab='', las=1, ylim=c(0, M+1), xlim=c(0, 5), axes = FALSE)
for (i in 1:M){
  points(spikesneurons[[i]],rep(i,length(spikesneurons[[i]])), pch='+', col=i)
}
axis(side=1, at=0:5)
axis(2, at = seq(1, M, by = 1), labels = 1:M)



#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################

### Estimation with jumps (Hawkes-diffusion model)

#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################



#########################################################################################################
# Conditional expectation estimation, function f
#########################################################################################################

adjM <- as.matrix(read.table('results/adj_true.txt', sep=''), M, M)
base <- t(as.matrix(read.table('results/baseline.txt', sep=''),  1, M))

paramhawkes = as.list(1:3)
paramhawkes[[1]]= base # vecteur xi taille M
paramhawkes[[2]]= adjM #matrice adj taille M
paramhawkes[[3]] = rep(5, M) #vecteur beta= decay size M, used for the simulation

intensity <- sapply(grid[-(n+1)], 
                    function(s) intensM(s, param = paramhawkes, times = spikesneurons))
#save(intensity, file = "intensitReduity.RData")
#intensM: list size M with each coordinate of the intensity 
#dim(intensity)
#plot(unlist(intensity[6,]), type='l')# example neuron 6 

###########################################################################
# Estimation of function $f$, needed for the estimation of $a$
# here the bandwidth is fixed bbut it can be chosen through cross-validation 
condiM = sapply(1:M, function(i){return(mNW(x = gridx, X = X[-length(X)], Y = unlist(intensity[i,]), h = 0.1, K = dnorm))})
#save(condiM, file = "condiMReduit.RData")

#########################################################################################################
# Estimation of g 
#########################################################################################################

kap <- 100
Nn <- 20

qq1 <- quantile(X, 0.05)
qq2 <- quantile(X, 0.95)
keep <- (gridx > qq1) & (gridx < qq2)
gridx2 <- gridx[keep]

Tquad <-  diff(X[2:length(X)])^2/Delta

collec <- collecestimcoeff(X = X, U = Tquad, q1 = q1, q2 = q2, Nn = Nn)
collecalphag <- collec$colleccoeffalpha
collecP <- collec$collecP
penaltyg<-penaltyg(Nn,n,Delta,kap)

collecestimg <- matrix(0, Nn, length(gridx))
for (k in 1:Nn){
  collecestimg[k, ] <- sapply(apply(projectionSm(gridx, q1, q2, k) * collecalphag[k,1:(2*k+1)], 2, sum), max, 0)
}
res_g <- adaptiveestim(colleccoeffalpha = collecalphag, collecmatP = collecP, U = Tquad, penalty=penaltyg)
mfinal_g = res_g$mhat
#paste('final selected dimension for g =', mfinal_g)


estimfinal_g <- collecestimg[mfinal_g, ]
#--- 
#par(mfrow=c(1,1))
#plot(gridx, estimfinal_g, xlab='', ylab='', xlim=c(qq1,qq2),col='red', type='l', main="Estimation of g", ylim=c(0, 500), las=1)

data_g = data.frame(x=gridx, y=estimfinal_g)
ggplot(data_g, aes(x=x, y=y))+ geom_line()+  ylim(0, 500) + xlim(qq1,qq2)+theme( panel.grid.minor = element_blank()) + theme_bw()

#########################################################################################################
# Estimation of sigma^2 
#########################################################################################################
kap = 100
Tquad <-  diff(X[2:length(X)])^2/Delta
seuil <- 1

#---
beta = 1/8
truncphi = sapply(diff(X[2:length(X)])/(seuil * Delta^beta), phifunc)
Tquadphi <-  Tquad  *truncphi
paste('nombre de valeurs seuillées:', sum(truncphi != 1))

collecsig <- collecestimcoeff(X = X, U = Tquadphi, q1 = q1, q2 = q2, Nn = Nn)
collecalphasig <- collecsig$colleccoeffalpha
collecPsig <- collecsig$collecP
penaltysig<-penaltysig(Nn,n,kap)

collecestimsig <- matrix(0, Nn, length(gridx))
for (k in 1:Nn){
  collecestimsig[k, ] <- sapply(apply(projectionSm(gridx, q1, q2, k) * collecalphasig[k,1:(2*k+1)], 2, sum), max, 0)
}
res_sig <- adaptiveestim(colleccoeffalpha = collecalphasig, collecmatP = collecPsig, U = Tquadphi, penalty=penaltysig)
mfinal_sig = res_sig$mhat
paste('final selected dimension for sig^2 =', mfinal_sig)

estimfinal_sig <- collecestimsig[mfinal_sig, ]

#---
#plot(gridx, estimfinal_sig,col='red', type='l', xlab='', ylab='', xlim=c(qq1,qq2),main="Estimation of sigma^2")

options(warn=-1)
data_sig = data.frame(x=gridx, y=estimfinal_sig)
ggplot(data_sig, aes(x=x, y=y))+ geom_line()+  xlim(qq1,qq2)+theme( panel.grid.minor = element_blank()) + theme_bw()



#########################################################################################################
#Estimation a
#########################################################################################################

sumcondiM = apply(condiM, 1, sum)

estim_a2 <- (estimfinal_g-estimfinal_sig)/sumcondiM
estim_a <- sqrt(estim_a2*(estim_a2>0))

#plot(gridx, estim_a, col = "red", type = 'l', xlab = '', ylab='', xlim = c(qq1,qq2), main="Estimation of a")

data_a = data.frame(x=gridx, y=estim_a)
ggplot(data_a, aes(x=x, y=y))+ geom_line()+  xlim(qq1,qq2)+theme(panel.grid.minor = element_blank()) + theme_bw()

#########################################################################################################
#  Drift estimation (b)
#########################################################################################################

sigmean = sqrt(mean(estimfinal_sig[keep]))

sigesti <- function(x){ return(0*x+ sigmean)}

gesti <- function(x){
  return(apply(projectionSm(x,q1,q2,mfinal_g) * collecalphag[mfinal_g,1:(2*mfinal_g+1)],2,sum))
}


# linear approximation 
coeff = lm(estim_a[keep]~ gridx[keep])$coefficients
aesti <- function(x){
  return(coeff[2]*x + coeff[1]) 
}

aestiX = aesti(X)


# Graph 
op = par(mfrow = c(1, 2), mar = c(2, 2, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0),
         cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)

plot(gridx2, collecestimsig[mfinal_sig, keep],col='red', type='l', xlab='', ylab='', xlim=c(qq1,qq2),main="estimation of sigma^2", ylim=c(000, 200))
lines(gridx2, sigesti(gridx2)^2, xlab='', ylab='',col='red', lty=2)

plot(gridx2, estim_a[keep], col = "red", type = 'l', xlab = '', ylab='', xlim = c(qq1,qq2),ylim=c(0, 4), main="estimation of a")
lines(gridx2, aesti(gridx2), xlab='', ylab='', col='red', lty=2)

#########################################################################################################

sigma02 <- max(estimfinal_sig[keep])

Y <-  diff(X[2:length(X)])/Delta
termT <- (1/Delta) * aestiX * isjumpN

U = Y - termT[2:(length(termT)-1)]

collecb = collecestimcoeff(X, U, q1, q2, Nn)
collecalphab = collecb$colleccoeffalpha

collecestimb = matrix(0,Nn,length(gridx))
for (k in 1:Nn){
  collecestimb[k, ] = apply(projectionSm(gridx,q1,q2,k) * collecalphab[k,1:(2*k+1)],2,sum)
}

rho=3
penaltyb <-sapply(1:Nn, penaltyb,n, Delta,rho, sigma02)
res_b = adaptiveestim(collecalphab, collecb$collecP, U, penalty=penaltyb)
mfinal_b = res_b$mhat
#mfinal_b 


estim_b = collecestimb[mfinal_b,]

coefflin_b = (lm(collecestimb[mfinal_b, ] ~ gridx))$coeff
#######

#par(mfrow=c(1,1))
#plot(gridx, estim_b, type='l', ylim=c(-1000, 1000), las=1,xlab='', ylab='')

# Graph
data_b = data.frame(x=gridx, y = estim_b, z = gridx*coefflin_b[2] + coefflin_b[1])
ggplot(data_b, aes(x=x, y=y))+ geom_line()+geom_line(aes(x=x, y=z), color='blue') + xlim(qq1,qq2)+ylim(-700, 700)+theme(panel.grid.minor = element_blank()) + theme_bw()

###### final graph ############################################################################

op = par(mfrow = c(1, 3), mar = c(2, 2, 2, 2), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0),
         cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)


plot(gridx2, collecestimsig[mfinal_sig, keep],col='red', type='l', xlab='', ylab='', xlim=c(qq1,qq2),main="estimation of sigma^2", ylim=c(000, 100))
lines(gridx2, sigesti(gridx2)^2, xlab='', ylab='',col='red', lty=2)

plot(gridx2, estim_a[keep], col = "red", type = 'l', xlab = '', ylab='', xlim = c(qq1,qq2),ylim=c(0, 6), main="estimation of a")
lines(gridx2, aesti(gridx2), xlab='', ylab='', col='red', lty=2)

plot(gridx2, collecestimb[mfinal_b,keep], col = "red", type = 'l', xlab = '', ylab='', xlim = c(qq1,qq2), main="estimation of b", ylim=c(-400, 400))
lines(gridx2, gridx2*coefflin_b[2] + coefflin_b[1], col='red', lty=2)


#################################################################################################
## Validation phase 
#################################################################################################


besti <- function(x){
  return(coefflin_b[2]*x+ coefflin_b[1])
}


########## 

jumpsimu = 0 
while(length(jumpsimu)<=1){
  jumpsimu <- simulateHawkes(lambda0= paramhawkes[[1]],alpha = paramhawkes[[2]] * paramhawkes[[3]][1],beta = paramhawkes[[3]],horizon= max(grid))
  jumpsimu = sort(unlist(jumpsimu))
}
isjumpSimu <- matrix(0,length(grid),length(jumpsimu))
for (k in 1:length(grid)-1){
  for (j in 1:length(jumpsimu)){
    isjumpSimu[k,j] <- (grid[k] < jumpsimu[j] & jumpsimu[j]<=grid[k+1])
  }
}
isjumpSimu <- apply(isjump1,1,sum)

Xvalid <- simu_jumpdiff(X0=runif(1, min=-55, max=-35), grid=grid, bfunc=besti, sigfunc=sigesti, afunc=aesti, Jumps=jumpsimu, isjumpN=isjumpSimu)

#--- example of simulated path
#par(mfrow=c(1,1))
#plot(NA,NA, type='l', las=1, xlab='', ylab='',xlim=c(0, max(grid)),ylim=c(-80, -20))
#lines(grid, Xvalid, type='l',col= 'blue')

data_Xvalid = data.frame(x=grid, y = Xvalid)
ggplot(data_Xvalid, aes(x=grid, y=Xvalid))+ geom_line(color='blue') + xlim(qq1,qq2)+ylim(-80, -20)+xlim(0, max(grid))+theme(panel.grid.minor = element_blank()) + theme_bw()


#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################

### Estimation without jumps

#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################


#########################################################################################################
# Estimation of sigma^2
#########################################################################################################

kap = 100
Tquad <-  diff(X[2:length(X)])^2/Delta

collecsig <- collecestimcoeff(X = X, U = Tquad, q1 = q1, q2 = q2, Nn = Nn)
collecalphasig <- collecsig$colleccoeffalpha
collecPsig <- collecsig$collecP
penaltysig<-penaltysig(Nn,n,kap)

collecestimsig <- matrix(0, Nn, length(gridx))
for (k in 1:Nn){
  collecestimsig[k, ] <- sapply(apply(projectionSm(gridx, q1, q2, k) * collecalphasig[k,1:(2*k+1)], 2, sum), max, 0)
}
res_sig <- adaptiveestim(colleccoeffalpha = collecalphasig, collecmatP = collecPsig, U = Tquad, penalty=penaltysig)
mfinal_sig = res_sig$mhat
paste('final selected dimension for sig^2 =', mfinal_sig)

estimfinal_sig_wtjump <- collecestimsig[mfinal_sig, ]

#---
plot(gridx, estimfinal_sig_wtjump,col='red', type='l', xlab='', ylab='', xlim=c(qq1,qq2),ylim=c(0, 400), las=1,main="estimation of sigma^2")

#########################################################################################################
#  Drift estimation (b)
#########################################################################################################

sigmean_wtjump = sqrt(mean(estimfinal_sig[keep]))

sigesti_wtjump <- function(x){ return(0*x+ sigmean_wtjump)}
sigma02 <- max(estimfinal_sig_wtjump[keep])

###############

Y <-  diff(X[2:length(X)])/Delta
collecb = collecestimcoeff(X, Y, q1, q2, Nn)
collecalphab = collecb$colleccoeffalpha

collecestimb = matrix(0,Nn,length(gridx))
for (k in 1:Nn){
  collecestimb[k, ] = apply(projectionSm(gridx,q1,q2,k) * collecalphab[k,1:(2*k+1)],2,sum)
}

rho=3
penaltyb <-sapply(1:Nn, penaltyb,n, Delta,rho, sigma02)
res_b = adaptiveestim(collecalphab, collecb$collecP, Y, penalty=penaltyb)
mfinal_b = res_b$mhat


# graph on [qq1, qq2]
ggplot() +aes(x=gridx2, y=collecestimb[mfinal_b, keep])+geom_line(col='blue', lwd=1)+ ylim(-500, 500)+ xlab('')+ylab('')+theme(panel.grid.minor = element_blank()) + theme_bw()



besti_wtjump <- function(x){
  return(apply(projectionSm( x, q1, q2, mfinal_b) * collecalphab[mfinal_b, 1:(2 * mfinal_b + 1)], 2, sum)) 
}


#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################

### Validation with Depth

#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################


library('curveDepth')

NrepMC = 50
Nsimu = 100
mrep = 100 #parameter of depthc.Tukey
depthXH = rep(0, NrepMC)
depthXwtH = rep(0, NrepMC)
rangXH = rep(0, NrepMC)
rangXwtH = rep(0, NrepMC)

for (nrepMC in 1:NrepMC){
  Xvalid <- matrix(0, Nsimu, n+2)
  Xvalid_wtjump <- matrix(0, Nsimu, n+2)
  for (i in 1:Nsimu){
    jumpsimu = 0 
    while(length(jumpsimu)<=1){
      jumpsimu <- simulateHawkes(lambda0= paramhawkes[[1]],alpha = paramhawkes[[2]] * paramhawkes[[3]][1],beta = paramhawkes[[3]],horizon= max(grid))
      jumpsimu = sort(unlist(jumpsimu))
    }
    isjumpSimu <- matrix(0,length(grid),length(jumpsimu))
    for (k in 1:length(grid)-1){
      for (j in 1:length(jumpsimu)){
        isjumpSimu[k,j] <- (grid[k] < jumpsimu[j] & jumpsimu[j]<=grid[k+1])
      }
    }
    isjumpSimu <- apply(isjump1,1,sum)
    
    Xvalid[i, ] <- simu_jumpdiff(X0=runif(1, min=-55, max=-35), grid=grid, bfunc=besti, sigfunc=sigesti, afunc=aesti, Jumps=jumpsimu, isjumpN=isjumpSimu)
    
    Xvalid_wtjump[i, ] <- simu_diff(X0=X0,grid=grid, bfunc=besti_wtjump, sigfunc= sigesti_wtjump)
  }
  
  Xvalidlist = as.list(Nsimu)
  for (i in 1:Nsimu){
    Xvalidlist[[i]] = as.list(1)
    Xvalidlist[[i]][[1]]= matrix(c(grid,Xvalid[i,]), ncol=2, nrow=length(grid), byrow=FALSE)
    names(Xvalidlist[[i]]) <- c('coords')
  }
  Xvalidlist_wtjump = as.list(Nsimu)
  for (i in 1:Nsimu){
    Xvalidlist_wtjump[[i]] = as.list(1)
    Xvalidlist_wtjump[[i]][[1]]= matrix(c(grid,Xvalid_wtjump[i,]), ncol=2, nrow=length(grid), byrow=FALSE)
    names(Xvalidlist_wtjump[[i]]) <- c('coords')
  }
  
  Xlist = as.list(1)
  Xlist[[1]] = as.list(1)
  Xlist[[1]][[1]] = matrix(c(grid,X), ncol=2, nrow=length(grid), byrow=FALSE)
  names(Xlist[[1]]) <- c('coords')
  
  depthXH[nrepMC] = depthc.Tukey(object = Xlist, data = Xvalidlist, m = mrep)#, nDirs = 100L, subs = TRUE, fracInt = 0.5,
  depthXwtH[nrepMC] = depthc.Tukey(object = Xlist, data = Xvalidlist_wtjump, m = mrep)#, nDirs = 100L, subs = TRUE, fracInt = 0.5,
  
  
  ################################
  # Rank for the first 'real' trajectory, in the simulated sample obtained with the estimated coefficients
  depthXHeachcurve = rep(0, Nsimu)
  for (i in 1:Nsimu){
    Xvalid_i = as.list(1)
    Xvalid_i[[1]] = as.list(1)
    Xvalid_i[[1]][[1]] = matrix(c(grid,Xvalid[i,]), ncol=2, nrow=length(grid), byrow=FALSE)
    names(Xvalid_i[[1]]) <- c('coords')
    
    depthXHeachcurve[i] = depthc.Tukey(object = Xvalid_i, data = Xvalidlist, m = mrep)#, nDirs = 100L, subs = TRUE, fracInt = 0.5,
  }
  
  depthXwtHeachcurve = rep(0, Nsimu)
  for (i in 1:Nsimu){
    Xvalidwtjump_i = as.list(1)
    Xvalidwtjump_i[[1]] = as.list(1)
    Xvalidwtjump_i[[1]][[1]] = matrix(c(grid,Xvalid_wtjump[i,]), ncol = 2, nrow = length(grid), byrow=FALSE)
    names(Xvalidwtjump_i[[1]]) <- c('coords')
    
    depthXwtHeachcurve[i] = depthc.Tukey(object=Xvalidwtjump_i, data = Xvalidlist_wtjump, m = mrep)#, nDirs = 100L, subs = TRUE, fracInt = 0.5,
  }
  
  rangXH[nrepMC]=which(sort(c(depthXHeachcurve,depthXH[nrepMC]))==depthXH[nrepMC])
  
  rangXwtH[nrepMC]=which(sort(c(depthXwtHeachcurve,depthXwtH[nrepMC]))==depthXwtH[nrepMC])
}

#mean(rangXH)
#mean(rangXwtH)

print(paste('Depth of the real trajectory in the bundle of diffusion paths:',mean(depthXwtH) , 
            'Depth of the real trajectory in the bundle of jump-diffusion paths:',mean(depthXH) ))


# Graph 
nsub = 5
sub = sample(1:Nsimu, size=nsub, replace=FALSE)

compar <- data.frame(grid=rep(grid, nsub*2),  paths= c(Xvalid[sub, ],Xvalid_wtjump[sub, ]),
                     model = c(rep("Jump-Diffusion", length(grid)*nsub), rep("Diffusion", length(grid)*nsub)))

plotfinal <- ggplot(compar, aes(x=grid, y=paths)) + xlab('')+ylab('')+
  geom_line(aes(color=model))+
  facet_wrap(~model , scales="free", ncol=2, nrow=1)+
  scale_color_manual(values=c('red', 'green'))+ylim(-60, -30)+
  theme_bw()+ theme(legend.position = "none")
plotfinal 




