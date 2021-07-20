#Hawkes : intensity functions for M Hawkes processes

intensM <- function(t, param, times) { 
  # Intensity function with phi = exponential decay
  # times: list of size M with the spike trains of the M neurons
  # param: list of size 3 with xi vector of size M, alpha matrix of size MxM, beta vector of size M
  M = length(param[[1]])
  intenst = as.list(1:M)
  
  for (i in 1:M){
    if (length(times[[i]])==0) {intenst[[i]] = param[[1]][i]}    
    if (length(times[[i]])>0){
      to_sumi <- times[[i]][times[[i]] < t]
      intenst[[i]] = param[[1]][i] + sum(param[[2]][i,i] * param[[3]][i] * exp(- param[[3]][i] * (t - to_sumi)))
      
      for (j in (1:M)[-i]){
        to_sumj <- times[[j]][times[[j]] < t]
        intenst[[i]] = intenst[[i]] + sum(param[[2]][i,j] * param[[3]][i] * exp(- param[[3]][i] * (t - to_sumj)))
      }
    }
  }
  return(intenst)
}



### sum of the intensities 
intens <- function(k, M, s,times, param,xi){
  # s: time at which the intensity is calculated
  # times: jump times, list of size M
  # param: list of size 3, with the parameters
  # k: number of the neuron up to which we sum up 
  # M: total number of neurons
  # xi: spontaneous firing rate 
  vect <- 0
  
  for (m in 1:k){
    for (l in 1:M){
      if(sum(which(times[[m]]<s)) ==0){ vect <- vect}
      if(sum(which(times[[m]]<s)) != 0){
        ind = which(times[[m]] < s)
        kern = function(t){return(param[[2]][m,l] * param[[3]][m] * exp( param[[3]][m] *(t)))}
        vect <- vect + sum(sapply((s - times[[m]])[ind], kern))
        #sum(sapply((s - times[[m]])[ind], functh, l, m, modelInt))
      }
    }
  }
  vect <- k * xi + vect 
  return(vect)
}
# it is the sum of the intensities up to k<M with the same xi



# HAWKES : simulation of the jump times from an exponential Hawkes process 
simuHawkesExpoM <- function(param, M){
  # param: list of size 3 
  # first list for xi, second list for adj, third list for decay 
  times <- as.list(rep(0, M))
  s <- 0 
  
  while (s < Tend){ 
    lambda_bar <- intens(M, M, xi=xi,s=s, times=times, modelInt)
    
    u <- runif(1, 0, 1)
    w <- -log(u)/lambda_bar # exponential, the interarrival to the next candidate point
    s <- s + w
    
    #if (s <= Tend){
    D <- runif(1, 0, 1)
    
    if( D * lambda_bar <= intens(M, M, xi=xi,s=s, times=times, modelInt) ){
      k <- 1
      if(D * lambda_bar > intens(k, M, xi=xi, s=s,times=times, modelInt)){
        k <- k+1; times[[k]] <- c(times[[k]], s)
      }
      else {times[[k]] <- c(times[[k]],s)}
    }
    #}
  }
  if (s > Tend) {
    if (k == 1) {times[[k]] <- times[[k]][1:(length(times[[k]])-1)]}
    else {times[[2]] <- times[[2]][1:(length(times[[2]])-1)]}
  }
  
  
  times[[1]] <- times[[1]][-1] # jump time of the 1st neuron, we remove the 1st time which is equal to 0
  times[[2]] <- times[[2]][-1] # jump time of the 2nd neuron
  return(times)
}


##########
# Choice for function a
a <- function(x,namea){
  if( namea == 'none'){
    return(rep(0, length(x)))
  }
  if( namea == 'constant'){
    return(rep(1,length(x)))
  }
  if (namea == 'lin'){
    test= function(x){if((x<=5)&(x>=-5)){return(x)}; if(x>5){return(5)}; if(x< -5){return(-5)}}
    return( sapply(x, test))
  }
  if (namea == 'lin2'){
    return(-0.1*x)
  }
  if (namea == 'lin3'){
    theta= mean(x)
    return(theta-x)
  }
}

##############################################################################
#############################################################################
# Estimation
#############################################################################
#############################################################################

# Projection on the trigonometric basis
projectionSm=function(x,q1,q2,m){
  Dm = 2*m+1
  c = sqrt(2)/sqrt(q2-q1)
  
  dimgrid=length(x)
  proj = matrix(0,Dm,dimgrid)
  proj[1,] = rep(1/sqrt(q2-q1),dimgrid)
  
  for (l in 1:m){
    proj[2*l, ] = c*cos(2*pi*l*(x-q1)/(q2-q1))
    proj[2*l+1,] = c*sin(2*pi*l*(x-q1)/(q2-q1))
  }
  return(proj)
}


##############################################################
# Estimations of the parameters
##############################################################

# estimators for the coefficients in the projection
alphachapeau = function(P,U){
  N = length(U)
  AA = P%*%t(P)
  B = P%*%matrix(U,N,1,byrow=FALSE)
  
  alpha = rep(0,dim(P)[1]) #of dimension Dm
  alpha = solve(AA)%*%B
  return(alpha)
}


collecestimcoeff = function(X, U, q1, q2, Nn){ 
  
  colleccoeffalpha = matrix(0, Nn, 2 * Nn + 1)
  N = length(U)
  collecP = list()
  
  
  for(k in 1:Nn){
    collecP[[k]] = projectionSm(X[2:(N+1)], q1, q2, k)
    M <- NA
    try(M <- alphachapeau(collecP[[k]], U),silent=TRUE)
    if(sum(is.na(M))==0){colleccoeffalpha[k,1:(2*k+1)] = M}
    else break
  }
  return(list(collecP = collecP, colleccoeffalpha = colleccoeffalpha))
}


##############################################################################
#############################################################################
# Adaptation 
#############################################################################
#############################################################################


# penalties for the different parameters
penaltyb <- function(m,n,Delta,rho, sigma02){
  penb <- rho * (2 * m + 1) * sigma02/(n*Delta)
  return(penb)
}

penaltyg<-function(Nn,n,Delta,kap){
  peng<-kap * (1:Nn)/(n*Delta)
  return(peng)
}

penaltysig<-function(Nn,n,kap){
  pensig<-kap * (1:Nn)/n
  return(pensig)
}

# Adaptation procedure 
adaptiveestim <- function(colleccoeffalpha,collecmatP,U,penalty){ 
  
  ind = length(collecmatP)  
  
  estimmhat <- as.list(1:ind)
  criteremhat <- rep(0,ind)
  
  for (l in 1:ind){
    estimmhat[[l]] <- apply(collecmatP[[l]]*colleccoeffalpha[l,1:(2*l+1)],2,sum)  
    criteremhat[l] <- mean((U-estimmhat[[l]])^2)   #calcul de tous les gamma_n(sig^2mhat)
  }
  
  pen <- penalty
  
  crit <- criteremhat + pen[1:ind]
  mhat <- which.min(crit)
  
  return(list(estimmhat=estimmhat, criteremhat=criteremhat, crit=crit, mhat=mhat))
}



####################################################################
################# Supplementary functions 
####################################################################
# Truncation function 
phifunc <- function(x){
  if(abs(x)<1){phi = 1}
  if(abs(x) >=2){phi = 0}
  if((abs(x)>=1)&(abs(x)<2)){ phi = exp((1/3)+1/(x^2-4))}
  return(phi)
}


# Nadaraya-Watson estimator
mNW <- function(x, X, Y, h, K = dnorm) {
  # Arguments
  # x: evaluation points
  # X: vector (size n) with the predictors
  # Y: vector (size n) with the response variable
  # h: bandwidth
  # K: kernel
  if (length(x) ==1){
    # Matrix of size n x length(x)
    Kx <- sapply(X, function(Xi) K((x - Xi) / h) / h)
    
    # Weights
    W <- Kx / sum(Kx) # Column recycling!
    
    # Means at x ("drop" to drop the matrix attributes)
    return(drop(W %*% Y))
    
  }
  if (length(x)>1){
    # Matrix of size n x length(x)
    Kx <- sapply(X, function(Xi) K((x - Xi) / h) / h)
    
    # Weights
    W <- Kx / rowSums(Kx) # Column recycling!
    
    # Means at x ("drop" to drop the matrix attributes)
    return(drop(W %*% Y))
  }
}



############### simulation of the diffusion process with and without jumps 
# With jumps: dXt= b(Xt)dt + sigma(Xt)dWt + a(Xt-)dNt
simu_jumpdiff <- function(X0, grid, bfunc, sigfunc, afunc, Jumps, isjumpN){
  
  numjump <- length(Jumps)
  
  newgrid <- sort(c(grid, Jumps))
  
  W <- rnorm(length(grid)-1, 0, 1)
  X <- rep(X0, length(grid))    
  
  ####### simulation X_Ti dXt= b(Xt)dt + sig(Xt) dWt + a(Xt) sum_j dNjt
  
  X[2] <- X[1] + (grid[2] - grid[1]) * bfunc(X[1]) + sqrt(grid[2] - grid[1]) *  sigfunc(X[1]) * W[1] + afunc(X[1])*isjumpN[1]
  for (i in 2:(length(X)-1)){
    X[i+1] <- X[i] + (grid[i+1] - grid[i]) * bfunc(X[i]) + sqrt(grid[i+1] - grid[i]) *  sigfunc(X[i]) * W[i] + afunc(X[i])* isjumpN[i] 
  } 
  
  return(X = X)
}




# Without jumps
simu_diff <- function(X0, grid, bfunc, sigfunc){
  
  
  W <- rnorm(length(grid)-1, 0, 1)
  X <- rep(X0, length(grid))    # X_0 jusqua X_n+1 pour avoir les Y
  
  ####### simulation X_Ti dXt= b(Xt)dt + sig(Xt) dWt 
  
  X[2] <- X[1] + (grid[2] - grid[1]) * bfunc(X[1]) + sqrt(grid[2] - grid[1]) *  sigfunc(X[1]) * W[1] 
  for (i in 2:(length(X)-1)){
    X[i+1] <- X[i] + (grid[i+1] - grid[i]) * bfunc(X[i]) + sqrt(grid[i+1] - grid[i]) *  sigfunc(X[i]) * W[i] 
  } 
  return(X = X)
}
# or use sde package
