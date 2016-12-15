################################################################################
# Author : Asma Bahamyirou
# Goal: The IRLS-BMD algorithm for solving the Tweedie model grouped elastic net
################################################################################
z_plus <- function(z){z*(z > 0)}
################################################################################
#Tweedie_irls_bmd(X, y, lam=0.5, rho=1.5, tau=0.1, tol=1e-7, maxit=2e4, b=NULL, b0=NULL)

Tweedie_irls_bmd <- function(X, y, lam=.2, rho=1.5, tau=0.5, tol=1e-7, maxit=2e4, b=NULL, b0=NULL)
{
  
  p <- ncol(X)
  n <-length(y)
  w <- rep( sqrt(3), 3)
  
  ##(1) use zeros starting value if unspecified
  warm.start<- (!is.null(b)) & (!is.null(b0))
  if(!warm.start)
  {
    b <- rep(0, p)
    b0 <- 0
  }
  ##########################################################
  # compute the value of the objective function for initial value
  Lq_0 <- 0
  # working weight and response
  vi_tilde_0 <-  rep(NA, n)
  yi_tilde_0 <-  rep(NA, n)
  for(i in 1:n)
  {
    # The working weight
    vi_tilde_0[i] <- (1/n)*( (rho-1)*y[i]*exp(-(rho-1)*(b0+crossprod(b,X[i,])) ) +
                               (2-rho)*exp( (2-rho)*(b0+crossprod(b,X[i,])) )      )
    # The working response
    yi_tilde_0[i] <- b0 + crossprod(b,X[i,]) + ((1/n)/vi_tilde_0[i])*( y[i]*exp(-(rho-1)*(b0+crossprod(b,X[i,])) )
                                                                       - exp( (2-rho)*(b0+crossprod(b,X[i,])) ) )
  }
  for(i in 1:n)
  {
    Lq_0  <- Lq_0 + vi_tilde_0[i]*( yi_tilde_0[i]-b0 - crossprod(b,X[i,]) )^2
  }
  Lq_0  <- 0.5* Lq_0 
  pen1 <- tau*w[1]*sqrt( sum(b[1]^2+b[2]^2+b[3]^2) ) +0.5*(1-tau)*(sum(b[1]^2+b[2]^2+b[3]^2))
  pen2 <- tau*w[2]*sqrt( sum(b[4]^2+b[5]^2+b[6]^2) ) +0.5*(1-tau)*(sum(b[4]^2+b[5]^2+b[6]^2))
  pen3 <- tau*w[3]*sqrt( sum(b[7]^2+b[8]^2+b[9]^2) ) +0.5*(1-tau)*(sum(b[7]^2+b[8]^2+b[9]^2))
  Penalty <- lam*(  pen1+pen2+pen3   )
  
  Pq_0 <-  Lq_0 + Penalty
  
  ###########################################################################################
  
  iterating=TRUE
  iter=0
  while(iterating)
  {
    iter=iter+1
    
    ##(2) Update the penalized WLS objective function (8)
    
    # working weight and response
    vi_tilde <-  rep(NA, n)
    yi_tilde <-  rep(NA, n)
    for(i in 1:n)
    {
      # The working weight
      vi_tilde[i] <- (1/n)*( (rho-1)*y[i]*exp(-(rho-1)*(b0+crossprod(b,X[i,])) ) +
                               (2-rho)*exp( (2-rho)*(b0+crossprod(b,X[i,])) )      )
      # The working response
      yi_tilde[i] <- b0 + crossprod(b,X[i,]) + ((1/n)/vi_tilde[i])*( y[i]*exp(-(rho-1)*(b0+crossprod(b,X[i,])) )
                                                                     - exp( (2-rho)*(b0+crossprod(b,X[i,])) ) )
    }
    
    #Compute the Hessian matrix and max eigenvalue by (10)
    gam_j <- rep(NA, 3)
    for(j in 1:3)
    {
      if(j==1) {k=1:3}
      if(j==2) {k=4:6}
      if(j==3) {k=7:9}
      
      #the Hessian matrix
      Hj<- matrix(0, nrow=3, ncol=3)
      for(i in 1:n)
      {
        Hj <- Hj + vi_tilde[i]*X[i,k]%*%t(X[i,k])
      }
      # max eigen value
      gam_j[j] <-  max( eigen(Hj)$values )
      
    }
    
    #Compute gam_0
    gam_0 <- sum(vi_tilde)
    
    
    ##(3) Apply BMD Algo to obtain the minimizer of the PWLS function (8)
    
    ## We use the same initial value
    bb  <- b
    bb0 <- b0
    
    ## Updating step for beta and beta_0
    
    ## Updating beta
    bj_new <- rep(NA, 3)
    for(j in 1:3)
    {
      # compute Mu_j
      if(j==1) {k=1:3}
      if(j==2) {k=4:6}
      if(j==3) {k=7:9}
      Mu <- rep(0, 3)
      for(i in 1:n)
      {
        Mu  <- Mu + vi_tilde[i]*(yi_tilde[i]-bb0 - crossprod(bb,X[i,]))%*%X[i,k]
      }
      
      Mu_j <- - Mu
      
      #compute beta_j
      denom <- sqrt( tcrossprod(gam_j[j]*bb[k]-Mu_j,gam_j[j]*bb[k]-Mu_j) )
      bj_new  <- (gam_j[j]*bb[k]-Mu_j)*as.numeric(z_plus( 1-(lam*tau*w[j])/denom )/(gam_j[j] + lam*(1-tau)))
      
      # Updating
      bb[k] <- bj_new
    }
    
    ## Updating beta_0 
    # compute U_0
    Mu_0 <- 0
    for(i in 1:n)
    {
      Mu_0  <- Mu_0 + vi_tilde[i]*( yi_tilde[i]-bb0 - crossprod(bb,X[i,]) )
    }
    Mu_0 <- - Mu_0
    
    # compute beta_0 
    b0_new <- bb0 -(1/gam_0)*Mu_0
    
    # Updating
    bb0 <- b0_new
    
    ## Update beta and beta_0
    b <- bb
    b0  <- bb0
    
    
    
    ## convergence check
    Lq_new <- 0
    for(i in 1:n)
    {
      Lq_new  <- Lq_new + vi_tilde[i]*( yi_tilde[i]-b0 - crossprod(b,X[i,]) )^2
    }
    Lq_new  <- 0.5* Lq_new 
    pen1 <- tau*w[1]*sqrt( sum(b[1]^2+b[2]^2+b[3]^2) ) +0.5*(1-tau)*(sum(b[1]^2+b[2]^2+b[3]^2))
    pen2 <- tau*w[2]*sqrt( sum(b[4]^2+b[5]^2+b[6]^2) ) +0.5*(1-tau)*(sum(b[4]^2+b[5]^2+b[6]^2))
    pen3 <- tau*w[3]*sqrt( sum(b[7]^2+b[8]^2+b[9]^2) ) +0.5*(1-tau)*(sum(b[7]^2+b[8]^2+b[9]^2))
    Penalty <- lam*(  pen1+pen2+pen3   )
    
    Pq_new <-  Lq_new + Penalty
    
    diff <- Pq_new-Pq_0
    Pq_0 <- Pq_new
    
    if( ( abs(diff)< tol ) | (iter >= maxit)) 
      iterating=FALSE
    
    
  }
  
  
  return(list(b0=b0, b=b, total.iterations=iter))      
  
}



################################################################################
# Data Simulation case 1
################################################################################

# n=500
#  p=8
# Function p

p1 <- function(x){x}
p2 <- function(x){(3*x^2-1)/6}
p3 <- function(x){(5*x^3-3*x)/10}


#######
simu <- function(n,p,theta)
{
  n=n; p=p; theta=theta
  Sigma <- matrix(NA, nrow=p, ncol=p)
  for(j in 1:p)
  {
    for( k in 1:p)
    {
      if(j==k) { Sigma[j,k] <- 1}
      else
      { Sigma[j,k] <- theta }
    }
  }
  ## generate T
  T <-  matrix(rnorm(n*p, sd = Sigma), nrow=n, ncol=p)
  
  # Generate the response mean mu
  p1Tj <- matrix(NA, nrow=n, ncol=3)
  p2Tj <- matrix(NA, nrow=n, ncol=3)
  p3Tj <- matrix(NA, nrow=n, ncol=3)
  for(j in 1:3)
  {
    p1Tj[,j] = p1(T[,j])
    p2Tj[,j] = p2(T[,j])
    p3Tj[,j] = p3(T[,j])
  }
  log_mu <- 0.3+(-1)^2*(0.5*p1Tj[,1]+0.2*p2Tj[,1]+0.5*p3Tj[,1])+
    (-1)^3*(0.5*p1Tj[,2]+0.2*p2Tj[,2]+0.5*p3Tj[,2])+
    (-1)^4*(0.5*p1Tj[,3]+0.2*p2Tj[,3]+0.5*p3Tj[,3])
  
  mu <- exp(log_mu)
  
  #generate the response Y-values
  library(mgcv)
  y  <- rTweedie(mu,p=1.5,phi=1)
  
  # Final data look like this
  data_X <- cbind(p1Tj[,1],p2Tj[,1],p3Tj[,1] , p1Tj[,2],p2Tj[,2],p3Tj[,2] ,p1Tj[,3],p2Tj[,3],p3Tj[,3])
  data_X <- as.matrix(data_X)
  #colnames(data_X) <- c("p1T1","p2T1","p3T1","p1T2","p2T2","p3T2","p1T3","p2T3","p3T3")
  
  data=list(Y=y,X=data_X)
  data
}

################################################################################
# Data Simulation case 2
################################################################################

simu2 <- function(n,p,theta)
{
  n=n; p=n; theta=theta
  Sigma <- matrix(NA, nrow=p, ncol=p)
  for(j in 1:p)
  {
    for( k in 1:p)
    {
      if(j==k) { Sigma[j,k] <- 1}
      else
      { Sigma[j,k] <- theta }
    }
  }
  ## generate T
  Z <-  matrix(rnorm(n*p, sd = Sigma), nrow=n, ncol=p)
  
  # generate T_j
  
  T1 <- Z[,1]+rnorm(n,sd=0.001)
  T2 <- Z[,1]+rnorm(n,sd=0.001)
  T3 <- Z[,1]+rnorm(n,sd=0.001)
  T <- cbind(T1,T2,T3)
  
  # Generate the response mean mu
  p1Tj <- matrix(NA, nrow=n, ncol=3)
  p2Tj <- matrix(NA, nrow=n, ncol=3)
  p3Tj <- matrix(NA, nrow=n, ncol=3)
  for(j in 1:3)
  {
    p1Tj[,j] = p1(T[,j])
    p2Tj[,j] = p2(T[,j])
    p3Tj[,j] = p3(T[,j])
  }
  log_mu <- 0.3+(-1)^2*(0.5*p1Tj[,1]+0.2*p2Tj[,1]+0.5*p3Tj[,1])+
    (-1)^3*(0.5*p1Tj[,2]+0.2*p2Tj[,2]+0.5*p3Tj[,2])+
    (-1)^4*(0.5*p1Tj[,3]+0.2*p2Tj[,3]+0.5*p3Tj[,3])
  
  mu <- exp(log_mu)
  
  #generate the response Y-values
  library(mgcv)
  y  <- rTweedie(mu,p=1.5,phi=1)
  
  # Final data look like this
  data_X <- cbind(p1Tj[,1],p2Tj[,1],p3Tj[,1] , p1Tj[,2],p2Tj[,2],p3Tj[,2] ,p1Tj[,3],p2Tj[,3],p3Tj[,3])
  data_X <- as.matrix(data_X)
  #colnames(data_X) <- c("p1T1","p2T1","p3T1","p1T2","p2T2","p3T2","p1T3","p2T3","p3T3")
  
  data=list(Y=y,X=data_X)
  data
}

# X=simu2(1000,8,0.5)$X ; y=simu2(1000,8,0.5)$Y



################################################################################
# Count function for Tweedie
###############################################################################

count <-  function(X=X,y=y,lam=0.5, rho=1.5, tau=0.1, tol=1e-7, maxit=2e4)
{
  sol <- Tweedie_irls_bmd(X, y, lam=lam, rho=rho, tau=tau, tol=1e-7, maxit=2e4, b=NULL, b0=NULL)$b
  b1  <- sol[1:3]
  b2  <- sol[4:6]
  b3  <- sol[7:9]
  k=rep(0,3)
  # Count the number of active covarites in each block
  for(i in 1:3)
  {
    if(b1[i]!=0){k[1]=k[1]+1}
    if(b2[i]!=0){k[2]=k[2]+1}
    if(b3[i]!=0){k[3]=k[3]+1}
  }
  # for each block where we observe at least one active covariate, set k to 1
  for(i in 1:3)
  {
    if(k[1]!=0){k[1]=1}
    if(k[2]!=0){k[2]=1}
    if(k[3]!=0){k[3]=1}
  }
  # Count the number of active block
  count <- sum(k)
  count
  
}
# count(X=X,y=y,lam=1, rho=1.5, tau=.005, tol=1e-7, maxit=2e4)
# set.seed(66)

#X=simu2(500,8,0.5)$X ; y=simu2(500,8,0.5)$Y
# Tweedie_irls_bmd(X, y, lam=1, rho=1.5, tau=.005, tol=1e-7, maxit=2e4, b=NULL, b0=NULL)

# X=simu2(500,8,0)$X ; y=simu2(500,8,0)$Y
# Tweedie_irls_bmd(X, y, lam=.05, rho=1.5, tau=.005, tol=1e-7, maxit=2e4, b=NULL, b0=NULL)

################################################################################
# Monte carlo simulations for Count function
###############################################################################


boot <- function(sim,n,p,theta,lam=0.5, rho=1.5, tau=0.1, tol=1e-7, maxit=2e4,B=B)
{
  
  
  
  #### monte carlo
  N=rep(0,B)
  for(b in 1:B)
  {
    set.seed(b)
    # apply count function
    data=sim(n,p,theta)
    X=data$X ; y=data$Y
    # apply count function
    N[b] <- count(X=X,y=y,lam=lam, rho=rho, tau=tau, tol=1e-7, maxit=2e4)
  }
  
  res <- sum(N)/B
  res
}
##############################################
#### results ####
################################

# boot(simu,500,8,0.5,lam=0.5, rho=1.5, tau=0.1, tol=1e-7, maxit=2e4,B=100)
# 3
#  boot(simu,500,8,0,lam=0.5, rho=1.5, tau=0.1, tol=1e-7, maxit=2e4,B=100)
# 2.96
#  boot(simu2,500,8,0.5,lam=1, rho=1.5, tau=.005, tol=1e-7, maxit=2e4,B=100)
# 3
# boot(simu2,500,8,0,lam=0.05, rho=1.5, tau=.005, tol=1e-7, maxit=2e4,B=100)
# 2.07





################################################################################
# gglasso for comparison with Tweedie model
################################################################################

## Generate data set
#X=simu(500,8,0)$X ; y=simu(500,8,0)$Y


###############################
#  
library(gglasso)
# load gglasso library

#group <- rep(1:3,each=3)
#m2 <- gglasso(x=X,y=y,group=group,loss="ls")
#cv <- cv.gglasso(x=X,y=y, group=group, loss="ls",
pred.loss="L2", lambda.factor=0.05, nfolds=5)

###########
#####
#cv <- cv.gglasso(x=X, y=y, group=group, loss="ls",
pred.loss="L2", lambda.factor=0.05, nfolds=5)
#coef(cv, s = "lambda.min")
###############################################################################

count_gglasso <-  function(X=X,y=y)
{
  # load gglasso library
  library(gglasso)
  
  group <- rep(1:3,each=3)
  #m2 <- gglasso(x=X,y=y,group=group,loss="ls")
  cv <- cv.gglasso(x=X, y=y, group=group, loss="ls",
                   pred.loss="L2", lambda.factor=0.05, nfolds=5)
  sol <- as.numeric( coef(cv, s = "lambda.min")) # get beta estimate for selected lambda
  b1  <- sol[2:4]
  b2  <-  sol[5:7]
  b3  <-  sol[8:10]
  k=rep(0,3)
  # Count the number of active covarites in each block
  for(i in 1:3)
  {
    if(b1[i]!=0){k[1]=k[1]+1}
    if(b2[i]!=0){k[2]=k[2]+1}
    if(b3[i]!=0){k[3]=k[3]+1}
  }
  # for each block where we observe at least one active covariate, set k to 1
  for(i in 1:3)
  {
    if(k[1]!=0){k[1]=1}
    if(k[2]!=0){k[2]=1}
    if(k[3]!=0){k[3]=1}
  }
  # Count the number of active block
  count <- sum(k)
  count
  
}

########### 
boot_gglasso <- function(sim,n,p,theta,B=B)
{
  #### monte carlo
  N=rep(0,B)
  for(b in 1:B)
  {
    set.seed(b)
    # apply count function
    data=sim(n,p,theta)
    X=data$X ; y=data$Y
    N[b] <- count_gglasso(X=X,y=y)
  }
  
  res <- sum(N)/B
  res
}

#boot_gglasso(simu,500,8,0.5,B=100) #1.11

# boot_gglasso(simu,500,8,0,B=100) #1.31

#boot_gglasso(simu2,500,8,0.5,B=100) #

#boot_gglasso(simu2,500,8,0,B=100) #


################################################################################
# Lasso for comparison with Tweedie model
################################################################################



# count_Lasso(X=X,y=y)
count_Lasso <-  function(X=X,y=y)
{
  library(glmnet)
  cvfit = cv.glmnet(X, y)
  sol <- as.numeric(coef(cvfit, s = "lambda.min")) # get beta estimate for selected lambda
  b1  <- sol[2:4]
  b2  <-  sol[5:7]
  b3  <-  sol[8:10]
  k=rep(0,3)
  # Count the number of active covarites in each block
  for(i in 1:3)
  {
    if(b1[i]!=0){k[1]=k[1]+1}
    if(b2[i]!=0){k[2]=k[2]+1}
    if(b3[i]!=0){k[3]=k[3]+1}
  }
  # for each block where we observe at least one active covariate, set k to 1
  for(i in 1:3)
  {
    if(k[1]!=0){k[1]=1}
    if(k[2]!=0){k[2]=1}
    if(k[3]!=0){k[3]=1}
  }
  # Count the number of active block
  count <- sum(k)
  count
  
}

########### 
boot_Lasso <- function(sim,n,p,theta,B=B)
{
  #### monte carlo
  N=rep(0,B)
  for(b in 1:B)
  {
    set.seed(b)
    # apply count function
    data=sim(n,p,theta)
    X=data$X ; y=data$Y
    N[b] <- count_Lasso(X=X,y=y)
  }
  
  res <- sum(N)/B
  res
}

# boot_Lasso(simu,500,8,0.5,B=100) #1.22
# boot_Lasso(simu,500,8,0,B=100) #1.79

# boot_Lasso(simu2,500,8,0.5,B=100) #2.27
#boot_Lasso(simu2,500,8,0,B=100) #0.56



