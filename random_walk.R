PolytopeVolume <- function(A, B, epsilon){
  # This functon implements randomzied algotihm for accurately approximating
  # for the polytope's volume in high dimensions.
  
  #-- INPUT VALUES--
  # Polytope, {x| Ax <= B}
  # A: nrow = m, ncol = n
  # B: nrow = m, ncol = 1
  # epsilon: object approximation, set epsilon <- 1
  
  #-- RETURN VALUES--
  # Volume: the computed volume
  
  # Packages needed for dependencies
  #install.packages("lpSolve")
  #install.packages("OptimalDesign")
  #install.packages("extraDistr")
  
  library(extraDistr)
  library(lpSolve)
  library(OptimalDesign)
  
  m <- nrow(A)
  n <- ncol(A)
  tr <- 1.414 #rounding threshold
  
  N <- floor((400*n*log(n,2))/(epsilon^2))
  W <- floor(10+n/10)
  
  #Changing the form from "Bl <= Ax = Bu" to "Ax <= b"
  #AMinus <- -1*A
  #AA <- rbind(A, AMinus)
  #BlMinus <- -1*Bl
  #BB <- rbind(Bu, BlMinus)
  
  
  # Comouting ChebyshevBall of polytope P = { x | (CA)x <= Cb}
  # Return Values: radius and center of the ball
  # e.g.   BALL <- ChebyshevBall(AA, BB)
  #        Radius <- BALL[[1,2]]
  #        Center <- BALL[[1,1]]
  ChebyshevBall <- function(CA, Cb, Cm, Cn){
    # This problem can be formulated a linear optimization problem.
    #Cm <- nrow(CA)
    # Cn <- ncol(CA)
    Cf.obj <- matrix(0, nrow = 2*Cn+1, ncol = 1)
    Cnnew <- 2*Cn
    for(i in 1:Cnnew)
      Cf.obj[i] <- 0
    Cf.obj[2*Cn+1] <- 1
    CC <- matrix(0, nrow = Cm, ncol = 1)
    for(i in 1:Cm){
      CC[i] <- sqrt(sum((CA[i, ])^2))
    }
    CD1 <- cbind(CA, CC)
    CAMinus <- -1*CA
    CD <- cbind(CA, CAMinus)
    CD <- cbind(CD, CC)
    Cf.con <- CD
    Cf.rhs <- t(Cb)
    Cf.dir <- matrix(nrow = Cm, ncol = 1)
    for(i in 1:Cm){
      Cf.dir[i] <- "<="
    }
    
    Cans <- lp ("max", Cf.obj, Cf.con, Cf.dir, Cf.rhs)$solution
    
    # Organizing the output format
    Cx <- matrix(nrow = Cn, ncol = 1)
    for(i in 1:Cn){
      Cx[i] <- Cans[i]-Cans[i+Cn]
    }
    CAns <- matrix(list(), nrow = 1, ncol = 2)
    CAns[[1,1]] <- Cx
    CAns[[1,2]] <- Cans[2*Cn+1]
    return(CAns)
  }
  
  BALL <- ChebyshevBall(A, B, m, n)
  Radius <- BALL[[1,2]]
  Center <- BALL[[1,1]]
  
  p0 <- Center
  
  
  # D-optimal Design Problem
  # Input Values: A = {a_1,..., a_m} from polytope
  # Output Values: p_1, p_2, ... p_m
  Rounding <- function(RA, Rm){
    if(rcond(RA%*%t(RA)) < (10^-15)) return (c(-1))
    else{
      PRINT <- capture.output(Rout <- od.AA(t(RA), N=1, crit="D",tab = NULL)$w.best)
      return(t(Rout))
    }
  }
  
  # Finding the smallest ellipsoid that cover S of O(n) random points in P
  # Input values: A = {a_1, a_2, ..., a_m} from polytope
  # Output Values: radius and center of this ellipsoid E = {x| (x-c_{epsilon}^TE(x-c_{epsilon}) <= 1}
  # e.g.     EEE <- Ellipsoid(S)
  #          E <- EEE[[1,2]]
  #          CEplison <- EEE[[1,1]]
  EA2 <- matrix(0, nrow = 1, ncol = 2*N)
  Ep <- matrix(0, nrow = N, ncol = 1)
  EI <- diag(n)
  EAns <- matrix(list(), nrow = 1, ncol = 2)
  EP <- matrix(0, nrow = N, ncol = N)
  for(i in 1:N){
    EA2[i] <- 1
    EA2[N+i] <- -1
  }
  Ellipsoid <- function(EA, Em, En){
    #Em <- nrow(EA)
    #En <- ncol(EA)
    EAMinus <- -1*EA
    EA1 <- cbind(EA, EAMinus)

    EA2 <- rbind(EA1, EA2)
    Ep1 <- Rounding(EA2, Em+1)
    if(is.null(nrow(Ep1)) == TRUE) return (c(-1))
    
    for(i in 1:En)
      Ep[i] <- Ep1[i]+Ep1[En+i]
    
    for( Ei in 1: En){
      EP[Ei,Ei] <- Ep[Ei]
    }
    Ec <- EA%*%Ep
    EQ1 <- EA%*%EP%*%t(EA)-Ec%*%t(Ec)
    
    EQinverse <- solve(EQ1, EI)
    EQ <- EQinverse/Em
    
    EAns[[1,1]] <- Ec    #c_{epsilon}
    EAns[[1,2]] <- EQ    #Q
    
    return(EAns)
  }
  
  
  # Hit-and-Run random walk Method: Generating a uniform distribution of points
  # Input Vaule: wp: the initial point
  #                  Polytope {x| wBl <= wAx <= wBu}, ball {(x-wc)^T(x-wc)<=wr^2}
  #                  ww is repeat times
  # Output Value: point generated 
  Walk <- function(wp, wA, wB, wFLAG, wc, wr, ww, wm, wn){
    #AlphapMax, AlphapMin, AlphaeMax, AlphaeMin
    #wm <- nrow(wA)
    #wn <- ncol(wA)
    for(wk in 1:ww){
      ii <- rdunif(1,1,wn)
      we <- matrix(0, nrow = wn, ncol = 1)
      we[ii] <- 1
      AlphapMax <- Inf
      AlphapMin <- -Inf
      for(wi in 1:wm){
        wlastp <- wp
        temp <- (wB[wi] - wA[wi, ]%*%wlastp)/wA[wi,ii]
        if(wA[wi, ii] < 0){
          if(temp > AlphapMin){
            AlphapMin <- temp        
          }
        }
        if(wA[wi,ii] > 0){
          if(temp < AlphapMax){
            AlphapMax <- temp        
          }
        }
      }
      if(wFLAG == 0){
        Alpha <- runif(1,AlphapMin,AlphapMax)
        wp <- wlastp +Alpha*we
      }
      if(wFLAG == 1){
        wD <- 1
        wE <- 2*(wlastp[ii]-wc[ii])
        wF <- t(wlastp)%*%wlastp-2*(t(wlastp)%*%wc)+t(wc)%*%wc-t(wr)%*%wr
        AlphaeMax <- (-1*wE+sqrt(wE*wE-4*wD*wF))/(2*wD)
        AlphaeMin <- (-1*wE-sqrt(wE*wE-4*wD*wF))/(2*wD)
        AlphaMax <- min(AlphapMax, AlphaeMax)
        AlphaMin <- max(AlphapMin, AlphaeMin)
        Alpha <- runif(1, AlphaMin, AlphaMax)
        wp <- wlastp + Alpha*we
      }
    }
    return (wp)
  }
  
  TagB <- B
  TagA <- A
  Newp <- p0
  VolL <- 1
  Snew <- matrix(0, nrow = n, ncol = N)
  S <- matrix(0, nrow = n, ncol = N)
  I <- diag(n)
  # Rounding the polytope
  repeat{
    # Sampling a set S of O(n) random points in P.
    
    for(kk in 1:N){
      
      if(kk == 1){
        S[ ,kk] <- Walk(Newp, TagA, TagB, 0, 0, 0, W, m, n)
        p0 <- S[ ,kk]
      }
      else{
        S[ ,kk] <- Walk(p0, TagA, TagB, 0, 0, 0, W, m, n)
        p0 <- S[ ,kk]
      }
    }
    
    # Computing min ellipsoid E of S, with p.s.d. E
    EEE <- Ellipsoid(S, n, N)
    if(is.null(nrow(EEE)) == TRUE) {next}
    E <- EEE[[1,2]]
    CEplison <- EEE[[1,1]]
    
    print(E)
    print(CEplison)
    
    Eeign <- eigen(E)$val
    
    # Set as E_{min} and E_{max} the min and max axes
    Emax <- Eeign[1]
    Emin <- Eeign[1]
    for(jj in 2:n){
      if(Eeign[jj]>Emax) Emax <- Eeign[jj]
      if(Eeign[jj]<Emin) Emin <- Eeign[jj]
    }
    
    # Computing the Cholesky decomposetion L^TL of E
    L <- chol(E)
    
    # Tranform P and p w.r.t L
    VolL <- VolL*(det(L))
    
    LTinverse <- solve(L,I)
    TagB <- TagB - TagA%*%CEplison
    TagA <- TagA%*%LTinverse
    Newp <- L%*%(p0-CEplison)
    
    print(VolL)
    
    if(sqrt(Emax/Emin) < tr) break
  }
  
  # Multiphase Monte Carlo
  
  # Computing the Chebyshev ball of the transformed polytope
  BALL <- ChebyshevBall(TagA, TagB, m, n)
  Radius <- BALL[[1,2]]
  Center <- BALL[[1,1]]
  
  for(jj in 1: N){
    S[ ,jj] <- L%*%(S[, jj]-CEplison)
  }
  
  
  # Set rho the largest distance from c to any point in S
  rho <- 0
  for(i in 1: N){
    if(sqrt(t(S[ ,i]-Center)%*%(S[ ,i]-Center)) >= rho)
      rho <- sqrt(t(S[ ,i]-Center)%*%(S[ ,i]-Center))
  }
  
  # Init
  Step <- 1
  alpha <- round(n*log(Radius,2),1)
  beta <-  round(n*log(rho,2),1)
  Fi <- beta
  
  Vol <- 2*(pi^(n/2))*((2^alpha)/(n*gamma(n/2)))
  
  #Vol <- (pi^(n/2)*Radius^(n))/(gamma(n/2+1))
  
  wcount <- N
  
  while(Fi > alpha){
    
    countprev <- wcount
    wcount <- 0
    LargeR <- 2^((Fi)/n)
    SmallR <- 2^((Fi-Step)/n)
    
    # Remove from S the points not in P_{small}
    for(Fj in 1:countprev){
      if(sqrt(t(S[ ,Fj]-Center)%*%(S[ ,Fj]-Center)) <= SmallR){
        wcount <- wcount + 1
        Snew[ ,wcount] <- S[, Fj]
      }
    }
    
    p0 <- Snew[ , 1] # Set p to be an arbitrary point from S
    if(N-countprev > 0){
      for(Fj in 1: (N-countprev)){
        p0 <- Walk(p0, A, B, 1, Center, LargeR, W, m, n)
        if(sqrt(t(p0-Center)%*%(p0-Center)) <= SmallR){
          wcount <- wcount+1
          Snew[ , wcount] <- p0
        }
      }
    }
    
    Vol <- Vol * (N/wcount)
    Fi <- Fi - Step
    S <- Snew
    print(Vol)
    
    
  }
  
  return(Vol/VolL)
  
}
