PolytopeVolume <- function(A, Bl, Bu, epsilon){
  # This functon implements randomzied algotihm for accurately approximating
  # for the polytope's volume in high dimensions.
  
  #-- INPUT VALUES--
  # Polytope, {x| Bl <= Ax <= Bu}
  # A: nrow = m, ncol = n
  # Bl: nrow = m, ncol = 1
  # Bu: nrow = m, ncol = 1
  # epsilon: object approximation, set epsilon <- 1
  
  #-- RETURN VALUES--
  # Volume: the computed volume
  
  # Packages needed for dependencies
  #install.packages("lpSolve")
  #install.packages("OptimalDesign")
  #install.packages("extraDistr")
  
  m <- nrow(A)
  n <- ncol(A)
  tr <- 1+10^(-1.2) #rounding threshold
  
  N <- floor((400*n*log(n,2))/(epsilon^2))
  W <- floor(10+n/10)
  
  #Changing the form from "Bl <= Ax = Bu" to "Ax <= b"
  AMinus <- -1*A
  AA <- rbind(A, AMinus)
  BlMinus <- -1*Bl
  BB <- rbind(Bu, BlMinus)
  
  # Comouting ChebyshevBall of polytope P = { x | (CA)x <= Cb}
  # Return Values: radius and center of the ball
  # e.g.   BALL <- ChebyshevBall(AA, BB)
  #        Radius <- BALL[[1,2]]
  #        Center <- BALL[[1,1]]
  ChebyshevBall <- function(CA, Cb){
    # This problem can be formulated a linear optimization problem.
    Cm <- nrow(CA)
    Cn <- ncol(CA)
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
    
    library(lpSolve)
    Cans <- lp ("max", Cf.obj, Cf.con, Cf.dir, Cf.rhs)$solution
    
    # return the output
    Cx <- matrix(nrow = Cn, ncol = 1)
    for(i in 1:Cn){
      Cx[i] <- Cans[i]-Cans[i+Cn]
    }
    CAns <- matrix(list(), nrow = 1, ncol = 2)
    CAns[[1,1]] <- Cx
    CAns[[1,2]] <- Cans[2*Cn+1]
    
    return(CAns)
  }
  
  BALL <- ChebyshevBall(AA, BB)
  
  Radius <- BALL[[1,2]]
  Center <- BALL[[1,1]]
  
  p0 <- Center
  
  
  # Input sets: A = {a_1, a_2,..., a_m} points from P
  Rounding <- function(RA){
    Rcol <- ncol(RA) #number of the sample points
    Rrow <- nrow(RA) # dimension of the points
    
    RP <- matrix(nrow = Rcol)
    for(Ri in 1:Rcol){
      RP[Ri, 1] <- 1/Rcol
    }
    
    Rk <- 1
    Pk1 <- matrix(0, nrow = Rcol, ncol = Rcol)
    
    repeat{
      Rj <- 0
      Rjans <- 0
      
      for(Rbb in 1:Rcol){
        Pk1[Rbb, Rbb] <- RP[Rbb, Rk]
      }
      
      Pm <- RA%*%Pk1%*%t(RA)
      
      PI <- diag(Rrow)
      Pminverse <- solve(Pm, PI)
      for(Ri in 1 : Rcol){
        Rq <- RA[ ,Ri]
        if(t(Rq)%*%Pminverse%*%Rq>Rjans){
          Rjans <- t(Rq)%*%Pminverse%*%Rq
          Rj <- Ri
        }
      }
      qj <- RA[, Rj]
      kappa <- t(qj)%*%Pminverse%*%qj
      
      
      beta <- c((kappa-Rrow)/(Rrow*(kappa-1)))
      
      Re <- matrix(0, nrow = Rcol, ncol = 1)
      Re[Rj] <- 1
      P2 <- (1-beta)*RP[ ,Rk]+beta*Re
      
      Rk <- Rk+1
      RP <- cbind(RP,P2)
      
      Rflag <- 0
      for(Rj in 1:Rcol){
        if(t(RA[ ,Rj])%*%Pminverse%*%RA[ ,Rj] > 2*Rrow){
          Rflag <- 1
          break
        }
      }
      
      if(Rflag == 0) 
        break
      
    }
    return(RP[ ,Rk])
  }
  
  
  Ellipsoid <- function(EA){
    Em <- nrow(EA)
    En <- ncol(EA)
    EAMinus <- -1*EA
    EA1 <- cbind(EA, EAMinus)
    EA2 <- matrix(0, nrow = Em+1, ncol = 2*En)
    
    for(i in 1:En){
      EA2[Em+1,i] <- 1
      EA2[Em+1, En+i] <- -1
    }
    
    for( i in 1:Em){
      for(j in 1:(2*En)){
        EA2[i,j] <- EA1[i,j]
      }
    }
    Ep1 <- Rounding(EA2)
    
    print(Ep1)
    
    if(is.null(nrow(Ep1)) == TRUE) return (c(-1))
    Ep <- matrix(0, nrow = En, ncol = 1)
    for(i in 1:En)
      Ep[i] <- Ep1[i]+Ep1[En+i]
    
    A11 <- matrix(0, nrow = Em, ncol = Em)
    for(i in 1:En){
      A11 <- A11 + Ep[i]*(EA[, i] %*% t(EA[, i]))
    }
    
    A12 <- matrix(0, nrow = Em, ncol = 1)
    for(i in 1:En){
      A12 <- A12 + Ep[i]*EA[, i]
    }
    
    A21 <- matrix(0, nrow = 1, ncol = Em)
    for(i in 1 :En){
      A21 <- A21 + Ep[i]*t(EA[ ,i])
    }
    A22 <- matrix(1, nrow = 1, ncol = 1)
    
    A1112 <- cbind(A11, A12)
    A2122 <- cbind(A21, A22)
    ATot <- rbind(A1112, A2122)
    
    I1 <- diag(Em+1)
    Ainv <- solve(ATot, I1)
    
    Ainv <- (1/sqrt(1+Em)) * Ainv
    
    Ainv11 <- matrix(0, nrow = Em, ncol = Em)
    Ainv12 <- matrix(0, nrow = Em, ncol = 1)
    Ainv21 <- matrix(0, nrow = 1, ncol = Em)
    Ainv22 <- matrix(0, nrow = 1, ncol = 1)
    for(i in 1:Em)
      for(j in 1:Em)
        Ainv11[i,j] <- Ainv[i,j]
    for(i in 1:Em)
      Ainv12[i] <- Ainv[Em+1, i]
    for(i in 1:Em)
      Ainv21[i] <- Ainv[i, Em+1]
    Ainv22 <- Ainv[Em+1,Em+1]
    
    I2 <- diag(Em)
    Ainv11inv <- solve(Ainv11, I2)
    
    EE <- -1*Ainv11inv%*%Ainv12
    
    EAns <- matrix(list(), nrow = 1, ncol = 2)
    EAns[[1,1]] <- EE   #c_{epsilon}
    EAns[[1,2]] <- Ainv11 #E
    
    return(EAns)
  }
  
  
  Boundary <- function(BPnew, BBl, BBu, BA, BFLAG, Bc, Br){
    BAPnew <- BA%*%BPnew
    BApnewBl <- (BA%*%BPnew >= BBl)
    BApnewBu <- (BA%*%BPnew <= BBu)
    Bm <- nrow(A)
    BFLAG1 <- 0
    for(Bj in 1:Bm){
      if(BApnewBl[Bj] == FALSE){
        BFLAG1 <-1
        break
      }
      if(BApnewBu[Bj] == FALSE){
        BFLAG1 <-1
        break
      }
    }
    if(BFLAG == 0){
      if(BFLAG1 == 1) {
        return (FALSE)
      }else {return (TRUE)
      }
    }
    if(BFLAG == 1){
      if(BFLAG1 == 1) return (FALSE)
      else{
        if(sqrt(sum((BPnew-Bc)^2)) <= Br){
          return(TRUE)
        }else{
          return(FALSE)
        }
      }
    }
  }
  
  Walk <- function(wp, wBl, wBu, wA, wFLAG, wc, wr, ww){
    library(extraDistr)
    Left <- -2^36
    wm <- nrow(wA)
    wn <- ncol(wA)
    for(wk in 1:ww){
      wlastp <- wp
      ii <- rdunif(1,1,wn)
      we <- matrix(0, nrow = wn, ncol = 1)
      we[ii] <- 1
      Left <- -2^36
      Right <- 0
      while(Right-Left>10^-10){
        Mid <- 0.5*(Left+Right)
        wPnew <- wp+Mid*we
        if(Boundary(wPnew,wBl, wBu, wA, wFLAG, wc, wr) == TRUE){
          Right <- Mid
        }
        else{
          Left <- Mid
        }
      }
      lambdamin <- Left
      Left <- 0
      Right <- 2^36
      while(Right-Left>10^-10){
        Mid <- 0.5*(Left+Right)
        wPnew <- wp+Mid*we
        if(Boundary(wPnew,wBl, wBu, wA, wFLAG, wc, wr) == TRUE){
          Left <- Mid
        }
        else{
          Right <- Mid
        }
      }
      lambdamax <- Left
      lambda <- runif(1,lambdamin,lambdamax)
      wp <- wlastp +lambda*we
    }
    return(wp)
  }
  
  TagBl <- Bl
  TagBu <- Bu
  TagA <- A
  Newp <- p0
  VolL <- 1
  
  repeat{
    S <- Walk(Newp, TagBl, TagBu, TagA, 0, 0, 0, W)
    p0 <- S
    for(kk in 2:N){
      p0 <- Walk(p0, TagBl, TagBu, TagA, 0, 0, 0, W)
      S <- cbind(S, p0)
    }
    
    EEE <- Ellipsoid(S)
    if(is.null(nrow(EEE)) == TRUE) {next}
    E <- EEE[[1,2]]
    CEplison <- EEE[[1,1]]
    Eeign <- eigen(E)$val
    Emax <- Eeign[1]
    Emin <- Eeign[1]
    for(jj in 2:n){
      if(Eeign[jj]>Emax) Emax <- Eeign[jj]
      if(Eeign[jj]<Emin) Emin <- Eeign[jj]
    }
    
    L <- chol(E)
    VolL <- VolL*(det(L))
    I <- diag(n)
    LTinverse <- solve(L,I)
    TagBl <- TagBl - TagA%*%CEplison
    TagBu <- TagBu - TagA%*%CEplison
    TagA <- TagA%*%LTinverse
    Newp <- L%*%(p0-CEplison)
    Scol <- ncol(S)
    for(jj in 1: Scol){
      S[ ,jj] <- L%*%(S[, jj]-CEplison)
    }
    
    if(sqrt(Emax/Emin) < tr) break
  }
  
  TagAMiuns <- -1*TagA
  TagAA <- rbind(TagA, TagAMiuns)
  TagBlMiuns <- -1*TagBl
  TagBB <- rbind(TagBu, TagBlMiuns)
  
  BALL <- ChebyshevBall(TagAA, TagBB)
  Radius <- BALL[[1,2]]
  Center <- BALL[[1,1]]
  
  
  
  rho <- 0
  for(i in 1: Scol){
    if(sqrt(t(S[ ,i]-Center)%*%(S[ ,i]-Center)) >= rho)
      rho <- sqrt(t(S[ ,i]-Center)%*%(S[ ,i]-Center))
  }
  
  
  alpha <- round(n*log(Radius,2),1)
  beta <-  round(n*log(rho,2),1)
  
  
  Fi <- beta
  Snew <- matrix(0, nrow = n, ncol = 80000)
  Vol <- 2*(pi^(n/2))*((2^alpha)/(n*gamma(n/2)))
  
  wcount <- ncol(S)
  
  while(Fi > alpha){
    
    countprev <- wcount
    wcount <- 0
    LargeR <- 2^((Fi)/n)
    SmallR <- 2^((Fi-0.1)/n)
    
    for(Fj in 1:countprev){
      if(Boundary(S[ , Fj], TagBl, TagBu, TagA, 1, Center, SmallR)==TRUE){
        wcount <- wcount + 1
        Snew[ ,wcount] <- S[, Fj]
      }
    }
    
    p0 <- Snew[ , 1]
    if(N-countprev > 0){
      for(Fj in 1: (N-countprev)){
        p0 <- Walk(p0, TagBl, TagBu, TagA, 1, Center, LargeR, W)
        if(sqrt(t(p0-Center)%*%(p0-Center)) <= SmallR){
          wcount <- wcount+1
          Snew[ , wcount] <- p0
        }
      }
    }
    
    Vol <- Vol * (N/wcount)
    Fi <- Fi - 0.1
    S <- Snew
  }
  
  return(Vol/VolL)
  
}
