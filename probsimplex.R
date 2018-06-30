ProbSimplex <- function(A, Bl, Bu, epsilon){
  # This functon computes the volumes of a polytope intersected with the probability simplex
  
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
  

  ChebyshevBall <- function(CA, Cb){
    # Comouting ChebyshevBall of polytope P = { x | (CA)x <= Cb} intersected with probability simplex
    # Return Values: radius and center of the ball
    # e.g.   BALL <- ChebyshevBall(AA, BB)
    #        Radius <- BALL[[1,2]]
    #        Center <- BALL[[1,1]]
    # This problem can be formulated a linear optimization problem.
    Cm <- nrow(CA); Cn <- ncol(CA)  # CA is a Cm-by-Cn matrix
    
    # object function: max r 
    Cf.obj <- matrix(0, nrow = Cn+1, ncol = 1)
    for(i in 1:Cn)
      Cf.obj[i] <- 0
    Cf.obj[Cn+1] <- 1
    
    # Constrain 
    CC <- matrix(0, nrow = Cm, ncol = 1)
    for(i in 1:Cm){
      CC[i] <- sqrt(sum((CA[i, ])^2))
    }
    
    CD1 <- cbind(CA, CC)
    CD2 <- matrix(0, nrow = 1, ncol = Cn+1)
    for(i in 1:Cn){
      CD2[i] <- 1
    }
    CD2[Cn+1] <- 0; CD <- rbind(CD1, CD2); Cf.con <- CD
    
    Cb2 <- matrix(0, nrow = 1, ncol = Cm+1)
    for(i in 1: Cm)
      Cb2[i] <- Cb[i]
    Cb2[Cm+1] <- 1; Cf.rhs <- t(Cb2)
    
    Cf.dir <- matrix(nrow = Cm+1, ncol = 1)
    for(i in 1:Cm){
      Cf.dir[i] <- "<="
    }
    Cf.dir[Cm+1] <- "=="
    
    print(Cf.obj)
    print(Cf.con)
    print(Cf.dir)
    print(Cf.rhs)
    
    library(lpSolve)
    Cans <- lp ("max", Cf.obj, Cf.con, Cf.dir, Cf.rhs)$solution
    
    # return the output
    Cx <- matrix(nrow = Cn, ncol = 1)
    for(i in 1:Cn){
      Cx[i] <- Cans[i]
    }
    CAns <- matrix(list(), nrow = 1, ncol = 2)
    CAns[[1,1]] <- Cx
    CAns[[1,2]] <- Cans[Cn+1]
    
    return(CAns)
  }
  
  BALL <- ChebyshevBall(AA, BB)
  
  Radius <- BALL[[1,2]]
  Center <- BALL[[1,1]]
  
  
  
}