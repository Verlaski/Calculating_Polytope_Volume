DOptimal <- function(EA){
  
  # transfroming this point set to the centrally symmetric set
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
  
  RA <- EA2
  
  Rcol <- ncol(RA) #number of the sample points
  Rrow <- nrow(RA) # dimension of the points
  
  # Start points: RP[ ,1] <- 1/m e
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
      if(t(RA[ ,Rj])%*%Pminverse%*%RA[ ,Rj] > 1.01*Rrow){
        Rflag <- 1
        break
      }
    }
    
    if(Rflag == 0) 
    break
    
  }
  return(RP[ ,Rk])
}
