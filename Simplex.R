ProbSimplex <- function(n){
  #t1 <- Sys.time()
  #Calculating a n-dimensional probability simplex
  A <- matrix(0, nrow = n, ncol =n-1)
  B <- matrix(0, nrow = n, ncol = 1)
  for(i in 1:n-1)
    A[i,i] <- -1
  for(i in 1:n-1)
    A[n,i] <- 1
  B[n] <- 1
  V <- PolytopeVolume(A, B,1)
  Ans <- matrix(list(), nrow = 1, ncol = 2)
  ansV <- sqrt(n)*V
  #t2 <- Sys.time()
  #df <- t2-t1
  Ans[[1,1]] <- ansV
  #Ans[[1,2]] <- df
  return(Ans[[1,1]])
}