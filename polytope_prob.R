PolytopeProbVolume <- function(A, b){
  # This function computes the volume of the polytope P = {x|Ax<=b} intersected with the 
  # probability simplex
  # Input Value: initial polytope P ={x| Ax<=b}
  # Output Value: the volume of the polytope intersected with the probability simplex
  
  #install.packages("pracma")
  library("pracma")
  m <- nrow(A)
  n <- ncol(A)
  Ans <- matrix(0, nrow = m+1, ncol = n+m+1)
  
  # add slack variables x_{n+1}, ..., x_{n+m}
  # Ans = [ A I b ]
  #         1 0 1
  for(i in 1:m)
    for(j in 1:n)
      Ans[i,j] <- A[i,j]
  for(i in 1:m)
      Ans[i, n+i] <- 1
  for(i in 1:n)
      Ans[m+1, i] <- 1
  for(i in 1:m)
      Ans[i, n+m+1] <- b[i]
  Ans[m+1, n+m+1] <- 1
  
  # use Gauss-Jordan elimination to transform into its unique reduced row echelon form
  rAns <- rref(Ans)
  
  # New polytope {x| newAx<=newB}, which is a full-dimensional polytope in R^{n-1}
  newA <- matrix(0, nrow = m+1+n-1, ncol = n-1)
  newB <- matrix(0, nrow = m+1+n-1, ncol = 1)
  for(i in 1:(m+1))
    for(j in (m+2):(m+n))
        newA[i, (j-m-1)] <- rAns[i,j]
  for(i in 1:(m+1))
        newB[i] <- rAns[i, m+n+1]
  
  # ensure that x_i >0
  for(i in 1:(n-1))
        newA[m+1+i, i] <- -1
  for(i in 1:(n-1))
        newB[m+1+i] <- 0
  
  Vol <- PolytopeVolume(newA, newB, 1)
  Vol <- Vol * sqrt(n)
  return(Vol)
}