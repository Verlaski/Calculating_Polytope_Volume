PolytopeVolume <- function (A, Bl, Bu){
# This functon is a algorithm to determine the volume of a polytope
# The idea is summing the volumes of simplices which from polytope
  
#-- INPUT VALUES--
# the polytope, {x| Bl <= Ax <= Bu}
# A: nrow = m, ncol = n
# Bl: nrow = m, ncol = 1
# Bu: nrow = m, ncol = 1

#-- RETURN VALUES--
# Volume: the computed volume
  
#-- Variables --
# m: the number of the inequalities
# n: the dimension of the polytope
# Vertex: the coordinates of the vertexs of the polytope
# vertexNumber: the number of the vertex
# F: logical matrix, F(i,j) = 1 if the set F_i contain jth vertex
# F_i =\psi(e_l^{n-1})
  
# This function needs installing package "combinat"
install.packages("combinat") 

#initialize C, D, m, n.  
AMinus <- -1*A
C <- rbind(A, AMinus)
BuMinus <- -1*Bu
D <- rbind(Bl, BuMinus)
m <- nrow(C)
n <- ncol(C)

# Storing the coordinates of the vertexs, updating F.
library(combinat)
Comb <- combn(c(1:m), n)
CombColumn <- ncol(Comb)
Vertex <- matrix(nrow = n,ncol = 1000)

vertexNumber <- 0

F <- matrix(0, nrow = 100, ncol = 100)

# Store the coordinates of a vertex and Update F_i
  for (i in 1:CombColumn){
      Aa <- C[Comb[1,i],]
      Bb <- D[Comb[1,i],]
      for(j in 2:n){
        Aa <- rbind(Aa, C[Comb[j,i],])
        Bb <- rbind(Bb, D[Comb[j,i],])
      }
      if(det(Aa) != 0){
        x <- solve(Aa, Bb)
        FLAG <- 0
        CxD <- C%*%x >= D
        
        for(j in 1:m){
          if(CxD[j] == FALSE){
            FLAG <- 1; break
          }
        }
        if(vertexNumber > 0){
          for(j in 1: vertexNumber){
            FLAG1 <- 0
            CxD <- x == Vertex[ ,j]
            for (k in 1:n){
              if(CxD[k] == FALSE){
                FLAG1 <- 1; break
              }
            }
            if(FLAG1 == 0){
              FLAG <- 2; break
            }
          }
        }
        if(FLAG == 0){
          vertexNumber <- vertexNumber + 1
          Vertex[, vertexNumber] <- x
          for (j in 1:n){
            F[Comb[j,i], vertexNumber] <- 1
          }
        }
        if(FLAG == 2){
          for (j in 1:n){
            F[Comb[j,i], vertexNumber] <- 1
          }
        }
      }
  }

VolumeCal <- function(d, last, S){
  # This function is a recursive depth-first procedure to generate all possible simplexs
  # initialize L to contain the empty set
  
  #--- Variable --
  # d: dimension
  # last: the vertices of the cell e_j'^{d+1}
  # last[i] == 1 if the set last contains the ith vertex 
  # S: the simplices currently being computed
  # S[i] == 1 if the ith vertex is the simplex is computed
  # SLocal: local variable, SLocal <- S \cup \eta{e_j^{d}}
  # VolumeSimplex: Volume of each simplex
  
  # initialize l, L, LRow, SLocal, VolumeSimplex
  L = matrix(0, nrow = 100, ncol = vertexNumber)
  l = matrix(0, nrow = 1, ncol = vertexNumber)
  LRow <- 0; SLocal <- S; VolumeSimplex <- 0 
  
  for(k in 1: m){
    # l = F_k \cup last, Select a candidate for e_j^d
    for(j in 1 : vertexNumber){ 
      if(F[k,j]==1 & last[j]==1)
        l[j] <- 1
      else
        l[j] <- 0
    }
    
    #Check if the candidate has already been generated
    FLAG3 <- 0
    if(LRow > 0){
      for(i in 1:LRow){
        FLAG2 <- 0
        for(j in 1: vertexNumber){
          if(L[i,j] != l[j]){
            FLAG2 <- 1; break
          }
        }
        if(FLAG2 == 0){
          FLAG3 <- 1; break
        }
      }
    }
    
    # Check if l \notin L and if \eta(l)\notin S
    EtaL <- -1
    if(FLAG3 == 0){ # if l \notin L
      LRow <- LRow + 1; L[LRow,] <- l
      for (i in 1:vertexNumber){
        if(l[i] == 1){
          EtaL <- i; break
        }
      }
      if(EtaL != -1){
        if(S[EtaL] == 0){ # if \eta(l)\notin S
          SLocal <- S; SLocal[EtaL] <- 1
          if(d > 0){
            VolumeSimplex <- VolumeSimplex + VolumeCal(d-1, l, SLocal)
          }
        }
      }
    }
    
    # Compute the volume of simplex
    if (d == 0){
      count <- 0
      Simplex <- matrix(0, nrow = n, ncol = n+1)
      for(i in 1: vertexNumber){
        if (SLocal[i] == 1){
          count <- count + 1
          Simplex[ ,count] <- Vertex[, i]
        }
      }
      SVector <- matrix(0, nrow = n, ncol = n)
      for(i in 1:n){
        for(j in 1:n){
          SVector[i,j] <- Simplex[i,j]-Simplex[i,n+1]
        }
      }
      VolumeD0 <- abs(det(SVector))
      for(i in 1:n){
        VolumeD0 <- VolumeD0/i
      }
      return(VolumeD0)
      
    }
    
  }
  return(VolumeSimplex)
  
}


Last <- matrix(1, nrow = 1, ncol = vertexNumber)  
S0 <- matrix(0, nrow = 1, ncol = vertexNumber)  
return (VolumeCal(n, Last, S0))

}


