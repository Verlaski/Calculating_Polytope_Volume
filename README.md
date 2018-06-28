# Calculating_Polytope_Volume
Motivation: given the polytope in the form of linear equalities, computing its volume exactly or approximately.

## Triangulation Method ([Cohen and Hickey, 1979](http://delivery.acm.org/10.1145/330000/322141/p401-cohen.pdf?ip=115.154.21.23&id=322141&acc=ACTIVE%20SERVICE&key=BF85BBA5741FDC6E%2EC4BFCDFF40C7237A%2E4D4702B0C3E38B35%2E4D4702B0C3E38B35&__acm__=1530123849_9da6556e22cd5f6647ed865f8d622d16))
Main idea: dividing polytope into simplices and sum up the volumes of the resulting simplices. 

## Random-Walk Method ([Ioannis Z. Emiris and Vissarion Fisikopoulos, 2014](https://www.cs.bgu.ac.il/~eurocg14/papers/paper_35.pdf))
Main procedures:
  *  Hit-and-run random walk, which generates a uniform distribution of points in convex body
  *  Rounding and sandwiching polytope, which is finding the smallest enclosed ellipsoid
  *  Multiphase Monte Carlo: approximating polytope volume by telescopic product
