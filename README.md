# Calculating_Polytope_Volume
Motivation: given the polytope in the form of linear equalities, computing its volume exactly or approximately.

## Triangulation Method ([Cohen and Hickey, 1979](http://delivery.acm.org/10.1145/330000/322141/p401-cohen.pdf?ip=115.154.21.23&id=322141&acc=ACTIVE%20SERVICE&key=BF85BBA5741FDC6E%2EC4BFCDFF40C7237A%2E4D4702B0C3E38B35%2E4D4702B0C3E38B35&__acm__=1530312975_ac85e063731769438044dc35e9cc35f4))
Main idea: dividing polytope into simplices and sum up the volumes of the resulting simplices. 

## Iterative Method (Lasserre, 1983)

## Random-Walk Method ([Ioannis Z. Emiris and Vissarion Fisikopoulos, 2018](http://delivery.acm.org/10.1145/3200000/3194656/a38-emiris.pdf?ip=115.154.21.23&id=3194656&acc=ACTIVE%20SERVICE&key=BF85BBA5741FDC6E%2EC4BFCDFF40C7237A%2E4D4702B0C3E38B35%2E4D4702B0C3E38B35&__acm__=1530312798_1acf8ff1fe0b3b313e61f9912037b75a))
Main procedures:
  *  Hit-and-run random walk, which generates a uniform distribution of points in convex body
  *  Rounding and sandwiching polytope, which is finding the smallest enclosed ellipsoid
  *  Multiphase Monte Carlo: approximating polytope volume by telescopic product

