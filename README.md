# Calculating_Polytope_Volume
Calculating the volume of the polytopes by several methods.

Triangulation Method based on Cohen and Hickey (1979).
Main idea: dividing polytope into simplices and sum up the volumes of the resulting simplices. 

Random-Walk Methods for approximating polytope volume base on Ioannis Z. Emiris and Vissarion Fisikopoulos(2014).
The main procedures are:
  1)  Hit-and-run random walk, which generates a uniform distribution of points in convex body
  2)  Rounding and sandwiching polytope, which is finding the smallest enclosed ellipsoid
  3) Multiphase Monte Carlo: approximating polytope volume by telescopic product
