# FastSweepReinitalisation
Danny van der Haven, dannyvdhaven@gmail.com,
last updated: 2022/01/24

This MATLAB code uses the Fast Sweeping method [1] to 
reinitialise a level set or signed-distance function.

Given a discrete function F that gives the distance to 
a surface at F = 0, this algorithm detects the grid 
points adjacent to the interface and then recomputes
the function F over the entire domain.

[1] Adam Chacon and Alexander Vladimirsk,
 SIAM Journal on Scientific computing, 2012
