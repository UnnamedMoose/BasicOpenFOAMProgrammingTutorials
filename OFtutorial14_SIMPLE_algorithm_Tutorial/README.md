# SIMPLE Algorithm tutorial on OpenFOAM
this tutorial is aimed to show how the matrix operations are performed and the
corresponding governing equations are solved.

As an example, the incompressible Navier-Stokes equations are coded in this
tutorial and SIMPLE algorithm is used for their solution steps.

The following are prerequisites that are required for understanding this
tutorial.
* Basic knowledge in Computational Fluid Dynamics, specially on incompressible flows
* An idea on SIMPLE Algorithm that is used to solve these incompressibel NS equations.
* Please go through this youtube video https://www.youtube.com/watch?v=ahdW5TKacok to understan the basic matrix notations and ideas, the same is implemented in this code. The code will be a lot easier to understand after watching that video.

It is worth noting that, the code is developed to show how the equations are solved only!. Hence the optimization and convergence check conditions are dropped to make the code simpler to understand.

A testCase that solves a 2d channel flow is included along with this code and it can be used to run and check this application after it is compiled.
