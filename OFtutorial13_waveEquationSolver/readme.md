this is a very basic solver that solves wave equation

the main motive of this tutorial is to induce a very basic idea on
how to solve custom equations using OpenFOAM C++.

the tutorial "OFtutorial10\_transportEquation" gives idea on solving a transport
equation on an already solved velocity field. Hence the equation is solved only
once. Whereas in here, the wave equation is solved in transient form.

It has included testCase/ where the solver is used to simulate 2D wave propagation.

this solver needs a new field named "h", which is the amplitude of the wave in
the domain.
