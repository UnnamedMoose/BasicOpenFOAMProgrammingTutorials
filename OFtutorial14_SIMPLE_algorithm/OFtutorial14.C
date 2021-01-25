/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 2021 AUTHOR,AFFILIATION
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    tutSimpleFoam

Description
    This tutorial is about explaining how the SIMPLE algorithm is implemented
    through coding in OpenFOAM solvers like "simpleFoam". This tutorial is
    developed on OpenFOAM v2006 version and likely to work on v19 and above
    versions

    SIMPLE algorithm stands for "Semi-Implicit Method for Pressure Linked
    Equations". it is one of the algorithms used in solving incompressible
    Navier-Stokes equations.

    this tutorial requires the following prerequisite
    1) basic knowledge in computational fluid dynamics, specifically to
        incompressible flow solutions.
    2) basic knowledge in Finite Volume Methods in discretizing partial
        differential equations.
    3) it is HIGHLY RECOMMENDED to go through this video on youtube

        https://www.youtube.com/watch?v=ahdW5TKacok

        it explains about Matrix notations/operations and the same is
        used in this tutorial.

\*---------------------------------------------------------------------------*/

// including basic functions and objects needed for OpenFOAM to run
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// incompressible flows are governed by 4 PDEs, 3 momentum equations and 1
// continuity equation. The field variables involved in computation are also
// 4, Ux,Uy,Uz and P. So we have 4 equations and 4 unknowns, but we dont have
// an explicit equation for pressure P, hence the continuity equation is
// modified to involve pressure and thus solved it as pressure correction
// equation. the details of derivation and matrix notations are explained well
// on the youtube video just shared above.

// main function declaration

int main(int argc, char *argv[])
{
    // basic header files that setup the time, mesh and case are included
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // declaring, reading and defining field values
    // the details on how the fields are defined/read/created can be seen
    // on the "createFields.H" file in the same directory
    #include "createFields.H"

    // the following factors are read from the fvSolution dictionary defined
    // in the "createFields.H" file. this file is present in the system/
    // directory of case folder

    // reading relaxation factor
    scalar alpha;
    fvSolution.lookup("alpha") >> alpha;
    // reading index of cell containing reference pressure
    scalar pRefCell;
    fvSolution.lookup("pRefCell") >> pRefCell;
    // reading reference pressure value
    scalar pRefValue;
    fvSolution.lookup("pRefValue") >> pRefValue;

    // the read values are printed on to screen for confirmation during run
    Info << nl << "following parameters read:" << endl;
    Info << tab << "relaxation factor \"alpha\" : " << alpha << endl;
    Info << tab << "index of cell containing reference pressure \"pRefCell\" : " << pRefCell << endl;
    Info << tab << "reference pressure value \"pRefValue\" : " << pRefValue << endl;

    // begining loop
    // the start and end times are read from the controlDict file in the
    // system directory. Since this is a steadyState solver, the timeStep
    // value are just used as iteration count, hence it is usually started on
    // a round number, it can be checked on the controlDict file of case directory
    while (runTime.loop())
    {

        Info << nl << "Iteration: " << runTime.timeName() << endl;

        // defining momentum equation
        // the momentum equation is defined in the format as
        // <convection term> - <diffusion term> == - <pressure gradient>
        fvVectorMatrix UEqn
        (
            fvm::div(phi,U) - fvm::laplacian(nu,U) == -fvc::grad(p)
        );

        // solving momentum equation
        UEqn.solve();

        // as can be seen from the youtube video, the matrix form of momentum
        // equation is M*U = Nab(P). And it is writen in the form of
        // A*U - H = Nab(P) for ease of inversion. Please see the video for
        // clear understanding. Those A and H matrices are received as field
        // values as shown below. Nab - stands for nabla operator

        // getting A and H matrices as field values
        volScalarField A = UEqn.A();
        volVectorField H = UEqn.H();

        // computing inverse of A matrix for ease of calculation; it is easy
        // as the A is a diagonal matrix.
        volScalarField A_inv = 1.0/A;
        // and interpolating it to surface field. this is done for the way
        // the laplacian operator works on OpenFOAM.
        surfaceScalarField A_inv_flux = fvc::interpolate(A_inv);
        // and computing HbyA field = H/A for ease of calculation
        volVectorField HbyA = A_inv * H;

        // framing pressure correction equation
        // this equation can be seen on the reference youtube video as
        // Nab(A^-1 Nab(p)) = Nab.(A^-1 * H)
        // the LHS can be defined using laplacian operator in OpenFOAM as below
        fvScalarMatrix pEqn
        (
            fvm::laplacian(A_inv_flux, p) == fvc::div(HbyA)
        );

        // setting reference pressure for equation
        pEqn.setReference(pRefCell,pRefValue);

        // solving pressure equation
        pEqn.solve();

        // under-relaxing pressure equation
        // two methods of relaxation was shown in the video, in that, the 2nd
        // one is implemented in here as to make things look simple
        p = alpha*p + (1 - alpha)*p_old;

        // updating velocity field with newly computed pressure field
        U = A_inv * H - A_inv * fvc::grad(p);
        // updating flux field with newly updated velocity field
        phi = fvc::interpolate(U) & mesh.Sf();

        // updating boundary conditions for both p and U fields
        U.correctBoundaryConditions();
        p.correctBoundaryConditions();

        // updating old_pressure field with new values
        p_old = p;

        // writing computed fields at the intervals insisted by controlDict
        runTime.write();

    }

    // printing execution time information at the end of simulation
    Info<< nl;
    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << endl;

    Info<< "End\n" << endl;

    return 0;
}

// P.S.: The main motive of this tutorial is to show how things work inside
// OpenFOAM solver. Hence the optimization or convergence check are ignored in
// this code. The user can check "icoFoam" solver after understanding this code
// there these optimizations are included and it will be easier to understand
// that.

// ************************************************************************* //
