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

\*---------------------------------------------------------------------------*/

// including basic functions and objects needed for OpenFOAM to run
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Incompressible flows are governed by 4 PDEs, 3 momentum equations and 1
// continuity equation. The field variables involved in computation are also
// 4, Ux,Uy,Uz and P. So we have 4 equations and 4 unknowns, but we dont have
// an explicit equation for pressure P, hence the continuity equation is
// modified to involve pressure and thus solved as a pressure correction
// equation.

int main(int argc, char *argv[])
{
    // basic header files that setup the time, mesh and case are included
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

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

    // Begin the outer loops.
    // The start and end iterations are read from the controlDict file in the
    // system directory. Since this is a steadyState solver, the timeStep
    // value are just used as iteration count, hence it is usually started on
    // a round number.
    while (runTime.loop())
    {

        Info << nl << "Iteration: " << runTime.timeName() << endl;

        // Define the momentum equations as:
        // <convection term> - <diffusion term> == - <pressure gradient>
        // Body forces, turbulence and the unsteady term are neglected (low
        // Reynolds number assumption, steady-state solution).
        fvVectorMatrix UEqn
        (
            fvm::div(phi,U) - fvm::laplacian(nu,U) == -fvc::grad(p)
        );

        // Solve momentum equation for the current values of pressure.
        UEqn.solve();

        // The matrix form of momentum equation is M*U = Nab(P). And it is writen
        // in the form of A*U - H = Nab(P) for ease of inversion. Those A and H
        // matrices are received as field values as shown below. Nab stands for
        // the nabla operator.

        // Getting A and H matrices as fields.
        volScalarField A = UEqn.A();
        volVectorField H = UEqn.H();

        // Computing inverse of A matrix for ease of calculation; it is easy
        // as A is a diagonal matrix.
        volScalarField A_inv = 1.0/A;
        // Interpolating it onto grid faces. This is done because of how
        // the laplacian operator works in OpenFOAM.
        surfaceScalarField A_inv_flux = fvc::interpolate(A_inv);
        // Computing HbyA field = H/A for ease of calculation
        volVectorField HbyA = A_inv * H;

        // Forming the pressure correction equation:
        // Nab(A^-1 Nab(p)) = Nab.(A^-1 * H)
        // The LHS can be defined using the laplacian operator in OpenFOAM as:
        fvScalarMatrix pEqn
        (
            fvm::laplacian(A_inv_flux, p) == fvc::div(HbyA)
        );

        // Setting reference pressure for the equation.
        pEqn.setReference(pRefCell, pRefValue);

        // Solving the pressure correction equation.
        pEqn.solve();

        // Under-relaxing the pressure equation using explicit relaxation:
        p = alpha*p + (1.0 - alpha)*p_old;

        // Updating the velocity field with newly computed pressure field.
        U = A_inv * H - A_inv * fvc::grad(p);
        // Updating the flux field with newly updated velocity field.
        phi = fvc::interpolate(U) & mesh.Sf();

        // Updating boundary conditions for both p and U fields.
        U.correctBoundaryConditions();
        p.correctBoundaryConditions();

        // Updating old_pressure field with new values
        p_old = p;

        // Writing computed fields at the intervals insisted by controlDict.
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
