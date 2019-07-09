/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
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
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
	// Set up the case, parse command line options and create the grid
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    
    // Create the scalar field and read BCs and the initial conditions
    // NOTE: beta is thus already subjects to the BCs specified in 0/beta
    Info << "Reading field beta" << nl << endl;
    volScalarField beta
    (
        IOobject
        (
            "beta",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    // Read the constant velocity field
    Info << "Reading field U" << nl << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    // Read transport properties and get the diffusion constant
    Info << "Reading transportProperties\n" << endl;
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    Info<< "Reading diffusivity\n" << endl;
    dimensionedScalar gamma (transportProperties.lookup("gamma"));
    
    // Create the flux field
    // NOTE: typically this is done by including createPhi.H from $FOAM_SRC/finiteVolume/cfdTools/incompressible
    Info << "Reading/calculating face flux field phi" << nl << endl;
	surfaceScalarField phi
	(
		IOobject
		(
		    "phi",
		    runTime.timeName(),
		    mesh,
		    IOobject::READ_IF_PRESENT,
		    IOobject::AUTO_WRITE
		),
		// Interpolates U onto the faces and does a dot product with the face area vectors
		// Yields a scalar representing rate of change of volume through each face, i.e. the flux
		// NOTE: the original implementation uses linearInterpolate(U); changed here to fvc::interpolate(U)
		// to show how the method searches system/fvSchemes for an interpolate(U) entry which allows
		// a different scheme to be chosen.
		fvc::interpolate(U) & mesh.Sf() // [(m s-1) * (m2) = (m3 s-1)] <=> flow rate
	);
	
	// Solve the steady scalar transport equation using the solver specified in the system/fvSolution dict.
	// Discretisation of the individual terms is specified in system/fvSchemes.
	// Boundary conditions form part of the beta field already, since it's been read from the file,
	// and thus do not need to be explicitly stated here - this keeps the syntax general.
	solve
    (
    	// Convective term - advection of beta due to the velocity field
		fvm::div(phi, beta) // [(m-1) * (m s-1) * (kg m-3) = (kg m-3 s-1)] <=> flux of beta
		// Diffusive term - diffusion of beta due to its own gradient and a proportionality constant gamma
		- fvm::laplacian(gamma, beta) // [(m2 s-1) * (m-2) * (kg m-3) = (kg m-3 s-1)] <=> flux of beta
		// NOTE: to apply an explicit source term, use the following:
		//== SourceTerm
		// NOTE: to make the source term implicit, use:
		//== Sp(SourceTerm)
		// NOTE: to use the OpenFOAM interface for applying arbitrary source terms from a dictionary, use:
		//== fvOptions(beta)
    );
    
    // Save the result under a different name - we don't do any time stepping so the result ends up
    // in the same folder as the initial conditions which we don't want to overwrite.
    volScalarField result
    (
        IOobject
        (
            "result",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        beta // copy all beta contents, including the BCs
    );
	result.write(); // force output

    Info << nl << "End" << nl << endl;

    return 0;
}

