/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Application
    tutSimpleFoam createFields

Description
	this file contains necessary field definitions for tutSimpleFoam solver

\*---------------------------------------------------------------------------*/

// creating/reading pressure field
// the pressure field is read from the solution folder set by the runTime,
// for example "0" is the initial folder where the boundary conditions are
// defined, and OpenFOAM by default starts on 0 folder
Info << nl << "Reading pressure field" << endl;
volScalarField p
(
 	IOobject
	(
	 	"p",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
	mesh
);

// creating another pressure field to store previous iteration pressure values
// this is done so that the relaxation factors can be applied for making
// the solution stable.
volScalarField p_old
(
 	IOobject
	(
	 	"p_old",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
    p
);

// creating/reading velocity field
// same as pressure field but velocity is a vector field.
Info << nl << "Reading velocity field" << endl;
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

// creating a flux field
// this field is used in the matrix operations like laplace and divergence
// it is computed by interpolating velocity (volume field) to faces (surface field)
// and multiplying with face area.
Info << nl << "Creating flux field" << endl;
surfaceScalarField phi
(
 	IOobject
	(
	 	"phi",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	fvc::interpolate(U) & mesh.Sf()
);

// declaring dictionary entry
// the transportProperties dictionary is read from constant folder which contains
// the kinematic viscosity "nu"
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

// reading kinematic viscosity from dictionary
dimensionedScalar nu
(
 	"nu",
	dimensionSet(0,2,-1,0,0,0,0),
	transportProperties
);

// declaring fvSolution dictionary entry
// this dictionary is explicity defined to read the reference pressure data
// otherwise, this dictionary is read in the background during case initialization
// for getting data like residual values and number of iterations for each field.
IOdictionary fvSolution
(
    IOobject
    (
        "fvSolution",
        runTime.system(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);
