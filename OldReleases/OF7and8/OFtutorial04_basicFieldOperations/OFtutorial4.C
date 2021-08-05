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

// This is a function declaration; this method will calculate some scalar value
// given the current time, location in space x, and a reference point x0. The
// function also accepts a scaling factor, scale.
// The actual implementation, or definition, is below.
scalar calculatePressure(scalar t, vector x, vector x0, scalar scale);

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

	// This reads a dictionary file.
	Info << "Reading transportProperties\n" << endl;

	IOdictionary transportProperties
	(
		IOobject
		(
		    "transportProperties", // name of the dictionary
		    runTime.constant(), // location in the case - this one is in constant
		    mesh, // needs the mesh object reference to do some voodoo - unimportant now
		    IOobject::MUST_READ_IF_MODIFIED, // the file will be re-read if it gets modified during time stepping
		    IOobject::NO_WRITE // read-only
		)
	);

	// Create a scalar constant for kinematic viscosity by reading the value from the dictionary.
	dimensionedScalar nu
	(
		"nu", // name of the variable
		dimViscosity, // dimensions
		// TIP: to check how this is defined, run:
		// grep -r dimViscosity $FOAM_SRC/OpenFOAM/
		// This returns:
		/*/opt/openfoam30/src/OpenFOAM/dimensionSet/dimensionSets.C:const dimensionSet dimViscosity(dimArea/dimTime);
		/opt/openfoam30/src/OpenFOAM/dimensionSet/dimensionSets.C:const dimensionSet dimDynamicViscosity(dimDensity*dimViscosity);
		/opt/openfoam30/src/OpenFOAM/dimensionSet/dimensionSets.H:extern const dimensionSet dimViscosity;*/
		// So, it becomes apparent we should check dimensionSets.C, which contain:
		/*const dimensionSet dimLength(0, 1, 0, 0, 0, 0, 0);
		const dimensionSet dimTime(0, 0, 1, 0, 0, 0, 0);
		const dimensionSet dimArea(sqr(dimLength));
		const dimensionSet dimViscosity(dimArea/dimTime);*/
		// This is what gets used here. But, an alternative would be to type in the units directly:
		// dimensionSet(0,2,-1,0,0,0,0),
		transportProperties.lookup("nu") // this takes the value from the dictionary and returns it, passing it to the object constructor as an argument
	);

	// These read the fields p and U from the time folders, as specified in system/controlDict (i.e. latestTime, startTime, etc.)
	Info<< "Reading field p\n" << endl;
	volScalarField p // note that pressure is a scalar field
	(
		IOobject
		(
		    "p", // name of the field
		    runTime.timeName(), // name of the current time, i.e. the time folder to read from
		    mesh,
		    IOobject::MUST_READ, // always gets imported, will throw an error if the field is missing
		    IOobject::AUTO_WRITE // will get saved automatically when the controlDict parameters will request it
		),
		mesh // initialises the field to match the size of the mesh with default (0) values
	);

	Info<< "Reading field U\n" << endl;
	volVectorField U // note that velocity is a vector field
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

	// Let us define a vector whose values will not change over the course of the program execution.
	const vector originVector(0.05,0.05,0.005);

	// Calculate the distance from the origin to the cell centre furthest away.
	// In Python, this is equivalent to:
	// np.sqrt(np.sum((x0-x)**2))
	// The .value() method is called to convert from a dimensionedScalar to a regular scalar.
	const scalar rFarCell = max( // find the maximum value from all distances
		// compute distance of each cell centre from x0; units of mesh.C() are those of length, as this field
		// describes position in the Cartesian reference frame.
		mag(dimensionedVector("x0",dimLength,originVector)-mesh.C())
		).value(); // convert to dim-less scalar

	// This part of the code performs time stepping for as long as is required by the simulation.
	Info<< "\nStarting time loop\n" << endl;

	// This will increment the current time automatically
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

		// Loop over all cells in the mesh and calculate the pressure value.
		for (label cellI=0; cellI<mesh.C().size(); cellI++)
		{
			// cellI describes a series of integers, each corresponding to an index of an individual cell in the grid.

			// Call the method and compute p.
			// Note how mesh.C() and p elements are accessed by using the [] operators, as in a regular C array.
			// .value() is also called to convert the time to a dim-less scalar
			p[cellI] = calculatePressure(runTime.time().value(), mesh.C()[cellI], originVector, rFarCell);

            // NOTE: it is also possbile to interact with boundary face values, but
            // this will be addressed in a separate tutorial.
		}

		// Calculate the gradient of p and substitute for U. Note that the units of grad(p) are not m/s,
		// so we multiply by a unit time variable to make them so. This is just for illustration purposes, of course.
		// The result will point either towards or away from the origin, depending on the sign of the
		// time varying "pressure".
		U = fvc::grad(p)*dimensionedScalar("tmp",dimTime,1.);

		// If requested by controlDict, save the fields to files.
		runTime.write();

		// NOTE: a more appropriate way to calculate things in OpenFOAM is through performing
		// operations on field objects and not iterating cell-by-cell, where possible.
		// How to do this has been shown above, where rFarCell is being computed.
		// The iterative approach has been presented for completeness and to illustrate certain
		// basic features, but is, generally, discouraged, unless absolutely necessary.
	}

	Info << "Finished! Best to visualise the results by plotting p iso-contours with range (-10,10) and applying a glyph filter to the U field in Paraview." << endl;

    Info<< "End\n" << endl;

    return 0;
}

// definition of the custom function
scalar calculatePressure(scalar t, vector x, vector x0, scalar scale)
{
	// Calculates the distance between the base point and x, which is given by the magnitude (hence the name mag)
	// of the vector between x0 and x. The value is scaled by the passed factor, with the intention of making
	// it vary between 0 and 1.
	scalar r (mag(x-x0)/scale);

	// Calculate the inverse of r and apply a limiter to avoid dividing by zero.
	scalar rR (1./(r+1e-12));

	// definition of a frequency
	scalar f (1.);

	// Return a sinusoidally varying pressure with maximum at x0.
	// Note how we call the OpenFOAM sin method by referring to the Foam namespace.
	// This is done to differentiate between the native C++ implementation of a method with the same name
	// and thus avoid an ambiguous expression.
	return Foam::sin(2.*Foam::constant::mathematical::pi*f*t)*rR;
}

// ************************************************************************* //
