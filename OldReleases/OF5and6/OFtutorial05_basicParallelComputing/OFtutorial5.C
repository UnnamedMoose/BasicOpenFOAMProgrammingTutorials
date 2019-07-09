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

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // For a case being run in parallel, the domain is decomposed into several
    // processor meshes. Each of them is run in a separate process and holds
    // instances of objects like mesh, U or p just as in a single-threaded (serial)
    // computation. These will have different sizes, of course, as they hold
    // fewer elements than the whole, undecomposed, mesh.
    // Pout is a stream to which each processor can write, unlike Info which only
    // gets used by the head process (processor0)
    Pout << "Hello from processor " << Pstream::myProcNo() << "! I am working on "
         << mesh.C().size() << " cells" << endl;

    // To exchange information between processes, special OpenMPI routines need
    // to be called.

    // This goes over each cell in the subdomain and integrates their volume.
    scalar meshVolume(0.);
    forAll(mesh.V(),cellI)
        meshVolume += mesh.V()[cellI];

    // Add the values from all processes together
    Pout << "Mesh volume on this processor: " << meshVolume << endl;
    reduce(meshVolume, sumOp<scalar>());
    Info << "Total mesh volume on all processors: " << meshVolume
        // Note how the reudction operation may be done in place without defning
        // a temporary variable, where appropriate.
         << " over " << returnReduce(mesh.C().size(), sumOp<label>()) << " cells" << endl;
    // During the reduction stage, different operations may be carried out, summation,
    // described by the sumOp template, being one of them.
    // Other very useful operations are minOp and maxOp.
    // Note how the type
    // of the variable must be added to make an instance of the template, here
    // this is done by adding <scalar> in front of the brackets.
    // Custom reduction operations are easy to implement but need fluency in
    // object-oriented programming in OpenFOAM, so we'll skip this for now.

    // Spreading a value across all processors is done using a scatter operation.
    Pstream::scatter(meshVolume);
    Pout << "Mesh volume on this processor is now " << meshVolume << endl;

    // It is often useful to check the distribution of something across all
    // processors. This may be done using a list, with each element of it
    // being written to by only one processor.
    List<label> nInternalFaces (Pstream::nProcs()), nBoundaries (Pstream::nProcs());
    nInternalFaces[Pstream::myProcNo()] = mesh.Cf().size();
    nBoundaries[Pstream::myProcNo()] = mesh.boundary().size();

    // The list may then be gathered on the head node as
    Pstream::gatherList(nInternalFaces);
    Pstream::gatherList(nBoundaries);
    // Scattering a list is also possbile
    Pstream::scatterList(nInternalFaces);
    Pstream::scatterList(nBoundaries);

    // It can also be useful to do things on the head node only
    // (in this case this is meaningless since we are using Info, which already
    // checks this and executes on the head node).
    // Note how the gathered lists hold information for all processors now.
    if (Pstream::master())
    {
        forAll(nInternalFaces,i)
            Info << "Processor " << i << " has " << nInternalFaces[i]
                 << " internal faces and " << nBoundaries[i] << " boundary patches" << endl;
    }

    // As the mesh is decomposed, interfaces between processors are turned
    // into patches, meaning each subdomain sees a processor boundary as a
    // boundary condition.
    forAll(mesh.boundary(),patchI)
        Pout << "Patch " << patchI << " named " << mesh.boundary()[patchI].name() << endl;

    // When looking for processor patches, it is useful to check their type,
    // similarly to how one can check if a patch is of empty type
    forAll(mesh.boundary(),patchI)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchI];
        if (isA<processorPolyPatch>(pp))
            Pout << "Patch " << patchI << " named " << mesh.boundary()[patchI].name()
                 << " is definitely a processor boundary!" << endl;
    }

    // ---
    // this is an example implementation of the code from tutoral 2 which
    // has been adjusted to run in parallel. Each difference is highlighted
    // as a NOTE.

    // It is conventional in OpenFOAM to move large parts of code to separate
    // .H files to make the code of the solver itself more readable. This is not
    // a standard C++ practice, as header files are normally associated with
    // declarations rather than definitions.
    // A very common include, apart from the setRootCase, createTime, and createMesh,
    // which are generic, is createFields, which is often unique for each solver.
    // Here we've moved all of the parts of the code dealing with setting up the fields
    // and transport constants into this include file.
    #include "createFields.H"

    // pre-calculate geometric information using field expressions rather than
    // cell-by-cell assignment.
	const dimensionedVector originVector("x0", dimLength, vector(0.05,0.05,0.005));
    volScalarField r (mag(mesh.C()-originVector));
    // NOTE: we need to get a global value; convert from dimensionedScalar to scalar
	const scalar rFarCell = returnReduce(max(r).value(), maxOp<scalar>());
    scalar f (1.);

	Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // assign values to the field;
        // sin function expects a dimensionless argument, hence need to convert
        // current time using .value().
        // r has dimensions of length, hence the small value being added to it
        // needs to match that.
        // Finally, the result has to match dimensions of pressure, which are
        // m^2 / s^-2/
		p = Foam::sin(2.*constant::mathematical::pi*f*runTime.time().value())
            / (r/rFarCell + dimensionedScalar("small", dimLength, 1e-12))
            * dimensionedScalar("tmp", dimensionSet(0, 3, -2, 0, 0), 1.);

        // NOTE: this is needed to update the values on the processor boundaries.
        // If this is not done, the gradient operator will get confused around the
        // processor patches.
        p.correctBoundaryConditions();

        // calculate velocity from gradient of pressure
		U = fvc::grad(p)*dimensionedScalar("tmp", dimTime, 1.);
		runTime.write();
	}

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
