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

// Include the code for the custom classes declared in .H and defined
// in .C files.
// NOTE: check how the Make/files changed to make sure the additional code gets
// compiled before the main utility.
#include "customClass.H"
#include "derivedClass.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // Create a custom class instance
    customClass customInstance;
    Info << "Default value " << customInstance.get() << endl;

    // Set a new value
    customInstance.set(10);
    Info << "New value " << customInstance.get() << endl;

    // Call a basic function
    customInstance.basicFunction();

    // Pass a reference to the mesh to the custom class and let it do its things
    customInstance.meshOpFunction(mesh);
    Info << "Yet another value " << customInstance.get() << endl;

    // Now, create an instance of a derived class which inherits from an IOdictionary object
    myDict myTransportProperties
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

    // Create a viscosity scalar using our new class
    dimensionedScalar nu
    (
        "nu",
        dimViscosity,
        myTransportProperties.lookup("nu")
    );
    Info << "Created a viscosity scalar: " << nu << endl;

    // List the contents of the dictionary using the derived class method
    // implemented specifically for this purpose
    myTransportProperties.printTokensInTheDict();

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
