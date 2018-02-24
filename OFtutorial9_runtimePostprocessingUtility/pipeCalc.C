/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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

#include "pipeCalc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(pipeCalc, 0);
}

// * * * * * * * * * * * * * * * * Protected members  * * * * * * * * * * * * * * //

// NOTE: this method first gets declared in functionObjectFile.H, from which this
// class is derived. This method gets called automatically when the base object
// function gets called too.
// The purpose of the function is to add the header to the output data file.
void Foam::pipeCalc::writeFileHeader(const label i)
{
        writeHeader(file(), "Flow rate through face zone");
        writeHeaderValue(file(), "Face zone name", faceZoneName_);
        writeCommented(file(), "Time [s] | Flow rate [m3s-1]");
        file() << endl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::pipeCalc::pipeCalc
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    // NOTE: call the base class constructor
    functionObjectFile(obr, name, type()),
    name_(name),
    obr_(obr),
    // NOTE: the instance only gets given a reference to the object registry
    // automaitcally (this keeps track of what objects, such as fields or
    // meshes) are registered in the simulation. The following line casts
    // (converts) the reference to a const fvMesh reference, which allows
    // the grid to be accessed but not modified.
    mesh_(refCast<const fvMesh>(obr_)),
    active_(true),
    UName_("U"),
    // NOTE: Read the face zone to integrate over. Get its name from the dict, find
    // it in the mesh, and get a reference to the list of its faces.
    faceZoneName_(dict.lookup("faceZoneName")),
    faceZoneLabel_( mesh_.faceZones().findZoneID(faceZoneName_) ),
    faces_( mesh_.faceZones()[faceZoneLabel_] )
{
    // Check if the available mesh is an fvMesh, otherwise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "pipeCalc::pipeCalc"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating " << name_ << nl
            << endl;
    }

    // NOTE: calls the separate .read() method to import the controls from the dict.
    // dict reference is passed automatically by the OpenFOAM runtime object manager.
    read(dict);

    if (active_)
    {
        // Any extra initialisation goes here, if necessary

        // NOTE: type() returns the typeName, as defined in the header file. Name
        // is the individual identifier of this instance, as specified in the dict
        Info << "Finished initialising " << type() << ": " << name_ << nl << endl;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pipeCalc::~pipeCalc()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pipeCalc::read(const dictionary& dict)
{
    if (active_)
    {
        UName_ = dict.lookupOrDefault<word>("UName", "U");
    }
}

void Foam::pipeCalc::execute()
{
    if (active_)
    {
        // This gets called before write, should put things on which other
        // function objects might depend on here (for instance field calculations)
    }
}

void Foam::pipeCalc::end()
{
    if (active_)
    {
        execute();
    }
}

void Foam::pipeCalc::timeSet()
{}

void Foam::pipeCalc::write()
{
    if (active_)
    {
        // NOTE: this is the main part of the function object and implements the
        // actual functionality.

        // Retrieve a reference to the velocity field
        const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);

        // itnerpolate onto the faces
		surfaceVectorField Uface = fvc::interpolate(U);

        // Go over each face zone face and compute the flow rate normal to the
        // face zone faces.
        // This assumes none of the face zone faces are on processor boundaries.
        // If they are, a seg fault will occur as faces_[faceI] will exceed the
        // internal field dimensions. A check could be made by using:
        // if ((faces_[faceI] > mesh_.owner().size()) || (mesh_.owner().size() == 0))
        // and the boundary values could be used instead. This is not done here
        // for simplicity.
        scalar flowRate(0.);

        forAll(faces_, faceI)
            // Flow rate = dot product of velocity and surface area vector; in Latex terms,
            // Q = \mathbf{U} \cdot \mathbf{\hat{n}} A
            flowRate += Uface[faces_[faceI]] & mesh_.Sf()[faces_[faceI]];

        // reduce for parallel running
        reduce(flowRate, sumOp<scalar>());

        Info << "Total flow rate " << flowRate << " through "
             << returnReduce(faces_.size(), sumOp<label>()) << " faces" << nl << endl;

        // Output to file - only execute on the master thread to avoid the file
        // getting written into from a few processors at the same time
        if (Pstream::master())
        {
            // Call the base class method which checks if the output file exists
            // and creates it, if necessary. That also calls the .writeFileHeader()
            // method of the derived class.
            functionObjectFile::write();

            // Add the entry for this time step that has just been computed.
            file() << obr_.time().value() << tab << flowRate << endl;
        }
    }
}

// ************************************************************************* //
