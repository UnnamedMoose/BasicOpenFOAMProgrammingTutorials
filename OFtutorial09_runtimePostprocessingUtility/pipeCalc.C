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
#include "addToRunTimeSelectionTable.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(pipeCalc, 0);
    addToRunTimeSelectionTable(functionObject, pipeCalc, dictionary);
}
}

// * * * * * * * * * * * * * * * * Protected members  * * * * * * * * * * * * * * //

// NOTE: this returns a list of file names which match indices of the enum
// defined in the header of this class. These names are used to create matching
// files by the logFiles object.
Foam::wordList Foam::functionObjects::pipeCalc::createFileNames
(
    const dictionary& dict
) const
{
    DynamicList<word> names(1);

    // use type of the utility as specified in the dict as the top-level dir name
    const word objectType(dict.lookup("type"));

    // Name for file(MAIN_FILE=0)
    names.append(objectType);

    return names;
}

// NOTE: this method first gets declared in logFiles.H, from which this
// class is derived. This method gets called automatically when the base object's
// write() function gets called too.
// The purpose of the function is to add the header to the output data file.
void Foam::functionObjects::pipeCalc::writeFileHeader(const label i)
{
    // Find the correct file to write to from the enum defined in the header.
    switch (fileID(i))
    {
        case MAIN_FILE:
        {
            writeHeader(file(i), "Flow rate through face zone");
            writeHeaderValue(file(i), "Face zone name", faceZoneName_);
            writeCommented(file(i), "Time [s] | Flow rate [m3s-1]");
            file() << endl;
            break; // exit the case structure
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled file index: " << i
                << abort(FatalError);
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::functionObjects::pipeCalc::pipeCalc
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    // NOTE: call the base class constructor
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),

    name_(name),
    active_(true),
    UName_("U"),
    // NOTE: Read the face zone to integrate over. Get its name from the dict, find
    // it in the mesh, and get a reference to the list of its faces.
    faceZoneName_(dict.lookup("faceZoneName")),
    faceZoneLabel_( mesh_.faceZones().findZoneID(faceZoneName_) ),
    faces_( mesh_.faceZones()[faceZoneLabel_] )
{
    // NOTE: calls the separate .read() method to import the controls from the dict.
    // dict reference is passed automatically by the OpenFOAM runtime object manager.
    read(dict);

    // built-in logFiles method for creating file streams.
    resetNames(createFileNames(dict));

    if (active_)
    {
        // Any extra initialisation goes here, if necessary

        // NOTE: type() returns the typeName, as defined in the header file. Name
        // is the individual identifier of this instance, as specified in the dict
        Info << "Finished initialising " << type() << ": " << name_ << nl << endl;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::pipeCalc::~pipeCalc()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::pipeCalc::read(const dictionary& dict)
{
    if (active_)
    {
        UName_ = dict.lookupOrDefault<word>("UName", "U");
    }
    return true;
}

bool Foam::functionObjects::pipeCalc::execute()
{
    if (active_)
    {
        // This gets called before write, should put things on which other
        // function objects might depend on here (for instance field calculations)
    }
    return true;
}

bool Foam::functionObjects::pipeCalc::end()
{
    if (active_)
    {
        execute();
    }
    return true;
}

void Foam::functionObjects::pipeCalc::timeSet()
{}

bool Foam::functionObjects::pipeCalc::write()
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
            logFiles::write();

            // Add the entry for this time step that has just been computed.
            file(MAIN_FILE) << obr_.time().value() << tab << flowRate << endl;
        }
    }
    return true;
}

// ************************************************************************* //
