/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

Class
    Foam::pipeCalc

Group
    grpUtilitiesFunctionObjects

Description
    This function object calculates and outputs various quantities relevant
    to flows through pipes and similar geometries.

SourceFiles
    pipeCalc.C
    IOpipeCalc.H

\*---------------------------------------------------------------------------*/

#ifndef pipeCalc_H
#define pipeCalc_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "OFstream.H"
#include "Switch.H"
#include "volFields.H"
#include "dictionary.H"
#include "fvc.H"
#include "functionObjectFile.H"
#include "writer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class polyMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                          Class pipeCalc Declaration
\*---------------------------------------------------------------------------*/

class pipeCalc
:
    // NOTE: derive from the base functionObjectFile in order to provide a ready
    // interface for dealing with output files
    public functionObjectFile
{
private:
    // Private data

        //- Name of this set of pipeCalc objects
        word name_;

        //- Reference to the database
        const objectRegistry& obr_;

        // reference to the mesh
        const fvMesh& mesh_;

        //- On/off switch
        bool active_;

        //- Name of velocity field, default is "U"
        word UName_;

        // NOTE: these fields hold the information about the faceZone
        // Name of the face zone
        word faceZoneName_;

        // index of the face zone in the mesh container
		label faceZoneLabel_;

        // indices of the faces of the CS in the mesh - keep as reference, should be faster for large sets
        const labelList& faces_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        pipeCalc(const pipeCalc&);

        //- Disallow default bitwise assignment
        void operator=(const pipeCalc&);

protected:
    // NOTE: this first gets declared in the baseline object and gets used to
    // write the header in the output file.
    virtual void writeFileHeader(const label i);

public:

    //- Runtime type information
    TypeName("pipeCalc");

    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        pipeCalc
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );

    //- Destructor
    virtual ~pipeCalc();

    // Member Functions

    //- Return name of the set of pipeCalc
    virtual const word& name() const { return name_; }

    //- Read the pipeCalc data
    virtual void read(const dictionary&);

    //- Execute, currently does nothing
    virtual void execute();

    //- Execute at the final time-loop, currently does nothing
    virtual void end();

    //- Called when time was set at the end of the Time::operator++
    virtual void timeSet();

    //- Calculate the pipeCalc and write
    virtual void write();

    //- Update for changes of mesh
    virtual void updateMesh(const mapPolyMesh&) {}

    //- Update for changes of mesh
    virtual void movePoints(const polyMesh&) {}
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
