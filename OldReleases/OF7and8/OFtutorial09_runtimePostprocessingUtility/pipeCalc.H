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

Description
    This function object calculates and outputs various quantities relevant
    to flows through pipes and similar geometries.

SourceFiles
    pipeCalc.C

\*---------------------------------------------------------------------------*/

#ifndef pipeCalc_H
#define pipeCalc_H

#include "volFieldsFwd.H"
#include "Switch.H"
#include "fvc.H"
#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class pipeCalc Declaration
\*---------------------------------------------------------------------------*/

class pipeCalc
:
    // NOTE: derive from the base fvMeshFunctionObject and logFiles in order to
    // provide a ready interface for dealing with output files
    public fvMeshFunctionObject,
    public logFiles
{
private:
    // Private data

        //- Name of this set of pipeCalc objects
        word name_;

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
    // NOTE: define a list of files this object writes to; the underlying logFiles
    // object will handle output to correct streams automatically.
    enum fileID
    {
        MAIN_FILE = 0
    };

    // NOTE: Create file names for each of the output files
    wordList createFileNames(const dictionary& dict) const;

    // NOTE: this first gets declared in the baseline object and gets used to
    // write the header in the output file.
    virtual void writeFileHeader(const label i);

public:

    //- Runtime type information
    TypeName("pipeCalc");

    // Constructors

        //- Construct for given Time and dictionary.
        pipeCalc
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

    //- Destructor
    virtual ~pipeCalc();

    // Member Functions

    //- Return name of the set of pipeCalc
    virtual const word& name() const { return name_; }

    //- Read the pipeCalc data
    virtual bool read(const dictionary&);

    //- Execute, currently does nothing
    virtual bool execute();

    //- Execute at the final time-loop, currently does nothing
    virtual bool end();

    //- Called when time was set at the end of the Time::operator++
    virtual void timeSet();

    // NOTE: Key method which implements the object's actual functionality
    virtual bool write();

    //- Update for changes of mesh
    virtual void updateMesh(const mapPolyMesh&) {}

    //- Update for changes of mesh
    virtual void movePoints(const polyMesh&) {}
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjectFile
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
