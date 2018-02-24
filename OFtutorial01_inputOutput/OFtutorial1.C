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
    // Initialise OF case
    #include "setRootCase.H"

    // These two create the time system (instance called runTime) and fvMesh (instance called mesh).
    #include "createTime.H"
    #include "createMesh.H"

    // ---
    // Get access to a custom dictionary
    dictionary customDict;
    const word dictName("customProperties");

    // Create and input-output object - this holds the path to the dict and its name
    IOobject dictIO
    (
        dictName, // name of the file
        mesh.time().constant(), // path to where the file is
        mesh, // reference to the mesh needed by the constructor
        IOobject::MUST_READ // indicate that reading this dictionary is compulsory
    );

    // Check the if the dictionary is present and follows the OF format
    if (!dictIO.typeHeaderOk<dictionary>(true))
        FatalErrorIn(args.executable()) << "Cannot open specified refinement dictionary "
            << dictName << exit(FatalError);

    // Initialise the dictionary object
    customDict = IOdictionary(dictIO);

    // ---
    // Read various pieces of information from the main part of the dictionary

    // Lookup which does not need to be told what type of variable we're looking for and
    // uses the standard C++ stringstream syntax
    word someWord;
    customDict.lookup("someWord") >> someWord;

    // This template method needs to know the type of the variable and can provide
    // a default value if the entry is not found in the dictionary
    scalar someScalar( customDict.lookupOrDefault<scalar>("someScalar", 1.0) );

    // A switch is a neat feature allowing boolean values to be read from a dict,
    // it supports the OpenFOAM yes/on/true/1 and no/off/false/0 values automatically.
    bool someBool ( customDict.lookupOrDefault<Switch>("someBool",true) );

    // Lists of values may also be read in the same way
    List<scalar> someList ( customDict.lookup("someList") );

    // This type of container is particularly interesting - it associates entries with
    // given key values (here of word type but can be anything); useful when
    // associating things by indices in a list is less handy
    HashTable<vector,word> someHashTable ( customDict.lookup("someHashTable") );

    // Summarise what's been read and print in the console
    Info << nl << "Read the following:" << nl << nl
         << "someWord " << someWord << nl << nl
         << "someScalar " << someScalar << nl << nl
         << "someList " << someList << nl << nl
         << "someHashTable " << someHashTable << nl << nl
         << "someBool " << someBool << nl << nl
         << endl;

    // ---
    // Create a custom directory and write an output file

    // Create the output path directory
    fileName outputDir = mesh.time().path()/"postProcessing";
    // Creathe the directory
    mkDir(outputDir);

    // File pointer to direct the output to
	autoPtr<OFstream> outputFilePtr;
    // Open the file in the newly created directory
    outputFilePtr.reset(new OFstream(outputDir/"customOutputFile.dat"));

    // Write stuff
    outputFilePtr() << "# This is a header" << endl;
    outputFilePtr() << "0 1 2 3 4 5" << endl;

    // Append to the imported hash table and wirte it too
    someHashTable.insert("newKey", vector(1., 0., 0.));
    outputFilePtr() << someHashTable << endl;

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
