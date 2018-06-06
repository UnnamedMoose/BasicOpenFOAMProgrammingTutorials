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

#include "cellModeller.H"

int main(int argc, char *argv[])
{
    // set up the case
    #include "setRootCase.H"

    // create the run time object
    Info<< "Create time\n" << endl;
    Time runTime
    (
        Time::controlDictName,
        args.rootPath(),
        args.caseName()
    );

    // disable post-processing etc.
    runTime.functionObjects().off();

    // ---
    // create cell types
    const cellModel& hex = *(cellModeller::lookup("hex"));
    const cellModel& prism = *(cellModeller::lookup("prism"));
    const cellModel& pyr = *(cellModeller::lookup("pyr"));
    const cellModel& tet = *(cellModeller::lookup("tet"));
    Info << "Created cell types" << endl;

    // create a point cloud between which the cells will stretch
    pointField points(12);
    points[ 0] = vector(0.0, 0, 0);
    points[ 1] = vector(0.5, 0, 0);
    points[ 2] = vector(1.0, 0, 0);
    points[ 3] = vector(0.0, 1, 0);
    points[ 4] = vector(0.5, 1, 0);
    points[ 5] = vector(1.0, 1, 0);
    points[ 6] = vector(0.0, 0, 0.1);
    points[ 7] = vector(0.5, 0, 0.1);
    points[ 8] = vector(1.0, 0, 0.1);
    points[ 9] = vector(0.0, 1, 0.1);
    points[10] = vector(0.5, 1, 0.1);
    points[11] = vector(1.0, 1, 0.1);
    Info << "Created points" << endl;

    // create the cells from a particular cell shape and list of points
    // TODO check if point ordering is important or not
    List<cellShape> cells;

    List<label> cellPoints(8);

    cellPoints[0] = 0;
    cellPoints[1] = 3;
    cellPoints[2] = 9;
    cellPoints[3] = 6;
    cellPoints[4] = 1;
    cellPoints[5] = 4;
    cellPoints[6] = 10;
    cellPoints[7] = 7;
    cells.append(cellShape(hex, cellPoints));
    cellPoints[0] = 1;
    cellPoints[1] = 4;
    cellPoints[2] = 10;
    cellPoints[3] = 7;
    cellPoints[4] = 2;
    cellPoints[5] = 5;
    cellPoints[6] = 11;
    cellPoints[7] = 8;
    cells.append(cellShape(hex, cellPoints));
    Info << "Created cells" << endl;

    // create patch definitions from lists of points defining individual faces
    // this is a list of List<List<label> > objects; each lowest-level list
    // contains indices of vertices making up the patch; these are grouped by
    // boundary index, i.e. patch number. The overall list holds all of those in
    // a single container.
    faceListList patchFaces;

    List<label> patchPoints(4);

    patchPoints[0] = 3;
    patchPoints[1] = 9;
    patchPoints[2] = 10;
    patchPoints[3] = 4;
    patchFaces.append(List<face>(1, face(patchPoints)));
    Info << "Created patches" << endl;

    // patch names
    List<word> boundaryPatchNames;
    boundaryPatchNames.append("side0");

    // types and physical types
    // TODO figure out how these are supposed to look really, now just use defaults
    wordList boundaryPatchTypes(patchFaces.size(), polyPatch::typeName);
    wordList boundaryPatchPhysicalTypes
    (
        patchFaces.size(),
        polyPatch::typeName
    );

    // default values for other fields
    word regionName = polyMesh::defaultRegion;
    word defaultFacesName = "defaultFaces";
    word defaultFacesType = emptyPolyPatch::typeName;
    Info << "Created physical and default values" << endl;

    // create the mesh
    polyMesh mesh
    (
        IOobject
        (
            regionName,
            runTime.constant(),
            runTime
        ),
        xferMove(points),
        cells,
        patchFaces,
        boundaryPatchNames,
        boundaryPatchTypes,
        defaultFacesName,
        defaultFacesType,
        boundaryPatchPhysicalTypes
    );
    Info << "Created the mesh object" << endl;

    // ---
    // Write the grid
    Info << nl << "Writing extruded mesh to time = " << runTime.timeName() << nl << endl;
    mesh.write();

/*
	// These two create the time system (instance called runTime) and fvMesh (instance called mesh).
    #include "createTime.H"
    #include "createMesh.H"

	// runTime and mesh are instances of objects (or classes).
	// If you are not familiar with what a class or object is, it is HIGHLY RECOMMENDED you visit this
	// website and only come back once you've read everything about classes, inheritance and polymorphism:
	// http://www.cplusplus.com/doc/tutorial/classes/
	// Note how the next lines call functions .timeName(), .C() and .Cf() implemented in the objects.
	// It is also important to realise that mesh.C() and .Cf() return vector fields denoting centres of each
    // cell and internal face.
	// Calling the mesh.C().size() method therefore yields the total size of the mesh.
	Info << "Hello there, the most recent time folder found is " << runTime.timeName() << nl
		 << "The mesh has " << mesh.C().size() << " cells and " << mesh.Cf().size()
         << " internal faces in it. Wubalubadubdub!" << nl << endl;

    // It's possible to iterate over every cell in a standard C++ for loop
    for (label cellI = 0; cellI < mesh.C().size(); cellI++)
        if (cellI%20 == 0) // only show every twentieth cell not to spam the screen too much
            Info << "Cell " << cellI << " with centre at " << mesh.C()[cellI] << endl;
    Info << endl; // spacer

    // Each cell is constructed of faces - these may either be internal or constitute a
    // boundary, or a patch in OpenFOAM terms; internal faces have an owner cell
    // and a neighbour.
    for (label faceI = 0; faceI < mesh.owner().size(); faceI++)
        if (faceI%40 == 0)
            Info << "Internal face " << faceI << " with centre at " << mesh.Cf()[faceI]
                 << " with owner cell " << mesh.owner()[faceI]
                 << " and neighbour " << mesh.neighbour()[faceI] << endl;
    Info << endl;

    // Boundary conditions may be accessed through the boundaryMesh object.
    // In reality, each boundary face is also included in the constant/polyMesh/faces
    // description. But, in that file, the internal faces are defined first.
    // In addition, the constant/polyMesh/boundary file defines the starting faceI
    // indices from which boundary face definitions start.
    // OpenFOAM also provides a macro definition for for loops over all entries
    // in a field or a list, which saves up on the amount of typing.
    forAll(mesh.boundaryMesh(), patchI)
        Info << "Patch " << patchI << ": " << mesh.boundary()[patchI].name() << " with "
             << mesh.boundary()[patchI].Cf().size() << " faces. Starts at total face "
             << mesh.boundary()[patchI].start() << endl;
    Info << endl;

    // Faces adjacent to boundaries may be accessed as follows.
    // Also, a useful thing to know about a face is its normal vector and face area.
    label patchFaceI(0);
    forAll(mesh.boundaryMesh(), patchI)
        Info << "Patch " << patchI << " has its face " << patchFaceI << " adjacent to cell "
             << mesh.boundary()[patchI].patch().faceCells()[patchFaceI]
             << ". It has normal vector " << mesh.boundary()[patchI].Sf()[patchFaceI]
             << " and surface area " << mag(mesh.boundary()[patchI].Sf()[patchFaceI])
             << endl;
    Info << endl;

    // For internal faces, method .Sf() can be called directly on the mesh instance.
    // Moreover, there is a shorthand method .magSf() which returns the surface area
    // as a scalar.
    // For internal faces, the normal vector points from the owner to the neighbour
    // and the owner has a smaller cellI index than the neighbour. For boundary faces,
    // the normals always point outside of the domain (they have "imaginary" neighbours
    // which do not exist).

    // It is possible to look at the points making up each face in more detail.
    // First, we define a few shorthands by getting references to the respective
    // objects in the mesh. These are defined as constants since we do not aim to
    // alter the mesh in any way.
    // NOTE: these lists refer to the physical definition of the mesh and thus
    // include boundary faces. Use can be made of the mesh.boundary()[patchI].Cf().size()
    // and mesh.boundary()[patchI].start() methods to check whether the face is internal
    // or lies on a boundary.
    const faceList& fcs = mesh.faces();
    const List<point>& pts = mesh.points();
    const List<point>& cents = mesh.faceCentres();

    forAll(fcs,faceI)
        if (faceI%80==0)
        {
            if (faceI<mesh.Cf().size())
                Info << "Internal face ";
            else
            {
                forAll(mesh.boundary(),patchI)
                    if ((mesh.boundary()[patchI].start()<= faceI) &&
                        (faceI < mesh.boundary()[patchI].start()+mesh.boundary()[patchI].Cf().size()))
                    {
                        Info << "Face on patch " << patchI << ", faceI ";
                        break; // exit the forAll loop prematurely
                    }
            }

            Info << faceI << " with centre at " << cents[faceI]
                 << " has " << fcs[faceI].size() << " vertices:";
            forAll(fcs[faceI],vertexI)
                // Note how fcs[faceI] holds the indices of points whose coordinates
                // are stored in the pts list.
                Info << " " << pts[fcs[faceI][vertexI]];
            Info << endl;
        }
    Info << endl;

    // In the original cavity tutorial, on which the test case is based,
    // the frontAndBack boundary is defined as and "empty" type. This is a special
    // BC case which may cause unexpected behaviour as its .Cf() field has size of 0.
    // Type of a patch may be checked to avoid running into this problem if there
    // is a substantial risk that an empty patch type will appear
    label patchID(0);
    const polyPatch& pp = mesh.boundaryMesh()[patchID];
    if (isA<emptyPolyPatch>(pp))
    {
        // patch patchID is of type "empty".
        Info << "You will not see this." << endl;
    }

    // Patches may also be retrieved from the mesh using their name. This could be
    // useful if the user were to refer to a particular patch from a dictionary
    // (like when you do when calculating forces on a particular patch).
    word patchName("movingWall");
    patchID = mesh.boundaryMesh().findPatchID(patchName);
    Info << "Retrieved patch " << patchName << " at index " << patchID << " using its name only." << nl << endl;
*/
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
