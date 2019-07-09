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
    pointField points(17);
    points[ 0] = vector(0.0, 0, 0); // hexes
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
    points[12] = vector(0.7, 0.5, -0.5); // tip of pyramid
    points[13] = vector(0.7, 1.0, -0.5); // last point of a tet
    points[14] = vector(0.7, 1.2, -0.5); // top face of a prism
    points[15] = vector(0.5, 1.2, 0.0);
    points[16] = vector(1.0, 1.2, 0.0);
    Info << "Created points" << endl;

    // create the cells from a particular cell shape and list of points
    // ordering is important
    List<cellShape> cells;

    List<label> cellPoints(8);
    cellPoints[0] = 0; // first face - x=0.0, clock-wise loop looking outwards
    cellPoints[1] = 3;
    cellPoints[2] = 9;
    cellPoints[3] = 6;
    cellPoints[4] = 1; // second face - x=0.5, points 4-7 share edges with 0-3, respectively
    cellPoints[5] = 4;
    cellPoints[6] = 10;
    cellPoints[7] = 7;
    cells.append(cellShape(hex, cellPoints));

    cellPoints[0] = 4;
    cellPoints[1] = 10;
    cellPoints[2] = 7;
    cellPoints[3] = 1;
    cellPoints[4] = 5;
    cellPoints[5] = 11;
    cellPoints[6] = 8;
    cellPoints[7] = 2;
    cells.append(cellShape(hex, cellPoints));

    cellPoints.resize(6);
    cellPoints[0] = 13;
    cellPoints[1] = 4;
    cellPoints[2] = 5;
    cellPoints[3] = 14;
    cellPoints[4] = 15;
    cellPoints[5] = 16;
    cells.append(cellShape(prism, cellPoints));

    // ordering of pyramid and tet points seems to be random
    cellPoints.resize(5);
    cellPoints[0] = 4;
    cellPoints[1] = 5;
    cellPoints[2] = 2;
    cellPoints[3] = 1;
    cellPoints[4] = 12;
    cells.append(cellShape(pyr, cellPoints));

    cellPoints.resize(4);
    cellPoints[0] = 12;
    cellPoints[1] = 4;
    cellPoints[2] = 5;
    cellPoints[3] = 13;
    cells.append(cellShape(tet, cellPoints));

    Info << "Created cells" << endl;

    // create patch definitions from lists of points defining individual faces
    // this is a list of List<List<label> > objects; each lowest-level list
    // contains indices of vertices making up the patch; these are grouped by
    // boundary index, i.e. patch number. The overall list holds all of those in
    // a single container.
    faceListList patchFaces;

    // create a single face
    List<label> patchPoints(4);
    patchPoints[0] = 3; // ordered clock-wise looking outside of fluid domain
    patchPoints[1] = 9;
    patchPoints[2] = 6;
    patchPoints[3] = 0;
    patchFaces.append(List<face>(1, face(patchPoints)));

    // create a list of faces constituting one patch
    List<face> faces(4);
    patchPoints[0] = 14;
    patchPoints[1] = 15;
    patchPoints[2] = 4;
    patchPoints[3] = 13;
    faces[0] = face(patchPoints);
    patchPoints[0] = 15;
    patchPoints[1] = 16;
    patchPoints[2] = 5;
    patchPoints[3] = 4;
    faces[1] = face(patchPoints);
    patchPoints[0] = 16;
    patchPoints[1] = 14;
    patchPoints[2] = 13;
    patchPoints[3] = 5;
    faces[2] = face(patchPoints);
    patchPoints.resize(3);
    patchPoints[0] = 15;
    patchPoints[1] = 14;
    patchPoints[2] = 16;
    faces[3] = face(patchPoints);
    patchFaces.append(faces);

    // add external faces of the tet
    faces.resize(2);
    patchPoints[0] = 4;
    patchPoints[1] = 12;
    patchPoints[2] = 13;
    faces[0] = face(patchPoints);
    patchPoints[0] = 13;
    patchPoints[1] = 12;
    patchPoints[2] = 5;
    faces[1] = face(patchPoints);
    patchFaces.append(faces);

    // add external faces of the pyramid
    faces.resize(3);
    patchPoints[0] = 4;
    patchPoints[1] = 1;
    patchPoints[2] = 12;
    faces[0] = face(patchPoints);
    patchPoints[0] = 12;
    patchPoints[1] = 1;
    patchPoints[2] = 2;
    faces[1] = face(patchPoints);
    patchPoints[0] = 12;
    patchPoints[1] = 2;
    patchPoints[2] = 5;
    faces[2] = face(patchPoints);
    patchFaces.append(faces);

    Info << "Created patches" << endl;

    // patch names
    List<word> boundaryPatchNames;
    boundaryPatchNames.append("hexSide0");
    boundaryPatchNames.append("prismFaces");
    boundaryPatchNames.append("tetFaces");
    boundaryPatchNames.append("pyramidFaces");

    // types and physical types
    wordList boundaryPatchTypes(patchFaces.size());
    boundaryPatchTypes[0] = "symmetryPlane";
    boundaryPatchTypes[1] = "patch";
    boundaryPatchTypes[2] = "wall";
    boundaryPatchTypes[3] = "symmetry";

    wordList boundaryPatchPhysicalTypes(patchFaces.size());
    boundaryPatchPhysicalTypes[0] = "symmetryPlane";
    boundaryPatchPhysicalTypes[1] = "patch";
    boundaryPatchPhysicalTypes[2] = "wall";
    boundaryPatchPhysicalTypes[3] = "symmetry";

    // name of the mesh region - use OF default
    word regionName = polyMesh::defaultRegion;

    // default values for undefined faces
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
        clone(points),
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

    // ---
    Info << nl << "To best visualise the results, load the mesh and extract all patches" << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
