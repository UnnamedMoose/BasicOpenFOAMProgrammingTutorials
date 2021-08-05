/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 AUTHOR,AFFILIATION
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

Application
    OFtutorial16_particleTracking

Description

this application computes the massless particle tracking over an already-solved
flow field. it helps in calculating the residence time i.e. time took by a
particle to leave the domain, and draw path of the particle.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

void saveTrajectoryToVtk(DynamicList<point>& path, fvMesh& mesh)
{
    // Used for writing a vtk file that contains a single polyline defined by
    // a list of points.
    // For VTK file format specification please refer to:
    // https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf

    Info << nl << "Writing particle's path as VTK file";

    // making a separate directory for storing VTK file
    fileName VTK_dir = mesh.time().path()/"VTK";
    mkDir(VTK_dir);

    // Create file pointer.
    autoPtr<OFstream> vtkFilePtr;
    vtkFilePtr.reset(new OFstream(VTK_dir/"particle_path.vtk"));

    // Write vtk header.
    vtkFilePtr() << "# vtk DataFile Version 2.0" << endl << "particle_path" << endl << "ASCII" << endl << "DATASET POLYDATA" << endl;
    vtkFilePtr() << nl;

    // Write points header.
    vtkFilePtr() << "POINTS " << path.size() << " DOUBLE" << endl;

    // Write point coordinates.
    forAll(path, ipt)
    {
      vtkFilePtr() << path[ipt].x() << " " << path[ipt].y() << " " << path[ipt].z() << endl;
    }
    vtkFilePtr() << nl;

    // Write lines header and list.
    vtkFilePtr() << "LINES 1 " << path.size() + 1 << endl; // <no of lines> <size of point indices list>
    vtkFilePtr() << path.size() << endl; // <no of points> <point list .....>
    forAll(path, ipt)
    {
      vtkFilePtr() << ipt << endl;
    }

    Info << tab << " Done writing VTK." << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // Select the latest available time step solution.
    instantList times = runTime.times(); // Get all the available times.
    runTime.setTime(times.last(), 0); // Set the latest/last one as the current time value.
    Info << nl << "Using flow field data from time = " << runTime.timeName() << endl;

    // Read the velocity field from latest solutiom time folder.
    Info << nl << "Reading velocity field" <<  endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    // Define initial position of the particle. Hard coded for simplicity - see
    // tutorial 01 for handling I/O if this were to be extended.
    point particle(0.01, 0.11875, 0);

    // Prepare the lists of particle positions, and a couple of scalar
    // variables to store its total travel time and distance.
    DynamicList<point> particlePositions;
    scalar timeTaken(0);
    scalar distanceTravelled(0);

    // Add the initial position to the particlePositions list.
    particlePositions.append(particle);

    // In order to start, we need to find the cell in which the massless
    // particle resides, so that the flow velocity of that cell can be used
    // to determine its next position.
    label cellId = mesh.findCell(particle);

    Info << nl << "particle initially recides inside cell with ID = "
        << cellId << endl;

    // declaring additional variables
    vector curPos = particle;
    vector newPos(0,0,0);
    label iterCount(1);
    const label maxIters(100);

    // starting particle tracking
    Info << nl << "Starting particle tracking" << endl;
    while (cellId != -1) // cellid will be -1 when particle leaves the domain.
    {
        // Get velocity for the cell in which the particle resides.
        vector velocity = U[cellId];

        // Get the current cell's characteristic length i.e. cubic root of its volume.
        scalar charLen = Foam::cbrt(mesh.V()[cellId]);

        // Calculate the characterstic timestep.
        scalar dt = charLen/mag(velocity);

        // Calculate new position of the particle based on its current position
        // and local flow velocity.
        newPos = curPos + velocity*dt;

        // Compute the distance travelled within this time step.
        scalar dist = mag(newPos - curPos);

        // updating on terminal
        Info << nl << "Lagrangian time step: " << iterCount << nl
            << tab << "current position = " << curPos << nl
            << tab << "new position = " << newPos << nl
            << tab << "local distance travelled = " << dist << nl
            << tab << "local time taken = " << dt << nl
            << tab << "currently reciding cell no. = " << cellId << endl;

        // Append the calculated values to the list and variables.
        distanceTravelled += dist;
        timeTaken += dt;
        particlePositions.append(newPos);
        curPos = newPos;
        iterCount++;

        // Find the new cell into which the particle has moved.
        // NOTE: this is the most expensive part of the code because it triggers
        //  a global mesh loop for each Lagrangian time step. The associated cost
        //  would grow rapidly with adding more particles so a better solution would
        //  be to check for crossing of the current particle path with the faces of
        //  the old owner cell and perform local updates. This is too involved for
        //  the purpose of this short demonstration, however.
        cellId = mesh.findCell(curPos);

        // Safeguard against infinite loops.
        if (iterCount > maxIters)
        {
            FatalErrorInFunction << "Maximum number of Lagrangian time steps exceeded." << abort(FatalError);
        }
    }

    if (cellId == -1)  // i.e. particle left the mesh domain
    {
        Info << nl << nl << "Particle left the domain! " << nl
            << "Total distance travelled = " << distanceTravelled << nl
            << "Total time taken = " << timeTaken <<  endl;
    }

    // Write particle's path as a polyline into a VTK file.
    saveTrajectoryToVtk(particlePositions, mesh);

    Info<< nl;
    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
