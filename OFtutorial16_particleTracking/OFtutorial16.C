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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // selecting latest available time step solution
    instantList times = runTime.times(); // getting all available times
    runTime.setTime(times.last(),0); // setting the latest/last one as current
    Info << nl << "setting the current simlation time = "
        << runTime.timeName() << endl;

    // reading velocity field from latest solutiom time folder
    Info << nl << "reading velocity field.. " <<  endl;
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

    // defining initial position of particle
    // in this tutorial only single particle is considered, and can later be
    // extended to number of particles easily by the learners.
    point particle(0.01,0.11875,0);

    // preparing the lists of particle's positions, and a couple of scalar
    // variables to store its total travel time and distance
    pointList particlePositions;
    scalar timeTaken(0);
    scalar distanceTravelled(0);

    // adding initial position to the particlePositions list
    particlePositions.append(particle);

    // now, in order to start, we need to find the cell in which our mass-less
    // particle resides, so that the flow velocity of that cell can be used
    // to determine its next position.
    // so finding id of a cell which encloses the particle's coordinates
    label cellId = mesh.findCell(particle);

    Info << nl << "particle initially recides inside cell with ID = "
        << cellId << endl;

    // declaring additional variables
    vector curPos = particle;
    vector newPos(0,0,0);
    label iterCount(1);

    // starting particle tracking
    Info << nl << "starting particle track" << endl;
    while(cellId != -1) // cellid will be -1 when particle leaves the mesh
    {
        // getting current reciding cell's velocity
        vector velocity = U[cellId];

        // getting current reciding cell's characteristic length i.e. cuberoot
        // of its volume
        scalar charLen = Foam::cbrt(mesh.V()[cellId]);

        // calculating characterstic timestep from above variables
        scalar dt = charLen/mag(velocity);

        // calculating the new position of particle based on its current
        // position and flow velocity
        newPos = curPos + velocity*dt;

        // calculating local distance travelled
        scalar dist = mag(newPos - curPos);

        // updating on terminal
        Info << nl << "iteration : " << iterCount << nl
            << tab << "current position = " << curPos << nl
            << tab << "new position = " << newPos << nl
            << tab << "local distance travelled = " << dist << nl
            << tab << "local time taken = " << dt << nl
            << tab << "currently reciding cell's id = " << cellId << endl;

        // appending the calculated values to the list and variables
        distanceTravelled += dist;
        timeTaken += dt;
        particlePositions.append(newPos);
        curPos = newPos;
        iterCount++;

        // finding the new cell's id, where the particle moved into
        cellId = mesh.findCell(curPos);
    }

    if (cellId == -1)  // i.e. particle left the mesh domain
    {
        Info << nl << nl << "Particle left the mesh domain! " << nl << endl;
        Info << nl << "Total distance travelled = " << distanceTravelled << nl
             << "Total time taken = " << timeTaken <<  endl;
    }

    // writing the particle's path as a poly line under VTK file
    #include "writeVTK.H"

    Info<< nl;
    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}
/*-----------------------------------------------------------------------------

  P.S:
  this code is made for tutorial purpose only as this method of tracking a
  massless particle can be expensive especially when the mesh size is quite
  high. the main CPU intensive task will be given by the following line

        cellId = mesh.findCell(curPos);

  as for each time, the function will check all the cells present in the
  domain. Hence this code can be further improved by the learners and make
  it optimized for larger domains.

  one suggestion is to reduce the search range by looking at the nearby cells
  to the current cell where the particle recides.

-----------------------------------------------------------------------------*/


// ************************************************************************* //
