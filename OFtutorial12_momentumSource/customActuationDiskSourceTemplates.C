/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "customActuationDiskSource.H"

// * * * * * * * * * * * * * * *  Member Functions * * * * * * * * * * * * * //

template<class RhoFieldType>
void Foam::fv::customActuationDiskSource::calculateMomentumSource
(
    vectorField& Usource,
    const labelList& cells,
    const scalarField& Vcells,
    const RhoFieldType& rho,
    const vectorField& U
) const
{
    // induced velocity coefficient
    scalar a = 1.0 - Cp_/Ct_;

    // find minimum velocity and density upstream of the actuator disk
    vector upU = vector(VGREAT, VGREAT, VGREAT);
    scalar upRho = VGREAT;
    if (upstreamCellId_ != -1)
    {
        upU =  U[upstreamCellId_];
        upRho = rho[upstreamCellId_];
    }
    reduce(upU, minOp<vector>());
    reduce(upRho, minOp<scalar>());

    // resolve the upstream velocity to get component normal to the disk
    upU *= (upU/mag(upU)) & (-1.*diskDir_);

    // calculate thrust
    // Units: [(kg m-3) * m2 * (m s-1) = (kg s-1)]
    scalar T = 2.0*upRho*diskArea()*mag(upU)*a*(1 - a);

    // distribute the total force around the disk proportionally to the volume of each cell
    forAll(cells, i)
    {
        // Units: [m3 * m-3 * (kg s-1) * (-) * (m s-1) = (kg m s-2) = (N)]
        Usource[cells[i]] += ((Vcells[cells[i]]/V()) * T * diskDir_ * mag(upU));
    }

    // show a quick summary
    Info << this->name() << ": added thrust = " << T << nl << endl;

    // ===
    // output a summary field highlighting active cells and showing thrust distribution
    if (mesh_.time().outputTime())
    {
        volVectorField momentumSourceField
        (
            IOobject
            (
                "momentumSourceField",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            vector::zero
        );
        forAll(cells_, i)
            momentumSourceField[cells[i]] = Usource[cells[i]];
        momentumSourceField.write();
    }
}


// ************************************************************************* //
