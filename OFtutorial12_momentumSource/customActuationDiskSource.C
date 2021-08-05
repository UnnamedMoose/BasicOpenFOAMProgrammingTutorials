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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(customActuationDiskSource, 0);
    addToRunTimeSelectionTable
    (
        fvModel,
        customActuationDiskSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::customActuationDiskSource::checkData() const
{
    if (magSqr(this->diskArea()) <= VSMALL)
    {
        FatalErrorInFunction
           << "diskArea is approximately zero"
           << exit(FatalIOError);
    }
    if ((Cp_ <= VSMALL) || (Ct_ <= VSMALL))
    {
        FatalErrorInFunction
           << "Cp and Ct must be greater than zero"
           << exit(FatalIOError);
    }
    if (mag(diskDir_) < VSMALL)
    {
        FatalErrorInFunction
           << "disk direction vector is approximately zero"
           << exit(FatalIOError);
    }
    if (returnReduce(upstreamCellId_, maxOp<label>()) == -1)
    {
        FatalErrorInFunction
           << "upstream location " << upstreamPoint_  << " not found in mesh"
           << exit(FatalIOError);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::customActuationDiskSource::customActuationDiskSource
(
    const word& name, // name as appears in fvOptions dictionary
    const word& modelType, // type of this instance
    const dictionary& dict,
    const fvMesh& mesh
)
:
    // call the base class constructor
    fvModel(name, modelType, dict, mesh),
    // read control parameters from dictionary
    Uname_(coeffs().lookup("Uname")),
    diskDir_(coeffs().lookup("diskDir")),
    diskCentre_(coeffs().lookup("diskCentre")),
    D_(readScalar(coeffs().lookup("D"))),
    t_(readScalar(coeffs().lookup("thickness"))),
    Cp_(readScalar(coeffs().lookup("Cp"))),
    Ct_(readScalar(coeffs().lookup("Ct"))),
    upstreamPoint_(coeffs().lookup("upstreamPoint")),
    // initialise other member fields
    upstreamCellId_(-1),
    V_(0.0)
{
    Info << tab << "- creating actuation disk zone: " << this->name() << endl;

    // Locate the upstream cell used to apply velocity magnitude augment
    upstreamCellId_ = mesh.findCell(upstreamPoint_);

    // ===
    // select cells forming part of the actuator disk

    // make sure direction vector is unit length
    diskDir_ /= mag(diskDir_);

    // Create a dimensioned origin position to use with default OF mesh types
    dimensionedVector x0 ("x0", dimLength, diskCentre_);

    // compute distance and normal vectors from origin of the disk to use for selection
    scalarField R (mag(mesh.C() - x0));
    vectorField rHat ((mesh.C() - x0) / mag(mesh.C() - x0));

    // go over each cell in the grid and comapre it against selection criteria
    for (label cellI = 0; cellI < mesh.C().size(); cellI++)
    {
        // determine distance from cell centre to disk axis and along the normal direction
        scalar dNormal = R[cellI] * (rHat[cellI] & diskDir_);
        scalar dRadial = sqrt(pow(R[cellI], 2.) - pow(dNormal, 2.));

        // see if the cell is within tolerance from the centre of the disk
        if ((mag(dNormal) < t_/2.) && (dRadial < D_/2.))
            cells_.append(cellI);
    }

    // ===
    // Set volume information
    V_ = 0.0;
    forAll(cells_, i)
        V_ += mesh.V()[cells_[i]];
    reduce(V_, sumOp<scalar>());

    // ===
    // validate inputs
    checkData();

    Info << tab
         << "- selected " << returnReduce(cells_.size(), sumOp<label>())
         << " cell(s) with volume " << V_ << endl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::customActuationDiskSource::addSupFields() const
{
    return wordList(1, Uname_);
}

void Foam::fv::customActuationDiskSource::addSup
(
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    const scalarField& cellsV = mesh().V();
    vectorField& Usource = eqn.source();
    const vectorField& U = eqn.psi();

    if (V() > VSMALL)
    {
        calculateMomentumSource
        (
            Usource,
            cells_,
            cellsV,
            geometricOneField(),
            U
        );
    }
}

void Foam::fv::customActuationDiskSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    const scalarField& cellsV = mesh().V();
    vectorField& Usource = eqn.source();
    const vectorField& U = eqn.psi();

    if (V() > VSMALL)
    {
        calculateMomentumSource
        (
            Usource,
            cells_,
            cellsV,
            rho,
            U
        );
    }
}

bool Foam::fv::customActuationDiskSource::read(const dictionary& dict)
{
    if (fvModel::read(dict))
    {
        coeffs().readIfPresent("Uname", Uname_);
        coeffs().readIfPresent("diskDir", diskDir_);
        coeffs().readIfPresent("diskCentre", diskCentre_);
        coeffs().readIfPresent("D", D_);
        coeffs().readIfPresent("thickness", t_);
        coeffs().readIfPresent("Cp", Cp_);
        coeffs().readIfPresent("Ct", Ct_);

        checkData();

        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
