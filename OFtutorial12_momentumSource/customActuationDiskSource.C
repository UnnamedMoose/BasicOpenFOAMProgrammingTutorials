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
        option,
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
    option(name, modelType, dict, mesh),
    // read control parameters from dictionary
    diskDir_(coeffs_.lookup("diskDir")),
    diskCentre_(coeffs_.lookup("diskCentre")),
    D_(readScalar(coeffs_.lookup("D"))),
    t_(readScalar(coeffs_.lookup("thickness"))),
    Cp_(readScalar(coeffs_.lookup("Cp"))),
    Ct_(readScalar(coeffs_.lookup("Ct"))),
    upstreamPoint_(coeffs_.lookup("upstreamPoint")),
    // initialise other member fields
    upstreamCellId_(-1),
    V_(0.0)
{
    Info << tab << "- creating actuation disk zone: " << this->name() << endl;

    // recover names of fields which are to be fed the sources into
    coeffs_.lookup("fields") >> fieldNames_;
    // prepare a list of boolean values for each field - used by the parent class
    // to determine whether the sources have been applied yet or not.
    applied_.setSize(fieldNames_.size(), false);

    // Locate the upstream cell used to apply velocity magnitude augment
    upstreamCellId_ = mesh.findCell(upstreamPoint_);

    // ===
    // select cells forming part of the actuator disk

    // make sure direction vector is unit length
    diskDir_ /= mag(diskDir_);

    // Create a dimensioned origin position to use with default OF mesh types
    dimensionedVector x0 ("x0", dimLength, diskCentre_);

    // compute distance and normal vectors from origin of the disk to use for selection
    scalarField R (mag(mesh_.C() - x0));
    vectorField rHat ((mesh_.C() - x0) / mag(mesh_.C() - x0));

    // go over each cell in the grid and comapre it against selection criteria
    for (label cellI = 0; cellI < mesh_.C().size(); cellI++)
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
        V_ += mesh_.V()[cells_[i]];
    reduce(V_, sumOp<scalar>());

    // ===
    // validate inputs
    checkData();

    Info << tab
         << "- selected " << returnReduce(cells_.size(), sumOp<label>())
         << " cell(s) with volume " << V_ << endl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::customActuationDiskSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    const scalarField& cellsV = mesh_.V();
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
    const label fieldi
)
{
    const scalarField& cellsV = mesh_.V();
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
    if (option::read(dict))
    {
        coeffs_.readIfPresent("diskDir", diskDir_);
        coeffs_.readIfPresent("diskCentre", diskCentre_);
        coeffs_.readIfPresent("D", D_);
        coeffs_.readIfPresent("thickness", t_);
        coeffs_.readIfPresent("Cp", Cp_);
        coeffs_.readIfPresent("Ct", Ct_);

        checkData();

        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
