/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::fv::customActuationDiskSource

Description
    Actuation disk source created directly from the mesh without using a cell set.
    Can be used to easily adjust the position of the disk without modifying the
    grid and topoSet dictionary.

    Constant values for momentum source for actuation disk
    \f[
        T = 2 \rho A U_{o}^2 a (1-a)
    \f]
    and
    \f[
        U_1 = (1 - a)U_{o}
    \f]

    where:
    \vartable
        A   | disk area
        U_o | upstream velocity resolved onto the disk normal direction
        a   | 1 - Cp/Ct
        U_1 | velocity at the disk
    \endvartable

Usage
    Example usage:
    \verbatim
    upstreamPoint (-1.5 0 0); // point for sampling U_o
    diskDir (-1.0 0.5 0); // disk normal direction (need not be unit vector)
    diskCentre (0 0 0); // origin of the disk
    D 1; // disk diameter
    thickness 0.1; // thickness of the disk; NOTE: must be chosen with care to
        // ensure enough cells are selected, use start-up print-out for guidance
    Cp 0.386; // power and thrust coefficients
    Ct 0.58;
    fields (U); // names of fields to be affected by the momentum source
    \endverbatim

SourceFiles
    customActuationDiskSource.C
    customActuationDiskSourceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef customActuationDiskSource_H
#define customActuationDiskSource_H

#include "fvCFD.H"
#include "fvModel.H"
#include "fvMesh.H"
#include "fvMatrix.H"
#include "geometricOneField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                     Class customActuationDiskSource Declaration
\*---------------------------------------------------------------------------*/

class customActuationDiskSource
:
    public fvModel
{

protected:

    // Protected data

        // ---
        // specified in the dictionary

        //- Name of the velocity field
        word Uname_;

        //- Disk area normal
        vector diskDir_;

        // centre of the disk
        point diskCentre_;

        // disk diameter
        scalar D_;

        // disk thicknesss
        scalar t_;

        // Power coefficient
        scalar Cp_;

        // Thrust coefficient
        scalar Ct_;

        //- Upstream point sample
        point upstreamPoint_;

        // ---
        // calculated

        // Upstream cell ID
        label upstreamCellId_;

        // Set of cells to apply the source to
        labelList cells_;

        // Sum of cell volumes
        scalar V_;

private:

    // Private Member Functions

        // Check validity of the input data from the dictionary
        void checkData() const;

        // Add momentum source to the UEqn
        template<class RhoFieldType>
        void calculateMomentumSource
        (
            vectorField& Usource,
            const labelList& cells,
            const scalarField& V,
            const RhoFieldType& rho,
            const vectorField& U
        ) const;

        //- Disallow default bitwise copy construct
        customActuationDiskSource(const customActuationDiskSource&);

        //- Disallow default bitwise assignment
        void operator=(const customActuationDiskSource&) = delete;

public:

    //- Runtime type information
    TypeName("customActuationDiskSource");

    // Constructors

        //- Construct from components
        customActuationDiskSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

    //- Destructor
    virtual ~customActuationDiskSource()
    {}

    // Member Functions

        // Access

            //- Return the list of fields for which the fvModel adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;

            // Return Cp
            scalar Cp() const
            {
                return Cp_;
            }

            // Return Ct
            scalar Ct() const
            {
                return Ct_;
            }

            // Normal disk direction
            const vector& diskDir() const
            {
                return diskDir_;
            }

            // Disk area
            inline scalar diskArea() const
            {
                return constant::mathematical::pi * pow(D_/2., 2.);
            }

            // Total volume of enclosed cells
            scalar V() const
            {
                return V_;
            }

         // explicit and implicit contributions to the field being modified

            //- Source term to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const word& fieldName
            ) const;

            //- Source term to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const word& fieldName
            ) const;

        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "customActuationDiskSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
