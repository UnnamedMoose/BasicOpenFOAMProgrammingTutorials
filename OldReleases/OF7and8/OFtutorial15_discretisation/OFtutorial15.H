/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::myCustomScheme

Description
    myCustomScheme differencing scheme class.

SourceFiles
    myCustomScheme.C

\*---------------------------------------------------------------------------*/

#ifndef myCustomScheme_H
#define myCustomScheme_H

#include "surfaceInterpolationScheme.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class myCustomScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class myCustomScheme
:
    public surfaceInterpolationScheme<Type>
{
    // Private Data

        const surfaceScalarField& faceFlux_;


public:

    //- Runtime type information
    TypeName("myCustomScheme");


    // Constructors

        //- Construct from faceFlux
        myCustomScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            faceFlux_(faceFlux)
        {}

        //- Construct from Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        myCustomScheme
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            faceFlux_
            (
                mesh.lookupObject<surfaceScalarField>
                (
                    word(is)
                )
            )
        {}

        //- Construct from faceFlux and Istream
        myCustomScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream&
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            faceFlux_(faceFlux)
        {}


    // Member Functions

        //- Return the interpolation weighting factors.
        // NOTE: key function which defines how much individual neighbour cells
        //  should contribute to the value at this face. It is inherited from the
        //  surfaceInterpolationScheme base class (in surfaceInterpolationScheme.H).
        //  As the prototype is defined as virtual (...) = 0, each derived class
        //  must implement this on its own due to there being no baseline functionality
        //  that might be inherited.
        // In the current example, the scheme is used to treat the div term of the
        //  scalar transport equation defined in scalarTransportFoam.C:
        //        fvScalarMatrix TEqn
        //        (
        //            fvm::ddt(T)
        //          + fvm::div(phi, T)      <=== THIS ONE
        //          - fvm::laplacian(DT, T)
        //         ==
        //            fvOptions(T)
        //        );
        // The computed weights are ultimately used in surfaceInterpolationScheme.C
        //  inside methods like 'interpolate' which use owner (P) and neighbout (N)
        //  values in order to determine face value:
        //      sfi[fi] = lambda[fi]*vfi[P[fi]] + y[fi]*vfi[N[fi]];
        virtual tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const
        {
            // Create the weights. Note the use of a temporary field to match function
            // return type and overall convention. Here we call straightforward linear
            // interpolation.
            tmp<surfaceScalarField> weights_linear = this->mesh().surfaceInterpolation::weights();
            
            // Also create upwind weights by picking values of cells which lie
            // in the positive flux direction 
            tmp<surfaceScalarField> weights_upwind = pos0(this->faceFlux_);
            
            // Store for later (become deallocated as soon as they are used due to how
            // the tmp structure works).
            surfaceScalarField wl = weights_linear();
            surfaceScalarField wu = weights_upwind();
            
            // Add the weights up in order to make a hybrid scheme.
            const scalar fBlend = 0.2;
            tmp<surfaceScalarField> retWeights = weights_linear*(1.0 - fBlend) + weights_upwind*fBlend;
            
            // For the purpose of the tutorial, instantine an actual scalar field.
            // This may be used just as any other one in OpenFOAM.
            surfaceScalarField w = retWeights();
            scalarField phi = this->faceFlux_();

            // - Print a detaield summary for one face for illustration purposes. -
            // Side note: The grid used in the test case has 20 cells so 19 internal
            //  faces (you may convince yourself of it by running checkMesh).
            label iFace = 10;
            
            // Get references to the important parts og grid information.
            const labelUList& owner = this->mesh().owner();
            const labelUList& neighbour = this->mesh().neighbour();
            const vectorField& Cf = this->mesh().faceCentres();
            const vectorField& C = this->mesh().cellCentres();
            const vectorField& Sf = this->mesh().faceAreas();
            
            Info << "Face " << iFace << " x=" << Cf[iFace] << " phi=" << phi[iFace]
                 << " (inferred velocity phi/A=" << (phi[iFace]/mag(Sf[iFace])) << ")"
                 << nl
                 << "x_owner=" << C[owner[iFace]] << " x_neighbour=" << C[neighbour[iFace]]
                 << nl
                 << "w_upwind=" << wu[iFace] << " w_linear=" << wl[iFace] << " w_total=" << w[iFace]
                 << nl << endl;
            
            // Recall that in OpenFOAM, positive flux is defined as flux from owner
            // to neighbour. In this case, velocity is along the x-axis and the
            // owner cell is also at a smaller x-coordinate (upstream of the face).
            // Hence weight of 1 assigns the owner cell value to the face, consistently
            // with the definition of an upwind scheme.
            
            // Because the grid is uniform, linear interpolation simply returns a weight of
            // 0.5, meaning half of each cell value is used.
            // It is worth looking at surfaceInterpolation.C::makeWeights() in order
            // to see in more detail how this is done, but, in short, relative
            // distance of owner and neighbour cell centre from the face centre
            // is used in order to determine the linear weight.

            // Return the originally computed weights.
            return retWeights;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const myCustomScheme&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

