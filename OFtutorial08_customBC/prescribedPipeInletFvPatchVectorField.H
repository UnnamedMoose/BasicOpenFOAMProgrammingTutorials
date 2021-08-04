/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::prescribedPipeInletFvPatchVectorField

Group
    grpInletBoundaryConditions grpWallBoundaryConditions

Description
    Prescribes an inlet profile for a circular pipe inlet patch. Flow direction is determined
    automatically from face normals, only velocity magnitude needs to be provided.

    \heading Patch usage

    \table
        Property    		| Description             								| Required  | Default value
        R           		| pipe radius 											| yes     	|
        flowSpeed   		| magnitude of flow velocity      						| yes       |
        deltaByR			| wall function thickness as a fraction of pipe radius  | yes 		|
        approximationType	| name of approximation profile							| no		| exponential
        centrePoint			| point on the centreline of the pipe					| no 		| (0 0 0)
        lambda				| pressure gradient coefficient 						| no 		| 0
    \endtable

    Example of the boundary condition specification:
    \verbatim
    inlet
	{
		// prescribes a parabolic inlet profile
	    type        prescribedPipeInlet;
	    // (optional, default exponential) type of profile used for approximation
	    approximationType exponential;
	    // pipe radius
	    R	$Rpipe;
	    // magnitude of flow velocity into the pipe
	    flowSpeed	$flowVelocity;
	    // (optional, defaults at origin) centrepoint of the pipe used to determine the radius of each face
	    centrepoint	(0 0 0);
	    // inlet BL thickness as a function of radius
	    deltaByR $deltaByR;
	    // (optional, default 0) pressure gradient coefficient (see Polhausen profile literature for explanation)
	    lambda 0.;
	    // dummy value overwritten by the BC, used for initialisation only
	    value		$internalField;
	}
	\endverbatim

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    prescribedPipeInletFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef prescribedPipeInletFvPatchVectorField_H
#define prescribedPipeInletFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "Switch.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class prescribedPipeInletFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

// NOTE: the patch is derived from a fixedValue BC since it does prescribe
// a fixed value, although non-uniform across all of the faces, unlike the base
// class does.
class prescribedPipeInletFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
private:
    // Private data
        // NOTE: these class fields are used to determine the inlet profile

    	// denotes the type of the profile used to perform the approximation
        // (optional, defaults to exponential)
    	word approximationType_;

    	// flow velocity magnitude
    	scalar flowSpeed_;

    	// bl thickness
    	scalar deltaByR_;

    	// centre of the pipe (optional, defaults to the origin)
    	vector centrepoint_;

    	// radius of the pipe
    	scalar R_;

    	// pressure gradient coefficient (optional, default to 0, not used by all approximations)
    	scalar lambda_;

public:

    //- Runtime type information
    // NOTE: this gets used by the runtimeSelector. In other words, this is the
    // name under which OpenFOAM knows this BC.
    TypeName("prescribedPipeInlet");

    // Constructors

        //- Construct from patch and internal field
        prescribedPipeInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct as copy setting internal field reference
        prescribedPipeInletFvPatchVectorField
        (
            const prescribedPipeInletFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new prescribedPipeInletFvPatchVectorField(*this, iF)
            );
        }

        //- Construct from patch, internal field and dictionary
        prescribedPipeInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&,
            const bool valueRequired=true
        );

        //- Construct by mapping given prescribedPipeInletFvPatchVectorField
        //  onto a new patch
        prescribedPipeInletFvPatchVectorField
        (
            const prescribedPipeInletFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&,
            const bool mappingRequired=true
        );

        //- Disallow copy without setting internal field reference
        prescribedPipeInletFvPatchVectorField(const prescribedPipeInletFvPatchVectorField&) = delete;

    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
