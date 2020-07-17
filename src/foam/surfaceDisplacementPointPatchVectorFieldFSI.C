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

\*---------------------------------------------------------------------------*/

#include "surfaceDisplacementPointPatchVectorFieldFSI.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "transformField.H"
#include "dynamicMotionSolverFvMesh.H"
#include "displacementMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

surfaceDisplacementPointPatchVectorFieldFSI::
surfaceDisplacementPointPatchVectorFieldFSI
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchVectorField(p, iF)
{}


surfaceDisplacementPointPatchVectorFieldFSI::
surfaceDisplacementPointPatchVectorFieldFSI
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchVectorField(p, iF, dict)
{}


surfaceDisplacementPointPatchVectorFieldFSI::
surfaceDisplacementPointPatchVectorFieldFSI
(
    const surfaceDisplacementPointPatchVectorFieldFSI& ppf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchVectorField(ppf, p, iF, mapper)
{}


surfaceDisplacementPointPatchVectorFieldFSI::
surfaceDisplacementPointPatchVectorFieldFSI
(
    const surfaceDisplacementPointPatchVectorFieldFSI& ppf
)
:
    fixedValuePointPatchVectorField(ppf)
{}


surfaceDisplacementPointPatchVectorFieldFSI::
surfaceDisplacementPointPatchVectorFieldFSI
(
    const surfaceDisplacementPointPatchVectorFieldFSI& ppf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchVectorField(ppf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void surfaceDisplacementPointPatchVectorFieldFSI::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const polyMesh& mesh = patch().boundaryMesh().mesh()();
    //const pointField&       points = mesh.points();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    const pointVectorField& incompingPointDisplacement = refCast<const pointVectorField>
    (
	    mesh.objectRegistry::lookupObject<pointVectorField>
	    (
	    "pointDisplacementNew"
	    )
    );

    vectorField patchNewPointDisplacement(*this);

    forAll(patches, ipatch)
    {
        const polyPatch& patch = patches[ipatch];
        const labelList& patchPoints = patch.meshPoints();
        if (patchPoints.size() == this->size())
        {
            // Loop over all nodes of boundary patch
            forAll(patchPoints, pointi)
            {
                const label& pointID = patch.meshPoints()[pointi];  // Node index
                // Do your calculations
                patchNewPointDisplacement[pointi] =
                    incompingPointDisplacement[pointID];
            }
            break;
        }
    }

    this->operator==(patchNewPointDisplacement);
    fixedValuePointPatchVectorField::updateCoeffs();

}


void surfaceDisplacementPointPatchVectorFieldFSI::write(Ostream& os) const
{
    fixedValuePointPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    fixedValuePointPatchVectorField,
    surfaceDisplacementPointPatchVectorFieldFSI
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
