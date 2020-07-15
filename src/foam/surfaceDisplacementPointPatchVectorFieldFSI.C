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

//template<>
//const char*
//NamedEnum<surfaceDisplacementPointPatchVectorFieldFSI::projectMode, 3>::
//names[] =
//{
//    "nearest",
//    "pointNormal",
//    "fixedNormal"
//};

//const NamedEnum<surfaceDisplacementPointPatchVectorFieldFSI::projectMode, 3>
//    surfaceDisplacementPointPatchVectorFieldFSI::projectModeNames_;


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
    const pointField&       points = mesh.points();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    const pointVectorField& incompingPointDisplacement = refCast<const pointVectorField>
    (
	    mesh.objectRegistry::lookupObject<pointVectorField>
	    (
	    "pointDisplacementNew"
	    )
    );

    vectorField patchNewPointDisplacement(*this);

/*
const Time& t = mesh.time();
patchNewPointDisplacement[0][0] = 0.00259899; patchNewPointDisplacement[0][1] = 0; patchNewPointDisplacement[0][2] =  -0.00401698;
patchNewPointDisplacement[1][0] = 0.00789971; patchNewPointDisplacement[1][1] =  0; patchNewPointDisplacement[1][2] =   -0.00816199;
patchNewPointDisplacement[2][0] = 0.00789971; patchNewPointDisplacement[2][1] =  0; patchNewPointDisplacement[2][2] =   -0.00816199;
patchNewPointDisplacement[3][0] = 0.00259899; patchNewPointDisplacement[3][1] =  0; patchNewPointDisplacement[3][2] =   -0.00401698;
patchNewPointDisplacement[4][0] = 0;          patchNewPointDisplacement[4][1] =  0; patchNewPointDisplacement[4][2] =   0;
patchNewPointDisplacement[5][0] = 0;          patchNewPointDisplacement[5][1] =  0; patchNewPointDisplacement[5][2] =   0;
patchNewPointDisplacement[6][0] = 0.00286153; patchNewPointDisplacement[6][1] =  0; patchNewPointDisplacement[6][2] =   0.00336792;
patchNewPointDisplacement[7][0] = 0;          patchNewPointDisplacement[7][1] =  0; patchNewPointDisplacement[7][2] =   0;
patchNewPointDisplacement[8][0] = 0;          patchNewPointDisplacement[8][1] =  0; patchNewPointDisplacement[8][2] =   0;
patchNewPointDisplacement[9][0] = 0.00286153; patchNewPointDisplacement[9][1] =  0; patchNewPointDisplacement[9][2] =   0.00336792;
patchNewPointDisplacement[10][0] = 0.00897453; patchNewPointDisplacement[10][1] =  0; patchNewPointDisplacement[10][2] =   0.00634141;
patchNewPointDisplacement[11][0] = 0.00897453; patchNewPointDisplacement[11][1] =  0; patchNewPointDisplacement[11][2] =   0.00634141;
patchNewPointDisplacement[12][0] = 0.0186637; patchNewPointDisplacement[12][1] =  0; patchNewPointDisplacement[12][2] =   0.00842845;
patchNewPointDisplacement[13][0] = 0.0186637; patchNewPointDisplacement[13][1] =  0; patchNewPointDisplacement[13][2] =   0.00842845;
patchNewPointDisplacement[14][0] = 0.0314006; patchNewPointDisplacement[14][1] =  0; patchNewPointDisplacement[14][2] =   0.00954133;
patchNewPointDisplacement[15][0] = 0.0314006; patchNewPointDisplacement[15][1] =  0; patchNewPointDisplacement[15][2] =   0.00954133;
patchNewPointDisplacement[16][0] = 0.0278702; patchNewPointDisplacement[16][1] =  0; patchNewPointDisplacement[16][2] =   -0.0167539;
patchNewPointDisplacement[17][0] = 0.0418323; patchNewPointDisplacement[17][1] =  0; patchNewPointDisplacement[17][2] =   -0.0214495;
patchNewPointDisplacement[18][0] = 0.0418323; patchNewPointDisplacement[18][1] =  0; patchNewPointDisplacement[18][2] =   -0.0214495;
patchNewPointDisplacement[19][0] = 0.0278702; patchNewPointDisplacement[19][1] =  0; patchNewPointDisplacement[19][2] =   -0.0167539;
patchNewPointDisplacement[20][0] = 0.0164762; patchNewPointDisplacement[20][1] =  0; patchNewPointDisplacement[20][2] =   -0.0123407;
patchNewPointDisplacement[21][0] = 0.0164762; patchNewPointDisplacement[21][1] =  0; patchNewPointDisplacement[21][2] =   -0.0123407;
patchNewPointDisplacement[22][0] = 0.0580029; patchNewPointDisplacement[22][1] =  0; patchNewPointDisplacement[22][2] =   -0.0263898;
patchNewPointDisplacement[23][0] = 0.0580029; patchNewPointDisplacement[23][1] =  0; patchNewPointDisplacement[23][2] =   -0.0263898;
patchNewPointDisplacement[24][0] = 0.0759499; patchNewPointDisplacement[24][1] =  0; patchNewPointDisplacement[24][2] =   -0.0314988;
patchNewPointDisplacement[25][0] = 0.0759499; patchNewPointDisplacement[25][1] =  0; patchNewPointDisplacement[25][2] =   -0.0314988;
patchNewPointDisplacement[26][0] = 0.0952359; patchNewPointDisplacement[26][1] =  0; patchNewPointDisplacement[26][2] =   -0.0367276;
patchNewPointDisplacement[27][0] = 0.0952359; patchNewPointDisplacement[27][1] =  0; patchNewPointDisplacement[27][2] =   -0.0367276;
patchNewPointDisplacement[28][0] = 0.0467775; patchNewPointDisplacement[28][1] =  0; patchNewPointDisplacement[28][2] =   0.00956404;
patchNewPointDisplacement[29][0] = 0.0467775; patchNewPointDisplacement[29][1] =  0; patchNewPointDisplacement[29][2] =   0.00956404;
patchNewPointDisplacement[30][0] = 0.0642932; patchNewPointDisplacement[30][1] =  0; patchNewPointDisplacement[30][2] =   0.00846738;
patchNewPointDisplacement[31][0] = 0.0642932; patchNewPointDisplacement[31][1] =  0; patchNewPointDisplacement[31][2] =   0.00846738;
patchNewPointDisplacement[32][0] = 0.0834096; patchNewPointDisplacement[32][1] =  0; patchNewPointDisplacement[32][2] =   0.00633716;
patchNewPointDisplacement[33][0] = 0.0834096; patchNewPointDisplacement[33][1] =  0; patchNewPointDisplacement[33][2] =   0.00633716;
patchNewPointDisplacement[34][0] = 0.103658; patchNewPointDisplacement[34][1] =  0; patchNewPointDisplacement[34][2] =   0.00336065;
patchNewPointDisplacement[35][0] = 0.103658; patchNewPointDisplacement[35][1] =  0; patchNewPointDisplacement[35][2] =   0.00336065;
patchNewPointDisplacement[36][0] = 0.124728; patchNewPointDisplacement[36][1] =  0; patchNewPointDisplacement[36][2] =   -0.000252435;
patchNewPointDisplacement[37][0] = 0.124728; patchNewPointDisplacement[37][1] =  0; patchNewPointDisplacement[37][2] =   -0.000252435;
patchNewPointDisplacement[38][0] = 0.18063; patchNewPointDisplacement[38][1] =  0; patchNewPointDisplacement[38][2] =   -0.0592771;
patchNewPointDisplacement[39][0] = 0.203513; patchNewPointDisplacement[39][1] =  0; patchNewPointDisplacement[39][2] =   -0.0654252;
patchNewPointDisplacement[40][0] = 0.203513; patchNewPointDisplacement[40][1] =  0; patchNewPointDisplacement[40][2] =   -0.0654252;
patchNewPointDisplacement[41][0] = 0.18063; patchNewPointDisplacement[41][1] =  0; patchNewPointDisplacement[41][2] =   -0.0592771;
patchNewPointDisplacement[42][0] = 0.115515; patchNewPointDisplacement[42][1] =  0; patchNewPointDisplacement[42][2] =   -0.0420828;
patchNewPointDisplacement[43][0] = 0.115515; patchNewPointDisplacement[43][1] =  0; patchNewPointDisplacement[43][2] =   -0.0420828;
patchNewPointDisplacement[44][0] = 0.136574; patchNewPointDisplacement[44][1] =  0; patchNewPointDisplacement[44][2] =   -0.0476056;
patchNewPointDisplacement[45][0] = 0.136574; patchNewPointDisplacement[45][1] =  0; patchNewPointDisplacement[45][2] =   -0.0476056;
patchNewPointDisplacement[46][0] = 0.1583;   patchNewPointDisplacement[46][1] =  0; patchNewPointDisplacement[46][2] =   -0.053333;
patchNewPointDisplacement[47][0] = 0.1583;   patchNewPointDisplacement[47][1] =  0; patchNewPointDisplacement[47][2] =   -0.053333;
patchNewPointDisplacement[48][0] = 0.226891; patchNewPointDisplacement[48][1] =  0; patchNewPointDisplacement[48][2] =   -0.071745;
patchNewPointDisplacement[49][0] = 0.226891; patchNewPointDisplacement[49][1] =  0; patchNewPointDisplacement[49][2] =   -0.071745;
patchNewPointDisplacement[50][0] = 0.250692; patchNewPointDisplacement[50][1] =  0; patchNewPointDisplacement[50][2] =   -0.0781904;
patchNewPointDisplacement[51][0] = 0.250692; patchNewPointDisplacement[51][1] =  0; patchNewPointDisplacement[51][2] =   -0.0781904;
patchNewPointDisplacement[52][0] = 0.274827; patchNewPointDisplacement[52][1] =  0; patchNewPointDisplacement[52][2] =   -0.0847105;
patchNewPointDisplacement[53][0] = 0.274827; patchNewPointDisplacement[53][1] =  0; patchNewPointDisplacement[53][2] =   -0.0847105;
patchNewPointDisplacement[54][0] = 0.299201; patchNewPointDisplacement[54][1] =  0; patchNewPointDisplacement[54][2] =   -0.0912577;
patchNewPointDisplacement[55][0] = 0.299201; patchNewPointDisplacement[55][1] =  0; patchNewPointDisplacement[55][2] =   -0.0912577;
patchNewPointDisplacement[56][0] = 0.323721; patchNewPointDisplacement[56][1] =  0; patchNewPointDisplacement[56][2] =   -0.0977962;
patchNewPointDisplacement[57][0] = 0.323721; patchNewPointDisplacement[57][1] =  0; patchNewPointDisplacement[57][2] =   -0.0977962;
patchNewPointDisplacement[58][0] = 0.146466; patchNewPointDisplacement[58][1] =  0; patchNewPointDisplacement[58][2] =   -0.00434417;
patchNewPointDisplacement[59][0] = 0.146466; patchNewPointDisplacement[59][1] =  0; patchNewPointDisplacement[59][2] =   -0.00434417;
patchNewPointDisplacement[60][0] = 0.168804; patchNewPointDisplacement[60][1] =  0; patchNewPointDisplacement[60][2] =   -0.00882705;
patchNewPointDisplacement[61][0] = 0.168804; patchNewPointDisplacement[61][1] =  0; patchNewPointDisplacement[61][2] =   -0.00882705;
patchNewPointDisplacement[62][0] = 0.348307; patchNewPointDisplacement[62][1] =  0; patchNewPointDisplacement[62][2] =   -0.104308;
patchNewPointDisplacement[63][0] = 0.348307; patchNewPointDisplacement[63][1] =  0; patchNewPointDisplacement[63][2] =   -0.104308;
patchNewPointDisplacement[64][0] = 0.372911; patchNewPointDisplacement[64][1] =  0; patchNewPointDisplacement[64][2] =   -0.110795;
patchNewPointDisplacement[65][0] = 0.372911; patchNewPointDisplacement[65][1] =  0; patchNewPointDisplacement[65][2] =   -0.110795;
patchNewPointDisplacement[66][0] = 0.379386; patchNewPointDisplacement[66][1] =  0; patchNewPointDisplacement[66][2] =   -0.0862009;
patchNewPointDisplacement[67][0] = 0.379386; patchNewPointDisplacement[67][1] =  0; patchNewPointDisplacement[67][2] =   -0.0862009;
patchNewPointDisplacement[68][0] = 0.385869; patchNewPointDisplacement[68][1] =  0; patchNewPointDisplacement[68][2] =   -0.0616025;
patchNewPointDisplacement[69][0] = 0.385869; patchNewPointDisplacement[69][1] =  0; patchNewPointDisplacement[69][2] =   -0.0616025;
patchNewPointDisplacement[70][0] = 0.191695; patchNewPointDisplacement[70][1] =  0; patchNewPointDisplacement[70][2] =   -0.0136619;
patchNewPointDisplacement[71][0] = 0.191695; patchNewPointDisplacement[71][1] =  0; patchNewPointDisplacement[71][2] =   -0.0136619;
patchNewPointDisplacement[72][0] = 0.215086; patchNewPointDisplacement[72][1] =  0; patchNewPointDisplacement[72][2] =   -0.0188303;
patchNewPointDisplacement[73][0] = 0.215086; patchNewPointDisplacement[73][1] =  0; patchNewPointDisplacement[73][2] =   -0.0188303;
patchNewPointDisplacement[74][0] = 0.238904; patchNewPointDisplacement[74][1] =  0; patchNewPointDisplacement[74][2] =   -0.0243164;
patchNewPointDisplacement[75][0] = 0.238904; patchNewPointDisplacement[75][1] =  0; patchNewPointDisplacement[75][2] =   -0.0243164;
patchNewPointDisplacement[76][0] = 0.263063; patchNewPointDisplacement[76][1] =  0; patchNewPointDisplacement[76][2] =   -0.0300938;
patchNewPointDisplacement[77][0] = 0.263063; patchNewPointDisplacement[77][1] =  0; patchNewPointDisplacement[77][2] =   -0.0300938;
patchNewPointDisplacement[78][0] = 0.287465; patchNewPointDisplacement[78][1] =  0; patchNewPointDisplacement[78][2] =   -0.0361215;
patchNewPointDisplacement[79][0] = 0.287465; patchNewPointDisplacement[79][1] =  0; patchNewPointDisplacement[79][2] =   -0.0361215;
patchNewPointDisplacement[80][0] = 0.312015; patchNewPointDisplacement[80][1] =  0; patchNewPointDisplacement[80][2] =   -0.0423461;
patchNewPointDisplacement[81][0] = 0.312015; patchNewPointDisplacement[81][1] =  0; patchNewPointDisplacement[81][2] =   -0.0423461;
patchNewPointDisplacement[82][0] = 0.336631; patchNewPointDisplacement[82][1] =  0; patchNewPointDisplacement[82][2] =   -0.0487061;
patchNewPointDisplacement[83][0] = 0.336631; patchNewPointDisplacement[83][1] =  0; patchNewPointDisplacement[83][2] =   -0.0487061;
patchNewPointDisplacement[84][0] = 0.361258; patchNewPointDisplacement[84][1] =  0; patchNewPointDisplacement[84][2] =   -0.0551419;
patchNewPointDisplacement[85][0] = 0.361258; patchNewPointDisplacement[85][1] =  0; patchNewPointDisplacement[85][2] =   -0.0551419;
*/
    int patchID{-1};
    forAll(patches, ipatch)
    {
        const polyPatch& patch = patches[ipatch];
        const labelList& patchPoints = patch.meshPoints();
        if (patchPoints.size() == this->size())
        {
            patchID = ipatch;
            

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

//Info << patchNewPointDisplacement << endl;
//Info << "In surfaceDisplacementPointPatchVectorFieldFSI" << endl;
//std::cin.get();

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
