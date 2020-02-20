#include "rocRhoCentral.H"

//^^^ DEFINITION OF CONSTRUCTORS ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
rhoCentral::rhoCentral()
    : posPtr(NULL),
      negPtr(NULL),
      amaxSfPtr(NULL),
      pThermoPtr(NULL),
      fluxScheme(""),
      inviscid(false)
{
    solverType = const_cast<char *>("rocRhoCentral");
};

rhoCentral::rhoCentral(int argc, char *argv[])
    : posPtr(NULL),
      negPtr(NULL),
      amaxSfPtr(NULL),
      pThermoPtr(NULL),
      fluxScheme(""),
      inviscid(false)
{
    solverType = const_cast<char *>("rocRhoCentral");
    initialize(argc, argv);
}
//===================================================================


//^^^ DEFINITION OF COM-RELATED MTHODS ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//^^^^^ LOAD MODULES ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
void rhoCentral::load(const char *name)
{
    //  Anouncing default communicator  ^^^^^^^^^^^^^^^^^^^
    MPI_Comm tmpComm;
    tmpComm = COM_get_default_communicator();  

    int tmpRank, tmpNProc;
    MPI_Comm_rank(tmpComm, &tmpRank);
    MPI_Comm_size(tmpComm, &tmpNProc);
    
    if (tmpRank == 0)
    {
        std::cout << "rocFoam.load: Loading rhoCentral with name "
                   << name << "." << std::endl;

        std::cout << "rocFoam.load: Rank = " << tmpRank
                  << ", NProc = " << tmpNProc
                  << ", COMM = " << tmpComm << std::endl;

        std::cout << std::endl;
    }


    //  Register module with COM ^^^^^^^^^^^^^^^^^^^^^^^^^^
    rhoCentral *comFoamPtr = new rhoCentral();

    //COM_new_window(name, MPI_COMM_NULL);
    COM_new_window(name, tmpComm);

    comFoamPtr->winNameVol = name;

    //MPI_Comm_dup(tmpComm, &(comFoamPtr->winComm));
    comFoamPtr->winComm = tmpComm;
    
    //Foam::PstreamGlobals::MPI_comFoam_to_openFoam = comFoamPtr->winComm;
    Foam::PstreamGlobals::MPI_COMM_FOAM = comFoamPtr->winComm;
    
    MPI_Comm_rank(comFoamPtr->winComm, &(comFoamPtr->winRank));
    MPI_Comm_size(comFoamPtr->winComm, &(comFoamPtr->winNProc));

    std::string globalName = name + string(".global");

    COM_new_dataitem(globalName.c_str(), 'w', COM_VOID, 1, "");

    COM_set_object(globalName.c_str(), 0, comFoamPtr);

    COM_window_init_done(name); 

    comFoamPtr->flowRegister();

    return;
}
//---------------------------------------------------------


//^^^^^ UNLOAD MODULES ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
void rhoCentral::unload(const std::string &name)
{
    Foam::Info << "rocFoam.unload: Unloading rocRhoCentral with name "
               << name << "." << Foam::endl;

    rhoCentral *comFoamPtr = NULL;
    std::string globalName(name+".global");

    COM_get_object(globalName.c_str(), 0, &comFoamPtr);

    //comFoamPtr->finalize();
    delete comFoamPtr;

    COM_delete_window(std::string(name));
}
//---------------------------------------------------------

//===================================================================


//^^^ DEFINITION OF OPENFOAM-RELATED MTHODS ^^^^^^^^^^^^^^^^^^^^^^^^^
int rhoCentral::initialize(int argc, char *argv[])
{
#define NO_CONTROL

    // Not quite sure where this line should be
    createArgs(argc, argv);

    //  postProcess.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    PostProcess(argc, argv);
    // ---------------------------------------------------

    //  setRootCaseLists.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    setRootCaseLists();
    // ---------------------------------------------------

    //  createTime.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    createTime();
    // ---------------------------------------------------

    //  createDynamicFvMesh.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^
    createDynamicFvMesh();
    // ---------------------------------------------------

    //  createFields.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    createFields();
    // ---------------------------------------------------

    //  createFieldRefs.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    createFieldRefs();
    // ---------------------------------------------------

    //  createTimeControls.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    createTimeControls();
    // ---------------------------------------------------

    //  readFluxScheme.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    readFluxScheme();
    // ---------------------------------------------------

    compressible::turbulenceModel &turbulence(*turbulencePtr);

    turbulence.validate();

    Foam::Info << "End of initialization of rhoCentral." << Foam::endl;

    // Setting comFoam variables that will be registered ^^
    //   with COM. This are needed for flow control when
    //   solver runs step-by-step.
    Foam::Time &runTime(*runTimePtr);
    winTime = runTime.value();
    winDeltaT = runTime.deltaTValue();
    winRun = static_cast<int>(runTime.run());
    //-----------------------------------------------------

    initializeStat = 0;
    return initializeStat;
}


int rhoCentral::loop()
{
    dynamicFvMesh &mesh(*meshPtr);
    Foam::Time &runTime(*runTimePtr);
    volScalarField &p(*pPtr);
    volVectorField &U(*UPtr);
    volVectorField &rhoU(*rhoUPtr);
    const volScalarField &T(*TPtr);
    const volScalarField &psi(*psiPtr);
    volScalarField &e(*ePtr);
    volScalarField &rho(*rhoPtr);
    volScalarField &rhoE(*rhoEPtr);
    surfaceScalarField &pos(*posPtr);
    surfaceScalarField &neg(*negPtr);
    surfaceScalarField &phi(*phiPtr);
    Foam::psiThermo &thermo(*pThermoPtr);
    compressible::turbulenceModel &turbulence(*turbulencePtr);

    dimensionedScalar v_zero("v_zero", dimVolume / dimTime, 0.0);

    // Courant numbers used to adjust the time-step
    // scalar CoNum = 0.0;
    // scalar meanCoNum = 0.0;

    Info << "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        //  readTimeControls.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        readTimeControls();
        // ---------------------------------------------------

        if (!LTS)
        {
            //  setDeltaT.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            setDeltaT();
            // ---------------------------------------------
            runTime++;

            // Do any mesh changes
            mesh.update();
        }

        // --- Directed interpolation of primitive fields onto faces

        surfaceScalarField rho_pos(interpolate(rho, pos));
        surfaceScalarField rho_neg(interpolate(rho, neg));

        surfaceVectorField rhoU_pos(interpolate(rhoU, pos, U.name()));
        surfaceVectorField rhoU_neg(interpolate(rhoU, neg, U.name()));

        volScalarField rPsi("rPsi", 1.0 / psi);
        surfaceScalarField rPsi_pos(interpolate(rPsi, pos, T.name()));
        surfaceScalarField rPsi_neg(interpolate(rPsi, neg, T.name()));

        surfaceScalarField e_pos(interpolate(e, pos, T.name()));
        surfaceScalarField e_neg(interpolate(e, neg, T.name()));

        surfaceVectorField U_pos("U_pos", rhoU_pos / rho_pos);
        surfaceVectorField U_neg("U_neg", rhoU_neg / rho_neg);

        surfaceScalarField p_pos("p_pos", rho_pos * rPsi_pos);
        surfaceScalarField p_neg("p_neg", rho_neg * rPsi_neg);

        surfaceScalarField phiv_pos("phiv_pos", U_pos & mesh.Sf());
        surfaceScalarField phiv_neg("phiv_neg", U_neg & mesh.Sf());

        // Make fluxes relative to mesh-motion
        if (mesh.moving())
        {
            phiv_pos -= mesh.phi();
            phiv_neg -= mesh.phi();
        }

        volScalarField c("c", sqrt(thermo.Cp() / thermo.Cv() * rPsi));
        surfaceScalarField cSf_pos
        (
            "cSf_pos",
            interpolate(c, pos, T.name()) * mesh.magSf()
        );

        surfaceScalarField cSf_neg
        (
            "cSf_neg",
            interpolate(c, neg, T.name()) * mesh.magSf()
        );

        surfaceScalarField ap
        (
            "ap",
            max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero)
        );
        
        surfaceScalarField am
        (
            "am",
            min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero)
        );

        surfaceScalarField a_pos
        (
            "a_pos",
            ap / (ap - am)
        );

        // surfaceScalarField amaxSf("amaxSf", max(mag(am), mag(ap)));
        if (amaxSfPtr == NULL)
        {
            amaxSfPtr = new surfaceScalarField("amaxSf", max(mag(am), mag(ap)));
        }
        
        surfaceScalarField &amaxSf(*amaxSfPtr);

        amaxSf = max(mag(am), mag(ap));

        surfaceScalarField aSf("aSf", am * a_pos);

        if (fluxScheme == "Tadmor")
        {
            aSf = -0.5 * amaxSf;
            a_pos = 0.5;
        }

        surfaceScalarField a_neg("a_neg", 1.0 - a_pos);

        phiv_pos *= a_pos;
        phiv_neg *= a_neg;

        surfaceScalarField aphiv_pos("aphiv_pos", phiv_pos - aSf);
        surfaceScalarField aphiv_neg("aphiv_neg", phiv_neg + aSf);

        // Reuse amaxSf for the maximum positive and negative fluxes
        // estimated by the central scheme
        amaxSf = max(mag(aphiv_pos), mag(aphiv_neg));

        //  centralCourantNo.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        centralCourantNo();
        // ---------------------------------------------------

        if (LTS)
        {
            // setRDeltaT.H
            // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            setRDeltaT();
            // -------------------------------
            runTime++;
        }

        Info << "Time = " << runTime.timeName() << nl << endl;

        phi = aphiv_pos * rho_pos + aphiv_neg * rho_neg;

        surfaceVectorField phiUp
        (
            (aphiv_pos * rhoU_pos + aphiv_neg * rhoU_neg) +
            (a_pos * p_pos + a_neg * p_neg) * mesh.Sf()
        );

        surfaceScalarField phiEp
        (
            "phiEp",
            aphiv_pos * (rho_pos * (e_pos + 0.5 * magSqr(U_pos)) + p_pos) +
            aphiv_neg * (rho_neg * (e_neg + 0.5 * magSqr(U_neg)) + p_neg) +
            aSf * p_pos - aSf * p_neg
        );

        // Make flux for pressure-work absolute
        if (mesh.moving())
        {
            phiEp += mesh.phi() * (a_pos * p_pos + a_neg * p_neg);
        }

        volScalarField muEff("muEff", turbulence.muEff());
        volTensorField tauMC("tauMC", muEff * dev2(Foam::T(fvc::grad(U))));

        // --- Solve density
        solve(fvm::ddt(rho) + fvc::div(phi));

        // --- Solve momentum
        solve(fvm::ddt(rhoU) + fvc::div(phiUp));

        U.ref() = rhoU() / rho();

        U.correctBoundaryConditions();
        rhoU.boundaryFieldRef() == rho.boundaryField() * U.boundaryField();

        if (!inviscid)
        {
            solve
            (
                fvm::ddt(rho, U) - fvc::ddt(rho, U)
              - fvm::laplacian(muEff, U)
              - fvc::div(tauMC)
            );
            rhoU = rho * U;
        }

        // --- Solve energy
        surfaceScalarField sigmaDotU
        (
            "sigmaDotU",
            (
                fvc::interpolate(muEff) * mesh.magSf() * fvc::snGrad(U) +
                fvc::dotInterpolate(mesh.Sf(), tauMC)
            ) & (a_pos * U_pos + a_neg * U_neg)
        );

        solve
        (
            fvm::ddt(rhoE)
          + fvc::div(phiEp)
          - fvc::div(sigmaDotU)
        );


        e = rhoE / rho - 0.5 * magSqr(U);
        e.correctBoundaryConditions();
        thermo.correct();
        rhoE.boundaryFieldRef() ==
            rho.boundaryField() *
                (e.boundaryField() + 0.5 * magSqr(U.boundaryField()));

        if (!inviscid)
        {
            solve
            (
                fvm::ddt(rho, e) - fvc::ddt(rho, e) -
                fvm::laplacian(turbulence.alphaEff(), e)
            );
            thermo.correct();
            rhoE = rho * (e + 0.5 * magSqr(U));
        }

        p.ref() = rho() / psi();

        p.correctBoundaryConditions();
        rho.boundaryFieldRef() == psi.boundaryField() * p.boundaryField();

        turbulence.correct();

        runTime.write();

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl
             << endl;
    }

    Info << "End\n" << endl;

    loopStat = 0;
    return loopStat;
}



int rhoCentral::centralizedIO()
{

    dynamicFvMesh &mesh(*meshPtr);
    Foam::Time &runTime(*runTimePtr);
    volScalarField &p(*pPtr);
    volVectorField &U(*UPtr);
    volVectorField &rhoU(*rhoUPtr);
    const volScalarField &T(*TPtr);
    const volScalarField &psi(*psiPtr);
    volScalarField &e(*ePtr);
    volScalarField &rho(*rhoPtr);
    volScalarField &rhoE(*rhoEPtr);
    surfaceScalarField &pos(*posPtr);
    surfaceScalarField &neg(*negPtr);
    surfaceScalarField &phi(*phiPtr);
    Foam::psiThermo &thermo(*pThermoPtr);
    compressible::turbulenceModel &turbulence(*turbulencePtr);



Foam::Info << "P size = " << p.size() << endl;
Foam::Info << "T size = " << T.size() << endl;
Foam::Info << "rho size = " << rho.size() << endl;
Foam::Info << "U size = " << U.size() << endl;







   //  Registering node data of this module to COM ^^^^^^^


    //dynamicFvMesh &mesh(*meshPtr);

    /*
    {
        const pointField &points = mesh.points();

        nComponents = 3;
        nNodes = mesh.nPoints();
        
        int nTotal = nComponents * nNodes;
        comXYZ = new double[nTotal];

        forAll(points, i)
        {
            for(int j=0; j< nComponents; j++)
            *(comXYZ + i * nComponents) = points[i][j];
        }
    }
    */

    {
        const pointField &points = mesh.points();
        int nPoints = mesh.nPoints();

        std::vector< std::vector<double> > vecPoints;
        std::vector< std::vector<double> > vecVel;
        std::vector<double> vecRho;
        std::vector<double> vecP;
        std::vector<double> vecT;
        
        std::vector<double> vecTmpDouble;
        forAll(points, i)
        {
            vecTmpDouble.clear();
            for(int j=0; j<nComponents; j++)
            {
                vecTmpDouble.push_back(points[i][j]);
            }
            vecPoints.push_back(vecTmpDouble);
            
            vecTmpDouble.clear();
            for(int j=0; j<nComponents; j++)
            {
                vecTmpDouble.push_back(U[i].component(j));
            }
            vecVel.push_back(vecTmpDouble);

            vecRho.push_back(rho[i]);
            vecP.push_back(p[i]);
            vecT.push_back(T[i]);
        }

/*        for (std::vector<std::vector<double>>::iterator vvIt=vecPoints.begin(); vvIt <= vecPoints.begin()+2; vvIt++)
        {
            Foam::Info << "Ponits Vector = ";
            for (std::vector<double>::iterator vIt=vvIt->begin(); vIt != vvIt->end(); vIt++)
            {
                Foam::Info << *vIt << " ";
            }
            Foam::Info << Foam::endl;
        }
        Foam::Info << Foam::endl;

        for (std::vector<std::vector<double>>::iterator vvIt=vecVel.begin(); vvIt <= vecVel.begin()+2; vvIt++)
        {
            Foam::Info << "Vel Vector    = ";
            for (std::vector<double>::iterator vIt=vvIt->begin(); vIt != vvIt->end(); vIt++)
            {
                Foam::Info << *vIt << " ";
            }
            Foam::Info << Foam::endl;
        }
        Foam::Info << Foam::endl;

        Foam::Info << "Rho = ";
        for (std::vector<double>::iterator vIt=vecRho.begin(); vIt <= vecRho.begin()+2; vIt++)
        {
            Foam::Info << *vIt << " ";
        }
        Foam::Info << Foam::endl;

        Foam::Info << "P = ";
        for (std::vector<double>::iterator vIt=vecP.begin(); vIt <= vecP.begin()+2; vIt++)
        {
            Foam::Info << *vIt << " ";
        }
        Foam::Info << Foam::endl;

        Foam::Info << "T = ";
        for (std::vector<double>::iterator vIt=vecT.begin(); vIt <= vecT.begin()+2; vIt++)
        {
            Foam::Info << *vIt << " ";
        }
        Foam::Info << Foam::endl;
*/
    }


    {
        const faceList& faces = mesh.faces();
        const labelList& faceOwner = mesh.faceOwner();
        const labelList& faceNeighb = mesh.faceNeighbour();
        int nFaces = faces.size();
        
        std::vector< std::vector<int> > vecFaces;
        std::vector<int> vecOwner;
        std::vector<int> vecNeighb;
        std::vector<int> vecTmpInt;
        
        forAll(faces, i)
        {
            const labelList &points = faces[i];
            
            vecOwner.push_back(faceOwner[i]);
            vecNeighb.push_back(faceNeighb[i]);

            vecTmpInt.clear();
            forAll(points, j)
            {
                vecTmpInt.push_back(points[j]);
            }
            
            vecFaces.push_back(vecTmpInt);
        }

        for (std::vector<std::vector<int>>::iterator vvIt=vecFaces.begin(); vvIt <= vecFaces.begin()+10; vvIt++) //vvIt != vecPoints.end(); vvIt++)
        {
            int index=std::distance(vecFaces.begin(), vvIt);
            Foam::Info << "Face points = ";
            for (std::vector<int>::iterator vIt=vvIt->begin(); vIt != vvIt->end(); vIt++)
            {
                Foam::Info << *vIt << " ";
            }
            
            Foam::Info << " | " << vecOwner[index] <<  " " << vecNeighb[index];
            
            Foam::Info << Foam::endl;
        }
    }

    {
        const polyBoundaryMesh &patches = mesh.boundaryMesh();
        const faceList& faces = mesh.faces();

        const int &nPatches = patches.size();

        std::vector< std::vector< std::vector<double> >> patchVel;
        std::vector< std::vector<double> > patchRho;
        std::vector< std::vector<double> > patchP;
        std::vector< std::vector<double> > patchT;

        std::vector<std::string> vecPatchName;
        std::vector<std::string> vecPatchType;
        // Note: I need to change this
        std::vector< wordList > vecPatchInGroup;

        std::vector<int> vecPatchStart;
        std::vector<int> vecPatchNFaces;
        std::vector<int> vecPatchNPoints;

        


        std::vector< std::vector<int> > patchGlobalPointIndex;
        std::vector< std::vector< std::vector<int> >> patchLocalConnectivity;


        forAll(patches, ipatch)
        {

//if (ipatch >= 1) break;

            const polyPatch &patch = patches[ipatch];

            const word &patchName = patch.name();
            const word &patchType = patch.type();
            const wordList &patchInGroup = patch.inGroups();

            const label& patchStart = patch.start();
            const int& patchSize = patch.size();

            vecPatchName   .push_back(patchName);
            vecPatchType   .push_back(patchType);
            vecPatchInGroup.push_back(patchInGroup);
            vecPatchStart  .push_back(patchStart);
            vecPatchNFaces .push_back(patchSize);

            //std::vector<int> vecTmpInt_globalFaceIndex;
            //std::vector<int> vecTmpInt_localFaceIndex;
            
            std::vector<int> vecTmpInt_patchGlobalPointIndex;

            //std::vector<int> VecTmpInt_globalConnectivity;
            std::vector<std::vector<int>> vecTmpInt_patchLocalConnectivity;

            Foam::Info << " Patch " << patchName << " patchStart = " << patchStart << " patchSize = " << patchSize << endl;

            //VecTmpInt_globalFaceIndex.clear();

            vecTmpInt_patchGlobalPointIndex.clear();
            vecTmpInt_patchLocalConnectivity.clear();
            for(int iface=0; iface<patchSize; iface++) //iface<5; iface++) 
            {
                const label& face = patchStart + iface; //patch.whichFace(patchStart + iFace);
                
                //vecTmpInt_localFaceIndex.push_back(iface);
                //vecTmpInt_globalFaceIndex.push_back(face);

                const labelList &points = faces[face];
                int nPoints = points.size();

//Foam::Info << "      face " << " local-to-global " << iface << endl << "       ";

                
                std::vector<int> vecTmpInt;
                //vecTmpInt.clear();
                forAll(points, ipoint)
                {
                    const int& point = points[ipoint];

//Foam::Info << point << " ";


                    //VecTmpInt_globalConnectivity.push_back( point );

                    std::vector<int>::iterator index = std::find
                                      (
                                        vecTmpInt_patchGlobalPointIndex.begin(),
                                        vecTmpInt_patchGlobalPointIndex.end(),
                                        point
                                      );

                    if  ( index == vecTmpInt_patchGlobalPointIndex.end() )
                    {
                        vecTmpInt_patchGlobalPointIndex.push_back(point);

                        index = vecTmpInt_patchGlobalPointIndex.end() - 1;
                    }

                    int indexVal = std::distance(vecTmpInt_patchGlobalPointIndex.begin(), index);
                    vecTmpInt.push_back( indexVal );
                }
                vecTmpInt_patchLocalConnectivity.push_back(vecTmpInt);


//Foam::Info << endl;

            }

            patchGlobalPointIndex.push_back(vecTmpInt_patchGlobalPointIndex);
            patchLocalConnectivity.push_back(vecTmpInt_patchLocalConnectivity);

            vecPatchNPoints.push_back(vecTmpInt_patchGlobalPointIndex.size()) ;
        }


Foam::Info << endl;

        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {

if (vecPatchName[ipatch] != string("outlet")) continue;

            Foam::Info << " Patch " << vecPatchName[ipatch]
                       << " patchStart = " << vecPatchStart[ipatch]
                       << " nFaces = " << vecPatchNFaces[ipatch]
                       << " nPoints = " << vecPatchNPoints[ipatch]
                       << endl;
        }



        //  Save Flow Quantities
        std::vector< std::vector<double> > vecVecTmpDouble;
        std::vector<double> vecTmpDouble;
        

        patchVel.clear();
        patchRho.clear();
        patchP.clear();
        patchT.clear();
        for(int ipatch=0; ipatch<nPatches ; ipatch++)
        {

            int faceStart = vecPatchStart[ipatch];
            int nFaces = vecPatchNFaces[ipatch];

            vecVecTmpDouble.clear();
            if (vecPatchType[ipatch] != string("empty"))
            {
                for(int iface=0; iface<nFaces; iface++)
                {
Foam::Info << ipatch << iface;

                    vecTmpDouble.clear();
                    for(int k=0; k<nComponents; k++)
                    {
                        vecTmpDouble.push_back( U.boundaryField()[ipatch][iface].component(k) );
                    }
                    vecVecTmpDouble.push_back(vecTmpDouble);

Foam::Info << " " << ipatch << iface << endl;

                }
            }
            patchVel.push_back(vecVecTmpDouble);



            //  Patch Face Density
            vecTmpDouble.clear();
            if (vecPatchType[ipatch] != string("empty"))
            {
                for(int iface=0; iface<nFaces; iface++)
                {
                    vecTmpDouble.push_back( rho.boundaryField()[ipatch][iface] );
                }
            }
            patchRho.push_back(vecTmpDouble);

            //  Patch Face Pressure
            vecTmpDouble.clear();
            if (vecPatchType[ipatch] != string("empty"))
            {
                for(int iface=0; iface<nFaces; iface++)
                {
                    vecTmpDouble.push_back( p.boundaryField()[ipatch][iface] );
                }
            }
            patchP.push_back(vecTmpDouble);


            //  Patch Face Temperature
            vecTmpDouble.clear();
            if (vecPatchType[ipatch] != string("empty"))
            {
                for(int iface=0; iface<nFaces; iface++)
                {
                    vecTmpDouble.push_back( T.boundaryField()[ipatch][iface] );
                }
            }
            patchT.push_back(vecTmpDouble);



        }


        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {
            Foam::Info << " Patch = " << vecPatchName[ipatch] << " patchID = " << ipatch << endl;

            int faceStart = vecPatchStart[ipatch];
            int nFaces =  patchVel[ipatch].size(); // vecPatchNFaces[ipatch];
            
            Foam::Info << " Global Face Indices = ";
            for(int iface=0; iface<nFaces; iface++)
            {
                Foam::Info << faceStart+iface << " ";
            }
            Foam::Info << endl;

            
            for(int iface=0; iface<nFaces; iface++)
            {
                Foam::Info  << "     "
                            << "faceID = " << iface
                            << " V = " << patchVel[ipatch][iface][0]
                            << ", " << patchVel[ipatch][iface][1] 
                            << ", " << patchVel[ipatch][iface][2]
                            << " Rho = " << patchRho[ipatch][iface]
                            << " P = " << patchP[ipatch][iface] 
                            << " T = " << patchT[ipatch][iface]
                            << endl;
            }
            Foam::Info << endl;
        }
        
Foam::Info << endl;






/*
        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {
            int npoints = vecPatchNPoints[ipatch];

            //  Patch Nodeal Velocity
            vecVecTmpDouble.clear();
            for(int ipoint=0; ipoint<npoints; ipoint++)
            {
                int index = patchGlobalPointIndex[ipatch][ipoint];
            
                vecTmpDouble.clear();
                for(int k=0; k<nComponents; k++)
                {
                    vecTmpDouble.push_back(U[index].component(k));
                }
                vecVecTmpDouble.push_back(vecTmpDouble);
            }
            patchVel.push_back(vecVecTmpDouble);

            //  Patch Nodeal Density
            vecTmpDouble.clear();
            for(int ipoint=0; ipoint<npoints; ipoint++)
            {
                int index = patchGlobalPointIndex[ipatch][ipoint];
                vecTmpDouble.push_back(rho[index]);
            }
            patchRho.push_back(vecTmpDouble);

            //  Patch Nodeal Pressure
            vecTmpDouble.clear();
            for(int ipoint=0; ipoint<npoints; ipoint++)
            {
                int index = patchGlobalPointIndex[ipatch][ipoint];
                vecTmpDouble.push_back(p[index]);
            }
            patchP.push_back(vecTmpDouble);

            //  Patch Nodeal Temperature
            vecTmpDouble.clear();
            for(int ipoint=0; ipoint<npoints; ipoint++)
            {
                int index = patchGlobalPointIndex[ipatch][ipoint];
                vecTmpDouble.push_back(T[index]);
            }
            patchT.push_back(vecTmpDouble);
        }

        for( std::vector<std::vector<int>>::iterator itVecVec=patchGlobalPointIndex.begin(); itVecVec != patchGlobalPointIndex.end(); itVecVec++)
        {
            int ipatch = std::distance( patchGlobalPointIndex.begin(), itVecVec);

if (vecPatchName[ipatch] != string("outlet")) continue;

            Foam::Info << " Patch " << vecPatchName[ipatch] << " local-to-global " << ipatch << endl;

            int npoints = vecPatchNPoints[ipatch];
            Foam::Info << " Global Node Indices = ";
            for(int ipoint=0; ipoint<npoints; ipoint++)
            {
                Foam::Info << patchGlobalPointIndex[ipatch][ipoint] << " ";
            }
            Foam::Info << endl;

            
            for( std::vector<int>::iterator itVec=itVecVec->begin(); itVec != itVecVec->end(); itVec++)
            {
                int ipoint = std::distance( itVecVec->begin(), itVec);
            
            
                Foam::Info  << "     "
                           << "PointID = " << *itVec
                           << " V = " << patchVel[ipatch][ipoint][1]
                              << ", " << patchVel[ipatch][ipoint][2] 
                              << ", " << patchVel[ipatch][ipoint][3]
                           << " Rho = " << patchRho[ipatch][ipoint]
                           << " P = " << patchP[ipatch][ipoint] 
                           << " T = " << patchT[ipatch][ipoint]  << endl;
            }
            Foam::Info << endl;
        }
        
        
*/
Foam::Info << endl;

/*
        for (
                std::vector<std::vector<std::vector<int>>>::iterator
                    itVecVecVec=patchLocalConnectivity.begin();
                itVecVecVec != patchLocalConnectivity.end();
                itVecVecVec++
            )
        {
            int index = std::distance( 
                                        patchLocalConnectivity.begin(),
                                        itVecVecVec
                                     );

            Foam::Info << " Patch " << vecPatchName[index] << " connectivity " << endl;

            for( std::vector<std::vector<int>>::iterator itVecVec=itVecVecVec->begin(); itVecVec != itVecVecVec->end(); itVecVec++)
            {
                int faceindex = std::distance( itVecVecVec->begin(), itVecVec);
                Foam::Info << "      face " << " local-to-global " << faceindex << endl << "       ";
                
                for( std::vector<int>::iterator itVec=itVecVec->begin(); itVec != itVecVec->end(); itVec++)
                {
                    Foam::Info << *itVec << " ";
                }
                Foam::Info << endl;
            }
        }
        
        Foam::Info << endl;
*/

    }

/*
            
            //  Patch Nodeal Velocity
            vecVecTmp.clear();
            for(int j=0; j<patchSize; j++)
            {
                int index = patchStart+j;

                vecTmp.clear();
                for(int k=0; k<nComponents; k++)
                {
                    vecTmp.push_back(U[index].component(k));
                }

                vecVecTmp.push_back(vecTmp);
            }
            vecPatchVel.push_back(vecVecTmp);



            //  Patch Nodeal Density
            vecTmp.clear();
            for(int j=0; j<patchSize; j++)
            {
                int index=patchStart+j;
                vecTmp.push_back(rho[index]);
            }
            vecPatchRho.push_back(vecTmp);

            //  Patch Nodeal Pressure
            vecTmp.clear();
            for(int j=0; j<patchSize; j++)
            {
                int index=patchStart+j;
                vecTmp.push_back(p[index]);
            }
            vecPatchP.push_back(vecTmp);

            //  Patch Nodeal Temperature
            vecTmp.clear();
            for(int j=0; j<patchSize; j++)
            {
                int index=patchStart+j;
                vecTmp.push_back(T[index]);
            }
            vecPatchT.push_back(vecTmp);


            //Foam::Info << i << " " << patchName << " " 
            //           << patchType << " " << patchInGroup 
            //           << " " << patchStart << " " << patchSize << endl;
        }
        
        Foam::Info <<  " SIZE OF T ==== " << vecPatchT.size() << endl;


        for(std::vector< std::vector<double> >::iterator vvI=vecPatchRho.begin(); vvI != vecPatchRho.end(); vvI++)
        {
            Foam::Info << " Rho of partch " << std::distance(vecPatchRho.begin(), vvI) << " = ";
            for(std::vector<double>::iterator vI=vvI->begin(); vI != vvI->begin()+5; vI++)
            {
                    Foam::Info << *vI << " ";
            }
            Foam::Info << endl;;
        }
*/



    //dataName = name+string(".nc");
    //COM_set_size( dataName.c_str(), 0, nNodes);
    //COM_set_array(dataName.c_str(), 0, comXYZ, 3);
        
    

    /*
    const faceList &faces = mesh.faces();
    int nFaces = faces.size();
    comFaces  = new int[nFaces];
    comOwner  = new int[nFaces];
    comNeighb = new int[nFaces];
    */




/*

    const cellList &cells = mesh.cells();
    const faceList &faces = mesh.faces();
    int nCells = cells.size();

    std::vector<int> vecPoints;
    std::vector< std::vector<int> > vecFaces;
    std::vector< std::vector< std::vector<int> > > vecCells;


Foam::Info << __LINE__ << endl;


    forAll(cells, i)
    {

        vecFaces.clear();

        const faceList &faces = mesh.faces();
        forAll(faces, j)
        {



            const labelList &points = faces[j];
            
            vecPoints.clear();
            forAll(points, k)
            {            
                vecPoints.push_back( points[k] );
            }
            vecFaces.push_back( vecPoints );

        }

//        vecCells.push_back( vecFaces );


    }

Foam::Info << __LINE__ << endl;

*/

//    std::vector< std::vector< std::vector<int> > >::iterator vvvIt = vecCells.begin();

/*    for (std::vector< std::vector<int> >::iterator vvIt=vvvIt->begin(); vvIt != vvvIt->end(); vvIt++)
    {
        for (std::vector<int>::iterator vIt=vvIt->begin(); vIt != vvIt->end(); vIt++)
        {
            Foam::Info << "Ponits Vector = " << *vIt;
        }
        Foam::Info << Foam::endl;
    }
*/


    //std::cout << "Cell Vector   = " << *vvvIt << endl;
    //Foam::Info << "Face Vector   = " << *vvIt << Foam::endl;
    //Foam::Info << "Ponits Vector = " << *vIt << Foam::endl;
    


return 0;
























}



int rhoCentral::step()
{
    dynamicFvMesh &mesh(*meshPtr);
    Foam::Time &runTime(*runTimePtr);
    volScalarField &p(*pPtr);
    volVectorField &U(*UPtr);
    volVectorField &rhoU(*rhoUPtr);
    const volScalarField &T(*TPtr);
    const volScalarField &psi(*psiPtr);
    volScalarField &e(*ePtr);
    volScalarField &rho(*rhoPtr);
    volScalarField &rhoE(*rhoEPtr);
    surfaceScalarField &pos(*posPtr);
    surfaceScalarField &neg(*negPtr);
    surfaceScalarField &phi(*phiPtr);
    Foam::psiThermo &thermo(*pThermoPtr);
    compressible::turbulenceModel &turbulence(*turbulencePtr);

    dimensionedScalar v_zero("v_zero", dimVolume / dimTime, 0.0);


centralizedIO();


    // Courant numbers used to adjust the time-step
    // scalar CoNum = 0.0;
    // scalar meanCoNum = 0.0;

    //Info << "\nStarting time loop\n" << endl;

    //while (runTime.run())
    {
        //  readTimeControls.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        readTimeControls();
        // ---------------------------------------------------

        if (!LTS)
        {
            //  setDeltaT.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            setDeltaT();
            // ---------------------------------------------
            runTime++;

            // Do any mesh changes
            mesh.update();
        }

        // --- Directed interpolation of primitive fields onto faces

        surfaceScalarField rho_pos(interpolate(rho, pos));
        surfaceScalarField rho_neg(interpolate(rho, neg));

        surfaceVectorField rhoU_pos(interpolate(rhoU, pos, U.name()));
        surfaceVectorField rhoU_neg(interpolate(rhoU, neg, U.name()));

        volScalarField rPsi("rPsi", 1.0 / psi);
        surfaceScalarField rPsi_pos(interpolate(rPsi, pos, T.name()));
        surfaceScalarField rPsi_neg(interpolate(rPsi, neg, T.name()));

        surfaceScalarField e_pos(interpolate(e, pos, T.name()));
        surfaceScalarField e_neg(interpolate(e, neg, T.name()));

        surfaceVectorField U_pos("U_pos", rhoU_pos / rho_pos);
        surfaceVectorField U_neg("U_neg", rhoU_neg / rho_neg);

        surfaceScalarField p_pos("p_pos", rho_pos * rPsi_pos);
        surfaceScalarField p_neg("p_neg", rho_neg * rPsi_neg);

        surfaceScalarField phiv_pos("phiv_pos", U_pos & mesh.Sf());
        surfaceScalarField phiv_neg("phiv_neg", U_neg & mesh.Sf());

        // Make fluxes relative to mesh-motion
        if (mesh.moving())
        {
            phiv_pos -= mesh.phi();
            phiv_neg -= mesh.phi();
        }

        volScalarField c("c", sqrt(thermo.Cp() / thermo.Cv() * rPsi));
        surfaceScalarField cSf_pos
        (
            "cSf_pos",
            interpolate(c, pos, T.name()) * mesh.magSf()
        );

        surfaceScalarField cSf_neg
        (
            "cSf_neg",
            interpolate(c, neg, T.name()) * mesh.magSf()
        );

        surfaceScalarField ap
        (
            "ap",
            max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero)
        );
        
        surfaceScalarField am
        (
            "am",
            min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero)
        );

        surfaceScalarField a_pos
        (
            "a_pos",
            ap / (ap - am)
        );

        // surfaceScalarField amaxSf("amaxSf", max(mag(am), mag(ap)));
        if (amaxSfPtr == NULL)
        {
            amaxSfPtr = new surfaceScalarField("amaxSf", max(mag(am), mag(ap)));
        }
        
        surfaceScalarField &amaxSf(*amaxSfPtr);

        amaxSf = max(mag(am), mag(ap));

        surfaceScalarField aSf("aSf", am * a_pos);

        if (fluxScheme == "Tadmor")
        {
            aSf = -0.5 * amaxSf;
            a_pos = 0.5;
        }

        surfaceScalarField a_neg("a_neg", 1.0 - a_pos);

        phiv_pos *= a_pos;
        phiv_neg *= a_neg;

        surfaceScalarField aphiv_pos("aphiv_pos", phiv_pos - aSf);
        surfaceScalarField aphiv_neg("aphiv_neg", phiv_neg + aSf);

        // Reuse amaxSf for the maximum positive and negative fluxes
        // estimated by the central scheme
        amaxSf = max(mag(aphiv_pos), mag(aphiv_neg));

        //  centralCourantNo.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        centralCourantNo();
        // ---------------------------------------------------

        if (LTS)
        {
            // setRDeltaT.H
            // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            setRDeltaT();
            // -------------------------------
            runTime++;
        }

        Info << "Time = " << runTime.timeName() << nl << endl;

        phi = aphiv_pos * rho_pos + aphiv_neg * rho_neg;

        surfaceVectorField phiUp
        (
            (aphiv_pos * rhoU_pos + aphiv_neg * rhoU_neg) +
            (a_pos * p_pos + a_neg * p_neg) * mesh.Sf()
        );

        surfaceScalarField phiEp
        (
            "phiEp",
            aphiv_pos * (rho_pos * (e_pos + 0.5 * magSqr(U_pos)) + p_pos) +
            aphiv_neg * (rho_neg * (e_neg + 0.5 * magSqr(U_neg)) + p_neg) +
            aSf * p_pos - aSf * p_neg
        );

        // Make flux for pressure-work absolute
        if (mesh.moving())
        {
            phiEp += mesh.phi() * (a_pos * p_pos + a_neg * p_neg);
        }

        volScalarField muEff("muEff", turbulence.muEff());
        volTensorField tauMC("tauMC", muEff * dev2(Foam::T(fvc::grad(U))));

        // --- Solve density
        solve(fvm::ddt(rho) + fvc::div(phi));

        // --- Solve momentum
        solve(fvm::ddt(rhoU) + fvc::div(phiUp));

        U.ref() = rhoU() / rho();

        U.correctBoundaryConditions();
        rhoU.boundaryFieldRef() == rho.boundaryField() * U.boundaryField();

        if (!inviscid)
        {
            solve
            (
                fvm::ddt(rho, U) - fvc::ddt(rho, U)
              - fvm::laplacian(muEff, U)
              - fvc::div(tauMC)
            );
            rhoU = rho * U;
        }

        // --- Solve energy
        surfaceScalarField sigmaDotU
        (
            "sigmaDotU",
            (
                fvc::interpolate(muEff) * mesh.magSf() * fvc::snGrad(U) +
                fvc::dotInterpolate(mesh.Sf(), tauMC)
            ) & (a_pos * U_pos + a_neg * U_neg)
        );

        solve
        (
            fvm::ddt(rhoE)
          + fvc::div(phiEp)
          - fvc::div(sigmaDotU)
        );


        e = rhoE / rho - 0.5 * magSqr(U);
        e.correctBoundaryConditions();
        thermo.correct();
        rhoE.boundaryFieldRef() ==
            rho.boundaryField() *
                (e.boundaryField() + 0.5 * magSqr(U.boundaryField()));

        if (!inviscid)
        {
            solve
            (
                fvm::ddt(rho, e) - fvc::ddt(rho, e) -
                fvm::laplacian(turbulence.alphaEff(), e)
            );
            thermo.correct();
            rhoE = rho * (e + 0.5 * magSqr(U));
        }

        p.ref() = rho() / psi();

        p.correctBoundaryConditions();
        rho.boundaryFieldRef() == psi.boundaryField() * p.boundaryField();

        turbulence.correct();

        runTime.write();

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl
             << endl;
    }

    //Info << "End\n" << endl;

    // Setting comFoam variables that will be registered ^^
    //   with COM. This are needed for flow control when
    //   solver runs step-by-step.
    winTime = runTime.value();
    winDeltaT = runTime.deltaTValue();
    winRun = static_cast<int>(runTime.run());
    //-----------------------------------------------------

    stepStat = 0;
    return stepStat;
}

int rhoCentral::createFields()
{
    Foam::Time &runTime(*runTimePtr);
    dynamicFvMesh &mesh(*meshPtr);

    //  createRDeltaT.H  ^^^^^^^^^^^^^
    createRDeltaT();
    // -------------------------------

    Info << "Reading thermophysical properties\n" << endl;

    pThermoPtr = autoPtr<psiThermo>(psiThermo::New(mesh));
    Foam::psiThermo &thermo(*pThermoPtr);

    ePtr = &thermo.he();
    volScalarField &e(*ePtr);

    Info << "Reading field U\n" << endl;

    UPtr = new volVectorField
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    volVectorField &U(*UPtr);

    rhoPtr = new volScalarField
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo.rho()
    );
    volScalarField &rho(*rhoPtr);

    rhoUPtr = new volVectorField
    (
        IOobject
        (
            "rhoU",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho * U
    );
    volVectorField &rhoU(*rhoUPtr);

    rhoEPtr = new volScalarField
    (
        IOobject
        (
            "rhoE",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho * (e + 0.5 * magSqr(U))
    );

    posPtr = new surfaceScalarField
    (
        IOobject
        (
            "pos",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimless, 1.0)
    );

    negPtr = new surfaceScalarField
    (
        IOobject
        (
            "neg",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimless, -1.0)
    );

    phiPtr = new surfaceScalarField
    (
        "phi",
        fvc::flux(rhoU)
    );
    surfaceScalarField &phi(*phiPtr);

    Info << "Creating turbulence model\n" << endl;

    turbulencePtr = autoPtr<compressible::turbulenceModel>
    (
        compressible::turbulenceModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );

    return 0;
}

int rhoCentral::createFieldRefs()
{
    Foam::psiThermo &thermo(*pThermoPtr);

    /*
    volScalarField& p = thermo.p();
    const volScalarField& T = thermo.T();
    const volScalarField& psi = thermo.psi(); */
    const volScalarField &mu = thermo.mu();

    pPtr = &thermo.p();
    TPtr = &thermo.T();
    psiPtr = &thermo.psi();

    // bool inviscid(true);
    inviscid = true;
    if (max(mu.primitiveField()) > 0.0)
    {
        inviscid = false;
    }

    return 0;
}

int rhoCentral::readFluxScheme()
{
    dynamicFvMesh &mesh(*meshPtr);

    // word fluxScheme("Kurganov");
    word fluxScheme("Kurganov");
    if (mesh.schemesDict().readIfPresent("fluxScheme", fluxScheme))
    {
        if ((fluxScheme == "Tadmor") || (fluxScheme == "Kurganov"))
        {
            Info << "fluxScheme: " << fluxScheme << endl;
        }
        else
        {
            FatalErrorInFunction
                << "fluxScheme: " << fluxScheme
                << " is not a valid choice. "
                << "Options are: Tadmor, Kurganov"
                << abort(FatalError);
        }
    }

    return 0;
}


int rhoCentral::centralCourantNo()
{
    Foam::Time &runTime(*runTimePtr);
    dynamicFvMesh &mesh(*meshPtr);
    surfaceScalarField &amaxSf(*amaxSfPtr);

    if (mesh.nInternalFaces())
    {
        scalarField sumAmaxSf(fvc::surfaceSum(amaxSf)().primitiveField());

        CoNum =
            0.5 * gMax(sumAmaxSf / mesh.V().field()) * runTime.deltaTValue();

        meanCoNum = 0.5 * (gSum(sumAmaxSf) / gSum(mesh.V().field())) *
                    runTime.deltaTValue();
    }

    Info << "Mean and max Courant Numbers = " << meanCoNum << " " << CoNum
         << endl;

    return 0;
}

int rhoCentral::setRDeltaT()
{
    Foam::Time &runTime(*runTimePtr);
    dynamicFvMesh &mesh(*meshPtr);
    surfaceScalarField &amaxSf(*amaxSfPtr);

    volScalarField &rDeltaT = trDeltaT.ref();

    scalar rDeltaTSmoothingCoeff
    (
        runTime.controlDict().lookupOrDefault<scalar>
        (
            "rDeltaTSmoothingCoeff",
            0.02
        )
    );

    // Set the reciprocal time-step from the local Courant number
    rDeltaT.ref() = max(1 / dimensionedScalar(dimTime, maxDeltaT),
                        fvc::surfaceSum(amaxSf)()() / ((2 * maxCo) * mesh.V()));

    // Update tho boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);

    Info << "Flow time scale min/max = "
         << gMin(1 / rDeltaT.primitiveField())
         << ", " << gMax(1 / rDeltaT.primitiveField()) << endl;

    return 0;
}

rhoCentral::~rhoCentral()
{
    finalizeStat = finalize();
}

int rhoCentral::finalize()
{
    // Delete thing that are allocated here
    if (posPtr != NULL) {delete posPtr; posPtr = NULL;}
    if (negPtr != NULL) {delete negPtr; negPtr = NULL;}
    if (amaxSfPtr != NULL) {delete amaxSfPtr; amaxSfPtr = NULL;}
    if (UPtr != NULL) {delete UPtr; UPtr=NULL;}
    if (rhoPtr != NULL) {delete rhoPtr; rhoPtr=NULL;}
    if (rhoUPtr != NULL) {delete rhoUPtr; rhoUPtr=NULL;}
    if (rhoEPtr != NULL) {delete rhoEPtr; rhoEPtr=NULL;}
    if (phiPtr != NULL) {delete phiPtr; phiPtr=NULL;}

    //delete argsPtr; Let it be the last thing to delete in the
    //                parrent class:rocFoam

    //delete pThermoPtr;
    //delete runTimePtr;
    //delete ePtr; a pointer
    //delete pPtr; a pointer
    //delete TPtr; a pointer
    //delete psiPtr; a pointer

    //delete meshPtr;
    //delete turbulencePtr;
    //delete trDeltaT;

    finalizeStat = 0;
    
    return finalizeStat;
}


double rhoCentral::errorEvaluate(int argc, char *argv[])
{
    createArgs(argc, argv);
    setRootCase();
    
    Foam::argList &args(*argsPtr);
    
    if (args.optionFound("list"))
    {
        functionObjectList::list();
        return 0;
    }
    
    createTime();
    createDynamicFvMesh();

    Foam::Time &runTime(*runTimePtr);
    
    Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);
    
    std::vector<volScalarField> rhoVec;
    std::vector<volVectorField> rhoUVec;
    std::vector<volScalarField> UMagVec;
    std::vector<volScalarField> rhoEVec;

    forAll(timeDirs, timei)
    {
        runTime.setTime(timeDirs[timei], timei);

        Info << "iTime = " << timei <<", " << "Time = " << runTime.timeName() << endl;

        FatalIOError.throwExceptions();

        createFields();

        volScalarField &rho(*rhoPtr);
        volVectorField &rhoU(*rhoUPtr);
        volScalarField &rhoE(*rhoEPtr);
        
        rhoVec.push_back(rho);
        rhoUVec.push_back(rhoU);
        UMagVec.push_back(magSqr(rhoU));
        rhoEVec.push_back(rhoE);
    }

    Info << "Field vectors created. " << endl;

    rhoVec[0]  = mag(rhoVec[2]  - rhoVec[1]);
    UMagVec[0] = mag(UMagVec[2]-UMagVec[1]);
    rhoEVec[0] = mag(rhoEVec[2] - rhoEVec[1]);
    
    Info << "Infinity norm = " << max(rhoVec[0]) << endl;
    Info << "Infinity norm = " << max(UMagVec[0]) << endl;
    Info << "Infinity norm = " << max(rhoEVec[0]) << endl;

    double maxError = std::max( max(rhoVec[0]).value(), max(UMagVec[0]).value() );
    testStat = std::max( maxError , max(rhoEVec[0]).value() );
    
    return testStat;
}

//===================================================================
