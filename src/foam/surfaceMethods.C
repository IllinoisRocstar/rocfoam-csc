#include "comFoam.H"

int comFoam::createSurfaceConnectivities()
{
    const dynamicFvMesh& mesh(*meshPtr);

    // Mesh and conmnectivities ^^^^^^^^^^^^^^^^^
    const faceList& faces = mesh.faces();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    //-------------------------------------------

    // Gather the number of patches in each processor
    int nPatches = patches.size();
    ca_nPatches = new int[ca_nProc]{};

    MPI_Allgather(&nPatches, 1, MPI_INT,
                   ca_nPatches, 1, MPI_INT, winComm);
    //-------------------------------------------

    // Taking care of patchStart & patchSize ^^^^^^^^^^^^^^
    int* patchStart_ = new int[nPatches]{};
    int* patchSize_  = new int[nPatches]{};
    forAll(patches, ipatch)
    {
        const polyPatch& patch = patches[ipatch];

        patchStart_[ipatch] = patch.start();
        patchSize_[ipatch]  = patch.size();
    }

    int* displace = new int[ca_nProc]{};
    int* recvCounts = new int[ca_nProc]{};
    int totalSize{0};
    for (int iproc=0; iproc<ca_nProc; iproc++)
    {
        displace[iproc] = totalSize;
        recvCounts[iproc] = ca_nPatches[iproc];
        totalSize += recvCounts[iproc];
    }


    int nPatchesTotal{0};
    for (int iproc=0; iproc<ca_nProc; iproc++)
    {
        nPatchesTotal += ca_nPatches[iproc];
    }
    
    if (ca_patchStart == nullptr)
        ca_patchStart = new int[nPatchesTotal]{};

    if (ca_patchSize == nullptr)
        ca_patchSize = new int[nPatchesTotal]{};

    int sendCount = nPatches;
    MPI_Allgatherv(patchStart_, sendCount, MPI_INT,
                   ca_patchStart, recvCounts, displace,
                   MPI_INT, winComm);

    MPI_Allgatherv(patchSize_, sendCount, MPI_INT,
                   ca_patchSize, recvCounts, displace,
                   MPI_INT, winComm);

    delete [] patchStart_; patchStart_ = nullptr;
    delete [] patchSize_;  patchSize_  = nullptr;
    delete [] displace;    displace    = nullptr;
    delete [] recvCounts;  recvCounts  = nullptr;
    //-----------------------------------------------------

    // Taking care of patchName ^^^^^^^^^^^^^^^^^^^^^^^^^^^
    //ca_patchInGroup = new wordList*[nPatches];
    
    if (patchNameStr == nullptr)
        patchNameStr = new std::string[nPatches]{};
    
    if (patchTypeStr == nullptr)
        patchTypeStr = new std::string[nPatches]{};

    forAll(patches, ipatch)
    {
        const polyPatch& patch = patches[ipatch];
        //const wordList& patchInGroup = patch.inGroups();

        patchNameStr[ipatch] = patch.name();
        patchTypeStr[ipatch] = patch.type();
    }

    int maxNameLength{0};
    forAll(patches, ipatch)
    {
        std::string strTmp = patchNameStr[ipatch];
        maxNameLength = std::max( maxNameLength, static_cast<int>(strTmp.length())+1 );
    }

    ca_maxNameLength = new int[ca_nProc];
    MPI_Allgather(&maxNameLength, 1, MPI_INT,
                   ca_maxNameLength, 1, MPI_INT, winComm);

    // Local patch Names; CHAR type
    char* patchName_ = new char[nPatches * maxNameLength]{' '};
    forAll(patches, ipatch)
    {
        int startIndex = ipatch * maxNameLength;
        int endIndex = startIndex + patchNameStr[ipatch].length();
        int count{0};
        for (int i=startIndex; i<endIndex; i++)
        {
            patchName_[i] = patchNameStr[ipatch][count];
            count++;
        }
        patchName_[endIndex] = '\0';

        for (int i=endIndex+1; i<startIndex+maxNameLength; i++)
        {
            patchName_[i] = ' ';
        }
    }

    // patchName displacements and counts
    displace = new int[ca_nProc]{};
    recvCounts = new int[ca_nProc]{};
    totalSize = 0;
    for (int iproc=0; iproc<ca_nProc; iproc++)
    {
        displace[iproc] = totalSize;
    
        recvCounts[iproc] = ca_maxNameLength[iproc] *
                               ca_nPatches[iproc];

        totalSize += recvCounts[iproc];
    }

    // global patchNames
    if (ca_patchName == nullptr)
        ca_patchName = new char[totalSize]{' '};
    sendCount = maxNameLength * nPatches;
    MPI_Allgatherv(patchName_, sendCount, MPI_CHAR,
                   ca_patchName, recvCounts, displace,
                   MPI_CHAR, winComm);

    delete [] patchName_; patchName_ = nullptr;
    delete [] displace;   displace   = nullptr;
    delete [] recvCounts; recvCounts = nullptr;
    //-----------------------------------------------------
    
    // Taking care of patchType ^^^^^^^^^^^^^^^^^^^^^^^^^^^
    int maxTypeLength{0};
    forAll(patches, ipatch)
    {
        std::string strTmp = patchTypeStr[ipatch];
        maxTypeLength = std::max( maxTypeLength, static_cast<int>(strTmp.length())+1 );
    }
    

    ca_maxTypeLength = new int[ca_nProc];
    MPI_Allgather(&maxTypeLength, 1, MPI_INT,
                   ca_maxTypeLength, 1, MPI_INT, winComm);

    // Local patch Type; CHAR type
    char* patchType_ = new char[nPatches * maxTypeLength]{' '};
    forAll(patches, ipatch)
    {
        int startIndex = ipatch * maxTypeLength;
        int endIndex = startIndex + patchTypeStr[ipatch].length();
        int count{0};
        for (int i=startIndex; i<endIndex; i++)
        {
            patchType_[i] = patchTypeStr[ipatch][count];
            count++;
        }
        patchType_[endIndex] = '\0';

        for (int i=endIndex+1; i<startIndex+maxTypeLength; i++)
        {
            patchType_[i] = ' ';
        }
    }

    // patchType displacements and counts
    
    displace = new int[ca_nProc]{};
    recvCounts = new int[ca_nProc]{};
    totalSize = 0;
    for (int iproc=0; iproc<ca_nProc; iproc++)
    {
        displace[iproc] = totalSize;
    
        recvCounts[iproc] = ca_maxTypeLength[iproc] *
                               ca_nPatches[iproc];

        totalSize += recvCounts[iproc];
    }
    
    // global patchNames
    if (ca_patchType == nullptr)
        ca_patchType = new char[totalSize]{' '};
    sendCount = maxTypeLength * nPatches;
    MPI_Allgatherv(patchType_, sendCount, MPI_CHAR,
                   ca_patchType, recvCounts, displace,
                   MPI_CHAR, winComm);
    delete [] patchType_; patchType_ = nullptr;
    delete [] displace;   displace   = nullptr;
    delete [] recvCounts; recvCounts = nullptr;
    //-----------------------------------------------------
    

    // Patch Connectivity Vectors ^^^^^^^^^^^^^^
    std::vector< std::map<int, std::vector<int> >> vecPatchFaceToFaceMap;
    std::vector< std::vector<int> >  vecPatchPointToPointMap;
    std::vector< std::map<int, std::vector< std::vector<int>>>>
        vecPatchFaceToPointConn;
    //-------------------------------------------

    // FaceToFace mapping vector ^^^^^^^^^^^^^^^^
    vecPatchFaceToFaceMap.clear();
    forAll(patches, ipatch)
    {
        const polyPatch& patch = patches[ipatch];
        const label& patchStart = patch.start();
        const int& patchSize = patch.size();

        std::map<int, std::vector<int>> mapTmpInt;
        for(int iface=0; iface<patchSize; iface++)
        {
            const label& faceID = patchStart + iface;
            const labelList& pointsList = faces[faceID];
                        
            int nPointsInFace = pointsList.size();

            mapTmpInt[nPointsInFace].push_back(iface);
        }
        vecPatchFaceToFaceMap.push_back(mapTmpInt);
    }
    //-------------------------------------------

    // PointToPoint mapping Vector ^^^^^^^^^^^^^^
    vecPatchPointToPointMap.clear();
    forAll(patches, ipatch)
    {
        const polyPatch& patch = patches[ipatch];

        const int& patchStart = patch.start();
        const int& patchSize = patch.size();

        std::vector<int> vecTmpInt;
        for(int iface=0; iface<patchSize; iface++)
        {
            const label& faceID = patchStart + iface;
            const labelList& pointsList = faces[faceID];

            forAll(pointsList, ipoint)
            {
                const int& pointID = pointsList[ipoint];
                auto index = std::find
                             (
                                vecTmpInt.begin(),
                                vecTmpInt.end(),
                                pointID
                             );
                if (index == vecTmpInt.end())
                {
                    vecTmpInt.push_back(pointID);
                }
            }
        }
        vecPatchPointToPointMap.push_back(vecTmpInt);
    }

    // FaceToPoint connectivity vectors ^^^^^^^^^
    vecPatchFaceToPointConn.clear();
    forAll(patches, ipatch)
    {
        const polyPatch& patch = patches[ipatch];
        const label& patchStart = patch.start();
        const int& patchSize = patch.size();
        
        const auto& pointToPointMap = vecPatchPointToPointMap[ipatch];
        

        std::map<int, std::vector< std::vector<int>>> mapVecTmpInt;

        for(int iface=0; iface<patchSize; iface++)
        {

            const label& faceID = patchStart + iface;
            const labelList& pointsList = faces[faceID];
            int nPointsInFace = pointsList.size();
            
            std::vector<int> vecTmpInt;
            forAll(pointsList, ipoint)
            {
                const int& pointID = pointsList[ipoint];
                auto indexPtr = std::find
                                (
                                    pointToPointMap.begin(),
                                    pointToPointMap.end(),
                                    pointID
                                );

                if (indexPtr == pointToPointMap.end())
                {
                    Foam::Info << "-------------Warnning-------------" << endl
                               << "Found an unregistered point: "
                               << "ipatch = " << ipatch
                               << ", iface = " << iface
                               << ", ipoint = "  << ipoint
                               << endl;
                   return -1;
                }

                int indexVal = std::distance
                                (
                                    pointToPointMap.begin(),
                                    indexPtr
                                );
                vecTmpInt.push_back(indexVal);
            }
            mapVecTmpInt[nPointsInFace].push_back(vecTmpInt);
        }
        vecPatchFaceToPointConn.push_back(mapVecTmpInt);
    }
    //-------------------------------------------

    // Create face mapping and connectivity arrays ^^^
    ca_patchFaceToPointConn_types = new int*[nPatches]{};
    ca_patchFaceToPointConn_map   = new int*[nPatches]{};
    ca_patchFaceToPointConn_size  = new int*[nPatches]{};
    ca_patchFaceToPointConn = new int**[nPatches]{};    
    
    ca_patchFaceToFaceMap = new int*[nPatches]{};
    ca_patchFaceToFaceMap_inverse = new int*[nPatches]{};
    forAll(patches, ipatch)
    {
        // Create faceToFace mapping arrays ^^^^^^^^^
        const polyPatch& patch = patches[ipatch];
        const int& nfacesTotal = patch.size();
        //if (nfacesTotal == 0)
        //    continue;
        
        ca_patchFaceToFaceMap[ipatch] = new int[nfacesTotal]{};
        ca_patchFaceToFaceMap_inverse[ipatch] = new int[nfacesTotal]{};

        const auto& mapFaceToFaceMap = vecPatchFaceToFaceMap[ipatch];

        int sortedFaceIndex = 0;
        for (auto it : mapFaceToFaceMap)
        {
            const auto& vecFaceToFaceMap = it.second;
            int nfaces = vecFaceToFaceMap.size();
            for(int iface=0; iface<nfaces; iface++)
            {
                ca_patchFaceToFaceMap[ipatch][sortedFaceIndex] =
                    vecFaceToFaceMap[iface];
                ca_patchFaceToFaceMap_inverse[ipatch][vecFaceToFaceMap[iface]] =
                    sortedFaceIndex;
                
                sortedFaceIndex++;
            }
        }
        //-------------------------------------------

        //  Create faceToPoint connectivity arrays ^^    
        int ntypes = mapFaceToFaceMap.size();
        ca_patchFaceToPointConn_types[ipatch] = new int(ntypes);
        ca_patchFaceToPointConn_map[ipatch]   = new int[ntypes]{};
        ca_patchFaceToPointConn_size[ipatch]  = new int[ntypes]{};
        ca_patchFaceToPointConn[ipatch] = new int*[ntypes]{};

        auto mapFaceToPointConn = vecPatchFaceToPointConn[ipatch];
        for (auto it=mapFaceToPointConn.begin(); it!=mapFaceToPointConn.end(); it++)
        {
            const auto& npoints = it->first;
            const auto& vecFaceToPointConn = it->second;
            int nfaces = vecFaceToPointConn.size();
            int itype = std::distance(mapFaceToPointConn.begin(), it);

            ca_patchFaceToPointConn_map[ipatch][itype]  = npoints;
            ca_patchFaceToPointConn_size[ipatch][itype] = nfaces;
    
            int nTypeConn = npoints * nfaces;
            ca_patchFaceToPointConn[ipatch][itype] = new int[nTypeConn]{};
            
            for(int iface=0; iface<nfaces; iface++)
            {
                for(int ipoint=0; ipoint<npoints; ipoint++)
                {
                    int index = ipoint+iface*npoints;
                    
                    ca_patchFaceToPointConn[ipatch][itype][index] =
                        vecFaceToPointConn[iface][ipoint] + 1; //CGNS starts form ID 1
                }
            }
        }
        //-------------------------------------------
    }
    //------------------------------------------------

    // Create pointToPoint mapping arrays ^^^^^^^
    ca_patchPointToPointMap_size  = new int*[nPatches]{};
    ca_patchPointToPointMap = new int*[nPatches]{};

    forAll(patches, ipatch)
    {
        int npoints = vecPatchPointToPointMap[ipatch].size();
        //if (npoints == 0)
        //    continue;

        ca_patchPointToPointMap_size[ipatch] = new int(npoints);
        ca_patchPointToPointMap[ipatch] = new int[npoints]{};
        
        for(int ipoint=0; ipoint<npoints; ipoint++)
        {
            ca_patchPointToPointMap[ipatch][ipoint] =
                vecPatchPointToPointMap[ipatch][ipoint];
        }
    }
    //-------------------------------------------

    return 0;
}

int comFoam::createSurfaceData()
{
    const dynamicFvMesh& mesh(*meshPtr);
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    int nPatches = patches.size();

    // Field-data
    ca_patchPoints = new double*[nPatches]{};
    // patchPointUpdated = new bool*[nPatches]{};

    ca_patchVel    = new double*[nPatches]{};
    ca_patchP      = new double*[nPatches]{};

    if (rhoPtr != nullptr)
        ca_patchRho = new double*[nPatches]{};

    if (TPtr != nullptr)
        ca_patchT   = new double*[nPatches]{};

    if (phiPtr != nullptr)
        ca_patchPhi = new double*[nPatches]{};
    
    autoPtr<surfaceVectorField>& rhoUf(rhoUfPtr);
    if (rhoUf.valid())
        ca_patchRhoUf = new double*[nPatches]{};

    // Turbulence data ^^^^^^^^^^^^^^^^^^^^^^^^^^
    const volVectorField& U(*UPtr);
#ifdef HAVE_OFE20
    IOdictionary turbProperties
    (
        IOobject
        (
            turbulenceModel::propertiesName,
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

     std::string simulationType = ITstream
                                  (
                                      turbProperties.lookup("simulationType")
                                  ).toString();
#elif defined(HAVE_OF7)
    IOdictionary turbProperties
    (
        IOobject
        (
            turbulenceModel::propertiesName,
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    word simulationType = turbProperties.lookup("simulationType");
#elif defined(HAVE_OF8)
    IOdictionary turbProperties
    ( 
        IOobject
        (
            momentumTransportModel::typeName,
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    word simulationType = turbProperties.lookup("simulationType");
#endif

    if (simulationType == "RAS")
    {
        const dictionary& subDict = turbProperties.subDict("RAS");
#ifdef HAVE_OFE20
        word RASModel = ITstream (subDict.lookup("RASModel")).toString();
#elif defined(HAVE_OF7)
        word RASModel = subDict.lookup("RASModel");
#elif defined(HAVE_OF8)
        word RASModel = subDict.lookup("model");
#endif

        if (RASModel == "kEpsilon")
        {
            ca_patchAlphaT = new double*[nPatches]{};
            ca_patchK = new double*[nPatches]{};
            ca_patchEpsilon = new double*[nPatches]{};
            ca_patchNuT = new double*[nPatches]{};
        }
        else if (RASModel == "kOmegaSST")
        {
            ca_patchAlphaT = new double*[nPatches]{};
            ca_patchK = new double*[nPatches]{};
            ca_patchOmega = new double*[nPatches]{};
            ca_patchNuT = new double*[nPatches]{};
        }
        else
        {
            FatalErrorInFunction
                << "Error: turbulence model not recongnized by the CSC Module."
                << nl << exit(FatalError);
        }
    }
    //-------------------------------------------

    // RocStar data ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ca_bcflag  = new int*[nPatches]{};
    ca_patchNf = new double*[nPatches]{};
    //ca_patchSf = new double*[nPatches];
    ca_patchTrac = new double*[nPatches]{};
    ca_patchDisp = new double*[nPatches]{};
    patchDispOld = new double*[nPatches]{};
    // patchPointUpdated = new bool*[nPatches]{};

    ca_patchMassFlux = new double*[nPatches]{};
    ca_patchFlameT   = new double*[nPatches]{};
    ca_patchMomentum = new double*[nPatches]{};
    //-------------------------------------------

    std::string dynamicSolverType = ca_dynamicSolverType;
    forAll(patches, ipatch)
    {
        if (ca_bcflag != nullptr)
        {
            ca_bcflag[ipatch] = new int(2);

            if
            (
                *ca_isDynamicFvMesh == 1 &&
                (
                    dynamicSolverType == "displacementLaplacian" ||
                    dynamicSolverType == "solidBodyDisplacementLaplacian"
                )
            )
            {

#if HAVE_OFE20
                const pointVectorField& pointDisplacement = mesh.lookupObject<displacementMotionSolver>
                        (
                            "dynamicMeshDict"
                        ).pointDisplacement();

#elif defined(HAVE_OF7) || defined(HAVE_OF8)

                const motionSolver& motion_ =
                    refCast<const dynamicMotionSolverFvMesh>(mesh).motion();

                const pointVectorField& pointDisplacement =
                        refCast<const displacementMotionSolver>(motion_).pointDisplacement();
#endif

                if (pointDisplacement.boundaryField()[ipatch].type() == movingWallTypeName)
                {
                    *ca_bcflag[ipatch] = 0;
                }

            }
        }

        // Points
        int procStartIndex{0};
        for (int iproc=0; iproc<ca_myRank; iproc++)
        {
            procStartIndex += ca_nPatches[iproc];
        }
        int index  = procStartIndex + ipatch;
        int nfaces = ca_patchSize[index];

        if (nfaces == 0)
            continue;

        int npoints = *ca_patchPointToPointMap_size[ipatch];
        int nTotal_ = npoints * nComponents;
        int nTotal  = nfaces * nComponents;
        
        ca_patchPoints[ipatch] = new double[nTotal_]{};

        // Field-data
        ca_patchVel[ipatch] = new double[nTotal]{};
        ca_patchP[ipatch]   = new double[nfaces]{};

        if (ca_patchRho != nullptr)
            ca_patchRho[ipatch] = new double[nfaces]{};

        if (ca_patchT != nullptr)
            ca_patchT[ipatch] = new double[nfaces]{};

        if (ca_patchPhi != nullptr)
            ca_patchPhi[ipatch] = new double[nfaces]{};

        if (ca_patchRhoUf != nullptr)
            ca_patchRhoUf[ipatch] = new double[nTotal]{};

        // Turbulence data ^^^^^^^^^^^^^^^^^^^^^^^^^^
        if (simulationType == "RAS")
        {
            if (ca_patchAlphaT != nullptr)
                ca_patchAlphaT[ipatch] = new double[nfaces]{};
            if (ca_patchK != nullptr)
                ca_patchK[ipatch] = new double[nfaces]{};
            if (ca_patchEpsilon != nullptr)
                ca_patchEpsilon[ipatch] = new double[nfaces]{};
            if (ca_patchOmega != nullptr)
                ca_patchOmega[ipatch] = new double[nfaces]{};
            if (ca_patchNuT != nullptr)
                ca_patchNuT[ipatch] = new double[nfaces]{};
        }
        //-------------------------------------------

        // RocStar data ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        if (ca_patchNf != nullptr)
            ca_patchNf[ipatch] = new double[nTotal]{};
        if (ca_patchSf != nullptr)
            ca_patchSf[ipatch] = new double[nTotal]{};
        if (ca_patchTrac != nullptr)
            ca_patchTrac[ipatch] = new double[nTotal]{};

        // input quantities to the fluid module
        if (ca_patchDisp != nullptr)
            ca_patchDisp[ipatch] = new double[nTotal_]{};

        if (patchDispOld != nullptr)
            patchDispOld[ipatch] = new double[nTotal_]{};

        // if (patchPointUpdated != nullptr)
        //     patchPointUpdated[ipatch] = new bool[npoints]{};

        if (ca_patchMassFlux != nullptr)
            ca_patchMassFlux[ipatch] = new double[nfaces]{};

        if (ca_patchFlameT != nullptr)
            ca_patchFlameT[ipatch] = new double[nfaces]{};

        if (ca_patchMomentum != nullptr)
            ca_patchMomentum[ipatch] = new double[nTotal]{};
        //-------------------------------------------------
    }

    return 0;
}

int comFoam::updateSurfaceData_outgoing()
{
    const dynamicFvMesh&    mesh(*meshPtr);
    const pointField&       points = mesh.points();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, ipatch)
    {
        int procStartIndex{0};
        for (int iproc=0; iproc<ca_myRank; iproc++)
        {
            procStartIndex += ca_nPatches[iproc];
        }
        int index  = procStartIndex + ipatch;
        int nfaces = ca_patchSize[index];
        if (nfaces == 0)
            continue;
        
        int npoints = *ca_patchPointToPointMap_size[ipatch];
        int localIndex = 0;
        for(int ipoint=0; ipoint<npoints; ipoint++)
        {
            int globalPointID = ca_patchPointToPointMap[ipatch][ipoint];

            for(int jcomp=0; jcomp<nComponents; jcomp++)
            {
                ca_patchPoints[ipatch][localIndex]
                    = points[globalPointID][jcomp];

                localIndex++;
            }
        }
    }

    // Face-centered data ^^^^^^^^^^^^^^^^^^^^^^^
    const volVectorField& U(*UPtr);
    const volScalarField& p(*pPtr);
    const volScalarField& T(*TPtr);
    const volScalarField& rho(*rhoPtr);
    const surfaceScalarField& phi(*phiPtr);
    autoPtr<surfaceVectorField>& rhoUf(rhoUfPtr);
    const surfaceScalarField& magSf = mesh.magSf();
    const surfaceVectorField& Sf = mesh.Sf();

#ifdef HAVE_OFE20
    const compressible::turbulenceModel& turbulence(*turbulencePtr);
#elif defined(HAVE_OF7)
    const compressible::turbulenceModel& turbulence(*turbulencePtr);
#elif defined(HAVE_OF8)
    const compressible::momentumTransportModel& turbulence(*turbulencePtr);
    //const fluidThermophysicalTransportModel& thermoTransModel(*thermophysicalTransportPtr);
#endif

    // Check the formulation bellow:
    // https://www.openfoam.com/documentation/guides/latest/doc/guide-turbulence-ras.html
    volTensorField gradU(fvc::grad(U));
    tmp<volScalarField> muEff(turbulence.muEff());
    volTensorField viscStress(muEff * (gradU + dev2(Foam::T(gradU))) - p*I);

    std::string dynamicSolverType = ca_dynamicSolverType;
    forAll(patches, ipatch)
    {
        double ca_patchVel_time{0};
        double ca_patchRhoUf_time{0};
        double ca_patchNf_time{0};
        double ca_patchSf_time{0};
        double ca_patchT_time{0};
        double ca_patchRho_time{0};
        double ca_patchPhi_time{0};
        double ca_patchTrac_time{0};
        double ca_patchAlphaT_time{0};
        double ca_patchEpsilon_time{0};
        double ca_patchOmega_time{0};
        double ca_patchK_time{0};
        double ca_patchNuT_time{0};

        // Skip zero-face patches ^^^^^^^^^^^^^^
        int procStartIndex{0};
        for (int iproc=0; iproc<ca_myRank; iproc++)
        {
            procStartIndex += ca_nPatches[iproc];
        }
        int index  = procStartIndex + ipatch;
        int nfacesTotal = ca_patchSize[index];
        if (nfacesTotal == 0)
            continue;
        //--------------------------------------

        // Interacting BCs should be set according to the BC name??!!
        if (ca_bcflag != nullptr)
        {
            *ca_bcflag[ipatch] = 2;
            
            if
            (
                *ca_isDynamicFvMesh == 1 &&
                (
                    dynamicSolverType == "displacementLaplacian" ||
                    dynamicSolverType == "solidBodyDisplacementLaplacian"
                )
            )
            {
#if HAVE_OFE20
                const pointVectorField& pointDisplacement = mesh.lookupObject<displacementMotionSolver>
                        (
                            "dynamicMeshDict"
                        ).pointDisplacement();

#elif defined(HAVE_OF7) || defined(HAVE_OF8)

                const motionSolver& motion_ =
                    refCast<const dynamicMotionSolverFvMesh>(mesh).motion();

                const pointVectorField& pointDisplacement =
                        refCast<const displacementMotionSolver>(motion_).pointDisplacement();
#endif

                if (pointDisplacement.boundaryField()[ipatch].type() == movingWallTypeName)
                {
                    *ca_bcflag[ipatch] = 0;
                }
            }
        }

        int ntypes = *ca_patchFaceToPointConn_types[ipatch];
        int faceIndex = 0;
        for(int itype=0; itype<ntypes; itype++)
        {
            int nfaces = ca_patchFaceToPointConn_size[ipatch][itype];
            
            for(int iface=0; iface<nfaces; iface++)
            {
                int localFaceID = ca_patchFaceToFaceMap[ipatch][faceIndex];
                if (std::string(patchTypeStr[ipatch]) == "empty")
                {
                    for(int jcomp=0; jcomp<nComponents; jcomp++)
                    {
                        int localComp = jcomp + faceIndex*nComponents;
                        ca_patchVel[ipatch][localComp] = 0;

                        if (ca_patchRhoUf != nullptr)
                            ca_patchRhoUf[ipatch][localComp] = 0;

                        if (ca_patchNf != nullptr)
                            ca_patchNf[ipatch][localComp] = 0;

                        if (ca_patchSf != nullptr)
                            ca_patchSf[ipatch][localComp] = 0;

                        if (ca_patchTrac != nullptr)
                            ca_patchTrac[ipatch][localComp] = 0;
                    }

                    if (ca_patchP != nullptr)
                        ca_patchP[ipatch][faceIndex] = 0;
                    if (ca_patchT != nullptr)
                        ca_patchT[ipatch][faceIndex] = 0;
                    if (ca_patchRho != nullptr)
                        ca_patchRho[ipatch][faceIndex] = 0;
                    if (ca_patchPhi != nullptr)
                        ca_patchPhi[ipatch][faceIndex] = 0;

                    // Turbulence data ^^^^^^^^^^^^^^^^^^^^
                    if (ca_patchAlphaT != nullptr)
                        ca_patchAlphaT[ipatch][faceIndex] = 0;
                    if (ca_patchK != nullptr)
                        ca_patchK[ipatch][faceIndex] = 0;
                    if (ca_patchEpsilon != nullptr)
                        ca_patchEpsilon[ipatch][faceIndex] = 0;
                    if (ca_patchOmega != nullptr)
                        ca_patchOmega[ipatch][faceIndex] = 0;
                    if (ca_patchNuT != nullptr)
                        ca_patchNuT[ipatch][faceIndex] = 0;
                    //-------------------------------------
                }
                else
                {
                    for(int jcomp=0; jcomp<nComponents; jcomp++)
                    {
                        int localComp = jcomp + faceIndex*nComponents;
                    
                        {
                            double timeIn = MPI_Wtime();

                            ca_patchVel[ipatch][localComp] = 
                                U.boundaryField()[ipatch][localFaceID].component(jcomp);

                            double timeOut = MPI_Wtime();
                            ca_patchVel_time += (timeOut - timeIn);
                        }

                        if (ca_patchRhoUf != nullptr)
                        {
                            double timeIn = MPI_Wtime();

                            ca_patchRhoUf[ipatch][localComp] =
                                rhoUf().boundaryField()[ipatch][localFaceID].component(jcomp);

                            double timeOut = MPI_Wtime();
                            ca_patchRhoUf_time += (timeOut - timeIn);
                        }

                        if (ca_patchNf != nullptr)
                        {
                            double timeIn = MPI_Wtime();

                            ca_patchNf[ipatch][localComp] =
                                Sf.boundaryField()[ipatch][localFaceID].component(jcomp)
                                    /magSf.boundaryField()[ipatch][localFaceID];

                            double timeOut = MPI_Wtime();
                            ca_patchNf_time += (timeOut - timeIn);
                        }
                        
                        if (ca_patchSf != nullptr)
                        {
                            double timeIn = MPI_Wtime();

                            ca_patchSf[ipatch][localComp] =
                                Sf.boundaryField()[ipatch][localFaceID].component(jcomp);

                            double timeOut = MPI_Wtime();
                            ca_patchSf_time += (timeOut - timeIn);
                        }
                    }
                    ca_patchP[ipatch][faceIndex] = p.boundaryField()[ipatch][localFaceID];

                    if (ca_patchT != nullptr)
                    {
                        double timeIn = MPI_Wtime();
                        
                        ca_patchT[ipatch][faceIndex] =
                            T.boundaryField()[ipatch][localFaceID];
                    
                        double timeOut = MPI_Wtime();
                        ca_patchT_time += (timeOut - timeIn);
                    }

                    if (ca_patchRho != nullptr)
                    {
                        double timeIn = MPI_Wtime();
                        
                        ca_patchRho[ipatch][faceIndex] =
                            rho.boundaryField()[ipatch][localFaceID];

                        double timeOut = MPI_Wtime();
                        ca_patchRho_time += (timeOut - timeIn);
                    }

                    if (ca_patchPhi != nullptr)
                    {
                        double timeIn = MPI_Wtime();

                        ca_patchPhi[ipatch][faceIndex] =
                            phi.boundaryField()[ipatch][localFaceID];

                        double timeOut = MPI_Wtime();
                        ca_patchPhi_time += (timeOut - timeIn);
                    }

                    if (ca_patchTrac != nullptr)
                    {
                        double timeIn = MPI_Wtime();

                        vector normal = Sf.boundaryField()[ipatch][localFaceID]
                                        /magSf.boundaryField()[ipatch][localFaceID];

                        vector tractionTmp =
                            viscStress.boundaryField()[ipatch][localFaceID]
                            & normal;

                        for(int jcomp=0; jcomp<nComponents; jcomp++)
                        {
                            int localComp = jcomp + faceIndex*nComponents;
                        
                            ca_patchTrac[ipatch][localComp]
                                = tractionTmp.component(jcomp);
                        }

                        double timeOut = MPI_Wtime();
                        ca_patchTrac_time += (timeOut - timeIn);
                    }

                    // Turblulence data ^^^^^^^^^
                    if (ca_patchAlphaT != nullptr)
                    {
                        double timeIn = MPI_Wtime();

#ifdef HAVE_OFE20
                        const tmp<volScalarField>& alphat = turbulence.alphat();
#elif defined(HAVE_OF7)
                        const tmp<volScalarField>& alphat = turbulence.alphat();
#elif defined(HAVE_OF8)
                        const tmp<volScalarField>& alphat = refCast<const volScalarField>
                            (
                                mesh.objectRegistry::lookupObject<volScalarField>
                                (
                                    "alphat"
                                )
                            );
#endif

                        ca_patchAlphaT[ipatch][faceIndex] =
                            alphat().boundaryField()[ipatch][localFaceID];

                        double timeOut = MPI_Wtime();
                        ca_patchAlphaT_time += (timeOut - timeIn);
                    }

                    if (ca_patchK != nullptr)
                    {
                        double timeIn = MPI_Wtime();
                        
                        const tmp<volScalarField>& k = turbulence.k();
                        ca_patchK[ipatch][faceIndex] =
                            k().boundaryField()[ipatch][localFaceID];

                        double timeOut = MPI_Wtime();
                        ca_patchK_time += (timeOut - timeIn);
                    }

                    if (ca_patchEpsilon != nullptr)
                    {
                        double timeIn = MPI_Wtime();
                        
                        const tmp<volScalarField>& epsilon = turbulence.epsilon();
                        ca_patchEpsilon[ipatch][faceIndex] =
                            epsilon().boundaryField()[ipatch][localFaceID];

                        double timeOut = MPI_Wtime();
                        ca_patchEpsilon_time += (timeOut - timeIn);
                    }

                    if (ca_patchOmega != nullptr)
                    {
                        double timeIn = MPI_Wtime();
                        
                        const tmp<volScalarField>& omega = 
                            mesh.objectRegistry::lookupObject<volScalarField>
                            (
                                "omega"
                            );
                        ca_patchOmega[ipatch][faceIndex] =
                            omega().boundaryField()[ipatch][localFaceID];

                        double timeOut = MPI_Wtime();
                        ca_patchOmega_time += (timeOut - timeIn);
                    }

                    if (ca_patchNuT != nullptr)
                    {
                        double timeIn = MPI_Wtime();
                        
                        const tmp<volScalarField>& nut = turbulence.nut();
                        ca_patchNuT[ipatch][faceIndex] =
                            nut().boundaryField()[ipatch][localFaceID];

                        double timeOut = MPI_Wtime();
                        ca_patchNuT_time += (timeOut - timeIn);
                    }
                    //---------------------------
                }

                faceIndex++;
            }
        }

        /*
        Info << "Elapsed time in updateSurfaceData_outgoing for patch["
             << ipatch << "] = " << ca_patchName[ipatch] << endl;
        Info << "  ca_patchVel_time =     " << ca_patchVel_time << " (s)" << endl;
        Info << "  ca_patchRhoUf_time =   " << ca_patchRhoUf_time << " (s)" << endl;
        Info << "  ca_patchNf_time asdasdsa=      " << ca_patchNf_time << " (s)" << endl;
        Info << "  ca_patchSf_time =      " << ca_patchSf_time << " (s)" << endl;
        Info << "  ca_patchT_time =       " << ca_patchT_time << " (s)" << endl;
        Info << "  ca_patchRho_time =     " << ca_patchRho_time << " (s)" << endl;
        Info << "  ca_patchPhi_time =     " << ca_patchPhi_time << " (s)" << endl;
        Info << "  ca_patchTrac_time =    " << ca_patchTrac_time << " (s)" << endl;
        Info << "  ca_patchAlphaT_time =  " << ca_patchAlphaT_time << " (s)" << endl;
        Info << "  ca_patchK_time =       " << ca_patchK_time << " (s)" << endl;
        Info << "  ca_patchEpsilon_time = " << ca_patchEpsilon_time << " (s)" << endl;
        Info << "  ca_patchOmega_time = " << ca_patchOmega_time << " (s)" << endl;
        Info << "  ca_patchNuT_time =     " << ca_patchNuT_time << " (s)" << endl;
        */
    }

//std::cin.get();
    
    return 0;
}

int comFoam::updateSurfaceData_incoming(const int& count)
{
    //ca_patchDisp: Displacement
    //ca_patchMassFlux: Mass flux (scalar)
    //ca_patchMomentum: Momentum flux (vector)

    if (ca_dynamicSolverType == nullptr)
        return 0;

    if (ca_isDynamicFvMesh == nullptr)
        return 0;

    const dynamicFvMesh& mesh(*meshPtr);
    std::string dynamicSolverType = ca_dynamicSolverType;
    if
    (
        *ca_isDynamicFvMesh == 1 &&
        (
            dynamicSolverType == "displacementLaplacian" ||
            dynamicSolverType == "solidBodyDisplacementLaplacian"
        )
    )
    {
        Info << "updateSurfaceData_incoming." << endl;
        
#if HAVE_OFE20
                const pointVectorField& pointDisplacement = mesh.lookupObject<displacementMotionSolver>
                        (
                            "dynamicMeshDict"
                        ).pointDisplacement();

#elif defined(HAVE_OF7) || defined(HAVE_OF8)

                const motionSolver& motion_ =
                    refCast<const dynamicMotionSolverFvMesh>(mesh).motion();

                const pointVectorField& pointDisplacement =
                        refCast<const displacementMotionSolver>(motion_).pointDisplacement();
#endif
        
        if (pointDisplacementNewPtr == nullptr)
        {
            pointDisplacementNewPtr = new pointVectorField
            (
                IOobject
                (
                    "pointDisplacementNew",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                pointMesh::New(mesh)
            );

            pointVectorField &pointDisplacementNew(*pointDisplacementNewPtr);
            //const pointField&    points = mesh.points();
            forAll(pointDisplacement, ipoint)
            {
                pointDisplacementNew[ipoint] = pointDisplacement[ipoint];
            }
        }
        pointVectorField &pointDisplacementNew(*pointDisplacementNewPtr);

        // This loop added to assure that each patchPoint is updated once;
        // otherwise mutual points might get updated multiple times.
        // Initilaization here, usage next loop
        /*
        const polyBoundaryMesh& patches = mesh.boundaryMesh();
        forAll(patches, ipatch)
        {
            const polyPatch& patch = patches[ipatch];
            const labelList& patchPoints = patch.meshPoints();
            forAll(patchPoints, ipoint)
            {
                patchPointUpdated[ipatch][ipoint] = false;
            }
        }*/

        forAll(pointDisplacement, ipoint)
        {
            pointUpdated[ipoint] = false;
        }

        const polyBoundaryMesh& patches = mesh.boundaryMesh();
        forAll(patches, ipatch)
        {
            int procStartIndex{0};
            for (int iproc=0; iproc<ca_myRank; iproc++)
            {
                procStartIndex += ca_nPatches[iproc];
            }
            int index  = procStartIndex + ipatch;
            int nfacesTotal = ca_patchSize[index];
            if (nfacesTotal == 0)
                continue;

            const polyPatch& patch = patches[ipatch];
            if (pointDisplacement.boundaryField()[ipatch].type() == movingWallTypeName)
            {
                // Loop over all nodes of boundary patch
                const labelList& patchPoints = patch.meshPoints();
                int ca_npoints = *ca_patchPointToPointMap_size[ipatch];
                compareWarningExit(ca_npoints, patchPoints.size(),
                              "ca_npoints", "patchPoints.size()");

                if (ca_npoints<=0 || patchPoints.size()<=0)
                {
                    FatalErrorInFunction
                        << "Error: ca_npoints = 0 "
                        << nl << exit(FatalError);
                }
                
                forAll(patchPoints, ipoint)
                {
                    int globalPointID = ca_patchPointToPointMap[ipatch][ipoint];

                    if (pointUpdated[globalPointID] == true)
                    {
                        /*std::stringstream output{};
                        output << "Patch " << ipatch 
                               << " Point " << ipoint
                               << " already updated.";
                        verbose_message(output.str(), true);*/

                        continue;
                    }

                    //const label& pointID = patch.meshPoints()[ipoint];  // Node index
                    if (count == 1)
                    {
                        for(int jcomp=0; jcomp<nComponents; jcomp++)
                        {
                            int localIndex = jcomp+ipoint*nComponents;
                            patchDispOld[ipatch][localIndex] = 0;
                        }
                    }

                    for(int jcomp=0; jcomp<nComponents; jcomp++)
                    {
                        int localIndex = jcomp+ipoint*nComponents;
                        pointDisplacementNew[globalPointID][jcomp]
                            += ca_patchDisp[ipatch][localIndex] - patchDispOld[ipatch][localIndex];
                    }

                    for(int jcomp=0; jcomp<nComponents; jcomp++)
                    {
                        int localIndex = jcomp+ipoint*nComponents;
                        patchDispOld[ipatch][localIndex] =
                            ca_patchDisp[ipatch][localIndex];
                    }

                    pointUpdated[globalPointID] = true;
                }
            }
        }
    }

    return 0;
}


int comFoam::registerSurfaceData(const char *name)
{
    const dynamicFvMesh& mesh(*meshPtr);
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    int nPatches = patches.size();

    std::string surfName = name+std::string("SURF");
    
    std::stringstream output{};
    output << "rocFoam.registerSurfaceData: "
               << "Registering flow data with name "
               << surfName;
    verbose_message(output.str(), true);

    // Genral patch data ^^^^^^^^^^^^^^^^^^^^^^^^
    std::string dataName = surfName+std::string(".nPatches");
    COM_new_dataitem( dataName, 'w', COM_INT, 1, "");
    COM_set_size( dataName, 0, ca_nProc);
    COM_set_array(dataName, 0, ca_nPatches);
    output = std::stringstream{};
    output << dataName.c_str() << " registered.";
    verbose_message(output.str(), true);

    dataName = surfName+std::string(".maxNameLength");
    COM_new_dataitem( dataName, 'w', COM_INT, 1, "");
    COM_set_size( dataName, 0, ca_nProc);
    COM_set_array(dataName, 0, ca_maxNameLength);
    output = std::stringstream{};
    output << dataName.c_str() << " registered.";
    verbose_message(output.str(), true);

    int patchNameTotalSize{0};
    for (int jproc=0; jproc<ca_nProc; jproc++)
    {
        patchNameTotalSize += ca_maxNameLength[jproc] *
                               ca_nPatches[jproc];
    }
    dataName = surfName+std::string(".patchName");
    COM_new_dataitem( dataName, 'w', COM_CHAR, 1, "");
    COM_set_size( dataName, 0, patchNameTotalSize);
    COM_set_array(dataName, 0, ca_patchName);
    output = std::stringstream{};
    output << dataName.c_str() << " registered.";
    verbose_message(output.str(), true);


    dataName = surfName+std::string(".maxTypeLength");
    COM_new_dataitem( dataName, 'w', COM_INT, 1, "");
    COM_set_size( dataName, 0, ca_nProc);
    COM_set_array(dataName, 0, ca_maxTypeLength);
    output = std::stringstream{};
    output << dataName.c_str() << " registered.";
    verbose_message(output.str(), true);

    int patchTypeTotalSize{0};
    for (int jproc=0; jproc<ca_nProc; jproc++)
    {
        patchTypeTotalSize += ca_maxTypeLength[jproc] *
                               ca_nPatches[jproc];
    }
    dataName = surfName+std::string(".patchType");
    COM_new_dataitem( dataName, 'w', COM_CHAR, 1, "");
    COM_set_size( dataName, 0, patchTypeTotalSize);
    COM_set_array(dataName, 0, ca_patchType);
    output = std::stringstream{};
    output << dataName.c_str() << " registered.";
    verbose_message(output.str(), true);
    
    //dataName = surfName+std::string(".patchInGroup");
    //COM_new_dataitem( dataName, 'p', COM_CHAR, 1, "");

    int nPatchesTotal{0};
    for (int iproc=0; iproc<ca_nProc; iproc++)
    {
        nPatchesTotal += ca_nPatches[iproc];
    }
    dataName = surfName+std::string(".patchStart");
    COM_new_dataitem( dataName, 'w', COM_INT, 1, "");
    COM_set_size( dataName, 0, nPatchesTotal);
    COM_set_array(dataName, 0, ca_patchStart);
    output = std::stringstream{};
    output << dataName.c_str() << " registered.";
    verbose_message(output.str(), true);

    dataName = surfName+std::string(".patchSize");
    COM_new_dataitem( dataName, 'w', COM_INT, 1, "");
    COM_set_size( dataName, 0, nPatchesTotal);
    COM_set_array(dataName, 0, ca_patchSize);
    output = std::stringstream{};
    output << dataName.c_str() << " registered.";
    verbose_message(output.str(), true);

    if (ca_bcflag != nullptr)
    {
        dataName = surfName+std::string(".bcflag");
        COM_new_dataitem( dataName, 'p', COM_INT, 1, "");
    }
    //-------------------------------------------

    // connefctivity-relatd stuff ^^^^^^^^^^^^^^^
    // Point-connectivity
    dataName = surfName+std::string(".patchPointToPointMap_size");
    COM_new_dataitem( dataName, 'p', COM_INT, 1, "");

    dataName = surfName+std::string(".patchPointToPointMap");
    COM_new_dataitem( dataName, 'n', COM_INT, 1, "");

    // face-connectivity
    dataName = surfName+std::string(".patchFaceToPointConn_types");
    COM_new_dataitem( dataName, 'p', COM_INT, 1, "");

    dataName = surfName+std::string(".patchFaceToPointConn_map");
    COM_new_dataitem( dataName, 'p', COM_INT, 1, "");

    dataName = surfName+std::string(".patchFaceToPointConn_size");
    COM_new_dataitem( dataName, 'p', COM_INT, 1, "");

    dataName = surfName+std::string(".patchFaceToFaceMap");
    COM_new_dataitem( dataName, 'e', COM_INT, 1, "");
   
    dataName = surfName+std::string(".patchFaceToFaceMap_inverse");
    COM_new_dataitem( dataName, 'e', COM_INT, 1, "");
    // ------------------------------------------

    // Element data registered with window ^^^^^^
    dataName = surfName+std::string(".vel");
    COM_new_dataitem( dataName, 'e', COM_DOUBLE, nComponents, "m/s");

    dataName = surfName+std::string(".pf"); //(".pres");
    COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "Pa");

    if (ca_patchT != nullptr)
    {
        dataName = surfName+std::string(".temp");
        COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "K");
    }

    if (ca_patchRho != nullptr)
    {
        dataName = surfName+std::string(".rhof_alp"); //(".rho");
        COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "kg/m^3");
    }

    if (ca_patchPhi != nullptr)
    {
        dataName = surfName+std::string(".phi");
        COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "kg/s");
    }

    if (ca_patchRhoUf != nullptr)
    {
        dataName = surfName+std::string(".rhoUf");
        COM_new_dataitem( dataName, 'e', COM_DOUBLE, nComponents, "kg/m^2*s");
    }

    // Turbulence data ^^^^^^^^^^^^^^^^^^^^^^^^^^
    if (ca_patchAlphaT != nullptr)
    {
        dataName = surfName+std::string(".alphaT");
        COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "kg/m/s");
    }

    if (ca_patchK != nullptr)
    {
        dataName = surfName+std::string(".k");
        COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "m^2/s^2");
    }

    if (ca_patchEpsilon != nullptr)
    {
        dataName = surfName+std::string(".epsilon");
        COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "m^2/s^3");
    }

    if (ca_patchOmega != nullptr)
    {
        dataName = surfName+std::string(".omega");
        COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "1/s");
    }

    if (ca_patchNuT != nullptr)
    {
        dataName = surfName+std::string(".nuT");
        COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "m^2/s");
    }
    //-------------------------------------------

    // RocStar data ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    if (ca_patchNf != nullptr)
    {
        dataName = surfName+std::string(".nf_alp"); //(".nf");
        COM_new_dataitem( dataName, 'e', COM_DOUBLE, nComponents, "");
    }

    if (ca_patchSf != nullptr)
    {
        dataName = surfName+std::string(".sf");
        COM_new_dataitem( dataName, 'e', COM_DOUBLE, nComponents, "m^2");
    }

    if (ca_patchTrac != nullptr)
    {
        dataName = surfName+std::string(".tf");
        COM_new_dataitem( dataName, 'e', COM_DOUBLE, nComponents, "Pa");
    }
    
    if (ca_patchDisp != nullptr)
    {
        dataName = surfName+std::string(".du_alp");
        COM_new_dataitem( dataName, 'n', COM_DOUBLE, nComponents, "m");
    }

    if (ca_patchMassFlux != nullptr)
    {
        dataName = surfName+std::string(".mdot_alp");
        COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "kg/(m^2s)");
    }

    if (ca_patchFlameT != nullptr)
    {
        dataName = surfName+std::string(".Tflm_alp");
        COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "K");
    }

    if (ca_patchMomentum != nullptr)
    {
        dataName = surfName+std::string(".rhofvf_alp");
        COM_new_dataitem( dataName, 'e', COM_DOUBLE, nComponents, "kg/(m^2s)");
    }
    // ------------------------------------------

    // paneID>2 reserved for patches
    int paneIDStart = 1;
    for (label iProc=Pstream::master(); iProc<Pstream::myProcNo(); iProc++)
    {
        paneIDStart += ca_nPatches[iProc];
    }
    //int paneIDEnd = paneIDStart+nPatches;
    

    for (int iproc=0; iproc<ca_nProc; iproc++)
    {
        MPI_Barrier(winComm);
        if(iproc!=ca_myRank)
            continue;

        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {
            int paneID = paneIDStart+ipatch;

            output = std::stringstream{};
            output << "procID = " << Pstream::myProcNo()
                   << ", paneID = " << paneID
                   << ", PatchID = " << ipatch << ","
                   << " ^^^^^^^^^^^^^^^";
            verbose_message(output.str(), true);

            int procStartIndex{0};
            for (int iproc_=0; iproc_<ca_myRank; iproc_++)
            {
                procStartIndex += ca_nPatches[iproc_];
            }
            int index  = procStartIndex + ipatch;
            int nfacesTotal = ca_patchSize[index];
            if (nfacesTotal == 0)
                continue;

            if (ca_bcflag != nullptr)
            {
                dataName = surfName+std::string(".bcflag");
                COM_set_size( dataName, paneID, 1);
                COM_set_array(dataName, paneID, ca_bcflag[ipatch]);
                output = std::stringstream{};
                output << dataName.c_str() << " registered.";
                verbose_message(output.str(), true);
            }
            //---------------------------------------
            
            // points
            int nPoints = *ca_patchPointToPointMap_size[ipatch];
            if (nPoints>0)
            {
                dataName = surfName+std::string(".nc");
                COM_set_size( dataName, paneID, nPoints);
                COM_set_array(dataName, paneID, ca_patchPoints[ipatch], nComponents);
                output = std::stringstream{};
                output << dataName.c_str() << " registered.";
                verbose_message(output.str(), true);
            }

            // point-mapping
            dataName = surfName+std::string(".patchPointToPointMap_size");
            COM_set_size( dataName, paneID, 1);
            COM_set_array(dataName, paneID, ca_patchPointToPointMap_size[ipatch]);
            output = std::stringstream{};
            output << dataName.c_str() << " registered.";
            verbose_message(output.str(), true);

            dataName = surfName+std::string(".patchPointToPointMap");
            COM_set_array(dataName, paneID, ca_patchPointToPointMap[ipatch], 1);
            output = std::stringstream{};
            output << dataName.c_str() << " registered.";
            verbose_message(output.str(), true);

            // face-connectivity
            dataName = surfName+std::string(".patchFaceToPointConn_types");
            COM_set_size(     dataName, paneID, 1);
            COM_set_array(    dataName, paneID, ca_patchFaceToPointConn_types[ipatch]);
            output = std::stringstream{};
            output << dataName.c_str() << " registered.";
            verbose_message(output.str(), true);

            int ntypes = *ca_patchFaceToPointConn_types[ipatch];

            dataName = surfName+std::string(".patchFaceToPointConn_map");
            COM_set_size(  dataName, paneID, ntypes);
            COM_set_array( dataName, paneID, ca_patchFaceToPointConn_map[ipatch]);
            output = std::stringstream{};
            output << dataName.c_str() << " registered.";
            verbose_message(output.str(), true);

            dataName = surfName+std::string(".patchFaceToPointConn_size");
            COM_set_size(     dataName, paneID, ntypes);
            COM_set_array(    dataName, paneID, ca_patchFaceToPointConn_size[ipatch]);
            output = std::stringstream{};
            output << dataName.c_str() << " registered.";
            verbose_message(output.str(), true);

            for(int itype=0; itype<ntypes; itype++)
            {
                int typeID = ca_patchFaceToPointConn_map[ipatch][itype];
                int nfaces = ca_patchFaceToPointConn_size[ipatch][itype];

                if (typeID == 3)
                { // Triangle
                    dataName = surfName+std::string(".:t3");
                }
                else if (typeID == 4)
                { // Quad
                    dataName = surfName+std::string(".:q4");
                }
                else
                { // Type not identified
                    FatalErrorInFunction
                        << "=================== ERROR ===================" << endl
                        << " Face typeID " << typeID << " with size = "
                        << nfaces << " not identified!"
                        << nl << exit(FatalError);
                }

                COM_set_size( dataName, paneID, nfaces);
                COM_set_array( dataName,
                               paneID,
                               ca_patchFaceToPointConn[ipatch][itype],
                               typeID
                             );
                output = std::stringstream{};
                output << dataName.c_str() << " registered.";
                verbose_message(output.str(), true);
            }

            dataName = surfName+std::string(".patchFaceToFaceMap");
            COM_set_array(dataName, paneID, ca_patchFaceToFaceMap[ipatch], 1);
            output = std::stringstream{};
            output << dataName.c_str() << " registered.";
            verbose_message(output.str(), true);

            dataName = surfName+std::string(".patchFaceToFaceMap_inverse");
            COM_set_array(dataName, paneID, ca_patchFaceToFaceMap_inverse[ipatch], 1);
            output = std::stringstream{};
            output << dataName.c_str() << " registered.";
            verbose_message(output.str(), true);
            // ------------------------------------------

            // Field variables
            dataName = surfName+std::string(".vel");
            COM_set_array(dataName, paneID, ca_patchVel[ipatch], nComponents);
            output = std::stringstream{};
            output << dataName.c_str() << " registered.";
            verbose_message(output.str(), true);

            dataName = surfName+std::string(".pf"); //(".pres");
            COM_set_array(dataName, paneID, ca_patchP[ipatch], 1);
            output = std::stringstream{};
            output << dataName.c_str() << " registered.";
            verbose_message(output.str(), true);

            if (ca_patchT != nullptr)
            {
                dataName = surfName+std::string(".temp");
                COM_set_array(dataName, paneID, ca_patchT[ipatch], 1);
                output = std::stringstream{};
                output << dataName.c_str() << " registered.";
                verbose_message(output.str(), true);
            }

            if (ca_patchRho != nullptr)
            {
                dataName = surfName+std::string(".rhof_alp"); //(".rho");
                COM_set_array(dataName, paneID, ca_patchRho[ipatch], 1);
                output = std::stringstream{};
                output << dataName.c_str() << " registered.";
                verbose_message(output.str(), true);
            }

            if (ca_patchPhi != nullptr)
            {
                dataName = surfName+std::string(".phi");
                COM_set_array(dataName, paneID, ca_patchPhi[ipatch], 1);
                output = std::stringstream{};
                output << dataName.c_str() << " registered.";
                verbose_message(output.str(), true);
            }

            if (ca_patchRhoUf != nullptr)
            {
                dataName = surfName+std::string(".rhoUf");
                COM_set_array(dataName, paneID, ca_patchRhoUf[ipatch], nComponents);
                output = std::stringstream{};
                output << dataName.c_str() << " registered.";
                verbose_message(output.str(), true);
            }
            // Turbulence data ^^^^^^^^^^^^^^^^^^^^^^
            if (ca_patchAlphaT != nullptr)
            {
                dataName = surfName+std::string(".alphaT");
                COM_set_array(dataName, paneID, ca_patchAlphaT[ipatch], 1);
                output = std::stringstream{};
                output << dataName.c_str() << " registered.";
                verbose_message(output.str(), true);
            }

            if (ca_patchK != nullptr)
            {
                dataName = surfName+std::string(".k");
                COM_set_array(dataName, paneID, ca_patchK[ipatch], 1);
                output = std::stringstream{};
                output << dataName.c_str() << " registered.";
                verbose_message(output.str(), true);
            }

            if (ca_patchEpsilon != nullptr)
            {
                dataName = surfName+std::string(".epsilon");
                COM_set_array(dataName, paneID, ca_patchEpsilon[ipatch], 1);
                output = std::stringstream{};
                output << dataName.c_str() << " registered.";
                verbose_message(output.str(), true);
            }

            if (ca_patchOmega != nullptr)
            {
                dataName = surfName+std::string(".omega");
                COM_set_array(dataName, paneID, ca_patchOmega[ipatch], 1);
                output = std::stringstream{};
                output << dataName.c_str() << " registered.";
                verbose_message(output.str(), true);
            }

            if (ca_patchNuT != nullptr)
            {
                dataName = surfName+std::string(".nuT");
                COM_set_array(dataName, paneID, ca_patchNuT[ipatch], 1);
                output = std::stringstream{};
                output << dataName.c_str() << " registered.";
                verbose_message(output.str(), true);
            }
            //---------------------------------------


            // RocStar data ^^^^^^^^^^^^^^^^^^^^^^^^^
            if (ca_patchNf != nullptr)
            {
                dataName = surfName+std::string(".nf_alp"); //(".nf");
                COM_set_array(dataName, paneID, ca_patchNf[ipatch], nComponents);
                output = std::stringstream{};
                output << dataName.c_str() << " registered.";
                verbose_message(output.str(), true);
            }

            if (ca_patchSf != nullptr)
            {
                dataName = surfName+std::string(".sf");
                COM_set_array(dataName, paneID, ca_patchSf[ipatch], nComponents);
                output = std::stringstream{};
                output << dataName.c_str() << " registered.";
                verbose_message(output.str(), true);
            }

            if (ca_patchTrac != nullptr)
            {
                dataName = surfName+std::string(".tf");
                COM_set_array(dataName, paneID, ca_patchTrac[ipatch], nComponents);
                output = std::stringstream{};
                output << dataName.c_str() << " registered.";
                verbose_message(output.str(), true);
            }

            if (ca_patchDisp != nullptr)
            {
                dataName = surfName+std::string(".du_alp");
                COM_set_array(dataName, paneID, ca_patchDisp[ipatch], nComponents);
                output = std::stringstream{};
                output << dataName.c_str() << " registered.";
                verbose_message(output.str(), true);
            }

            if (ca_patchMassFlux != nullptr)
            {
                dataName = surfName+std::string(".mdot_alp");
                COM_set_array(dataName, paneID, ca_patchMassFlux[ipatch], 1);
                output = std::stringstream{};
                output << dataName.c_str() << " registered.";
                verbose_message(output.str(), true);
            }

            if (ca_patchFlameT != nullptr)
            {
                dataName = surfName+std::string(".Tflm_alp");
                COM_set_array(dataName, paneID, ca_patchFlameT[ipatch], 1);
                output = std::stringstream{};
                output << dataName.c_str() << " registered.";
                verbose_message(output.str(), true);
            }

            if (ca_patchMomentum != nullptr)
            {
                dataName = surfName+std::string(".rhofvf_alp");
                COM_set_array(dataName, paneID, ca_patchMomentum[ipatch], nComponents);
                output = std::stringstream{};
                output << dataName.c_str() << " registered.";
                verbose_message(output.str(), true);
            }
            //---------------------------------------
            output = std::stringstream{};
            output << "----------------------------------------------------"
                    << std::endl;
            verbose_message(output.str(), true);

            /*
            // VTK output: gas-phase grid data ^^^^^^^^^^^^^^^^
            std::string content;
            content  = "# vtk DataFile Version 3.0\n";
            content += "UNSTRUCTURED_GRID example\n";
            content += "ASCII\n";
            content += "DATASET UNSTRUCTURED_GRID\n";
            
            int npoints = ca_patchPointToPointMap_size[ipatch];
            
            content += "POINTS "+std::to_string(npoints)+" float\n";
            
            int localIndex = 0;
            for(int ipoint=0; ipoint<npoints; ipoint++)
            {
                for(int jcomp=0; jcomp<nComponents; jcomp++)
                {
                    content += std::to_string(ca_patchPoints[ipatch][localIndex]);
                                // points[globalPointID][jcomp];
                    if (jcomp<nComponents-1)
                    {
                        content += " ";
                    }
                    else
                    {
                            content += "\n";
                    }
                    localIndex++;
                }
            }
            
            int size{0};
            for (int itype=0; itype<ntypes; itype++)
            {
                int npoints = ca_patchFaceToPointConn_map[ipatch][itype];
                int nfaces = ca_patchFaceToPointConn_size[ipatch][itype];
                for(int iface=0; iface<nfaces; iface++)
                {
                    size++;
                    for(int ipoint=0; ipoint<npoints; ipoint++)
                    {
                        size++;
                    }
                }
            }

            int nfacesTotal = ca_patchSize[ipatch];
            content += "CELLS "+std::to_string(nfacesTotal)
                    +" "+std::to_string(size)+"\n";
            for (int itype=0; itype<ntypes; itype++)
            {
                int npoints = ca_patchFaceToPointConn_map[ipatch][itype];
                int nfaces = ca_patchFaceToPointConn_size[ipatch][itype];
        
                for(int iface=0; iface<nfaces; iface++)
                {

                    content += std::to_string(npoints);

                    for(int ipoint=0; ipoint<npoints; ipoint++)
                    {
                        content +=" ";
                    
                        int index = ipoint+iface*npoints;
                        
                        int ID = ca_patchFaceToPointConn[ipatch][itype][index] - 1;

                        content += std::to_string(ID);
                            //vecFaceToPointConn[iface][ipoint];
                    }
                    content +="\n";
                }
            }

            content += "CELL_TYPES "+std::to_string(nfacesTotal)+"\n";
            for (int itype=0; itype<ntypes; itype++)
            {
                int npoints = ca_patchFaceToPointConn_map[ipatch][itype];
                int nfaces = ca_patchFaceToPointConn_size[ipatch][itype];
                for(int iface=0; iface<nfaces; iface++)
                {

                    if (npoints == 3)
                    {
                        content += "5\n";
                    }
                    else if (npoints == 4)
                    {
                        content += "9\n";
                    }
                    else
                    {
                        content += "XXXX\n";
                    }
                }
            }

            std::ofstream outFile;
            std::string fileName;
            fileName = "SURFACE/patch"+std::to_string(ipatch)+".vtk";
            outFile.open(fileName, std::ios::out);
            if (!outFile.is_open())
            {
                std::cout << "Writing to file " << fileName
                        << " not successfull" << std::endl;
                exit(-1);
            }
            outFile << content;
            outFile.close();
            */
        }
    }
    COM_window_init_done(surfName); 

    return 0;
}

int comFoam::reconstSurfaceData(const char *name)
{
    std::string surfName = name+std::string("SURF");

    std::stringstream output{};
    output << "rocFoam.reconstCaSurfaceData, procID = "
           << ca_myRank
           << ", Retreiving surface data form window "
           << surfName << ".";
    verbose_message(output.str(), true);

    std::string regNames;
    int numDataItems=0;
    
    COM_get_dataitems(surfName.c_str(), &numDataItems, regNames);

    std::vector<std::string> dataItemNames;
    dataItemNames.clear();
    std::istringstream Istr(regNames);
    for (int i=0; i<numDataItems; ++i)
    {
        std::string nameTmp;
        Istr >> nameTmp;
        dataItemNames.push_back(nameTmp);

        output = std::stringstream{};
        output << "  DataItem[" << i << "] = " << nameTmp;
        verbose_message(output.str(), true);
    }
    output = std::stringstream{};
    output << "  Number of items = " << numDataItems << std::endl;
    verbose_message(output.str(), true);

    // Surface data ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    std::string dataName = std::string("nPatches");
    nameExists(dataItemNames, dataName);
    std::string regName = surfName+std::string(".")+dataName;
    int nComp;
    COM_get_array(regName.c_str(), 0, &ca_nPatches);
    COM_get_size(regName.c_str(), 0, &nComp);
    for(int icomp=0; icomp<nComp; icomp++)
    {
        output = std::stringstream{};
        output << "  " << dataName.c_str() << "[" << icomp << "] = "
               << ca_nPatches[icomp];
        verbose_message(output.str(), true);
    }

    // ca_myrank has already been set in YYY::load method
    int nPatches = ca_nPatches[ca_myRank];

    // Primary allocations ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    
    // pathcName ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    dataName = std::string("maxNameLength");
    nameExists(dataItemNames, dataName);
    regName = surfName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), 0, &ca_maxNameLength);
    COM_get_size(regName.c_str(), 0, &nComp);    
    for(int icomp=0; icomp<nComp; icomp++)
    {
        output = std::stringstream{};
        output << "  " << dataName.c_str() << "[" << icomp << "] = "
               << ca_maxNameLength[icomp];
        verbose_message(output.str(), true);
    }

    dataName = std::string("patchName");
    nameExists(dataItemNames, dataName);
    regName = surfName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), 0, &ca_patchName);
    COM_get_size(regName.c_str(), 0, &nComp);

    // Double-check the size of ca_patchName is correct
    int patchNameTotalSize{0};
    for (int jproc=0; jproc<ca_nProc; jproc++)
    {
        patchNameTotalSize += ca_maxNameLength[jproc] *
                               ca_nPatches[jproc];
    }    
    compareWarningExit(patchNameTotalSize, nComp,
                      "patchNameTotalSize", "nComp");
    
    int procStartIndex{0};
    for (int iproc=0; iproc<ca_myRank; iproc++)
    {
        procStartIndex += ca_maxNameLength[iproc] *
                               ca_nPatches[iproc];
    }
    int maxNameLength = ca_maxNameLength[ca_myRank];
    patchNameStr = new std::string[nPatches]{};
    for (int i=0; i<nPatches; i++)
    {
        int startIndex = procStartIndex + i*maxNameLength;
        char* charTmp = &ca_patchName[startIndex];
        patchNameStr[i] = charTmp;
        
        output = std::stringstream{};
        output << "    " << dataName.c_str() << "[" << i << "] = "
               << charTmp << ", " << patchNameStr[i];
        verbose_message(output.str(), true);
    }
    //-----------------------------------------------------

    // pathcType ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    dataName = std::string("maxTypeLength");
    nameExists(dataItemNames, dataName);
    regName = surfName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), 0, &ca_maxTypeLength);
    COM_get_size(regName.c_str(), 0, &nComp);    
    for(int icomp=0; icomp<nComp; icomp++)
    {
        output = std::stringstream{};
        output << "  " << dataName.c_str() << "[" << icomp << "] = "
                  << ca_maxTypeLength[icomp];
        verbose_message(output.str(), true);
    }

    dataName = std::string("patchType");
    nameExists(dataItemNames, dataName);
    regName = surfName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), 0, &ca_patchType);
    COM_get_size(regName.c_str(), 0, &nComp);

    // Double-check if the size of ca_patchType is correct
    int patchTypeTotalSize{0};
    for (int jproc=0; jproc<ca_nProc; jproc++)
    {
        patchTypeTotalSize += ca_maxTypeLength[jproc] *
                               ca_nPatches[jproc];
    }
    compareWarningExit(patchTypeTotalSize, nComp,
                      "patchTypeTotalSize", "nComp");
    
    procStartIndex = 0;
    for (int iproc=0; iproc<ca_myRank; iproc++)
    {
        procStartIndex += ca_maxTypeLength[iproc] *
                               ca_nPatches[iproc];
    }
    int maxTypeLength = ca_maxTypeLength[ca_myRank];
    patchTypeStr = new std::string[nPatches]{};
    for (int i=0; i<nPatches; i++)
    {
        int startIndex = procStartIndex + i*maxTypeLength;
        char* charTmp = &ca_patchType[startIndex];
        patchTypeStr[i] = charTmp;
        
        output = std::stringstream{};
        output << "    " << dataName.c_str() << "[" << i << "] = "
                  << charTmp << ", " << patchTypeStr[i];
        verbose_message(output.str(), true);
    }
    //-----------------------------------------------------

    //ca_patchInGroup = new wordList*[nPatches];


    // pathcStart ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    dataName = std::string("patchStart");
    nameExists(dataItemNames, dataName);
    regName = surfName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), 0, &ca_patchStart);
    COM_get_size(regName.c_str(), 0, &nComp);

    // Double-check if the size of ca_patchType is correct
    int nPatchesTotal{0};
    for (int jproc=0; jproc<ca_nProc; jproc++)
    {
        nPatchesTotal += ca_nPatches[jproc];
    }
    compareWarningExit(nComp, nPatchesTotal,
                      "nComp", "nPatchesTotal");
    procStartIndex = 0;
    for (int iproc=0; iproc<ca_myRank; iproc++)
    {
        procStartIndex += ca_nPatches[iproc];
    }
    for(int i=0; i<nPatches; i++)
    {
        int index = procStartIndex + i;

        output = std::stringstream{};
        output << "  " << dataName.c_str() << "[" << i << "] = "
                  << ca_patchStart[index];
        verbose_message(output.str(), true);
    }
    //-----------------------------------------------------

    // pathcSize ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    dataName = std::string("patchSize");
    nameExists(dataItemNames, dataName);
    regName = surfName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), 0, &ca_patchSize);
    COM_get_size(regName.c_str(), 0, &nComp);
    compareWarningExit(nComp, nPatchesTotal,
                      "nComp", "nPatchesTotal");

    for(int i=0; i<nPatches; i++)
    {
        int index = procStartIndex + i;

        output = std::stringstream{};
        output << "  " << dataName.c_str() << "[" << i << "] = "
                  << ca_patchSize[index];
        verbose_message(output.str(), true);
    }
    //-----------------------------------------------------

    dataName = std::string("bcflag");
    if (nameExists(dataItemNames, dataName))
        ca_bcflag = new int*[nPatches]{};

    ca_patchPointToPointMap_size = new int*[nPatches]{};
    ca_patchPointToPointMap = new int*[nPatches]{};
    ca_patchFaceToFaceMap = new int*[nPatches]{};
    ca_patchFaceToFaceMap_inverse = new int*[nPatches]{};

    ca_patchFaceToPointConn_types = new int*[nPatches]{};
    ca_patchFaceToPointConn_map = new int*[nPatches]{};
    ca_patchFaceToPointConn_size = new int*[nPatches]{};
    ca_patchFaceToPointConn = new int**[nPatches]{};

    ca_patchPoints = new double*[nPatches]{};
    ca_patchVel    = new double*[nPatches]{};
    ca_patchP      = new double*[nPatches]{};

    dataName = std::string("temp");
    if (nameExists(dataItemNames, dataName))
        ca_patchT = new double*[nPatches]{};
    
    dataName = std::string("rhof_alp"); //("rho");
    if (nameExists(dataItemNames, dataName))
        ca_patchRho = new double*[nPatches]{};

    dataName = std::string("phi");
    if (nameExists(dataItemNames, dataName))
        ca_patchPhi = new double*[nPatches]{};

    dataName = std::string("rhoUf");
    if (nameExists(dataItemNames, dataName))
        ca_patchRhoUf = new double*[nPatches]{};

    // Turbulence data ^^^^^^^^^^^^^^^^^^^^^^^^^^
    dataName = std::string("alphaT");
    if (nameExists(dataItemNames, dataName))
        ca_patchAlphaT = new double*[nPatches]{};

    dataName = std::string("k"); //("rho");
    if (nameExists(dataItemNames, dataName))
        ca_patchK = new double*[nPatches]{};

    dataName = std::string("epsilon"); //("rho");
    if (nameExists(dataItemNames, dataName))
        ca_patchEpsilon = new double*[nPatches]{};

    dataName = std::string("omega"); //("rho");
    if (nameExists(dataItemNames, dataName))
        ca_patchOmega = new double*[nPatches]{};

    dataName = std::string("nuT"); //("rho");
    if (nameExists(dataItemNames, dataName))
        ca_patchNuT = new double*[nPatches]{};
    //-------------------------------------------
    
    // RocStar data ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    dataName = std::string("nf_alp"); //("nf");
    if (nameExists(dataItemNames, dataName))
        ca_patchNf = new double*[nPatches]{};

    dataName = std::string("sf");
    if (nameExists(dataItemNames, dataName))
        ca_patchSf = new double*[nPatches]{};

    dataName = std::string("tf");
    if (nameExists(dataItemNames, dataName))
        ca_patchTrac = new double*[nPatches]{};

    dataName = std::string("du_alp");
    if (nameExists(dataItemNames, dataName))
    {
        ca_patchDisp = new double*[nPatches]{};
        patchDispOld = new double*[nPatches]{};
        // patchPointUpdated = new bool*[nPatches]{};
    }

    dataName = std::string("mdot_alp");
    if (nameExists(dataItemNames, dataName))
        ca_patchMassFlux = new double*[nPatches]{};

    dataName = std::string("Tflm_alp");
    if (nameExists(dataItemNames, dataName))
        ca_patchFlameT = new double*[nPatches]{};

    dataName = std::string("rhofvf_alp");
    if (nameExists(dataItemNames, dataName))
        ca_patchMomentum = new double*[nPatches]{};
    //-------------------------------------------
    //-----------------------------------------------------

    //  List of panes in this window ^^^^^^^^^^^^^^^^^^^^^^
    int nPanes;
    int* paneList;
    COM_get_panes(surfName.c_str(), &nPanes, &paneList);
    output = std::stringstream{};
    output << "  Number of Panes = " << nPanes;
    verbose_message(output.str(), true);

    int paneIDStart = 1;
    for (int iProc=0; iProc<ca_myRank; iProc++)
    {
        paneIDStart += ca_nPatches[iProc];
    }
   
    for(int ipatch=0; ipatch<nPatches; ipatch++)
    {
        int paneID = paneIDStart+ipatch;

        output = std::stringstream{};
        output << "  Patch[" << ipatch
             << "], paneID = " << paneID
             << " ^^^^^^^^^^^^^^^^^^^^^^^";
        verbose_message(output.str(), true);

        procStartIndex = 0;
        for (int iproc=0; iproc<ca_myRank; iproc++)
        {
            procStartIndex += ca_nPatches[iproc];
        }
        int index  = procStartIndex + ipatch;
        int nfacesTotal = ca_patchSize[index];
        if (nfacesTotal == 0)
            continue;

        dataName = std::string("bcflag");
        if (nameExists(dataItemNames, dataName))
        {
            regName = surfName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_bcflag[ipatch]);

            output = std::stringstream{};
            output << "    " << dataName.c_str()
                 << " = " << *ca_bcflag[ipatch];
            verbose_message(output.str(), true);
        }

        dataName = std::string("patchFaceToPointConn_types");
        nameExists(dataItemNames, dataName);
        regName = surfName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_patchFaceToPointConn_types[ipatch]);
        output = std::stringstream{};
        output << "    " << dataName.c_str() << " = "
             << *ca_patchFaceToPointConn_types[ipatch];
        verbose_message(output.str(), true);

        dataName = std::string("patchFaceToPointConn_map");
        nameExists(dataItemNames, dataName);
        regName = surfName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_patchFaceToPointConn_map[ipatch]);
        COM_get_size(regName.c_str(), paneID, &nComp);
        for(int icomp=0; icomp<nComp; icomp++)
        {
            output = std::stringstream{};
            output << "    " << dataName.c_str() << "[" << icomp << "] = "
                 << ca_patchFaceToPointConn_map[ipatch][icomp];
            verbose_message(output.str(), true);
        }

        dataName = std::string("patchFaceToPointConn_size");
        nameExists(dataItemNames, dataName);
        regName = surfName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_patchFaceToPointConn_size[ipatch]);
        COM_get_size(regName.c_str(), paneID, &nComp);
        for(int icomp=0; icomp<nComp; icomp++)
        {
            output = std::stringstream{};
            output << "    " << dataName.c_str() << "[" << icomp << "] = "
                 << ca_patchFaceToPointConn_size[ipatch][icomp];
            verbose_message(output.str(), true);
        }

        // Point and connectivity stuff ^^^^^^^^^
        dataName = std::string("nc");
        //nameExists(dataItemNames, dataName);
        regName = surfName+std::string(".")+dataName;
        int nPoints;
            
        COM_get_array(regName.c_str(), paneID, &ca_patchPoints[ipatch], &nComp);
        COM_get_size(regName.c_str(), paneID, &nPoints);
        output = std::stringstream{};
        output << "    " << dataName.c_str() << " points = " << nPoints
             << ", components = " << nComp;
        verbose_message(output.str(), true);

        if (patchDispOld != nullptr)
            patchDispOld[ipatch] = new double[nPoints*nComp]{};

        // if (patchPointUpdated != nullptr)
        //     patchPointUpdated[ipatch] = new bool[nPoints]{};

        int nConn;
        int numElem;
        std::string connNames;
        COM_get_connectivities(surfName.c_str(), paneID, &nConn, connNames);
        std::istringstream connISS(connNames);

        // Secondary allocation ^^^^^^^^^^^^^^^^^
        ca_patchFaceToPointConn[ipatch] = new int*[nConn]{};
        //---------------------------------------

        for (int icon=0; icon<nConn; ++icon)
        {
            std::string connName;
            connISS >> connName;
            //connNames.push_back(connName);

            dataName = surfName+std::string(".")+connName;
            //nameExists(dataItemNames, dataName);

            COM_get_array(dataName.c_str(), paneID, &ca_patchFaceToPointConn[ipatch][icon], &nComp);
            COM_get_size(dataName.c_str(), paneID, &numElem);
            output = std::stringstream{};
            output << "    Connectivity[" << icon << "] = " << connName
                 << ", elements = " << numElem
                 << ", components =" << nComp;
            verbose_message(output.str(), true);
        }

        // Mapping data ^^^^^^^^^^^^^^^^^^^^^^^^^
        dataName = std::string("patchPointToPointMap_size");
        nameExists(dataItemNames, dataName);
        regName = surfName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_patchPointToPointMap_size[ipatch], &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);
        output = std::stringstream{};
        output << "    " << dataName.c_str() << " elements = " << numElem
             << ", components = " << nComp;
        verbose_message(output.str(), true);

        dataName = std::string("patchPointToPointMap");
        nameExists(dataItemNames, dataName);
        regName = surfName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_patchPointToPointMap[ipatch], &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);
        output = std::stringstream{};
        output << "    " << dataName.c_str() << " elements = " << numElem
             << ", components = " << nComp;
        verbose_message(output.str(), true);

        dataName = std::string("patchFaceToFaceMap");
        nameExists(dataItemNames, dataName);
        regName = surfName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_patchFaceToFaceMap[ipatch], &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);
        output = std::stringstream{};
        output << "    " << dataName.c_str() << " elements = " << numElem
             << ", components = " << nComp;
        verbose_message(output.str(), true);

        dataName = std::string("patchFaceToFaceMap_inverse");
        nameExists(dataItemNames, dataName);
        regName = surfName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_patchFaceToFaceMap_inverse[ipatch], &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);
        output = std::stringstream{};
        output << "    " << dataName.c_str() << " elements = " << numElem
             << ", components = " << nComp;
        verbose_message(output.str(), true);


        // Field data ^^^^^^^^^^^^^^^^^^^^^^^^^^^
        dataName = std::string("vel");
        nameExists(dataItemNames, dataName);
        regName = surfName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_patchVel[ipatch], &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);
        output = std::stringstream{};
        output << "    " << dataName.c_str() << " elements = " << numElem
             << ", components = " << nComp;
        verbose_message(output.str(), true);
    
        dataName = std::string("pf"); //("pres");
        nameExists(dataItemNames, dataName);
        regName = surfName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_patchP[ipatch], &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);
        output = std::stringstream{};
        output << "    " << dataName.c_str() << " elements = " << numElem
             << ", components = " << nComp;
        verbose_message(output.str(), true);

        dataName = std::string("temp");
        if (nameExists(dataItemNames, dataName))
        {
            regName = surfName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_patchT[ipatch], &nComp);
            COM_get_size(regName.c_str(), paneID, &numElem);
            output = std::stringstream{};
            output << "    " << dataName.c_str() << " elements = " << numElem
                << ", components = " << nComp;
            verbose_message(output.str(), true);
        }

        dataName = std::string("rhof_alp"); //("rho");
        if (nameExists(dataItemNames, dataName))
        {
            regName = surfName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_patchRho[ipatch], &nComp);
            COM_get_size(regName.c_str(), paneID, &numElem);
            output = std::stringstream{};
            output << "    " << dataName.c_str() << " elements = " << numElem
                << ", components = " << nComp;
            verbose_message(output.str(), true);
        }

        dataName = std::string("phi");
        if (nameExists(dataItemNames, dataName))
        {
            regName = surfName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_patchPhi[ipatch], &nComp);
            COM_get_size(regName.c_str(), paneID, &numElem);
            output = std::stringstream{};
            output << "    " << dataName.c_str() << " elements = " << numElem
                    << ", components = " << nComp;
            verbose_message(output.str(), true);
        }

        dataName = std::string("rhoUf");
        if (nameExists(dataItemNames, dataName))
        {
            regName = surfName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_patchRhoUf[ipatch], &nComp);
            COM_get_size(regName.c_str(), paneID, &numElem);
            output = std::stringstream{};
            output << "    " << dataName.c_str() << " elements = " << numElem
                    << ", components = " << nComp;
            verbose_message(output.str(), true);
        }

        // Turbulence Data ^^^^^^^^^^^^^^^^^^^^^^
        dataName = std::string("alphaT");
        if (nameExists(dataItemNames, dataName))
        {
            regName = surfName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_patchAlphaT[ipatch], &nComp);
            COM_get_size(regName.c_str(), paneID, &numElem);
            output = std::stringstream{};
            output << "    " << dataName.c_str() << " elements = " << numElem
                << ", components = " << nComp;
            verbose_message(output.str(), true);
        }

        dataName = std::string("k");
        if (nameExists(dataItemNames, dataName))
        {
            regName = surfName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_patchK[ipatch], &nComp);
            COM_get_size(regName.c_str(), paneID, &numElem);
            output = std::stringstream{};
            output << "    " << dataName.c_str() << " elements = " << numElem
                << ", components = " << nComp;
            verbose_message(output.str(), true);
        }

        dataName = std::string("epsilon");
        if (nameExists(dataItemNames, dataName))
        {
            regName = surfName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_patchEpsilon[ipatch], &nComp);
            COM_get_size(regName.c_str(), paneID, &numElem);
            output = std::stringstream{};
            output << "    " << dataName.c_str() << " elements = " << numElem
                << ", components = " << nComp;
            verbose_message(output.str(), true);
        }

        dataName = std::string("omega");
        if (nameExists(dataItemNames, dataName))
        {
            regName = surfName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_patchOmega[ipatch], &nComp);
            COM_get_size(regName.c_str(), paneID, &numElem);
            output = std::stringstream{};
            output << "    " << dataName.c_str() << " elements = " << numElem
                << ", components = " << nComp;
            verbose_message(output.str(), true);
        }

        dataName = std::string("nuT");
        if (nameExists(dataItemNames, dataName))
        {
            regName = surfName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_patchNuT[ipatch], &nComp);
            COM_get_size(regName.c_str(), paneID, &numElem);
            output = std::stringstream{};
            output << "    " << dataName.c_str() << " elements = " << numElem
                << ", components = " << nComp;
            verbose_message(output.str(), true);
        }
        //---------------------------------------

        // RocStar Data ^^^^^^^^^^^^^^^^^^^^^^^^^
        dataName = std::string("nf_alp");
        if (nameExists(dataItemNames, dataName))
        {
            regName = surfName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_patchNf[ipatch], &nComp);
            COM_get_size(regName.c_str(), paneID, &numElem);
            output = std::stringstream{};
            output << "    " << dataName.c_str() << " elements = " << numElem
                << ", components = " << nComp;
            verbose_message(output.str(), true);
        }

        dataName = std::string("sf");
        if (nameExists(dataItemNames, dataName))
        {
            regName = surfName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_patchSf[ipatch], &nComp);
            COM_get_size(regName.c_str(), paneID, &numElem);
            output = std::stringstream{};
            output << "    " << dataName.c_str() << " elements = " << numElem
                << ", components = " << nComp;
            verbose_message(output.str(), true);
        }

        dataName = std::string("tf");
        if (nameExists(dataItemNames, dataName))
        {
            regName = surfName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_patchTrac[ipatch], &nComp);
            COM_get_size(regName.c_str(), paneID, &numElem);
            output = std::stringstream{};
            output << "    " << dataName.c_str() << " elements = " << numElem
                << ", components = " << nComp;
            verbose_message(output.str(), true);
        }

        dataName = std::string("du_alp");
        if (nameExists(dataItemNames, dataName))
        {
            regName = surfName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_patchDisp[ipatch], &nComp);
            COM_get_size(regName.c_str(), paneID, &numElem);
            output = std::stringstream{};
            output << "    " << dataName.c_str() << " elements = " << numElem
                << ", components = " << nComp;
            verbose_message(output.str(), true);
        }


        dataName = std::string("mdot_alp");
        if (nameExists(dataItemNames, dataName))
        {
            regName = surfName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_patchMassFlux[ipatch], &nComp);
            COM_get_size(regName.c_str(), paneID, &numElem);
            output = std::stringstream{};
            output << "    " << dataName.c_str() << " elements = " << numElem
                << ", components = " << nComp;
            verbose_message(output.str(), true);
        }

        dataName = std::string("Tflm_alp");
        if (nameExists(dataItemNames, dataName))
        {
            regName = surfName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_patchFlameT[ipatch], &nComp);
            COM_get_size(regName.c_str(), paneID, &numElem);
            output = std::stringstream{};
            output << "    " << dataName.c_str() << " elements = " << numElem
                << ", components = " << nComp;
            verbose_message(output.str(), true);
        }

        dataName = std::string("rhofvf_alp");
        if (nameExists(dataItemNames, dataName))
        {
            regName = surfName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_patchMomentum[ipatch], &nComp);
            COM_get_size(regName.c_str(), paneID, &numElem);
            output = std::stringstream{};
            output << "    " << dataName.c_str() << " elements = " << numElem
                << ", components = " << nComp;
            verbose_message(output.str(), true);
        }
        //---------------------------------------

        output = std::stringstream{};
        output << "  --------------------------------------------------";
        verbose_message(output.str(), true);

        /*
        if (false)
        {
            // VTK output: gas-phase grid data ^^^^^^^^^^^^^^^^
            
            int ntypes = *ca_patchFaceToPointConn_types[ipatch];
            std::string content;
            content  = "# vtk DataFile Version 3.0\n";
            content += "UNSTRUCTURED_GRID example\n";
            content += "ASCII\n";
            content += "DATASET UNSTRUCTURED_GRID\n";
            
            int npoints = ca_patchPointToPointMap_size[ipatch];
            
            content += "POINTS "+std::to_string(npoints)+" float\n";
            
            int localIndex = 0;
            for(int ipoint=0; ipoint<npoints; ipoint++)
            {
                for(int jcomp=0; jcomp<nComponents; jcomp++)
                {
                    content += std::to_string(ca_patchPoints[ipatch][localIndex]);
                              // points[globalPointID][jcomp];
                    if (jcomp<nComponents-1)
                    {
                        content += " ";
                    }
                    else
                    {
                            content += "\n";
                    }
                    localIndex++;
                }
            }
            
            int size{0};
            for (int itype=0; itype<ntypes; itype++)
            {
                int npoints = ca_patchFaceToPointConn_map[ipatch][itype];
                int nfaces = ca_patchFaceToPointConn_size[ipatch][itype];
                for(int iface=0; iface<nfaces; iface++)
                {
                    size++;
                    for(int ipoint=0; ipoint<npoints; ipoint++)
                    {
                        size++;
                    }
                }
            }

            int nfacesTotal = ca_patchSize[ipatch];
            content += "CELLS "+std::to_string(nfacesTotal)
                    +" "+std::to_string(size)+"\n";
            for (int itype=0; itype<ntypes; itype++)
            {
                int npoints = ca_patchFaceToPointConn_map[ipatch][itype];
                int nfaces = ca_patchFaceToPointConn_size[ipatch][itype];
        
                for(int iface=0; iface<nfaces; iface++)
                {

                    content += std::to_string(npoints);

                    for(int ipoint=0; ipoint<npoints; ipoint++)
                    {
                        content +=" ";
                    
                        int index = ipoint+iface*npoints;
                        
                        int ID = ca_patchFaceToPointConn[ipatch][itype][index] - 1;

                        content += std::to_string(ID);
                            //vecFaceToPointConn[iface][ipoint];
                    }
                    content +="\n";
                }
            }

            content += "CELL_TYPES "+std::to_string(nfacesTotal)+"\n";
            for (int itype=0; itype<ntypes; itype++)
            {
                int npoints = ca_patchFaceToPointConn_map[ipatch][itype];
                int nfaces = ca_patchFaceToPointConn_size[ipatch][itype];
                for(int iface=0; iface<nfaces; iface++)
                {

                    if (npoints == 3)
                    {
                        content += "5\n";
                    }
                    else if (npoints == 4)
                    {
                        content += "9\n";
                    }
                    else
                    {
                        content += "XXXX\n";
                    }
                }
            }

            std::ofstream outFile;
            std::string fileName;
            fileName = "SURFACE_REC/patch"+std::to_string(ipatch)+".vtk";
            outFile.open(fileName, std::ios::out);
            if (!outFile.is_open())
            {
                std::cout << "Writing to file " << fileName
                     << " not successfull" << std::endl;
                exit(-1);
            }
            outFile << content;
            outFile.close();
        }
        */
    }


    output = std::stringstream{};
    output << "----------------------------------------------------";
    verbose_message(output.str(), true);

    COM_free_buffer(&paneList);

    return 0;
}

int comFoam::deleteSurfaceData()
{
    const dynamicFvMesh& mesh(*meshPtr);
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    int nPatches = patches.size();

    if (ca_patchInGroup != nullptr)
    {
        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {

            if (ca_patchInGroup[ipatch] != nullptr)
            {
                delete [] ca_patchInGroup[ipatch];
                ca_patchInGroup[ipatch] = nullptr;
            }
        }

        delete [] ca_patchInGroup;
        ca_patchInGroup = nullptr;
    }

    if (ca_patchStart != nullptr)
    {
        delete [] ca_patchStart;
        ca_patchStart = nullptr;
    }

    if (ca_patchSize != nullptr)
    {
        delete [] ca_patchSize;
        ca_patchSize = nullptr;
    }

    if (ca_patchName != nullptr)
    {
        delete [] ca_patchName;
        ca_patchName = nullptr;
    }

    if (ca_maxNameLength != nullptr)
    {
        delete [] ca_maxNameLength;
        ca_maxNameLength = nullptr;
    }

    if (ca_patchType != nullptr)
    {
        delete [] ca_patchType;
        ca_patchType = nullptr;
    }

    if (ca_maxTypeLength != nullptr)
    {
        delete [] ca_maxTypeLength;
        ca_maxTypeLength = nullptr;
    }

    if (patchNameStr != nullptr)
    {
        delete [] patchNameStr;
        patchNameStr = nullptr;
    }
    
    if (patchTypeStr != nullptr)
    {
        delete [] patchTypeStr;
        patchTypeStr = nullptr;
    }

    if (ca_patchFaceToFaceMap != nullptr)
    {
        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {
            if (ca_patchFaceToFaceMap[ipatch] != nullptr)
            {
                delete [] ca_patchFaceToFaceMap[ipatch];
                ca_patchFaceToFaceMap[ipatch] = nullptr;
            }
        }

        delete [] ca_patchFaceToFaceMap;
        ca_patchFaceToFaceMap = nullptr;
    }

    if (ca_patchFaceToFaceMap_inverse != nullptr)
    {
        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {
            if (ca_patchFaceToFaceMap_inverse[ipatch] != nullptr)
            {
                delete [] ca_patchFaceToFaceMap_inverse[ipatch];
                ca_patchFaceToFaceMap_inverse[ipatch] = nullptr;
            }
        }

        delete [] ca_patchFaceToFaceMap_inverse;
        ca_patchFaceToFaceMap_inverse = nullptr;
    }

    if (ca_patchPoints != nullptr)
    {
        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {
            if (ca_patchPoints[ipatch] != nullptr)
            {
                delete [] ca_patchPoints[ipatch];
                ca_patchPoints[ipatch] = nullptr;
            }
        }
        delete [] ca_patchPoints;
        ca_patchPoints = nullptr;
    }

    if (ca_patchVel != nullptr)
    {
        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {
            if (ca_patchVel[ipatch] != nullptr)
            {
                delete [] ca_patchVel[ipatch];
                ca_patchVel[ipatch] = nullptr;
            }
        }
        delete [] ca_patchVel;
        ca_patchVel = nullptr;
    }

    if (ca_patchRho != nullptr)
    {
        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {
            if (ca_patchRho[ipatch] != nullptr)
            {
                delete [] ca_patchRho[ipatch];
                ca_patchRho[ipatch] = nullptr;
            }
        }
        delete [] ca_patchRho;
        ca_patchRho = nullptr;
    }

    if (ca_patchPhi != nullptr)
    {
        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {
            if (ca_patchPhi[ipatch] != nullptr)
            {
                delete [] ca_patchPhi[ipatch];
                ca_patchPhi[ipatch] = nullptr;
            }
        }
        delete [] ca_patchPhi;
        ca_patchPhi = nullptr;
    }

    if (ca_patchRhoUf != nullptr)
    {
        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {
            if (ca_patchRhoUf[ipatch] != nullptr)
            {
                delete [] ca_patchRhoUf[ipatch];
                ca_patchRhoUf[ipatch] = nullptr;
            }
        }
        delete [] ca_patchRhoUf;
        ca_patchRhoUf = nullptr;
    }

    if (ca_patchP != nullptr)
    {
        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {
            if (ca_patchP[ipatch] != nullptr)
            {
                delete [] ca_patchP[ipatch];
                ca_patchP[ipatch] = nullptr;
            }
        }
        delete [] ca_patchP;
        ca_patchP = nullptr;
    }

    if (ca_patchT != nullptr)
    {
        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {
            if (ca_patchT[ipatch] != nullptr)
            {
                delete [] ca_patchT[ipatch];
                ca_patchT[ipatch] = nullptr;
            }
        }
        delete [] ca_patchT;
        ca_patchT = nullptr;
    }

    // Turbulence data ^^^^^^^^^^^^^^^^^^^^^^^^^^
    if (ca_patchAlphaT != nullptr)
    {
        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {
            if (ca_patchAlphaT[ipatch] != nullptr)
            {
                delete [] ca_patchAlphaT[ipatch];
                ca_patchAlphaT[ipatch] = nullptr;
            }
        }
        delete [] ca_patchAlphaT;
        ca_patchAlphaT = nullptr;
    }

    if (ca_patchK != nullptr)
    {
        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {
            if (ca_patchK[ipatch] != nullptr)
            {
                delete [] ca_patchK[ipatch];
                ca_patchK[ipatch] = nullptr;
            }
        }
        delete [] ca_patchK;
        ca_patchK = nullptr;
    }

    if (ca_patchEpsilon != nullptr)
    {
        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {
            if (ca_patchEpsilon[ipatch] != nullptr)
            {
                delete [] ca_patchEpsilon[ipatch];
                ca_patchEpsilon[ipatch] = nullptr;
            }
        }
        delete [] ca_patchEpsilon;
        ca_patchEpsilon = nullptr;
    }

    if (ca_patchOmega != nullptr)
    {
        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {
            if (ca_patchOmega[ipatch] != nullptr)
            {
                delete [] ca_patchOmega[ipatch];
                ca_patchOmega[ipatch] = nullptr;
            }
        }
        delete [] ca_patchOmega;
        ca_patchOmega = nullptr;
    }

    if (ca_patchNuT != nullptr)
    {
        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {
            if (ca_patchNuT[ipatch] != nullptr)
            {
                delete [] ca_patchNuT[ipatch];
                ca_patchNuT[ipatch] = nullptr;
            }
        }
        delete [] ca_patchNuT;
        ca_patchNuT = nullptr;
    }
    //-------------------------------------------

    // RocStar data ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    if (ca_bcflag != nullptr)
    {
        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {
            if (ca_bcflag[ipatch] != nullptr)
            {
                delete ca_bcflag[ipatch];
                ca_bcflag[ipatch] = nullptr;
            }
        }
        delete [] ca_bcflag;
        ca_bcflag = nullptr;
    }

    if (ca_patchNf != nullptr)
    {
        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {
            if (ca_patchNf[ipatch] != nullptr)
            {
                delete [] ca_patchNf[ipatch];
                ca_patchNf[ipatch] = nullptr;
            }
        }
        delete [] ca_patchNf;
        ca_patchNf = nullptr;
    }

    if (ca_patchSf != nullptr)
    {
        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {
            if (ca_patchSf[ipatch] != nullptr)
            {
                delete [] ca_patchSf[ipatch];
                ca_patchSf[ipatch] = nullptr;
            }
        }
        delete [] ca_patchSf;
        ca_patchSf = nullptr;
    }

    if (ca_patchTrac != nullptr)
    {
        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {
            if (ca_patchTrac[ipatch] != nullptr)
            {
                delete [] ca_patchTrac[ipatch];
                ca_patchTrac[ipatch] = nullptr;
            }
        }
        delete [] ca_patchTrac;
        ca_patchTrac = nullptr;
    }

    if (ca_patchDisp != nullptr)
    {
        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {
            if (ca_patchDisp[ipatch] != nullptr)
            {
                delete [] ca_patchDisp[ipatch];
                ca_patchDisp[ipatch] = nullptr;
            }
        }
        delete [] ca_patchDisp;
        ca_patchDisp = nullptr;
    }

    if (patchDispOld != nullptr)
    {
        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {
            if (patchDispOld[ipatch] != nullptr)
            {
                delete [] patchDispOld[ipatch];
                patchDispOld[ipatch] = nullptr;
            }
        }
        delete [] patchDispOld;
        patchDispOld = nullptr;
    }

    // if (patchPointUpdated != nullptr)
    // {
    //     for(int ipatch=0; ipatch<nPatches; ipatch++)
    //     {
    //         if (patchPointUpdated[ipatch] != nullptr)
    //         {
    //             delete [] patchPointUpdated[ipatch];
    //             patchPointUpdated[ipatch] = nullptr;
    //         }
    //     }
    //     delete [] patchPointUpdated;
    //     patchPointUpdated = nullptr;
    // }

    if (ca_patchMassFlux != nullptr)
    {
        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {
            if (ca_patchMassFlux[ipatch] != nullptr)
            {
                delete [] ca_patchMassFlux[ipatch];
                ca_patchMassFlux[ipatch] = nullptr;
            }
        }
        delete [] ca_patchMassFlux;
        ca_patchMassFlux = nullptr;
    }

    if (ca_patchFlameT != nullptr)
    {
        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {
            if (ca_patchFlameT[ipatch] != nullptr)
            {
                delete [] ca_patchFlameT[ipatch];
                ca_patchFlameT[ipatch] = nullptr;
            }
        }
        delete [] ca_patchFlameT;
        ca_patchFlameT = nullptr;
    }

    if (ca_patchMomentum != nullptr)
    {
        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {
            if (ca_patchMomentum[ipatch] != nullptr)
            {
                delete [] ca_patchMomentum[ipatch];
                ca_patchMomentum[ipatch] = nullptr;
            }
        }
        delete [] ca_patchMomentum;
        ca_patchMomentum = nullptr;
    }
    //-------------------------------------------
    
    //  Delete faceToPoint connectivity arrays ^^    
    if (ca_patchFaceToPointConn != nullptr)
    {
        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {
            if (ca_patchFaceToPointConn[ipatch] != nullptr)
            {
                int ntypes = *ca_patchFaceToPointConn_types[ipatch];
                for(int itype=0; itype<ntypes; itype++)
                {
                    if (ca_patchFaceToPointConn[ipatch][itype] != nullptr)
                    {
                        delete [] ca_patchFaceToPointConn[ipatch][itype];
                        ca_patchFaceToPointConn[ipatch][itype] = nullptr;
                    }
                }

                delete [] ca_patchFaceToPointConn[ipatch];
                ca_patchFaceToPointConn[ipatch] = nullptr;
            }
        }
        delete [] ca_patchFaceToPointConn;
        ca_patchFaceToPointConn = nullptr;
    }

    if (ca_patchFaceToPointConn_types != nullptr)
    {
        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {
            if (ca_patchFaceToPointConn_types[ipatch] != nullptr)
            {
                delete [] ca_patchFaceToPointConn_types[ipatch];
                ca_patchFaceToPointConn_types[ipatch] = nullptr;
            }
        }

        delete [] ca_patchFaceToPointConn_types;
        ca_patchFaceToPointConn_types = nullptr;
    }

    if (ca_patchFaceToPointConn_map != nullptr)
    {
        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {
            if (ca_patchFaceToPointConn_map[ipatch] != nullptr)
            {
                delete [] ca_patchFaceToPointConn_map[ipatch];
                ca_patchFaceToPointConn_map[ipatch] = nullptr;
            }
        }
        delete [] ca_patchFaceToPointConn_map;
        ca_patchFaceToPointConn_map = nullptr;
    }

    if (ca_patchFaceToPointConn_size != nullptr)
    {
        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {
            if (ca_patchFaceToPointConn_size[ipatch] != nullptr)
            {
                delete [] ca_patchFaceToPointConn_size[ipatch];
                ca_patchFaceToPointConn_size[ipatch] = nullptr;
            }
        }
        delete [] ca_patchFaceToPointConn_size;
        ca_patchFaceToPointConn_size = nullptr;
    }
    //-------------------------------------------

    // Delete pointToPoint mapping arrays ^^^^^^^
    if (ca_patchPointToPointMap != nullptr)
    {
        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {
            if (ca_patchPointToPointMap[ipatch] != nullptr)
            {
                delete [] ca_patchPointToPointMap[ipatch];
                ca_patchPointToPointMap[ipatch] = nullptr;
            }
        }
        delete [] ca_patchPointToPointMap;
        ca_patchPointToPointMap = nullptr;
    }

    if (ca_patchPointToPointMap_size != nullptr)
    {
        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {
            if (ca_patchPointToPointMap_size[ipatch] != nullptr)
            {
                delete [] ca_patchPointToPointMap_size[ipatch];
                ca_patchPointToPointMap_size[ipatch] = nullptr;
            }
        }

        delete [] ca_patchPointToPointMap_size;
        ca_patchPointToPointMap_size = nullptr;
    }
    //-------------------------------------------

    if (ca_nPatches != nullptr)
    {
        delete [] ca_nPatches;
        ca_nPatches = nullptr;
    }

    return 0;
}

void comFoam::compareWarningExit(
                    const int& val1,
                    const int& val2,
                    const string& name1,
                    const string& name2)
{
    if (val1 != val2)
    {
        FatalErrorInFunction
            << "Error: " << name1
            << " not equal to " << name2
            << val1 << " vs " << val2
            << nl << exit(FatalError);
    }
    
    return;
}

