int comFoam::createSurfaceConnectivities()
{
    const dynamicFvMesh& mesh(*meshPtr);

    // Mesh and conmnectivities ^^^^^^^^^^^^^^^^^
    const faceList& faces = mesh.faces();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    //-------------------------------------------

    // Gather the number of patches in each processor
    int nPatches = patches.size();
    ca_nPatches = new int[Pstream::nProcs()];
    for(int i=0; i<Pstream::nProcs(); i++)
    {
        ca_nPatches[i] = 0;
    }
    
    //send local cell addressing to master process
    if (Pstream::master())
    {
        ca_nPatches[Pstream::myProcNo()] = nPatches;
        for(label jSlave=Pstream::firstSlave(); jSlave<=Pstream::lastSlave(); jSlave++)
        {
	        IPstream fromSlave(Pstream::commsTypes::scheduled, jSlave);
            fromSlave >> ca_nPatches[jSlave];
        }
    }
    else
    {
        OPstream toMaster(Pstream::commsTypes::scheduled, Pstream::masterNo());
        toMaster << nPatches;
    }

    //redistribute cell addressing to slave processes
    if (Pstream::master())
    {
        for(label jSlave=Pstream::firstSlave(); jSlave<=Pstream::lastSlave(); jSlave++)
        {
	        OPstream toSlave(Pstream::commsTypes::scheduled, jSlave);
	        for(int i=0; i<Pstream::nProcs(); i++)
	        {
	            toSlave << ca_nPatches[i];
	        }
        }
    }
    else
    {
        IPstream fromMaster(
                            Pstream::commsTypes::scheduled,
                            Pstream::masterNo());
        for(int i=0; i<Pstream::nProcs(); i++)
        {
            fromMaster >> ca_nPatches[i];
        }
    }
    //-------------------------------------------

    //  Create patch general data  arrays ^^^^^^^
    //ca_nPatches = new int(patches.size());
    ca_patchName    = new char*[nPatches];
    if (patchNameStr == nullptr)
        patchNameStr = new std::string[nPatches];

    ca_patchType    = new char*[nPatches];
    //ca_patchInGroup = new wordList*[nPatches];
    ca_patchStart   = new int*[nPatches];
    ca_patchSize    = new int*[nPatches];

    forAll(patches, ipatch)
    {
        const polyPatch& patch = patches[ipatch];

        const word& patchName = patch.name();
        const word& patchType = patch.type();
        //const wordList& patchInGroup = patch.inGroups();
        const int& patchStart = patch.start();
        const int& patchSize = patch.size();

        std::string tmpStr = patchName;
        ca_patchName[ipatch] = new char [tmpStr.length()+1];
        std::strcpy(ca_patchName[ipatch], tmpStr.c_str());
        patchNameStr[ipatch] = ca_patchName[ipatch];

        tmpStr = patchType;
        ca_patchType[ipatch] = new char [tmpStr.length()+1];
        std::strcpy(ca_patchType[ipatch], tmpStr.c_str());

        ca_patchStart[ipatch] = new int(patchStart);
        ca_patchSize[ipatch]  = new int(patchSize);
    }
    //-------------------------------------------

    // Temporary Vectors ^^^^^^^^^^^^^^^^^^^^^^^^
    std::vector<int> vecTmpInt;

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

        const label& patchStart = patch.start();
        const int& patchSize = patch.size();

        vecTmpInt.clear();
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

        std::map<int, std::vector< std::vector<int>>> mapVecTmpInt;
        for(int iface=0; iface<patchSize; iface++)
        {
            const label& faceID = patchStart + iface;
            const labelList& pointsList = faces[faceID];

            int nPointsInFace = pointsList.size();
            
            vecTmpInt.clear();
            forAll(pointsList, ipoint)
            {
                const int& pointID = pointsList[ipoint];
                
                const auto& pointToPointMap = vecPatchPointToPointMap[ipatch];
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
    ca_patchFaceToPointConn_types = new int*[nPatches];
    ca_patchFaceToPointConn_map   = new int*[nPatches];
    ca_patchFaceToPointConn_size  = new int*[nPatches];
    ca_patchFaceToPointConn = new int**[nPatches];    
    
    ca_patchFaceToFaceMap = new int*[nPatches];
    ca_patchFaceToFaceMap_inverse = new int*[nPatches];
    forAll(patches, ipatch)
    {
        // Create faceToFace mapping arrays ^^^^^^^^^
        int nfacesTotal = *ca_patchSize[ipatch];
        ca_patchFaceToFaceMap[ipatch] = new int[nfacesTotal];
        ca_patchFaceToFaceMap_inverse[ipatch] = new int[nfacesTotal];

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
        ca_patchFaceToPointConn_map[ipatch]   = new int[ntypes];
        ca_patchFaceToPointConn_size[ipatch]  = new int[ntypes];
        ca_patchFaceToPointConn[ipatch] = new int*[ntypes];

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
            ca_patchFaceToPointConn[ipatch][itype] = new int[nTypeConn];
            
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
    ca_patchPointToPointMap_size  = new int*[nPatches];
    ca_patchPointToPointMap = new int*[nPatches];

    forAll(patches, ipatch)
    {
        int npoints = vecPatchPointToPointMap[ipatch].size();
        ca_patchPointToPointMap_size[ipatch] = new int(npoints);
    
        ca_patchPointToPointMap[ipatch] = new int[npoints];
        
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
    ca_patchPoints = new double*[nPatches];
    ca_patchVel    = new double*[nPatches];
    ca_patchP      = new double*[nPatches];
    if (rhoPtr != nullptr)
        ca_patchRho = new double*[nPatches];
    if (TPtr != nullptr)
        ca_patchT   = new double*[nPatches];
    if (phiPtr != nullptr)
        ca_patchPhi = new double*[nPatches];
    
    autoPtr<surfaceVectorField>& rhoUf(rhoUfPtr);
    if (rhoUf.valid())
        ca_patchRhoUf = new double*[nPatches];

    // Turbulence data ^^^^^^^^^^^^^^^^^^^^^^^^^^
    const volVectorField& U(*UPtr);
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
    if (simulationType == "RAS")
    {
        const dictionary& subDict = turbProperties.subDict("RAS");
        word RASModel = subDict.lookup("RASModel");

        if (RASModel == "kEpsilon")
        {
            ca_patchAlphaT = new double*[nPatches];
            ca_patchEpsilon = new double*[nPatches];
            ca_patchK = new double*[nPatches];
            ca_patchNuT = new double*[nPatches];
        }
    }
    //-------------------------------------------

    // RocStar data ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ca_bcflag  = new int*[nPatches];
    ca_patchNf = new double*[nPatches];
    //ca_patchSf = new double*[nPatches];
    ca_patchTrac = new double*[nPatches];
    ca_patchDisp = new double*[nPatches];
    ca_patchMassFlux = new double*[nPatches];
    ca_patchFlameT   = new double*[nPatches];
    ca_patchMomentum = new double*[nPatches];
    //-------------------------------------------
 
    forAll(patches, ipatch)
    {
        // Points
        int npoints = *ca_patchPointToPointMap_size[ipatch];
        int nTotal = npoints * nComponents;
        ca_patchPoints[ipatch] = new double[nTotal];
        
        // Field-data
        int nfaces = *ca_patchSize[ipatch];
        nTotal = nfaces * nComponents;

        ca_patchVel[ipatch] = new double[nTotal]{0};
        ca_patchP[ipatch]   = new double[nfaces]{0};
        if (ca_patchRho != nullptr)
            ca_patchRho[ipatch] = new double[nfaces]{0};

        if (ca_patchT != nullptr)
            ca_patchT[ipatch] = new double[nfaces]{0};

        if (ca_patchPhi != nullptr)
            ca_patchPhi[ipatch] = new double[nfaces]{0};

        if (ca_patchRhoUf != nullptr)
            ca_patchRhoUf[ipatch] = new double[nTotal]{0};

        // Turbulence data ^^^^^^^^^^^^^^^^^^^^^^^^^^
        if (simulationType == "RAS")
        {
            ca_patchAlphaT[ipatch] = new double[nfaces]{0};
            ca_patchEpsilon[ipatch] = new double[nfaces]{0};
            ca_patchK[ipatch] = new double[nfaces]{0};
            ca_patchNuT[ipatch] = new double[nfaces]{0};
        }
        //-------------------------------------------

        // RocStar data ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        if (ca_bcflag != nullptr)
            ca_bcflag[ipatch] = new int(2);
        if (ca_patchNf != nullptr)
            ca_patchNf[ipatch] = new double[nTotal]{0};
        if (ca_patchSf != nullptr)
            ca_patchSf[ipatch] = new double[nTotal]{0};
        if (ca_patchTrac != nullptr)
            ca_patchTrac[ipatch] = new double[nTotal]{0};

        // input quantities to the fluid module
        if (ca_patchDisp != nullptr)
        {
            int nTotal_ = npoints * nComponents;
            ca_patchDisp[ipatch] = new double[nTotal_]{0};
        }
        if (ca_patchMassFlux != nullptr)
            ca_patchMassFlux[ipatch] = new double[nfaces]{0};

        if (ca_patchFlameT != nullptr)
            ca_patchFlameT[ipatch] = new double[nfaces]{0};

        if (ca_patchMomentum != nullptr)
            ca_patchMomentum[ipatch] = new double[nTotal]{0};
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
    const compressible::turbulenceModel& turbulence(*turbulencePtr);

    // Check the formulation bellow:
    // https://www.openfoam.com/documentation/guides/latest/doc/guide-turbulence-ras.html
    volTensorField gradU(fvc::grad(U));
    tmp<volScalarField> muEff(turbulence.muEff());
    volTensorField viscStress(muEff * (gradU + dev2(Foam::T(gradU))) - p*I);

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
        double ca_patchK_time{0};
        double ca_patchNuT_time{0};

        // Interacting BCs should be set according to the BC name??!!
        if (ca_bcflag != nullptr)
        {
            if (patchNameStr[ipatch] == movingWallName)
            {
                *ca_bcflag[ipatch] = 0;
            }
            else
            {
                *ca_bcflag[ipatch] = 2;
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
                if (std::string(ca_patchType[ipatch]) == "empty")
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
                    if (ca_patchEpsilon != nullptr)
                        ca_patchEpsilon[ipatch][faceIndex] = 0;
                    if (ca_patchK != nullptr)
                        ca_patchK[ipatch][faceIndex] = 0;
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

                        const tmp<volScalarField>& alphat = turbulence.alphat();

                        ca_patchAlphaT[ipatch][faceIndex] =
                            alphat().boundaryField()[ipatch][localFaceID];

                        double timeOut = MPI_Wtime();
                        ca_patchAlphaT_time += (timeOut - timeIn);
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

                    if (ca_patchK != nullptr)
                    {
                        double timeIn = MPI_Wtime();
                        
                        const tmp<volScalarField>& k = turbulence.k();
                        ca_patchK[ipatch][faceIndex] =
                            k().boundaryField()[ipatch][localFaceID];

                        double timeOut = MPI_Wtime();
                        ca_patchK_time += (timeOut - timeIn);
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
        Info << "  ca_patchEpsilon_time = " << ca_patchEpsilon_time << " (s)" << endl;
        Info << "  ca_patchK_time =       " << ca_patchK_time << " (s)" << endl;
        Info << "  ca_patchNuT_time =     " << ca_patchNuT_time << " (s)" << endl;
        */
    }


    
    return 0;
}

int comFoam::updateSurfaceData_incoming()
{
    //ca_patchDisp: Displacement
    //ca_patchMassFlux: Mass flux (scalar)
    //ca_patchMomentum: Momentum flux (vector)

    dynamicFvMesh& mesh(*meshPtr);
    const Time& t = mesh.time();
    pointVectorField &pointDisplacementNew(*pointDisplacementNewPtr);

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    int patchFSIid = mesh.boundaryMesh().findPatchID(movingWallName);
    //const fvPatch& patchWallFaces = mesh.boundary()[patchWallID];

    Info << "rocFoam.updateSurfaceData_incoming:"
         << " Assigning pointDisplacement to"
         << " patch[" << patchFSIid << "]=" << movingWallName << " patch."
         << endl;



    /*
    const motionSolver& motion_ =
        refCast<const dynamicMotionSolverFvMesh>(mesh).motion();
    const pointField& pointDisplacement =
            refCast<const displacementMotionSolver>(motion_).pointDisplacement();
    */

ca_patchDisp[patchFSIid][0] = 0.00259899; ca_patchDisp[patchFSIid][1] = 0;   ca_patchDisp[patchFSIid][2] =  -0.00401698;
ca_patchDisp[patchFSIid][3] = 0.00789971; ca_patchDisp[patchFSIid][4] =  0;  ca_patchDisp[patchFSIid][5] =   -0.00816199;
ca_patchDisp[patchFSIid][6] = 0.00789971; ca_patchDisp[patchFSIid][7] =  0;  ca_patchDisp[patchFSIid][8] =   -0.00816199;
ca_patchDisp[patchFSIid][9] = 0.00259899; ca_patchDisp[patchFSIid][10] =  0; ca_patchDisp[patchFSIid][11] =   -0.00401698;
ca_patchDisp[patchFSIid][12] = 0;          ca_patchDisp[patchFSIid][13] =  0; ca_patchDisp[patchFSIid][14] =   0;
ca_patchDisp[patchFSIid][15] = 0;          ca_patchDisp[patchFSIid][16] =  0; ca_patchDisp[patchFSIid][17] =   0;
ca_patchDisp[patchFSIid][18] = 0.00286153; ca_patchDisp[patchFSIid][19] =  0; ca_patchDisp[patchFSIid][20] =   0.00336792;
ca_patchDisp[patchFSIid][21] = 0;          ca_patchDisp[patchFSIid][22] =  0; ca_patchDisp[patchFSIid][23] =   0;
ca_patchDisp[patchFSIid][24] = 0;          ca_patchDisp[patchFSIid][25] =  0; ca_patchDisp[patchFSIid][26] =   0;
ca_patchDisp[patchFSIid][27] = 0.00286153; ca_patchDisp[patchFSIid][28] =  0; ca_patchDisp[patchFSIid][29] =   0.00336792;
ca_patchDisp[patchFSIid][30] = 0.00897453; ca_patchDisp[patchFSIid][31] =  0; ca_patchDisp[patchFSIid][32] =   0.00634141;
ca_patchDisp[patchFSIid][33] = 0.00897453; ca_patchDisp[patchFSIid][34] =  0; ca_patchDisp[patchFSIid][35] =   0.00634141;
ca_patchDisp[patchFSIid][36] = 0.0186637; ca_patchDisp[patchFSIid][37] =  0; ca_patchDisp[patchFSIid][38] =   0.00842845;
ca_patchDisp[patchFSIid][39] = 0.0186637; ca_patchDisp[patchFSIid][40] =  0; ca_patchDisp[patchFSIid][41] =   0.00842845;
ca_patchDisp[patchFSIid][42] = 0.0314006; ca_patchDisp[patchFSIid][43] =  0; ca_patchDisp[patchFSIid][44] =   0.00954133;
ca_patchDisp[patchFSIid][45] = 0.0314006; ca_patchDisp[patchFSIid][46] =  0; ca_patchDisp[patchFSIid][47] =   0.00954133;
ca_patchDisp[patchFSIid][48] = 0.0278702; ca_patchDisp[patchFSIid][49] =  0; ca_patchDisp[patchFSIid][50] =   -0.0167539;
ca_patchDisp[patchFSIid][51] = 0.0418323; ca_patchDisp[patchFSIid][52] =  0; ca_patchDisp[patchFSIid][53] =   -0.0214495;
ca_patchDisp[patchFSIid][54] = 0.0418323; ca_patchDisp[patchFSIid][55] =  0; ca_patchDisp[patchFSIid][56] =   -0.0214495;
ca_patchDisp[patchFSIid][57] = 0.0278702; ca_patchDisp[patchFSIid][58] =  0; ca_patchDisp[patchFSIid][59] =   -0.0167539;
ca_patchDisp[patchFSIid][60] = 0.0164762; ca_patchDisp[patchFSIid][61] =  0; ca_patchDisp[patchFSIid][62] =   -0.0123407;
ca_patchDisp[patchFSIid][63] = 0.0164762; ca_patchDisp[patchFSIid][64] =  0; ca_patchDisp[patchFSIid][65] =   -0.0123407;
ca_patchDisp[patchFSIid][66] = 0.0580029; ca_patchDisp[patchFSIid][67] =  0; ca_patchDisp[patchFSIid][68] =   -0.0263898;
ca_patchDisp[patchFSIid][69] = 0.0580029; ca_patchDisp[patchFSIid][70] =  0; ca_patchDisp[patchFSIid][71] =   -0.0263898;
ca_patchDisp[patchFSIid][72] = 0.0759499; ca_patchDisp[patchFSIid][73] =  0; ca_patchDisp[patchFSIid][74] =   -0.0314988;
ca_patchDisp[patchFSIid][75] = 0.0759499; ca_patchDisp[patchFSIid][76] =  0; ca_patchDisp[patchFSIid][77] =   -0.0314988;
ca_patchDisp[patchFSIid][78] = 0.0952359; ca_patchDisp[patchFSIid][79] =  0; ca_patchDisp[patchFSIid][80] =   -0.0367276;
ca_patchDisp[patchFSIid][81] = 0.0952359; ca_patchDisp[patchFSIid][82] =  0; ca_patchDisp[patchFSIid][83] =   -0.0367276;
ca_patchDisp[patchFSIid][84] = 0.0467775; ca_patchDisp[patchFSIid][85] =  0; ca_patchDisp[patchFSIid][86] =   0.00956404;
ca_patchDisp[patchFSIid][87] = 0.0467775; ca_patchDisp[patchFSIid][88] =  0; ca_patchDisp[patchFSIid][89] =   0.00956404;
ca_patchDisp[patchFSIid][90] = 0.0642932; ca_patchDisp[patchFSIid][91] =  0; ca_patchDisp[patchFSIid][92] =   0.00846738;
ca_patchDisp[patchFSIid][93] = 0.0642932; ca_patchDisp[patchFSIid][94] =  0; ca_patchDisp[patchFSIid][95] =   0.00846738;
ca_patchDisp[patchFSIid][96] = 0.0834096; ca_patchDisp[patchFSIid][97] =  0; ca_patchDisp[patchFSIid][98] =   0.00633716;
ca_patchDisp[patchFSIid][99] = 0.0834096; ca_patchDisp[patchFSIid][100] =  0; ca_patchDisp[patchFSIid][101] =   0.00633716;
ca_patchDisp[patchFSIid][102] = 0.103658; ca_patchDisp[patchFSIid][103] =  0; ca_patchDisp[patchFSIid][104] =   0.00336065;
ca_patchDisp[patchFSIid][105] = 0.103658; ca_patchDisp[patchFSIid][106] =  0; ca_patchDisp[patchFSIid][107] =   0.00336065;
ca_patchDisp[patchFSIid][108] = 0.124728; ca_patchDisp[patchFSIid][109] =  0; ca_patchDisp[patchFSIid][110] =   -0.000252435;
ca_patchDisp[patchFSIid][111] = 0.124728; ca_patchDisp[patchFSIid][112] =  0; ca_patchDisp[patchFSIid][113] =   -0.000252435;
ca_patchDisp[patchFSIid][114] = 0.18063;  ca_patchDisp[patchFSIid][115] =  0; ca_patchDisp[patchFSIid][116] =   -0.0592771;
ca_patchDisp[patchFSIid][117] = 0.203513; ca_patchDisp[patchFSIid][118] =  0; ca_patchDisp[patchFSIid][119] =   -0.0654252;
ca_patchDisp[patchFSIid][120] = 0.203513; ca_patchDisp[patchFSIid][121] =  0; ca_patchDisp[patchFSIid][122] =   -0.0654252;
ca_patchDisp[patchFSIid][123] = 0.18063;  ca_patchDisp[patchFSIid][124] =  0; ca_patchDisp[patchFSIid][125] =   -0.0592771;
ca_patchDisp[patchFSIid][126] = 0.115515; ca_patchDisp[patchFSIid][127] =  0; ca_patchDisp[patchFSIid][128] =   -0.0420828;
ca_patchDisp[patchFSIid][129] = 0.115515; ca_patchDisp[patchFSIid][130] =  0; ca_patchDisp[patchFSIid][131] =   -0.0420828;
ca_patchDisp[patchFSIid][132] = 0.136574; ca_patchDisp[patchFSIid][133] =  0; ca_patchDisp[patchFSIid][134] =   -0.0476056;
ca_patchDisp[patchFSIid][135] = 0.136574; ca_patchDisp[patchFSIid][136] =  0; ca_patchDisp[patchFSIid][137] =   -0.0476056;
ca_patchDisp[patchFSIid][138] = 0.1583;   ca_patchDisp[patchFSIid][139] =  0; ca_patchDisp[patchFSIid][140] =   -0.053333;
ca_patchDisp[patchFSIid][141] = 0.1583;   ca_patchDisp[patchFSIid][142] =  0; ca_patchDisp[patchFSIid][143] =   -0.053333;
ca_patchDisp[patchFSIid][144] = 0.226891; ca_patchDisp[patchFSIid][145] =  0; ca_patchDisp[patchFSIid][146] =   -0.071745;
ca_patchDisp[patchFSIid][147] = 0.226891; ca_patchDisp[patchFSIid][148] =  0; ca_patchDisp[patchFSIid][149] =   -0.071745;
ca_patchDisp[patchFSIid][150] = 0.250692; ca_patchDisp[patchFSIid][151] =  0; ca_patchDisp[patchFSIid][152] =   -0.0781904;
ca_patchDisp[patchFSIid][153] = 0.250692; ca_patchDisp[patchFSIid][154] =  0; ca_patchDisp[patchFSIid][155] =   -0.0781904;
ca_patchDisp[patchFSIid][156] = 0.274827; ca_patchDisp[patchFSIid][157] =  0; ca_patchDisp[patchFSIid][158] =   -0.0847105;
ca_patchDisp[patchFSIid][159] = 0.274827; ca_patchDisp[patchFSIid][160] =  0; ca_patchDisp[patchFSIid][161] =   -0.0847105;
ca_patchDisp[patchFSIid][162] = 0.299201; ca_patchDisp[patchFSIid][163] =  0; ca_patchDisp[patchFSIid][164] =   -0.0912577;
ca_patchDisp[patchFSIid][165] = 0.299201; ca_patchDisp[patchFSIid][166] =  0; ca_patchDisp[patchFSIid][167] =   -0.0912577;
ca_patchDisp[patchFSIid][168] = 0.323721; ca_patchDisp[patchFSIid][169] =  0; ca_patchDisp[patchFSIid][170] =   -0.0977962;
ca_patchDisp[patchFSIid][171] = 0.323721; ca_patchDisp[patchFSIid][172] =  0; ca_patchDisp[patchFSIid][173] =   -0.0977962;
ca_patchDisp[patchFSIid][174] = 0.146466; ca_patchDisp[patchFSIid][175] =  0; ca_patchDisp[patchFSIid][176] =   -0.00434417;
ca_patchDisp[patchFSIid][177] = 0.146466; ca_patchDisp[patchFSIid][178] =  0; ca_patchDisp[patchFSIid][179] =   -0.00434417;
ca_patchDisp[patchFSIid][180] = 0.168804; ca_patchDisp[patchFSIid][181] =  0; ca_patchDisp[patchFSIid][182] =   -0.00882705;
ca_patchDisp[patchFSIid][183] = 0.168804; ca_patchDisp[patchFSIid][184] =  0; ca_patchDisp[patchFSIid][185] =   -0.00882705;
ca_patchDisp[patchFSIid][186] = 0.348307; ca_patchDisp[patchFSIid][187] =  0; ca_patchDisp[patchFSIid][188] =   -0.104308;
ca_patchDisp[patchFSIid][189] = 0.348307; ca_patchDisp[patchFSIid][190] =  0; ca_patchDisp[patchFSIid][191] =   -0.104308;
ca_patchDisp[patchFSIid][192] = 0.372911; ca_patchDisp[patchFSIid][193] =  0; ca_patchDisp[patchFSIid][194] =   -0.110795;
ca_patchDisp[patchFSIid][195] = 0.372911; ca_patchDisp[patchFSIid][196] =  0; ca_patchDisp[patchFSIid][197] =   -0.110795;
ca_patchDisp[patchFSIid][198] = 0.379386; ca_patchDisp[patchFSIid][199] =  0; ca_patchDisp[patchFSIid][200] =   -0.0862009;
ca_patchDisp[patchFSIid][201] = 0.379386; ca_patchDisp[patchFSIid][202] =  0; ca_patchDisp[patchFSIid][203] =   -0.0862009;
ca_patchDisp[patchFSIid][204] = 0.385869; ca_patchDisp[patchFSIid][205] =  0; ca_patchDisp[patchFSIid][206] =   -0.0616025;
ca_patchDisp[patchFSIid][207] = 0.385869; ca_patchDisp[patchFSIid][208] =  0; ca_patchDisp[patchFSIid][209] =   -0.0616025;
ca_patchDisp[patchFSIid][210] = 0.191695; ca_patchDisp[patchFSIid][211] =  0; ca_patchDisp[patchFSIid][212] =   -0.0136619;
ca_patchDisp[patchFSIid][213] = 0.191695; ca_patchDisp[patchFSIid][214] =  0; ca_patchDisp[patchFSIid][215] =   -0.0136619;
ca_patchDisp[patchFSIid][216] = 0.215086; ca_patchDisp[patchFSIid][217] =  0; ca_patchDisp[patchFSIid][218] =   -0.0188303;
ca_patchDisp[patchFSIid][219] = 0.215086; ca_patchDisp[patchFSIid][220] =  0; ca_patchDisp[patchFSIid][221] =   -0.0188303;
ca_patchDisp[patchFSIid][222] = 0.238904; ca_patchDisp[patchFSIid][223] =  0; ca_patchDisp[patchFSIid][224] =   -0.0243164;
ca_patchDisp[patchFSIid][225] = 0.238904; ca_patchDisp[patchFSIid][226] =  0; ca_patchDisp[patchFSIid][227] =   -0.0243164;
ca_patchDisp[patchFSIid][228] = 0.263063; ca_patchDisp[patchFSIid][229] =  0; ca_patchDisp[patchFSIid][230] =   -0.0300938;
ca_patchDisp[patchFSIid][231] = 0.263063; ca_patchDisp[patchFSIid][232] =  0; ca_patchDisp[patchFSIid][233] =   -0.0300938;
ca_patchDisp[patchFSIid][234] = 0.287465; ca_patchDisp[patchFSIid][235] =  0; ca_patchDisp[patchFSIid][236] =   -0.0361215;
ca_patchDisp[patchFSIid][237] = 0.287465; ca_patchDisp[patchFSIid][238] =  0; ca_patchDisp[patchFSIid][239] =   -0.0361215;
ca_patchDisp[patchFSIid][240] = 0.312015; ca_patchDisp[patchFSIid][241] =  0; ca_patchDisp[patchFSIid][242] =   -0.0423461;
ca_patchDisp[patchFSIid][243] = 0.312015; ca_patchDisp[patchFSIid][244] =  0; ca_patchDisp[patchFSIid][245] =   -0.0423461;
ca_patchDisp[patchFSIid][246] = 0.336631; ca_patchDisp[patchFSIid][247] =  0; ca_patchDisp[patchFSIid][248] =   -0.0487061;
ca_patchDisp[patchFSIid][249] = 0.336631; ca_patchDisp[patchFSIid][250] =  0; ca_patchDisp[patchFSIid][251] =   -0.0487061;
ca_patchDisp[patchFSIid][252] = 0.361258; ca_patchDisp[patchFSIid][253] =  0; ca_patchDisp[patchFSIid][254] =   -0.0551419;
ca_patchDisp[patchFSIid][255] = 0.361258; ca_patchDisp[patchFSIid][256] =  0; ca_patchDisp[patchFSIid][257] =   -0.0551419;



    int patchID{-1};
    forAll(patches, ipatch)
    {
        const polyPatch& patch = patches[ipatch];
        int localIndex = 0;
        if (ipatch == patchFSIid)
        {
            // Loop over all nodes of boundary patch
            const labelList& patchPoints = patch.meshPoints();
            int ca_npoints = *ca_patchPointToPointMap_size[ipatch];
            
            if (ca_npoints != patchPoints.size())
            {
                std::cout << "Warning: patchPoints.size() != ca_npoints "
                          << patchPoints.size() << " vs " << ca_npoints
                          << std::endl;
            }
            
            forAll(patchPoints, ipoint)
            {
                int globalPointID = ca_patchPointToPointMap[ipatch][ipoint];
            
                const label& pointID = patch.meshPoints()[ipoint];  // Node index

                for(int jcomp=0; jcomp<nComponents; jcomp++)
                {

                    pointDisplacementNew[pointID][jcomp]
                        = ca_patchDisp[ipatch][localIndex] *(t.value()-8.9)/10;
                
                    localIndex++;
                }

//////////////////////////////////
//////////////////////////////////
//////////////////////////////////
//////////////////////////////////
// FIX THIS


            }
            break;
        }
    }


//Info << pointDisplacementNew << endl;
//Info << "In updateSurfaceData_incoming" << endl;
//std::cin.get();


    /*
    forAll(patches, ipatch)
    {
        if (*ca_bcflag[ipatch] == 2)
            continue;

        int npoints = *ca_patchPointToPointMap_size[ipatch];

        int localIndex = 0;
        for(int ipoint=0; ipoint<npoints; ipoint++)
        {
            int globalPointID = ca_patchPointToPointMap[ipatch][ipoint];

            for(int jcomp=0; jcomp<nComponents; jcomp++)
            {

                newPointDisplacement[globalPointID][jcomp]
                    = ca_patchDisp[ipatch][localIndex];

                localIndex++;
            }
        }
    }
    */



    return 0;
}


int comFoam::registerSurfaceData(const char *name)
{
    const dynamicFvMesh& mesh(*meshPtr);
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    int nPatches = patches.size();

    std::string surfName = name+std::string("SURF");
    
    Foam::Info << endl
               << "rocFoam.registerSurfaceData: "
               << "Registering flow data with name "
               << surfName
               << endl;

    // Genral patch data ^^^^^^^^^^^^^^^^^^^^^^^^
    std::string dataName = surfName+std::string(".nPatches");
    COM_new_dataitem( dataName, 'w', COM_INT, 1, "");
    COM_set_size( dataName, 0, Pstream::nProcs());
    COM_set_array(dataName, 0, ca_nPatches);
    Info << dataName.c_str() << " registered." << endl;

    dataName = surfName+std::string(".patchName");
    COM_new_dataitem( dataName, 'p', COM_CHAR, 1, "");

    dataName = surfName+std::string(".patchType");
    COM_new_dataitem( dataName, 'p', COM_CHAR, 1, "");
    
    //dataName = surfName+std::string(".patchInGroup");
    //COM_new_dataitem( dataName, 'p', COM_CHAR, 1, "");

    dataName = surfName+std::string(".patchStart");
    COM_new_dataitem( dataName, 'p', COM_INT, 1, "");

    dataName = surfName+std::string(".patchSize");
    COM_new_dataitem( dataName, 'p', COM_INT, 1, "");

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

    if (ca_patchEpsilon != nullptr)
    {
        dataName = surfName+std::string(".epsilon");
        COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "m^2/s^3");
    }

    if (ca_patchK != nullptr)
    {
        dataName = surfName+std::string(".k");
        COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "m^2/s^2");
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
    for(label iProc=Pstream::master(); iProc<Pstream::myProcNo(); iProc++)
    {
        paneIDStart += ca_nPatches[iProc];
    }
    //int paneIDEnd = paneIDStart+nPatches;
    
    for(int ipatch=0; ipatch<nPatches; ipatch++)
    {
        int paneID = paneIDStart+ipatch;

        Info << "procID = " << Pstream::myProcNo()
             << ", paneID = " << paneID
             << ", PatchID = " << ipatch << ","
             << " ^^^^^^^^^^^^^^^" << endl;


        // Genral patch data ^^^^^^^^^^^^^^^^^^^^^
        std::string charToStr = std::string(ca_patchName[ipatch]);
        int charSize = charToStr.size()+1;
        dataName = surfName+std::string(".patchName");
        COM_set_size( dataName, paneID, charSize);
        COM_set_array(dataName, paneID, ca_patchName[ipatch]);
        Info << "  " << dataName.c_str() << " registered." << endl;

        charToStr = std::string(ca_patchType[ipatch]);
        charSize = charToStr.size()+1;
        dataName = surfName+std::string(".patchType");
        COM_set_size( dataName, paneID, charSize);
        COM_set_array(dataName, paneID, ca_patchType[ipatch]);
        Info << "  " << dataName.c_str() << " registered." << endl;
        
        /*
        charPtr = const_cast<char*>(ca_patchInGroup[ipatch].c_str());
        charSize = ca_patchInGroup[ipatch].size();
        dataName = surfName+std::string(".patchInGroup");
        COM_set_size( dataName, paneID, charSize);
        COM_set_array(dataName, paneID, charPtr);
        Foam::Info << "   patchInGroup registered." << endl;
        */

        dataName = surfName+std::string(".patchStart");
        COM_set_size( dataName, paneID, 1);
        COM_set_array(dataName, paneID, ca_patchStart[ipatch]);
        Info << "  " << dataName.c_str() << " registered." << endl;

        dataName = surfName+std::string(".patchSize");
        COM_set_size( dataName, paneID, 1);
        COM_set_array(dataName, paneID, ca_patchSize[ipatch]);
        Info << "  " << dataName.c_str() << " registered." << endl;

        if (ca_bcflag != nullptr)
        {
            dataName = surfName+std::string(".bcflag");
            COM_set_size( dataName, paneID, 1);
            COM_set_array(dataName, paneID, ca_bcflag[ipatch]);
            Info << "  " << dataName.c_str() << " registered." << endl;
        }
        //---------------------------------------
        
        // points
        dataName = surfName+std::string(".patchPointToPointMap_size");
        COM_set_size( dataName, paneID, 1);
        COM_set_array(dataName, paneID, ca_patchPointToPointMap_size[ipatch]);
        Info << "  " << dataName.c_str() << " registered." << endl;

        int npoints = *ca_patchPointToPointMap_size[ipatch];
        dataName = surfName+std::string(".nc");
        COM_set_size( dataName, paneID, npoints);
        COM_set_array(dataName, paneID, ca_patchPoints[ipatch], nComponents);
        Info << "  " << dataName.c_str() << " registered." << endl;

        // point-mapping
        dataName = surfName+std::string(".patchPointToPointMap");
        COM_set_array(dataName, paneID, ca_patchPointToPointMap[ipatch], 1);
        Info << "  " << dataName.c_str() << " registered." << endl;

        // face-connectivity
        dataName = surfName+std::string(".patchFaceToPointConn_types");
        COM_set_size(     dataName, paneID, 1);
        COM_set_array(    dataName, paneID, ca_patchFaceToPointConn_types[ipatch]);
        Info << "  " << dataName.c_str() << " registered." << endl;

        int ntypes = *ca_patchFaceToPointConn_types[ipatch];

        dataName = surfName+std::string(".patchFaceToPointConn_map");
        COM_set_size(  dataName, paneID, ntypes);
        COM_set_array( dataName, paneID, ca_patchFaceToPointConn_map[ipatch]);
        Info << "  " << dataName.c_str() << " registered." << endl;

        dataName = surfName+std::string(".patchFaceToPointConn_size");
        COM_set_size(     dataName, paneID, ntypes);
        COM_set_array(    dataName, paneID, ca_patchFaceToPointConn_size[ipatch]);
        Info << "  " << dataName.c_str() << " registered." << endl;

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

                Foam::Info << "=================== WARNING ==================="
                           << " Face typeID " << typeID << " with size = "
                           << nfaces << " not identified!"
                           << endl;
                return -1;
            }

            COM_set_size( dataName, paneID, nfaces);
            COM_set_array( dataName,
                           paneID,
                           ca_patchFaceToPointConn[ipatch][itype],
                           typeID
                         );
            Info << "  " << dataName.c_str() << " registered." << endl;
        }

        dataName = surfName+std::string(".patchFaceToFaceMap");
        COM_set_array(dataName, paneID, ca_patchFaceToFaceMap[ipatch], 1);
        Info << "  " << dataName.c_str() << " registered." << endl;

        dataName = surfName+std::string(".patchFaceToFaceMap_inverse");
        COM_set_array(dataName, paneID, ca_patchFaceToFaceMap_inverse[ipatch], 1);
        Info << "  " << dataName.c_str() << " registered." << endl;
        // ------------------------------------------

        // Field variables
        dataName = surfName+std::string(".vel");
        COM_set_array(dataName, paneID, ca_patchVel[ipatch], nComponents);
        Info << "  " << dataName.c_str() << " registered." << endl;

        dataName = surfName+std::string(".pf"); //(".pres");
        COM_set_array(dataName, paneID, ca_patchP[ipatch], 1);
        Info << "  " << dataName.c_str() << " registered." << endl;

        if (ca_patchT != nullptr)
        {
            dataName = surfName+std::string(".temp");
            COM_set_array(dataName, paneID, ca_patchT[ipatch], 1);
            Info << "  " << dataName.c_str() << " registered." << endl;
        }

        if (ca_patchRho != nullptr)
        {
            dataName = surfName+std::string(".rhof_alp"); //(".rho");
            COM_set_array(dataName, paneID, ca_patchRho[ipatch], 1);
            Info << "  " << dataName.c_str() << " registered." << endl;
        }

        if (ca_patchPhi != nullptr)
        {
            dataName = surfName+std::string(".phi");
            COM_set_array(dataName, paneID, ca_patchPhi[ipatch], 1);
            Info << "  " << dataName.c_str() << " registered." << endl;
        }

        if (ca_patchRhoUf != nullptr)
        {
            dataName = surfName+std::string(".rhoUf");
            COM_set_array(dataName, paneID, ca_patchRhoUf[ipatch], nComponents);
            Info << "  " << dataName.c_str() << " registered." << endl;
        }
        // Turbulence data ^^^^^^^^^^^^^^^^^^^^^^
        if (ca_patchAlphaT != nullptr)
        {
            dataName = surfName+std::string(".alphaT");
            COM_set_array(dataName, paneID, ca_patchAlphaT[ipatch], 1);
            Info << "  " << dataName.c_str() << " registered." << endl;
        }

        if (ca_patchEpsilon != nullptr)
        {
            dataName = surfName+std::string(".epsilon");
            COM_set_array(dataName, paneID, ca_patchEpsilon[ipatch], 1);
            Info << "  " << dataName.c_str() << " registered." << endl;
        }

        if (ca_patchK != nullptr)
        {
            dataName = surfName+std::string(".k");
            COM_set_array(dataName, paneID, ca_patchK[ipatch], 1);
            Info << "  " << dataName.c_str() << " registered." << endl;
        }

        if (ca_patchNuT != nullptr)
        {
            dataName = surfName+std::string(".nuT");
            COM_set_array(dataName, paneID, ca_patchNuT[ipatch], 1);
            Info << "  " << dataName.c_str() << " registered." << endl;
        }
        //---------------------------------------


        // RocStar data ^^^^^^^^^^^^^^^^^^^^^^^^^
        if (ca_patchNf != nullptr)
        {
            dataName = surfName+std::string(".nf_alp"); //(".nf");
            COM_set_array(dataName, paneID, ca_patchNf[ipatch], nComponents);
            Info << "  " << dataName.c_str() << " registered." << endl;
        }

        if (ca_patchSf != nullptr)
        {
            dataName = surfName+std::string(".sf");
            COM_set_array(dataName, paneID, ca_patchSf[ipatch], nComponents);
            Info << "  " << dataName.c_str() << " registered." << endl;
        }

        if (ca_patchTrac != nullptr)
        {
            dataName = surfName+std::string(".tf");
            COM_set_array(dataName, paneID, ca_patchTrac[ipatch], nComponents);
            Info << "  " << dataName.c_str() << " registered." << endl;
        }

        if (ca_patchDisp != nullptr)
        {
            dataName = surfName+std::string(".du_alp");
            COM_set_array(dataName, paneID, ca_patchDisp[ipatch], nComponents);
            Info << "  " << dataName.c_str() << " registered." << endl;
        }

        if (ca_patchMassFlux != nullptr)
        {
            dataName = surfName+std::string(".mdot_alp");
            COM_set_array(dataName, paneID, ca_patchMassFlux[ipatch], 1);
            Info << "  " << dataName.c_str() << " registered." << endl;
        }

        if (ca_patchFlameT != nullptr)
        {
            dataName = surfName+std::string(".Tflm_alp");
            COM_set_array(dataName, paneID, ca_patchFlameT[ipatch], 1);
            Info << "  " << dataName.c_str() << " registered." << endl;
        }

        if (ca_patchMomentum != nullptr)
        {
            dataName = surfName+std::string(".rhofvf_alp");
            COM_set_array(dataName, paneID, ca_patchMomentum[ipatch], nComponents);
            Info << "  " << dataName.c_str() << " registered." << endl;
        }
        //---------------------------------------
        Info << "----------------------------------------------------"
             << endl << endl;

        if (false)
        {
            // VTK output: gas-phase grid data ^^^^^^^^^^^^^^^^
            std::string content;
            content  = "# vtk DataFile Version 3.0\n";
            content += "UNSTRUCTURED_GRID example\n";
            content += "ASCII\n";
            content += "DATASET UNSTRUCTURED_GRID\n";
            
            int npoints = *ca_patchPointToPointMap_size[ipatch];
            
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

            int nfacesTotal = *ca_patchSize[ipatch];
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
                exit(1);
            }
            outFile << content;
            outFile.close();
        }
    }

    COM_window_init_done(surfName); 

    return 0;
}

int comFoam::reconstSurfaceData(const char *name)
{
    std::string surfName = name+std::string("SURF");
    std::cout << "rocFoam.reconstCaSurfaceData, procID = "
              << Pstream::myProcNo()
              << ", Retreiving surface data form window "
              << surfName << "."
              << std::endl;

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
        std::cout << "  DataItem[" << i << "] = " << nameTmp << std::endl;
    }
    std::cout << "  Number of items = " << numDataItems
              << std::endl
              << std::endl;
    
    // Surface data ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    std::string dataName = std::string("nPatches");
    nameExists(dataItemNames, dataName);
    std::string regName = surfName+std::string(".")+dataName;
    int nComp;
    COM_get_array(regName.c_str(), 0, &ca_nPatches);
    COM_get_size(regName.c_str(), 0, &nComp);
    for(int icomp=0; icomp<nComp; icomp++)
    {
        std::cout << "  " << dataName.c_str() << "[" << icomp << "] = "
                  << ca_nPatches[icomp] << std::endl;
    }

    // ca_myrank has already been set in YYY::load method
    int nPatches = ca_nPatches[ca_myRank];

    // Primary allocation ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ca_patchName = new char*[nPatches];
    ca_patchType = new char*[nPatches];

    patchNameStr = new std::string[nPatches];
    patchTypeStr = new std::string[nPatches];

    //ca_patchInGroup = new wordList*[nPatches];
    ca_patchStart = new int*[nPatches];
    ca_patchSize  = new int*[nPatches];

    dataName = std::string("bcflag");
    if (nameExists(dataItemNames, dataName))
        ca_bcflag   = new int*[nPatches];

    ca_patchPointToPointMap_size = new int*[nPatches];
    ca_patchPointToPointMap = new int*[nPatches];
    ca_patchFaceToFaceMap = new int*[nPatches];
    ca_patchFaceToFaceMap_inverse = new int*[nPatches];

    ca_patchFaceToPointConn_types = new int*[nPatches];
    ca_patchFaceToPointConn_map = new int*[nPatches];
    ca_patchFaceToPointConn_size = new int*[nPatches];
    ca_patchFaceToPointConn = new int**[nPatches];

    ca_patchPoints = new double*[nPatches];
    ca_patchVel    = new double*[nPatches];
    ca_patchP      = new double*[nPatches];

    dataName = std::string("temp");
    if (nameExists(dataItemNames, dataName))
        ca_patchT = new double*[nPatches];
    
    dataName = std::string("rhof_alp"); //("rho");
    if (nameExists(dataItemNames, dataName))
        ca_patchRho = new double*[nPatches];

    dataName = std::string("phi");
    if (nameExists(dataItemNames, dataName))
        ca_patchPhi = new double*[nPatches];

    dataName = std::string("rhoUf");
    if (nameExists(dataItemNames, dataName))
        ca_patchRhoUf = new double*[nPatches];

    // Turbulence data ^^^^^^^^^^^^^^^^^^^^^^^^^^
    dataName = std::string("alphaT");
    if (nameExists(dataItemNames, dataName))
        ca_patchAlphaT = new double*[nPatches];

    dataName = std::string("epsilon"); //("rho");
    if (nameExists(dataItemNames, dataName))
        ca_patchEpsilon = new double*[nPatches];

    dataName = std::string("k"); //("rho");
    if (nameExists(dataItemNames, dataName))
        ca_patchK = new double*[nPatches];

    dataName = std::string("nuT"); //("rho");
    if (nameExists(dataItemNames, dataName))
        ca_patchNuT = new double*[nPatches];
    //-------------------------------------------
    
    // RocStar data ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    dataName = std::string("nf_alp"); //("nf");
    if (nameExists(dataItemNames, dataName))
        ca_patchNf = new double*[nPatches];

    dataName = std::string("sf");
    if (nameExists(dataItemNames, dataName))
        ca_patchSf = new double*[nPatches];

    dataName = std::string("tf");
    if (nameExists(dataItemNames, dataName))
        ca_patchTrac = new double*[nPatches];

    dataName = std::string("du_alp");
    if (nameExists(dataItemNames, dataName))
        ca_patchDisp = new double*[nPatches];

    dataName = std::string("mdot_alp");
    if (nameExists(dataItemNames, dataName))
        ca_patchMassFlux = new double*[nPatches];

    dataName = std::string("Tflm_alp");
    if (nameExists(dataItemNames, dataName))
        ca_patchFlameT = new double*[nPatches];

    dataName = std::string("rhofvf_alp");
    if (nameExists(dataItemNames, dataName))
        ca_patchMomentum = new double*[nPatches];
    //-------------------------------------------
    //-----------------------------------------------------

    //  List of panes in this window ^^^^^^^^^^^^^^^^^^^^^^
    int nPanes;
    int* paneList;
    COM_get_panes(surfName.c_str(), &nPanes, &paneList);
    std::cout << "  Number of Panes = "
              << nPanes << std::endl;

    for (int ipane=0; ipane<nPanes; ++ipane)
    {
        int paneID = paneList[ipane];

        std::cout << "  Pane[" << ipane
             << "], paneID = " << paneID
             << " ^^^^^^^^^^^^^^^^^^^^^^^" << std::endl;

        dataName = std::string("patchName");
        nameExists(dataItemNames, dataName);
        regName = surfName+std::string(".")+dataName;
        int nComp;

        COM_get_array(regName.c_str(), paneID, &ca_patchName[ipane]);
        COM_get_size(regName.c_str(), paneID, &nComp);

        patchNameStr[ipane] = ca_patchName[ipane];
        //patchNameStr[ipane].resize(nComp);
        std::cout << "    " << dataName.c_str()
             << " = " << patchNameStr[ipane].c_str() << std::endl;

        dataName = std::string("patchType");
        nameExists(dataItemNames, dataName);
        regName = surfName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_patchType[ipane]);
        COM_get_size(regName.c_str(), paneID, &nComp);

        patchTypeStr[ipane] = ca_patchType[ipane];
        //patchTypeStr[ipane].resize(nComp);
        std::cout << "    " << dataName.c_str()
             << " = " << patchTypeStr[ipane].c_str() << std::endl;

        dataName = std::string("patchStart");
        nameExists(dataItemNames, dataName);
        regName = surfName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_patchStart[ipane]);
        std::cout << "    " << dataName.c_str()
             << " = " << *ca_patchStart[ipane] << std::endl;

        dataName = std::string("patchSize");
        nameExists(dataItemNames, dataName);
        regName = surfName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_patchSize[ipane]);
        std::cout << "    " << dataName.c_str()
             << " = " << *ca_patchSize[ipane] << std::endl;

        dataName = std::string("bcflag");
        if (nameExists(dataItemNames, dataName))
        {
            regName = surfName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_bcflag[ipane]);
            std::cout << "    " << dataName.c_str()
                 << " = " << *ca_bcflag[ipane] << std::endl;
        }

        dataName = std::string("patchPointToPointMap_size");
        nameExists(dataItemNames, dataName);
        regName = surfName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_patchPointToPointMap_size[ipane]);
        std::cout << "    " << dataName.c_str()
             << " = " << *ca_patchPointToPointMap_size[ipane] << std::endl;

        dataName = std::string("patchFaceToPointConn_types");
        nameExists(dataItemNames, dataName);
        regName = surfName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_patchFaceToPointConn_types[ipane]);
        std::cout << "    " << dataName.c_str() << " = "
             << *ca_patchFaceToPointConn_types[ipane] << std::endl;

        dataName = std::string("patchFaceToPointConn_map");
        nameExists(dataItemNames, dataName);
        regName = surfName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_patchFaceToPointConn_map[ipane]);
        COM_get_size(regName.c_str(), paneID, &nComp);
        //std::cout << "    " << dataName.c_str() << " size = " << nComp << std::endl;
        for(int icomp=0; icomp<nComp; icomp++)
        {
            std::cout << "    " << dataName.c_str() << "[" << icomp << "] = "
                 << ca_patchFaceToPointConn_map[ipane][icomp] << std::endl;
        }

        dataName = std::string("patchFaceToPointConn_size");
        nameExists(dataItemNames, dataName);
        regName = surfName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_patchFaceToPointConn_size[ipane]);
        COM_get_size(regName.c_str(), paneID, &nComp);
        //std::cout << "    " << dataName.c_str() << " size = " << nComp << std::endl;
        for(int icomp=0; icomp<nComp; icomp++)
        {
            std::cout << "    " << dataName.c_str() << "[" << icomp << "] = "
                 << ca_patchFaceToPointConn_size[ipane][icomp] << std::endl;
        }

        // Point and connectivity stuff ^^^^^^^^^
        dataName = std::string("nc");
        //nameExists(dataItemNames, dataName);
        regName = surfName+std::string(".")+dataName;
        int nPoints;
            
        COM_get_array(regName.c_str(), paneID, &ca_patchPoints[ipane], &nComp);
        COM_get_size(regName.c_str(), paneID, &nPoints);
        std::cout << "    " << dataName.c_str() << " points = " << nPoints
             << ", components = " << nComp << std::endl;

        int nConn;
        int numElem;
        std::string connNames;
        COM_get_connectivities(surfName.c_str(), paneID, &nConn, connNames);
        std::istringstream connISS(connNames);

        // Secondary allocation ^^^^^^^^^^^^^^^^^
        ca_patchFaceToPointConn[ipane] = new int*[nConn];
        //---------------------------------------

        for (int icon=0; icon<nConn; ++icon)
        {
            std::string connName;
            connISS >> connName;
            //connNames.push_back(connName);

            dataName = surfName+std::string(".")+connName;
            //nameExists(dataItemNames, dataName);

            COM_get_array(dataName.c_str(), paneID, &ca_patchFaceToPointConn[ipane][icon], &nComp);
            COM_get_size(dataName.c_str(), paneID, &numElem);
            std::cout << "    Connectivity[" << icon << "] = " << connName
                 << ", elements = " << numElem
                 << ", components =" << nComp << std::endl;
        }
        
        // Mapping data ^^^^^^^^^^^^^^^^^^^^^^^^^
        dataName = std::string("patchPointToPointMap");
        nameExists(dataItemNames, dataName);
        regName = surfName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_patchPointToPointMap[ipane], &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);
        std::cout << "    " << dataName.c_str() << " elements = " << numElem
             << ", components = " << nComp << std::endl;

        dataName = std::string("patchFaceToFaceMap");
        nameExists(dataItemNames, dataName);
        regName = surfName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_patchFaceToFaceMap[ipane], &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);
        std::cout << "    " << dataName.c_str() << " elements = " << numElem
             << ", components = " << nComp << std::endl;

        dataName = std::string("patchFaceToFaceMap_inverse");
        nameExists(dataItemNames, dataName);
        regName = surfName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_patchFaceToFaceMap_inverse[ipane], &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);
        std::cout << "    " << dataName.c_str() << " elements = " << numElem
             << ", components = " << nComp << std::endl;


        // Field data ^^^^^^^^^^^^^^^^^^^^^^^^^^^
        dataName = std::string("vel");
        nameExists(dataItemNames, dataName);
        regName = surfName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_patchVel[ipane], &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);
        std::cout << "    " << dataName.c_str() << " elements = " << numElem
             << ", components = " << nComp << std::endl;
    
        dataName = std::string("pf"); //("pres");
        nameExists(dataItemNames, dataName);
        regName = surfName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_patchP[ipane], &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);
        std::cout << "    " << dataName.c_str() << " elements = " << numElem
             << ", components = " << nComp << std::endl;

        dataName = std::string("temp");
        if (nameExists(dataItemNames, dataName))
        {
            regName = surfName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_patchT[ipane], &nComp);
            COM_get_size(regName.c_str(), paneID, &numElem);
            std::cout << "    " << dataName.c_str() << " elements = " << numElem
                 << ", components = " << nComp << std::endl;
        }

        dataName = std::string("rhof_alp"); //("rho");
        if (nameExists(dataItemNames, dataName))
        {
            regName = surfName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_patchRho[ipane], &nComp);
            COM_get_size(regName.c_str(), paneID, &numElem);
            std::cout << "    " << dataName.c_str() << " elements = " << numElem
                 << ", components = " << nComp << std::endl;
        }

        dataName = std::string("phi");
        if (nameExists(dataItemNames, dataName))
        {
            regName = surfName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_patchPhi[ipane], &nComp);
            COM_get_size(regName.c_str(), paneID, &numElem);
            std::cout << "    " << dataName.c_str() << " elements = " << numElem
                 << ", components = " << nComp << std::endl;
        }

        dataName = std::string("rhoUf");
        if (nameExists(dataItemNames, dataName))
        {
            regName = surfName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_patchRhoUf[ipane], &nComp);
            COM_get_size(regName.c_str(), paneID, &numElem);
            std::cout << "    " << dataName.c_str() << " elements = " << numElem
                 << ", components = " << nComp << std::endl;
        }

        // Turbulence Data ^^^^^^^^^^^^^^^^^^^^^^
        dataName = std::string("alphaT");
        if (nameExists(dataItemNames, dataName))
        {
            regName = surfName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_patchAlphaT[ipane], &nComp);
            COM_get_size(regName.c_str(), paneID, &numElem);
            std::cout << "    " << dataName.c_str() << " elements = " << numElem
                 << ", components = " << nComp << std::endl;
        }

        dataName = std::string("epsilon");
        if (nameExists(dataItemNames, dataName))
        {
            regName = surfName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_patchEpsilon[ipane], &nComp);
            COM_get_size(regName.c_str(), paneID, &numElem);
            std::cout << "    " << dataName.c_str() << " elements = " << numElem
                 << ", components = " << nComp << std::endl;
        }

        dataName = std::string("k");
        if (nameExists(dataItemNames, dataName))
        {
            regName = surfName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_patchK[ipane], &nComp);
            COM_get_size(regName.c_str(), paneID, &numElem);
            std::cout << "    " << dataName.c_str() << " elements = " << numElem
                 << ", components = " << nComp << std::endl;
        }

        dataName = std::string("nuT");
        if (nameExists(dataItemNames, dataName))
        {
            regName = surfName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_patchNuT[ipane], &nComp);
            COM_get_size(regName.c_str(), paneID, &numElem);
            std::cout << "    " << dataName.c_str() << " elements = " << numElem
                 << ", components = " << nComp << std::endl;
        }
        //---------------------------------------

        // RocStar Data ^^^^^^^^^^^^^^^^^^^^^^^^^
        dataName = std::string("nf_alp");
        if (nameExists(dataItemNames, dataName))
        {
            regName = surfName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_patchNf[ipane], &nComp);
            COM_get_size(regName.c_str(), paneID, &numElem);
            std::cout << "    " << dataName.c_str() << " elements = " << numElem
                 << ", components = " << nComp << std::endl;
        }

        dataName = std::string("sf");
        if (nameExists(dataItemNames, dataName))
        {
            regName = surfName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_patchSf[ipane], &nComp);
            COM_get_size(regName.c_str(), paneID, &numElem);
            std::cout << "    " << dataName.c_str() << " elements = " << numElem
                 << ", components = " << nComp << std::endl;
        }

        dataName = std::string("tf");
        if (nameExists(dataItemNames, dataName))
        {
            regName = surfName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_patchTrac[ipane], &nComp);
            COM_get_size(regName.c_str(), paneID, &numElem);
            std::cout << "    " << dataName.c_str() << " elements = " << numElem
                 << ", components = " << nComp << std::endl;
        }

        dataName = std::string("du_alp");
        if (nameExists(dataItemNames, dataName))
        {
            regName = surfName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_patchDisp[ipane], &nComp);
            COM_get_size(regName.c_str(), paneID, &numElem);
            std::cout << "    " << dataName.c_str() << " elements = " << numElem
                 << ", components = " << nComp << std::endl;
        }


        dataName = std::string("mdot_alp");
        if (nameExists(dataItemNames, dataName))
        {
            regName = surfName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_patchMassFlux[ipane], &nComp);
            COM_get_size(regName.c_str(), paneID, &numElem);
            std::cout << "    " << dataName.c_str() << " elements = " << numElem
                 << ", components = " << nComp << std::endl;
        }

        dataName = std::string("Tflm_alp");
        if (nameExists(dataItemNames, dataName))
        {
            regName = surfName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_patchFlameT[ipane], &nComp);
            COM_get_size(regName.c_str(), paneID, &numElem);
            std::cout << "    " << dataName.c_str() << " elements = " << numElem
                 << ", components = " << nComp << std::endl;
        }

        dataName = std::string("rhofvf_alp");
        if (nameExists(dataItemNames, dataName))
        {
            regName = surfName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_patchMomentum[ipane], &nComp);
            COM_get_size(regName.c_str(), paneID, &numElem);
            std::cout << "    " << dataName.c_str() << " elements = " << numElem
                 << ", components = " << nComp << std::endl;
        }
        //---------------------------------------

        std::cout << "  --------------------------------------------------"
             << std::endl;


        if (false)
        {
            // VTK output: gas-phase grid data ^^^^^^^^^^^^^^^^
            
            int ntypes = *ca_patchFaceToPointConn_types[ipane];
            std::string content;
            content  = "# vtk DataFile Version 3.0\n";
            content += "UNSTRUCTURED_GRID example\n";
            content += "ASCII\n";
            content += "DATASET UNSTRUCTURED_GRID\n";
            
            int npoints = *ca_patchPointToPointMap_size[ipane];
            
            content += "POINTS "+std::to_string(npoints)+" float\n";
            
            int localIndex = 0;
            for(int ipoint=0; ipoint<npoints; ipoint++)
            {
                for(int jcomp=0; jcomp<nComponents; jcomp++)
                {
                    content += std::to_string(ca_patchPoints[ipane][localIndex]);
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
                int npoints = ca_patchFaceToPointConn_map[ipane][itype];
                int nfaces = ca_patchFaceToPointConn_size[ipane][itype];
                for(int iface=0; iface<nfaces; iface++)
                {
                    size++;
                    for(int ipoint=0; ipoint<npoints; ipoint++)
                    {
                        size++;
                    }
                }
            }

            int nfacesTotal = *ca_patchSize[ipane];
            content += "CELLS "+std::to_string(nfacesTotal)
                    +" "+std::to_string(size)+"\n";
            for (int itype=0; itype<ntypes; itype++)
            {
                int npoints = ca_patchFaceToPointConn_map[ipane][itype];
                int nfaces = ca_patchFaceToPointConn_size[ipane][itype];
        
                for(int iface=0; iface<nfaces; iface++)
                {

                    content += std::to_string(npoints);

                    for(int ipoint=0; ipoint<npoints; ipoint++)
                    {
                        content +=" ";
                    
                        int index = ipoint+iface*npoints;
                        
                        int ID = ca_patchFaceToPointConn[ipane][itype][index] - 1;

                        content += std::to_string(ID);
                            //vecFaceToPointConn[iface][ipoint];
                    }
                    content +="\n";
                }
            }

            content += "CELL_TYPES "+std::to_string(nfacesTotal)+"\n";
            for (int itype=0; itype<ntypes; itype++)
            {
                int npoints = ca_patchFaceToPointConn_map[ipane][itype];
                int nfaces = ca_patchFaceToPointConn_size[ipane][itype];
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
            fileName = "SURFACE_REC/patch"+std::to_string(ipane)+".vtk";
            outFile.open(fileName, std::ios::out);
            if (!outFile.is_open())
            {
                std::cout << "Writing to file " << fileName
                     << " not successfull" << std::endl;
                exit(1);
            }
            outFile << content;
            outFile.close();
        }



    }

    std::cout << "----------------------------------------------------"
         << std::endl;

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
        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {

            if (ca_patchStart[ipatch] != nullptr)
            {
                delete ca_patchStart[ipatch];
                ca_patchStart[ipatch] = nullptr;
            }
        }

        delete [] ca_patchStart;
        ca_patchStart = nullptr;
    }


    if (ca_patchSize != nullptr)
    {
        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {

            if (ca_patchSize[ipatch] != nullptr)
            {
                delete ca_patchSize[ipatch];
                ca_patchSize[ipatch] = nullptr;
            }
        }

        delete [] ca_patchSize;
        ca_patchSize = nullptr;
    }


    if (ca_patchName != nullptr)
    {
        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {

            if (&ca_patchName[ipatch] != nullptr)
            {
                delete [] ca_patchName[ipatch];
                ca_patchName[ipatch] = nullptr;
            }
        }

        delete [] ca_patchName;
        ca_patchName = nullptr;
    }

    if (ca_patchType != nullptr)
    {
        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {

            if (&ca_patchType[ipatch] != nullptr)
            {
                delete [] ca_patchType[ipatch];
                ca_patchType[ipatch] = nullptr;
            }
        }

        delete [] ca_patchType;
        ca_patchType = nullptr;
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
                delete ca_patchPointToPointMap_size[ipatch];
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

