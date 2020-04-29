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

        tmpStr = patchType;
        ca_patchType[ipatch] = new char [tmpStr.length()+1];
        std::strcpy(ca_patchType[ipatch], tmpStr.c_str());

        //ca_patchName[ipatch] = const_cast<char*>(patchName.c_str());
        //ca_patchType[ipatch] = const_cast<char*>(patchType.c_str());
        //ca_patchInGroup[ipatch] = new std::string(patchInGroup);
        ca_patchStart[ipatch] = new int(patchStart);
        ca_patchSize[ipatch]  = new int(patchSize);
    }
    //-------------------------------------------

    // Temporary Vectors ^^^^^^^^^^^^^^^^^^^^^^^^
    std::vector<int> vecTmpInt;

    // Patch Connectivity Vectors ^^^^^^^^^^^^^^
    std::vector< std::vector< std::vector<int> >> vecPatchFaceToFaceMap;
    std::vector< std::vector<int> >  vecPatchPointToPointMap;
    std::vector< std::vector< std::vector< std::vector<int> >>>
        vecPatchFaceToPointConn;
    //-------------------------------------------

    // FaceToFace mapping vector ^^^^^^^^^^^^^^^^
    vecPatchFaceToFaceMap.clear();
    vecPatchFaceToFaceMap.resize(nPatches);
    forAll(patches, ipatch)
    {
        vecPatchFaceToFaceMap[ipatch].resize(faceToPointTypeSize);
        for (int itype=0; itype<faceToPointTypeSize; itype++)
        {
            vecPatchFaceToFaceMap[ipatch][itype].clear();
        }
    }

    forAll(patches, ipatch)
    {
        const polyPatch& patch = patches[ipatch];

        const label& patchStart = patch.start();
        const int& patchSize = patch.size();

        for(int iface=0; iface<patchSize; iface++)
        {
            const label& faceID = patchStart + iface;
            const labelList& pointsList = faces[faceID];
                        
            int nPointsInFace = pointsList.size();
            vecPatchFaceToFaceMap[ipatch][nPointsInFace-1]
                .push_back(iface);
        }
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

                std::vector<int>::iterator index = std::find
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
    vecPatchFaceToPointConn.resize(nPatches);
    for (int ipatch=0; ipatch<nPatches; ipatch++)
    {
        vecPatchFaceToPointConn[ipatch].resize(faceToPointTypeSize);
        for (int itype=0; itype<faceToPointTypeSize; itype++)
        {
            vecPatchFaceToPointConn[ipatch][itype].clear();
        }
    }

    forAll(patches, ipatch)
    {
        const polyPatch& patch = patches[ipatch];

        const label& patchStart = patch.start();
        const int& patchSize = patch.size();

        for(int iface=0; iface<patchSize; iface++)
        {
            const label& faceID = patchStart + iface;
            const labelList& pointsList = faces[faceID];

            int nPointsInFace = pointsList.size();
            
            vecTmpInt.clear();
            forAll(pointsList, ipoint)
            {
                const int& pointID = pointsList[ipoint];
                
                std::vector<int>::iterator indexPtr = std::find
                                  (
                                    vecPatchPointToPointMap[ipatch].begin(),
                                    vecPatchPointToPointMap[ipatch].end(),
                                    pointID
                                  );
                if (indexPtr == vecPatchPointToPointMap[ipatch].end())
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
                                    vecPatchPointToPointMap[ipatch].begin(),
                                    indexPtr
                                );

                vecTmpInt.push_back(indexVal);
            }
            vecPatchFaceToPointConn[ipatch][nPointsInFace-1].push_back(vecTmpInt);
        }
    }
    //-------------------------------------------
    
    // Create face mapping and connectivity arrays ^^^
    ca_patchFaceToPointConn_types = new int*[nPatches];
    ca_patchFaceToPointConn_map   = new int*[nPatches];
    ca_patchFaceToPointConn_size  = new int*[nPatches];
    ca_patchFaceToPointConn = new int**[nPatches];    
    
    ca_patchFaceToFaceMap = new int*[nPatches];
    forAll(patches, ipatch)
    {
        // Create faceToFace mapping arrays ^^^^^^^^^
        int nfacesTotal = *ca_patchSize[ipatch];
        ca_patchFaceToFaceMap[ipatch] = new int[nfacesTotal];

        int typeCount = 0;
        int sortedFaceIndex = 0;
        int ntypesTotal = vecPatchFaceToFaceMap[ipatch].size();
        for (int itype=0; itype<ntypesTotal; itype++)
        {
            int nfaces = vecPatchFaceToFaceMap[ipatch][itype].size();
            if (nfaces>0)
            {
                for(int iface=0; iface<nfaces; iface++)
                {
                    ca_patchFaceToFaceMap[ipatch][sortedFaceIndex] =
                        vecPatchFaceToFaceMap[ipatch][itype][iface];
                    
                    sortedFaceIndex++;
                }
                typeCount++;
            }
        }
        //-------------------------------------------

        //  Create faceToPoint connectivity arrays ^^    
        ca_patchFaceToPointConn_types[ipatch] = new int(typeCount);
        ca_patchFaceToPointConn_map[ipatch]   = new int[typeCount];
        ca_patchFaceToPointConn_size[ipatch]  = new int[typeCount];
        ca_patchFaceToPointConn[ipatch] = new int*[typeCount];

        typeCount = 0;
        for (int itype=0; itype<ntypesTotal; itype++)
        {
            int nfaces = vecPatchFaceToPointConn[ipatch][itype].size();
            if (nfaces>0)
            {
                int npoints = itype+1;
                ca_patchFaceToPointConn_map[ipatch][typeCount]  = npoints;
                ca_patchFaceToPointConn_size[ipatch][typeCount] = nfaces;
        
                int nTypeConn = npoints * nfaces;
                ca_patchFaceToPointConn[ipatch][typeCount] = new int[nTypeConn];
                
                for(int iface=0; iface<nfaces; iface++)
                {
                    for(int ipoint=0; ipoint<npoints; ipoint++)
                    {
                        int index = ipoint+iface*npoints;
                        
                        ca_patchFaceToPointConn[ipatch][typeCount][index] =
                            vecPatchFaceToPointConn[ipatch][itype][iface][ipoint];
                    }
                }

                typeCount++;
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
    ca_patchRho    = new double*[nPatches];
    ca_patchP      = new double*[nPatches];
    ca_patchT      = new double*[nPatches];
 
    forAll(patches, ipatch)
    {
        // Points
        int npoints = *ca_patchPointToPointMap_size[ipatch];
        int nTotal = npoints * nComponents;
        ca_patchPoints[ipatch] = new double[nTotal];
        
        // Field-data
        int nfaces = *ca_patchSize[ipatch];
        nTotal = nfaces * nComponents;
        ca_patchVel[ipatch] = new double[nTotal];
        ca_patchRho[ipatch] = new double[nfaces];
        ca_patchP[ipatch]   = new double[nfaces];
        ca_patchT[ipatch]   = new double[nfaces];
    }
    
    return 0;
}

int comFoam::updateSurfaceData()
{
    const dynamicFvMesh& mesh(*meshPtr);
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

    // Cell-centered data ^^^^^^^^^^^^^^^^^^^^^^^
    const volVectorField& U(*UPtr);
    const volScalarField& p(*pPtr);
    const volScalarField& T(*TPtr);
    const volScalarField& rho(*rhoPtr);

    forAll(patches, ipatch)
    {
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
                    }
                    ca_patchP[ipatch][faceIndex] = 0;
                    ca_patchT[ipatch][faceIndex] = 0;
                    ca_patchRho[ipatch][faceIndex] = 0;
                }
                else
                {
                    for(int jcomp=0; jcomp<nComponents; jcomp++)
                    {
                        int localComp = jcomp + faceIndex*nComponents;
                    
                        ca_patchVel[ipatch][localComp] = 
                            U.boundaryField()[ipatch][localFaceID].component(jcomp);
                    }
                    ca_patchP[ipatch][faceIndex] = p.boundaryField()[ipatch][localFaceID];
                    ca_patchT[ipatch][faceIndex] = T.boundaryField()[ipatch][localFaceID];
                    ca_patchRho[ipatch][faceIndex] = rho.boundaryField()[ipatch][localFaceID];
                }
                
                faceIndex++;
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

    Foam::Info << "rocFoam.registerSurfaceData: "
               << "Registering flow data with name "
               << name
               << endl;

    std::string surfName = name+std::string("SURF");

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
    // ------------------------------------------

    // Element data registered with window ^^^^^^
    dataName = surfName+std::string(".vel");
    COM_new_dataitem( dataName, 'e', COM_DOUBLE, nComponents, "m/s");

    dataName = surfName+std::string(".pres");
    COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "Pa");

    dataName = surfName+std::string(".temp");
    COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "K");

    dataName = surfName+std::string(".rho");
    COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "kg/m^3");
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
        // ------------------------------------------

        // Field variables
        dataName = surfName+std::string(".vel");
        COM_set_array(dataName, paneID, ca_patchVel[ipatch], nComponents);
        Info << "  " << dataName.c_str() << " registered." << endl;

        dataName = surfName+std::string(".pres");
        COM_set_array(dataName, paneID, ca_patchP[ipatch], 1);
        Info << "  " << dataName.c_str() << " registered." << endl;

        dataName = surfName+std::string(".temp");
        COM_set_array(dataName, paneID, ca_patchT[ipatch], 1);
        Info << "  " << dataName.c_str() << " registered." << endl;

        dataName = surfName+std::string(".rho");
        COM_set_array(dataName, paneID, ca_patchRho[ipatch], 1);
        Info << "  " << dataName.c_str() << " registered." << endl;

        Info << "----------------------------------------------------"
             << endl << endl;
    }

    COM_window_init_done(surfName); 

    return 0;
}


int comFoam::reconstCaSurfaceData(const char *name)
{
    std::string volName = name+std::string("SURF");
    std::cout << "rocFoam.reconstCaSurfaceData, proID = "
              << Pstream::myProcNo()
              << ", Retreiving surface data form window "
              << volName << "."
              << std::endl;

    std::string regNames;
    int numDataItems=0;
    
    COM_get_dataitems(volName.c_str(), &numDataItems, regNames);
    std::cout << "  numDataItems = " << numDataItems << std::endl;

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
    std::cout << std::endl;

    
    // Surface data ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    std::string dataName = std::string("nPatches");
    nameExists(dataItemNames, dataName);
    std::string regName = volName+std::string(".")+dataName;
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

    // Primary allocation ^^^^^^^^^^^^^^^^^^^^^^^
    ca_patchName    = new char*[nPatches];
    ca_patchType    = new char*[nPatches];

    patchNameStr    = new std::string[nPatches];
    patchTypeStr    = new std::string[nPatches];

    //ca_patchInGroup = new wordList*[nPatches];
    ca_patchStart   = new int*[nPatches];
    ca_patchSize    = new int*[nPatches];

    ca_patchPointToPointMap_size = new int*[nPatches];
    ca_patchPointToPointMap = new int*[nPatches];
    ca_patchFaceToFaceMap = new int*[nPatches];

    ca_patchFaceToPointConn_types = new int*[nPatches];
    ca_patchFaceToPointConn_map = new int*[nPatches];
    ca_patchFaceToPointConn_size = new int*[nPatches];
    ca_patchFaceToPointConn = new int**[nPatches];

    ca_patchPoints = new double*[nPatches];
    ca_patchVel = new double*[nPatches];
    ca_patchP = new double*[nPatches];
    ca_patchT = new double*[nPatches];
    ca_patchRho = new double*[nPatches];
    //-------------------------------------------

    //  List of panes in this window ^^^^^^^^^^^^
    int nPanes;
    int* paneList;
    COM_get_panes(volName.c_str(), &nPanes, &paneList);
    std::cout << "  Number of Panes = "
              << nPanes << std::endl;

    //int paneIDStart = 0;
    //for(label iProc=Pstream::master(); iProc<Pstream::myProcNo(); iProc++)
    //{
    //    paneIDStart += ca_nPatches[iProc];
    //}
    //int paneIDEnd = paneIDStart+nPatches;
    //for (int ipane=paneIDStart; ipane<paneIDEnd; ++ipane)
    for (int ipane=0; ipane<nPanes; ++ipane)
    {
        int paneID = paneList[ipane];

        std::cout << "  Pane[" << ipane
             << "], paneID = " << paneID
             << " ^^^^^^^^^^^^^^^^^^^^^^^" << std::endl;

        dataName = std::string("patchName");
        nameExists(dataItemNames, dataName);
        regName = volName+std::string(".")+dataName;
        int nComp;

        COM_get_array(regName.c_str(), paneID, &ca_patchName[ipane]);
        COM_get_size(regName.c_str(), paneID, &nComp);

        patchNameStr[ipane] = ca_patchName[ipane];
        //patchNameStr[ipane].resize(nComp);
        std::cout << "    " << dataName.c_str()
             << " = " << patchNameStr[ipane].c_str() << std::endl;

        dataName = std::string("patchType");
        nameExists(dataItemNames, dataName);
        regName = volName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_patchType[ipane]);
        COM_get_size(regName.c_str(), paneID, &nComp);

        patchTypeStr[ipane] = ca_patchType[ipane];
        //patchTypeStr[ipane].resize(nComp);
        std::cout << "    " << dataName.c_str()
             << " = " << patchTypeStr[ipane].c_str() << std::endl;

        dataName = std::string("patchStart");
        nameExists(dataItemNames, dataName);
        regName = volName+std::string(".")+dataName;

        COM_get_array(regName.c_str(), paneID, &ca_patchStart[ipane]);
        std::cout << "    " << dataName.c_str()
             << " = " << *ca_patchStart[ipane] << std::endl;

        dataName = std::string("patchSize");
        nameExists(dataItemNames, dataName);
        regName = volName+std::string(".")+dataName;

        COM_get_array(regName.c_str(), paneID, &ca_patchSize[ipane]);
        std::cout << "    " << dataName.c_str()
             << " = " << *ca_patchSize[ipane] << std::endl;

        dataName = std::string("patchPointToPointMap_size");
        nameExists(dataItemNames, dataName);
        regName = volName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_patchPointToPointMap_size[ipane]);
        std::cout << "    " << dataName.c_str()
             << " = " << *ca_patchPointToPointMap_size[ipane] << std::endl;

        dataName = std::string("patchFaceToPointConn_types");
        nameExists(dataItemNames, dataName);
        regName = volName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_patchFaceToPointConn_types[ipane]);
        std::cout << "    " << dataName.c_str() << " = "
             << *ca_patchFaceToPointConn_types[ipane] << std::endl;

        dataName = std::string("patchFaceToPointConn_map");
        nameExists(dataItemNames, dataName);
        regName = volName+std::string(".")+dataName;
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
        regName = volName+std::string(".")+dataName;
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
        regName = volName+std::string(".")+dataName;
        int nPoints;
            
        COM_get_array(regName.c_str(), paneID, &ca_patchPoints[ipane], &nComp);
        COM_get_size(regName.c_str(), paneID, &nPoints);
        std::cout << "    " << dataName.c_str() << " points = " << nPoints
             << ", components = " << nComp << std::endl;
/*        for(int ipoint=0; ipoint<nPoints; ipoint++)*/
/*        {*/
/*          std::cout << "Node " << ipoint << " ca_patchPoints = ";*/
/*          for(int icomp=0; icomp<nComp; icomp++)*/
/*          {*/
/*              std::cout << *(ca_patchPoints+ipoint*nComp+icomp) << " ";*/
/*          }*/
/*          std::cout << std::endl;*/
/*        }*/

        int nConn;
        int numElem;
        std::string connNames;
        COM_get_connectivities(volName.c_str(), paneID, &nConn, connNames);
        std::istringstream connISS(connNames);

        // Secondary allocation ^^^^^^^^^^^^^^^^^
        ca_patchFaceToPointConn[ipane] = new int*[nConn];
        //---------------------------------------

        for (int icon=0; icon<nConn; ++icon)
        {
            std::string connName;
            connISS >> connName;
            //connNames.push_back(connName);

            dataName = volName+std::string(".")+connName;
            //nameExists(dataItemNames, dataName);

            COM_get_array(dataName.c_str(), paneID, &ca_patchFaceToPointConn[ipane][icon], &nComp);
            COM_get_size(dataName.c_str(), paneID, &numElem);
            
            std::cout << "    Connectivity[" << icon << "] = " << connName
                 << ", elements = " << numElem
                 << ", components =" << nComp << std::endl;
/*            for(int icell=0; icell<numElem; icell++)*/
/*            {*/
/*                std::cout << "Cell " << icell << " velocity = ";*/
/*                for(int icomp=0; icomp<nComp; icomp++)*/
/*                {*/
/*                    std::cout << *(cellVel+icell*nComp+icomp) << " ";*/
/*                }*/
/*                std::cout << std::endl;*/
/*            }*/
        }
        //---------------------------------------

        // Mapping data ^^^^^^^^^^^^^^^^^^^^^^^^^
        dataName = std::string("patchPointToPointMap");
        nameExists(dataItemNames, dataName);
        regName = volName+std::string(".")+dataName;
        
        COM_get_array(regName.c_str(), paneID, &ca_patchPointToPointMap[ipane], &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);
        std::cout << "    " << dataName.c_str() << " elements = " << numElem
             << ", components = " << nComp << std::endl;
/*        for(int icell=0; icell<numElem; icell++)*/
/*        {*/
/*            std::cout << "Cell " << icell << " velocity = ";*/
/*            for(int icomp=0; icomp<nComp; icomp++)*/
/*            {*/
/*                std::cout << *(cellVel+icell*nComp+icomp) << " ";*/
/*            }*/
/*            std::cout << std::endl;*/
/*        }*/

        dataName = std::string("patchFaceToFaceMap");
        nameExists(dataItemNames, dataName);
        regName = volName+std::string(".")+dataName;
            
        COM_get_array(regName.c_str(), paneID, &ca_patchFaceToFaceMap[ipane], &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);
        std::cout << "    " << dataName.c_str() << " elements = " << numElem
             << ", components = " << nComp << std::endl;
/*        for(int icell=0; icell<numElem; icell++)*/
/*        {*/
/*            std::cout << "Cell " << icell << " velocity = ";*/
/*            for(int icomp=0; icomp<nComp; icomp++)*/
/*            {*/
/*                std::cout << *(cellVel+icell*nComp+icomp) << " ";*/
/*            }*/
/*            std::cout << std::endl;*/
/*        }*/
        //---------------------------------------


        // Field data ^^^^^^^^^^^^^^^^^^^^^^^^^^^
        dataName = std::string("vel");
        nameExists(dataItemNames, dataName);
        regName = volName+std::string(".")+dataName;

        COM_get_array(regName.c_str(), paneID, &ca_patchVel[ipane], &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);
        std::cout << "    " << dataName.c_str() << " elements = " << numElem
             << ", components = " << nComp << std::endl;
/*        for(int icell=0; icell<numCells; icell++)*/
/*        {*/
/*            std::cout << "Cell " << icell << " velocity = ";*/
/*            for(int icomp=0; icomp<nComp; icomp++)*/
/*            {*/
/*                std::cout << *(cellVel+icell*nComp+icomp) << " ";*/
/*            }*/
/*            std::cout << std::endl;*/
/*        }*/
            
    
        dataName = std::string("pres");
        nameExists(dataItemNames, dataName);
        regName = volName+std::string(".")+dataName;

        COM_get_array(regName.c_str(), paneID, &ca_patchP[ipane], &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);
        std::cout << "    " << dataName.c_str() << " elements = " << numElem
             << ", components = " << nComp << std::endl;
/*        for(int icell=0; icell<numCells; icell++)*/
/*        {*/
/*            std::cout << "Cell " << icell << " velocity = ";*/
/*            for(int icomp=0; icomp<nComp; icomp++)*/
/*            {*/
/*                std::cout << *(cellVel+icell*nComp+icomp) << " ";*/
/*            }*/
/*            std::cout << std::endl;*/
/*        }*/

        dataName = std::string("temp");
        nameExists(dataItemNames, dataName);
        regName = volName+std::string(".")+dataName;

        COM_get_array(regName.c_str(), paneID, &ca_patchT[ipane], &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);
        std::cout << "    " << dataName.c_str() << " elements = " << numElem
             << ", components = " << nComp << std::endl;
/*        for(int icell=0; icell<numCells; icell++)*/
/*        {*/
/*            std::cout << "Cell " << icell << " velocity = ";*/
/*            for(int icomp=0; icomp<nComp; icomp++)*/
/*            {*/
/*                std::cout << *(cellVel+icell*nComp+icomp) << " ";*/
/*            }*/
/*            std::cout << std::endl;*/
/*        }*/

        dataName = std::string("rho");
        nameExists(dataItemNames, dataName);
        regName = volName+std::string(".")+dataName;

        COM_get_array(regName.c_str(), paneID, &ca_patchRho[ipane], &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);
        std::cout << "    " << dataName.c_str() << " elements = " << numElem
             << ", components = " << nComp << std::endl;
/*        for(int icell=0; icell<numCells; icell++)*/
/*        {*/
/*            std::cout << "Cell " << icell << " velocity = ";*/
/*            for(int icomp=0; icomp<nComp; icomp++)*/
/*            {*/
/*                std::cout << *(cellVel+icell*nComp+icomp) << " ";*/
/*            }*/
/*            std::cout << std::endl;*/
/*        }*/

        std::cout << "  --------------------------------------------------"
             << std::endl;
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

    if (ca_patchInGroup != NULL)
    {
        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {

            if (ca_patchInGroup[ipatch] != NULL)
            {
                delete [] ca_patchInGroup[ipatch];
                ca_patchInGroup[ipatch] = NULL;
            }
        }

        delete [] ca_patchInGroup;
        ca_patchInGroup = NULL;
    }

    if (ca_patchStart != NULL)
    {
        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {

            if (ca_patchStart[ipatch] != NULL)
            {
                delete [] ca_patchStart[ipatch];
                ca_patchStart[ipatch] = NULL;
            }
        }

        delete [] ca_patchStart;
        ca_patchStart = NULL;
    }


    if (ca_patchSize != NULL)
    {
        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {

            if (ca_patchSize[ipatch] != NULL)
            {
                delete [] ca_patchSize[ipatch];
                ca_patchSize[ipatch] = NULL;
            }
        }

        delete [] ca_patchSize;
        ca_patchSize = NULL;
    }


    if (ca_patchName != NULL)
    {
        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {

            if (&ca_patchName[ipatch] != NULL)
            {
                delete [] ca_patchName[ipatch];
                ca_patchName[ipatch] = NULL;
            }
        }

        delete [] ca_patchName;
        ca_patchName = NULL;
    }

    if (ca_patchType != NULL)
    {
        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {

            if (&ca_patchType[ipatch] != NULL)
            {
                delete [] ca_patchType[ipatch];
                ca_patchType[ipatch] = NULL;
            }
        }

        delete [] ca_patchType;
        ca_patchType = NULL;
    }

    if (patchNameStr != NULL)
    {
        delete [] patchNameStr;
        patchNameStr = NULL;
    }
    
    if (patchTypeStr != NULL)
    {
        delete [] patchTypeStr;
        patchTypeStr = NULL;
    }


    for(int ipatch=0; ipatch<nPatches; ipatch++)
    {
        if (ca_patchFaceToFaceMap[ipatch] != NULL)
        {
            delete [] ca_patchFaceToFaceMap[ipatch];
            ca_patchFaceToFaceMap[ipatch] = NULL;
        }

        if (ca_patchPoints[ipatch] != NULL)
        {
            delete [] ca_patchPoints[ipatch];
            ca_patchPoints[ipatch] = NULL;
        }

        if (ca_patchVel[ipatch] != NULL)
        {
            delete [] ca_patchVel[ipatch];
            ca_patchVel[ipatch] = NULL;
        }

        if (ca_patchRho[ipatch] != NULL)
        {
            delete [] ca_patchRho[ipatch];
            ca_patchRho[ipatch] = NULL;
        }

        if (ca_patchP[ipatch] != NULL)
        {
            delete [] ca_patchP[ipatch];
            ca_patchP[ipatch] = NULL;
        }

        if (ca_patchT[ipatch] != NULL)
        {
            delete [] ca_patchT[ipatch];
            ca_patchT[ipatch] = NULL;
        }
    }

    if (ca_patchFaceToFaceMap != NULL)
    {
        delete [] ca_patchFaceToFaceMap;
        ca_patchFaceToFaceMap = NULL;
    }

    if (ca_patchPoints != NULL)
    {
        delete [] ca_patchPoints;
        ca_patchPoints = NULL;
    }

    if (ca_patchVel != NULL)
    {
        delete [] ca_patchVel;
        ca_patchVel = NULL;
    }

    if (ca_patchRho != NULL)
    {
        delete [] ca_patchRho;
        ca_patchRho = NULL;
    }

    if (ca_patchP != NULL)
    {
        delete [] ca_patchP;
        ca_patchP = NULL;
    }

    if (ca_patchT != NULL)
    {
        delete [] ca_patchT;
        ca_patchT = NULL;
    }

    //  Delete faceToPoint connectivity arrays ^^    
    for(int ipatch=0; ipatch<nPatches; ipatch++)
    {
        int ntypes = *ca_patchFaceToPointConn_types[ipatch];
        for(int itype=0; itype<ntypes; itype++)
        {
            if (ca_patchFaceToPointConn[ipatch][itype] != NULL)
            {
                delete [] ca_patchFaceToPointConn[ipatch][itype];
                ca_patchFaceToPointConn[ipatch][itype] = NULL;
            }
        }

        if (ca_patchFaceToPointConn[ipatch] != NULL)
        {
            delete [] ca_patchFaceToPointConn[ipatch];
            ca_patchFaceToPointConn[ipatch] = NULL;
        }

        if (ca_patchFaceToPointConn_map[ipatch] != NULL)
        {
            delete [] ca_patchFaceToPointConn_map[ipatch];
            ca_patchFaceToPointConn_map[ipatch] = NULL;
        }

        if (ca_patchFaceToPointConn_size[ipatch] != NULL)
        {
            delete [] ca_patchFaceToPointConn_size[ipatch];
            ca_patchFaceToPointConn_size[ipatch] = NULL;
        }

        if (ca_patchFaceToPointConn_types[ipatch] != NULL)
        {
            delete [] ca_patchFaceToPointConn_types[ipatch];
            ca_patchFaceToPointConn_types[ipatch] = NULL;
        }
    }

    if (ca_patchFaceToPointConn != NULL)
    {
        delete [] ca_patchFaceToPointConn;
        ca_patchFaceToPointConn = NULL;
    }

    if (ca_patchFaceToPointConn != NULL)
    {
        delete [] ca_patchFaceToPointConn;
        ca_patchFaceToPointConn = NULL;
    }

    if (ca_patchFaceToPointConn_map != NULL)
    {
        delete [] ca_patchFaceToPointConn_map;
        ca_patchFaceToPointConn_map = NULL;
    }

    if (ca_patchFaceToPointConn_size != NULL)
    {
        delete [] ca_patchFaceToPointConn_size;
        ca_patchFaceToPointConn_size = NULL;
    }

    if (ca_patchFaceToPointConn_types != NULL)
    {
        delete [] ca_patchFaceToPointConn_types;
        ca_patchFaceToPointConn_types = NULL;
    }
    //-------------------------------------------

    // Delete pointToPoint mapping arrays ^^^^^^^
    for(int ipatch=0; ipatch<nPatches; ipatch++)
    {
        if (ca_patchPointToPointMap[ipatch]!= NULL)
        {
            delete [] ca_patchPointToPointMap[ipatch];
            ca_patchPointToPointMap[ipatch] = NULL;
        }

        if (ca_patchPointToPointMap_size[ipatch] != NULL)
        {
            delete [] ca_patchPointToPointMap_size[ipatch];
            ca_patchPointToPointMap_size[ipatch] = NULL;
        }
    }

    if (ca_patchPointToPointMap != NULL)
    {
        delete [] ca_patchPointToPointMap;
        ca_patchPointToPointMap = NULL;
    }

    if (ca_patchPointToPointMap_size != NULL)
    {
        delete [] ca_patchPointToPointMap_size;
        ca_patchPointToPointMap_size = NULL;
    }
    //-------------------------------------------

    if (ca_nPatches != NULL)
    {
        delete [] ca_nPatches;
        ca_nPatches = NULL;
    }

    return 0;
}

