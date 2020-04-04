int comFoam::createSurfaceConnectivities()
{
    const dynamicFvMesh& mesh(*meshPtr);

    // Mesh and conmnectivities ^^^^^^^^^^^^^^^^^
    const faceList& faces = mesh.faces();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    //-------------------------------------------

    //  Create patch general data  arrays ^^^^^^^
    ca_nPatches = new int(patches.size());
    
    ca_patchName    = new char*[*ca_nPatches];
    ca_patchType    = new char*[*ca_nPatches];
    //ca_patchInGroup = new wordList*[*ca_nPatches];
    ca_patchStart   = new int*[*ca_nPatches];
    ca_patchSize    = new int*[*ca_nPatches];


    forAll(patches, ipatch)
    {
        const polyPatch& patch = patches[ipatch];

        const word& patchName = patch.name();
        const word& patchType = patch.type();
        const wordList& patchInGroup = patch.inGroups();
        const int& patchStart = patch.start();
        const int& patchSize = patch.size();

        ca_patchName[ipatch] = const_cast<char*>(patchName.c_str());
        ca_patchType[ipatch] = const_cast<char*>(patchType.c_str());
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
    vecPatchFaceToFaceMap.resize(*ca_nPatches);
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
    vecPatchFaceToPointConn.resize(*ca_nPatches);
    for (int ipatch=0; ipatch<*ca_nPatches; ipatch++)
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
    ca_patchFaceToPointConn_types = new int*[*ca_nPatches];
    ca_patchFaceToPointConn_map   = new int*[*ca_nPatches];
    ca_patchFaceToPointConn_size  = new int*[*ca_nPatches];
    ca_patchFaceToPointConn = new int**[*ca_nPatches];    
    
    ca_patchFaceToFaceMap = new int*[*ca_nPatches];
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
    ca_patchPointToPointMap_size  = new int*[*ca_nPatches];
    ca_patchPointToPointMap = new int*[*ca_nPatches];

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

    // Field-data
    ca_patchPoints = new double*[*ca_nPatches];
    ca_patchVel    = new double*[*ca_nPatches];
    ca_patchRho    = new double*[*ca_nPatches];
    ca_patchP      = new double*[*ca_nPatches];
    ca_patchT      = new double*[*ca_nPatches];
 
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
                
                if (string(ca_patchType[ipatch]) == "empty")
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
    Foam::Info << "rocFoam.registerSurfaceData: "
               << "Registering flow data with name "
               << name
               << endl;

    std::string surfName = name+string("SURF");

    // Genral patch data ^^^^^^^^^^^^^^^^^^^^^^^^
    std::string dataName = surfName+string(".nPatches");
    COM_new_dataitem( dataName, 'w', COM_INT, 1, "");
    COM_set_size(     dataName, 0, 1);
    COM_set_array(    dataName, 0, ca_nPatches);
    Foam::Info << dataName << " registered." << endl;

    dataName = surfName+string(".patchName");
    COM_new_dataitem( dataName, 'p', COM_CHAR, 1, "");

    dataName = surfName+string(".patchType");
    COM_new_dataitem( dataName, 'p', COM_CHAR, 1, "");
    
    //dataName = surfName+string(".patchInGroup");
    //COM_new_dataitem( dataName, 'p', COM_CHAR, 1, "");

    dataName = surfName+string(".patchStart");
    COM_new_dataitem( dataName, 'p', COM_INT, 1, "");

    dataName = surfName+string(".patchSize");
    COM_new_dataitem( dataName, 'p', COM_INT, 1, "");
    //-------------------------------------------

    // connefctivity-relatd stuff ^^^^^^^^^^^^^^^
    // Point-connectivity
    dataName = surfName+string(".patchPointToPointMap_size");
    COM_new_dataitem( dataName, 'p', COM_INT, 1, "");

    dataName = surfName+string(".patchPointToPointMap");
    COM_new_dataitem( dataName, 'n', COM_INT, 1, "");

    // face-connectivity    
    dataName = surfName+string(".patchFaceToPointConn_types");
    COM_new_dataitem( dataName, 'p', COM_INT, 1, "");

    dataName = surfName+string(".patchFaceToPointConn_map");
    COM_new_dataitem( dataName, 'p', COM_INT, 1, "");

    dataName = surfName+string(".patchFaceToPointConn_size");
    COM_new_dataitem( dataName, 'p', COM_INT, 1, "");

    dataName = surfName+string(".patchFaceToFaceMap");
    COM_new_dataitem( dataName, 'e', COM_INT, 1, "");
    // ------------------------------------------

    // Element data registered with window ^^^^^^
    dataName = surfName+string(".vel");
    COM_new_dataitem( dataName, 'e', COM_DOUBLE, nComponents, "m/s");

    dataName = surfName+string(".pres");
    COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "Pa");

    dataName = surfName+string(".temp");
    COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "K");

    dataName = surfName+string(".rho");
    COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "kg/m^3");
    // ------------------------------------------

    // paneID>2 reserved for patches
    for(int ipatch=0; ipatch<*ca_nPatches; ipatch++)
    {
        int paneID = ipatch+3;

        Info << "Patch[" << ipatch
             <<"]^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;

        // Genral patch data ^^^^^^^^^^^^^^^^^^^^z
        std::string charToStr = string(ca_patchName[ipatch]);
        int charSize = charToStr.size();
        dataName = surfName+string(".patchName");
        COM_set_size( dataName, paneID, charSize);
        COM_set_array(dataName, paneID, ca_patchName[ipatch]);
        Foam::Info << "   patchName registered." << endl;

        charToStr = string(ca_patchType[ipatch]);
        charSize = charToStr.size();
        dataName = surfName+string(".patchType");
        COM_set_size( dataName, paneID, charSize);
        COM_set_array(dataName, paneID, ca_patchType[ipatch]);
        Foam::Info << "   patchType registered." << endl;
        
        /*
        charPtr = const_cast<char*>(ca_patchInGroup[ipatch].c_str());
        charSize = ca_patchInGroup[ipatch].size();
        dataName = surfName+string(".patchInGroup");
        COM_set_size( dataName, paneID, charSize);
        COM_set_array(dataName, paneID, charPtr);
        Foam::Info << "   patchInGroup registered." << endl;
        */

        dataName = surfName+string(".patchStart");
        COM_set_size( dataName, paneID, 1);
        COM_set_array(dataName, paneID, ca_patchStart[ipatch]);
        Foam::Info << "   patchStart registered." << endl;

        dataName = surfName+string(".patchSize");
        COM_set_size( dataName, paneID, 1);
        COM_set_array(dataName, paneID, ca_patchSize[ipatch]);
        Foam::Info << "   patchSize registered." << endl;
        //---------------------------------------
        
        // points
        dataName = surfName+string(".patchPointToPointMap_size");
        COM_set_size( dataName, paneID, 1);
        COM_set_array(dataName, paneID, ca_patchPointToPointMap_size[ipatch]);
        Foam::Info << "   patchPointToPointMap_size registered." << endl;

        int npoints = *ca_patchPointToPointMap_size[ipatch];
        dataName = surfName+string(".nc");
        COM_set_size( dataName, paneID, npoints);
        COM_set_array(dataName, paneID, ca_patchPoints[ipatch], nComponents);
        Foam::Info << "   nc registered." << endl;

        // point-mapping
        dataName = surfName+string(".patchPointToPointMap");
        COM_set_array(dataName, paneID, ca_patchPointToPointMap[ipatch], 1);
        Foam::Info << "   patchPointToPointMap registered." << endl;

        // face-connectivity
        dataName = surfName+string(".patchFaceToPointConn_types");
        COM_set_size(     dataName, paneID, 1);
        COM_set_array(    dataName, paneID, ca_patchFaceToPointConn_types[ipatch]);
        Foam::Info << "   patchFaceToPointConn_types registered." << endl;

        int ntypes = *ca_patchFaceToPointConn_types[ipatch];

        dataName = surfName+string(".patchFaceToPointConn_map");
        COM_set_size(  dataName, paneID, ntypes);
        COM_set_array( dataName, paneID, ca_patchFaceToPointConn_map[ipatch]);
        Foam::Info << "   patchFaceToPointConn_map registered." << endl;

        dataName = surfName+string(".patchFaceToPointConn_size");
        COM_set_size(     dataName, paneID, ntypes);
        COM_set_array(    dataName, paneID, ca_patchFaceToPointConn_size[ipatch]);
        Foam::Info << "   patchFaceToPointConn_size registered." << endl;

        for(int itype=0; itype<ntypes; itype++)
        {
            int typeID = ca_patchFaceToPointConn_map[ipatch][itype];
            int nfaces = ca_patchFaceToPointConn_size[ipatch][itype];

            if (typeID == 3)
            { // Triangle
                dataName = surfName+string(".:t3");
            }
            else if (typeID == 4)
            { // Quad
                dataName = surfName+string(".:q4");
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
            Foam::Info << "   " << dataName << " registered." << endl;
        }

        dataName = surfName+string(".patchFaceToFaceMap");
        COM_set_array(dataName, paneID, ca_patchFaceToFaceMap[ipatch], 1);
        Foam::Info << "   patchFaceToFaceMap registered." << endl;
        // ------------------------------------------

        // Field variables
        dataName = surfName+string(".vel");
        COM_set_array(dataName, paneID, ca_patchVel[ipatch], nComponents);
        Foam::Info << "   vel registered." << endl;

        dataName = surfName+string(".pres");
        COM_set_array(dataName, paneID, ca_patchP[ipatch], 1);
        Foam::Info << "   pres registered." << endl;

        dataName = surfName+string(".temp");
        COM_set_array(dataName, paneID, ca_patchT[ipatch], 1);
        Foam::Info << "   temp registered." << endl;

        dataName = surfName+string(".rho");
        COM_set_array(dataName, paneID, ca_patchRho[ipatch], 1);
        Foam::Info << "   rho registered." << endl;

        Info << "----------------------------------------------------"
             << endl << endl;
    }

    COM_window_init_done(surfName); 

    return 0;
}


int comFoam::reconstCaSurfaceData(const char *name)
{

    std::string volName = name+string("SURF");
    Info << "rocFoam.main: Retreiving surface data form window "
         << volName << "."
         << endl;

    std::string regNames;
    int numDataItems=0;
    
    COM_get_dataitems(volName.c_str(), &numDataItems, regNames);
    Info << "  numDataItems = " << numDataItems << endl;

    std::vector<std::string> dataItemNames;
    dataItemNames.clear();
    std::istringstream Istr(regNames);
    for (int i=0; i<numDataItems; ++i)
    {
        std::string nameTmp;
        Istr >> nameTmp;
        dataItemNames.push_back(nameTmp);
        Info << "  DataItem[" << i << "] = " << nameTmp << endl;
    }
    Info << endl;

    
    // Surface data ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    std::string dataName = std::string("nPatches");
    nameExists(dataItemNames, dataName);
    std::string regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), 0, &ca_nPatches);
    Info << "  " << dataName.c_str() << " = " << *ca_nPatches << endl;

    // Primary allocation ^^^^^^^^^^^^^^^^^^^^^^^
    ca_patchName    = new char*[*ca_nPatches];
    ca_patchType    = new char*[*ca_nPatches];
    //ca_patchInGroup = new wordList*[*ca_nPatches];
    ca_patchStart   = new int*[*ca_nPatches];
    ca_patchSize    = new int*[*ca_nPatches];

    ca_patchPointToPointMap_size = new int*[*ca_nPatches];
    ca_patchPointToPointMap = new int*[*ca_nPatches];
    ca_patchFaceToFaceMap = new int*[*ca_nPatches];

    ca_patchFaceToPointConn_types = new int*[*ca_nPatches];
    ca_patchFaceToPointConn_map = new int*[*ca_nPatches];
    ca_patchFaceToPointConn_size = new int*[*ca_nPatches];
    ca_patchFaceToPointConn = new int**[*ca_nPatches];

    ca_patchPoints = new double*[*ca_nPatches];
    ca_patchVel = new double*[*ca_nPatches];
    ca_patchP = new double*[*ca_nPatches];
    ca_patchT = new double*[*ca_nPatches];
    ca_patchRho = new double*[*ca_nPatches];
    //-------------------------------------------

    //  List of panes in this window ^^^^^^^^^^^^
    int nPanes;
    int* paneList;
    COM_get_panes(volName.c_str(), &nPanes, &paneList);
    Info << "  Number of Panes = " << nPanes << endl;

    for (int ipane=0; ipane<nPanes; ++ipane)
    {
        int paneID = paneList[ipane];

        Info << "  Pane[" << ipane << "], paneID = " << paneID
             << " ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;

        dataName = std::string("patchName");
        nameExists(dataItemNames, dataName);
        regName = volName+std::string(".")+dataName;
        int nComp;

        COM_get_array(regName.c_str(), paneID, &ca_patchName[ipane]);
        COM_get_size(regName.c_str(), paneID, &nComp);

        std::string strPatchName=ca_patchName[ipane];
        strPatchName.resize(nComp);
        Info << "    " << dataName.c_str() << " = " << strPatchName.c_str() << endl;

        dataName = std::string("patchType");
        nameExists(dataItemNames, dataName);
        regName = volName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_patchType[ipane]);
        COM_get_size(regName.c_str(), paneID, &nComp);

        std::string strPatchType=ca_patchType[ipane];
        strPatchType.resize(nComp);
        Info << "    " << dataName.c_str() << " = " << strPatchType.c_str() << endl;

        dataName = std::string("patchStart");
        nameExists(dataItemNames, dataName);
        regName = volName+std::string(".")+dataName;

        COM_get_array(regName.c_str(), paneID, &ca_patchStart[ipane]);
        Info << "    " << dataName.c_str() << " = " << *ca_patchStart[ipane] << endl;

        dataName = std::string("patchSize");
        nameExists(dataItemNames, dataName);
        regName = volName+std::string(".")+dataName;

        COM_get_array(regName.c_str(), paneID, &ca_patchSize[ipane]);
        Info << "    " << dataName.c_str() << " = " << *ca_patchSize[ipane] << endl;

        dataName = std::string("patchPointToPointMap_size");
        nameExists(dataItemNames, dataName);
        regName = volName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_patchPointToPointMap_size[ipane]);
        Info << "    " << dataName.c_str() << " = "
             << *ca_patchPointToPointMap_size[ipane] << endl;

        dataName = std::string("patchFaceToPointConn_types");
        nameExists(dataItemNames, dataName);
        regName = volName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_patchFaceToPointConn_types[ipane]);
        Info << "    " << dataName.c_str() << " = "
             << *ca_patchFaceToPointConn_types[ipane] << endl;

        dataName = std::string("patchFaceToPointConn_map");
        nameExists(dataItemNames, dataName);
        regName = volName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_patchFaceToPointConn_map[ipane]);
        COM_get_size(regName.c_str(), paneID, &nComp);
        //Info << "    " << dataName.c_str() << " size = " << nComp << endl;
        for(int icomp=0; icomp<nComp; icomp++)
        {
            Info << "    " << dataName.c_str() << "[" << icomp << "] = "
                 << ca_patchFaceToPointConn_map[ipane][icomp] << endl;
        }

        dataName = std::string("patchFaceToPointConn_size");
        nameExists(dataItemNames, dataName);
        regName = volName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_patchFaceToPointConn_size[ipane]);
        COM_get_size(regName.c_str(), paneID, &nComp);
        //Info << "    " << dataName.c_str() << " size = " << nComp << endl;
        for(int icomp=0; icomp<nComp; icomp++)
        {
            Info << "    " << dataName.c_str() << "[" << icomp << "] = "
                 << ca_patchFaceToPointConn_size[ipane][icomp] << endl;
        }

        // Point and connectivity stuff ^^^^^^^^^
        dataName = std::string("nc");
        //nameExists(dataItemNames, dataName);
        regName = volName+std::string(".")+dataName;
        int nPoints;
            
        COM_get_array(regName.c_str(), paneID, &ca_patchPoints[ipane], &nComp);
        COM_get_size(regName.c_str(), paneID, &nPoints);
        Info << "    " << dataName.c_str() << " points = " << nPoints
             << ", components = " << nComp << endl;
/*        for(int ipoint=0; ipoint<nPoints; ipoint++)*/
/*        {*/
/*          Info << "Node " << ipoint << " ca_patchPoints = ";*/
/*          for(int icomp=0; icomp<nComp; icomp++)*/
/*          {*/
/*              Info << *(ca_patchPoints+ipoint*nComp+icomp) << " ";*/
/*          }*/
/*          Info << endl;*/
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
            
            Info << "    Connectivity[" << icon << "] = " << connName
                 << ", elements = " << numElem
                 << ", components =" << nComp << endl;
/*            for(int icell=0; icell<numElem; icell++)*/
/*            {*/
/*                Info << "Cell " << icell << " velocity = ";*/
/*                for(int icomp=0; icomp<nComp; icomp++)*/
/*                {*/
/*                    Info << *(cellVel+icell*nComp+icomp) << " ";*/
/*                }*/
/*                Info << endl;*/
/*            }*/
        }
        //---------------------------------------

        // Mapping data ^^^^^^^^^^^^^^^^^^^^^^^^^
        dataName = std::string("patchPointToPointMap");
        nameExists(dataItemNames, dataName);
        regName = volName+std::string(".")+dataName;
        
        COM_get_array(regName.c_str(), paneID, &ca_patchPointToPointMap[ipane], &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);
        Info << "    " << dataName.c_str() << " elements = " << numElem
             << ", components = " << nComp << endl;
/*        for(int icell=0; icell<numElem; icell++)*/
/*        {*/
/*            Info << "Cell " << icell << " velocity = ";*/
/*            for(int icomp=0; icomp<nComp; icomp++)*/
/*            {*/
/*                Info << *(cellVel+icell*nComp+icomp) << " ";*/
/*            }*/
/*            Info << endl;*/
/*        }*/

        dataName = std::string("patchFaceToFaceMap");
        nameExists(dataItemNames, dataName);
        regName = volName+std::string(".")+dataName;
            
        COM_get_array(regName.c_str(), paneID, &ca_patchFaceToFaceMap[ipane], &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);
        Info << "    " << dataName.c_str() << " elements = " << numElem
             << ", components = " << nComp << endl;
/*        for(int icell=0; icell<numElem; icell++)*/
/*        {*/
/*            Info << "Cell " << icell << " velocity = ";*/
/*            for(int icomp=0; icomp<nComp; icomp++)*/
/*            {*/
/*                Info << *(cellVel+icell*nComp+icomp) << " ";*/
/*            }*/
/*            Info << endl;*/
/*        }*/
        //---------------------------------------


        // Field data ^^^^^^^^^^^^^^^^^^^^^^^^^^^
        dataName = std::string("vel");
        nameExists(dataItemNames, dataName);
        regName = volName+std::string(".")+dataName;

        COM_get_array(regName.c_str(), paneID, &ca_patchVel[ipane], &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);
        Info << "    " << dataName.c_str() << " elements = " << numElem
             << ", components = " << nComp << endl;
/*        for(int icell=0; icell<numCells; icell++)*/
/*        {*/
/*            Info << "Cell " << icell << " velocity = ";*/
/*            for(int icomp=0; icomp<nComp; icomp++)*/
/*            {*/
/*                Info << *(cellVel+icell*nComp+icomp) << " ";*/
/*            }*/
/*            Info << endl;*/
/*        }*/
            
    
        dataName = std::string("pres");
        nameExists(dataItemNames, dataName);
        regName = volName+std::string(".")+dataName;

        COM_get_array(regName.c_str(), paneID, &ca_patchP[ipane], &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);
        Info << "    " << dataName.c_str() << " elements = " << numElem
             << ", components = " << nComp << endl;
/*        for(int icell=0; icell<numCells; icell++)*/
/*        {*/
/*            Info << "Cell " << icell << " velocity = ";*/
/*            for(int icomp=0; icomp<nComp; icomp++)*/
/*            {*/
/*                Info << *(cellVel+icell*nComp+icomp) << " ";*/
/*            }*/
/*            Info << endl;*/
/*        }*/

        dataName = std::string("temp");
        nameExists(dataItemNames, dataName);
        regName = volName+std::string(".")+dataName;

        COM_get_array(regName.c_str(), paneID, &ca_patchT[ipane], &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);
        Info << "    " << dataName.c_str() << " elements = " << numElem
             << ", components = " << nComp << endl;
/*        for(int icell=0; icell<numCells; icell++)*/
/*        {*/
/*            Info << "Cell " << icell << " velocity = ";*/
/*            for(int icomp=0; icomp<nComp; icomp++)*/
/*            {*/
/*                Info << *(cellVel+icell*nComp+icomp) << " ";*/
/*            }*/
/*            Info << endl;*/
/*        }*/

        dataName = std::string("rho");
        nameExists(dataItemNames, dataName);
        regName = volName+std::string(".")+dataName;

        COM_get_array(regName.c_str(), paneID, &ca_patchRho[ipane], &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);
        Info << "    " << dataName.c_str() << " elements = " << numElem
             << ", components = " << nComp << endl;
/*        for(int icell=0; icell<numCells; icell++)*/
/*        {*/
/*            Info << "Cell " << icell << " velocity = ";*/
/*            for(int icomp=0; icomp<nComp; icomp++)*/
/*            {*/
/*                Info << *(cellVel+icell*nComp+icomp) << " ";*/
/*            }*/
/*            Info << endl;*/
/*        }*/

        Info << "  --------------------------------------------------"
             << endl;
    }

    Info << "----------------------------------------------------"
         << endl;

    COM_free_buffer(&paneList);
    return 0;
}

int comFoam::deleteSurfaceData()
{
    for(int ipatch=0; ipatch<*ca_nPatches; ipatch++)
    {
        if (ca_patchInGroup != NULL)
        {
            delete [] ca_patchInGroup[ipatch];
            ca_patchInGroup[ipatch] = NULL;
        }

        if (ca_patchStart != NULL)
        {
            delete [] ca_patchStart[ipatch];
            ca_patchStart[ipatch] = NULL;
        }    

        if (ca_patchSize != NULL)
        {
            delete [] ca_patchSize[ipatch];
            ca_patchSize[ipatch] = NULL;
        }    
    }

    if (ca_patchName != NULL)
    {
        delete [] ca_patchName;
        ca_patchName = NULL;
    }
    
    if (ca_patchType != NULL)
    {
        delete [] ca_patchType;
        ca_patchType = NULL;
    }

    if (ca_patchInGroup != NULL)
    {
        delete [] ca_patchInGroup;
        ca_patchInGroup = NULL;
    }

    if (ca_patchStart != NULL)
    {
        delete [] ca_patchStart;
        ca_patchStart = NULL;
    }    

    if (ca_patchSize != NULL)
    {
        delete [] ca_patchSize;
        ca_patchSize = NULL;
    }

    for(int ipatch=0; ipatch<*ca_nPatches; ipatch++)
    {
        if (ca_patchFaceToFaceMap != NULL)
        {
            delete [] ca_patchFaceToFaceMap[ipatch];
            ca_patchFaceToFaceMap[ipatch] = NULL;
        }

        if (ca_patchPoints != NULL)
        {
            delete [] ca_patchPoints[ipatch];
            ca_patchPoints[ipatch] = NULL;
        }

        if (ca_patchVel != NULL)
        {
            delete [] ca_patchVel[ipatch];
            ca_patchVel[ipatch] = NULL;
        }

        if (ca_patchRho != NULL)
        {
            delete [] ca_patchRho[ipatch];
            ca_patchRho[ipatch] = NULL;
        }

        if (ca_patchP != NULL)
        {
            delete [] ca_patchP[ipatch];
            ca_patchP[ipatch] = NULL;
        }

        if (ca_patchT != NULL)
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
    for(int ipatch=0; ipatch<*ca_nPatches; ipatch++)
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
    for(int ipatch=0; ipatch<*ca_nPatches; ipatch++)
    {
        if (ca_patchPointToPointMap != NULL)
        {
            delete [] ca_patchPointToPointMap[ipatch];
            ca_patchPointToPointMap[ipatch] = NULL;
        }

        if (ca_patchPointToPointMap_size != NULL)
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
