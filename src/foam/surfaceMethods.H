int comFoam::createSurfaceConnectivities()
{
    const dynamicFvMesh& mesh(*meshPtr);

    // Mesh and conmnectivities ^^^^^^^^^^^^^^^^^
    const faceList& faces = mesh.faces();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    //-------------------------------------------

    //  Create patch general data  arrays ^^^^^^^
    ca_nPatches = patches.size();
    ca_patchName    = new std::string[ca_nPatches];
    ca_patchType    = new std::string[ca_nPatches];
    ca_patchInGroup = new wordList[ca_nPatches];
    ca_patchStart   = new int[ca_nPatches];
    ca_patchSize    = new int[ca_nPatches];



    forAll(patches, ipatch)
    {
        const polyPatch& patch = patches[ipatch];

        const word& patchName = patch.name();
        const word& patchType = patch.type();
        const wordList& patchInGroup = patch.inGroups();
        const int& patchStart = patch.start();
        const int& patchSize = patch.size();

        ca_patchName[ipatch] = patchName;
        ca_patchType[ipatch] = patchType;
        ca_patchInGroup[ipatch] = patchInGroup;
        ca_patchStart[ipatch] = patchStart;
        ca_patchSize[ipatch] = patchSize;
    }
    //-------------------------------------------

//patchName2 = const_cast<char*>(ca_patchName[0].c_str());
/*std::string patchName2 = ca_patchName[0] + string(" ");*/

/*Info << "test = " << patchName2 << endl;*/
/*Info << "test.size() = " <<patchName2.size() << endl;*/
/*Info << "test.length() = " <<patchName2.length() << endl;*/


/*patchName2 = const_cast<char*>(ca_patchName[0].c_str());*/
/*char* a= const_cast<char*>(" ");*/
/*Info << "patchName2 = " << patchName2 << endl;*/
/*Info << "size(patchName2) = " <<sizeof(patchName2) << endl;*/
/*Info << "size(char) = " <<sizeof(char) << endl;*/
/*Info << "size(patchName2) = " <<sizeof(patchName2)/sizeof(char) << endl;*/
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
    vecPatchFaceToFaceMap.resize(ca_nPatches);
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
    vecPatchFaceToPointConn.resize(ca_nPatches);
    for (int ipatch=0; ipatch<ca_nPatches; ipatch++)
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
    ca_patchFaceToPointConn_types = new int[ca_nPatches];
    ca_patchFaceToPointConn_map   = new int*[ca_nPatches];
    ca_patchFaceToPointConn_size  = new int*[ca_nPatches];
    ca_patchFaceToPointConn = new int**[ca_nPatches];    
    
    ca_patchFaceToFaceMap = new int*[ca_nPatches];
    forAll(patches, ipatch)
    {
        // Create faceToFace mapping arrays ^^^^^^^^^
        int nfacesTotal = ca_patchSize[ipatch];
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
        ca_patchFaceToPointConn_types[ipatch] = typeCount;
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
    ca_patchPointToPointMap_size  = new int[ca_nPatches];
    ca_patchPointToPointMap = new int*[ca_nPatches];

    forAll(patches, ipatch)
    {
        int npoints = vecPatchPointToPointMap[ipatch].size();
        ca_patchPointToPointMap_size[ipatch] = npoints;
    
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

int comFoam::deleteSurfaceData()
{
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


    for(int ipatch=0; ipatch<ca_nPatches; ipatch++)
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
    for(int ipatch=0; ipatch<ca_nPatches; ipatch++)
    {
        int ntypes = ca_patchFaceToPointConn_types[ipatch];
        
        for(int itype=0; itype<ntypes; itype++)
        {
            if (ca_patchFaceToPointConn[ipatch][itype] != NULL)
            {
                delete [] ca_patchFaceToPointConn[ipatch][itype];
                ca_patchFaceToPointConn[ipatch][itype] = NULL;
            }
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
    for(int ipatch=0; ipatch<ca_nPatches; ipatch++)
    {
        if (ca_patchPointToPointMap[ipatch] != NULL)
        {
            delete [] ca_patchPointToPointMap[ipatch];
            ca_patchPointToPointMap[ipatch] = NULL;
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

    return 0;
}


int comFoam::createSurfaceData()
{
    const dynamicFvMesh& mesh(*meshPtr);
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Field-data
    ca_patchPoints = new double*[ca_nPatches];
    ca_patchVel    = new double*[ca_nPatches];
    ca_patchRho    = new double*[ca_nPatches];
    ca_patchP      = new double*[ca_nPatches];
    ca_patchT      = new double*[ca_nPatches];
 
    forAll(patches, ipatch)
    {
        // Points
        int npoints = ca_patchPointToPointMap_size[ipatch];
        int nTotal = npoints * nComponents;
        ca_patchPoints[ipatch] = new double[nTotal];
        
        // Field-data
        int nfaces = ca_patchSize[ipatch];
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
        int npoints = ca_patchPointToPointMap_size[ipatch];

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
        if (ca_patchType[ipatch] == string("empty")) continue;

        int ntypes = ca_patchFaceToPointConn_types[ipatch];
        
        int faceIndex = 0;
        for(int itype=0; itype<ntypes; itype++)
        {
            int nfaces = ca_patchFaceToPointConn_size[ipatch][itype];
            
            for(int iface=0; iface<nfaces; iface++)
            {
                int localFaceID = ca_patchFaceToFaceMap[ipatch][faceIndex];
                
                for(int jcomp=0; jcomp<nComponents; jcomp++)
                {
                    int localComp = jcomp + faceIndex*nComponents;
                
                    ca_patchVel[ipatch][localComp] = 
                        U.boundaryField()[ipatch][localFaceID].component(jcomp);
                }

                ca_patchP[ipatch][faceIndex] = p.boundaryField()[ipatch][localFaceID];
                ca_patchT[ipatch][faceIndex] = T.boundaryField()[ipatch][localFaceID];
                ca_patchRho[ipatch][faceIndex] = rho.boundaryField()[ipatch][localFaceID];
                
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

/*int size = patchName2.size();*/
/*dataName = surfName+string(".patchName2");*/
/*COM_new_dataitem( dataName, 'w', COM_CHAR, 1, "");*/
/*COM_set_size(     dataName, 0, patchName2.size());*/
/*COM_set_array(    dataName, 0, &patchName2);*/
/*Foam::Info << dataName << " registered." << endl;*/

    // Genral patch data ^^^^^^^^^^^^^^^^^^^^^^^^
    std::string dataName = surfName+string(".nPatches");
    COM_new_dataitem( dataName, 'w', COM_INT, 1, "");
    COM_set_size(     dataName, 0, 1);
    COM_set_array(    dataName, 0, &(ca_nPatches));
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


    for(int ipatch=0; ipatch<ca_nPatches; ipatch++)
    {
        int paneID = ipatch+1;

        Info << "Patch[" << ipatch
             <<"]^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;

        // Genral patch data ^^^^^^^^^^^^^^^^^^^^z
        char* charPtr = const_cast<char*>(ca_patchName[ipatch].c_str());
        int charSize = ca_patchName[ipatch].size();
        dataName = surfName+string(".patchName");
        COM_set_size( dataName, paneID, charSize);
        COM_set_array(dataName, paneID, charPtr);
        Foam::Info << "   patchName registered." << endl;

        charPtr = const_cast<char*>(ca_patchType[ipatch].c_str());
        charSize = ca_patchType[ipatch].size();
        dataName = surfName+string(".patchType");
        COM_set_size( dataName, paneID, charSize);
        COM_set_array(dataName, paneID, charPtr);
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
        COM_set_array(dataName, paneID, &ca_patchStart[ipatch]);
        Foam::Info << "   patchStart registered." << endl;

        dataName = surfName+string(".patchSize");
        COM_set_size( dataName, paneID, 1);
        COM_set_array(dataName, paneID, &ca_patchSize[ipatch]);
        Foam::Info << "   patchSize registered." << endl;
        //---------------------------------------
        
        // points
        dataName = surfName+string(".patchPointToPointMap_size");
        COM_set_size( dataName, paneID, 1);
        COM_set_array(dataName, paneID, &ca_patchPointToPointMap_size[ipatch]);
        Foam::Info << "   patchPointToPointMap_size registered." << endl;

        int npoints = ca_patchPointToPointMap_size[ipatch];
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
        COM_set_array(    dataName, paneID, &ca_patchFaceToPointConn_types[ipatch]);
        Foam::Info << "   patchFaceToPointConn_types registered." << endl;

        int ntypes = ca_patchFaceToPointConn_types[ipatch];

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

