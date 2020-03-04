int comFoam::createSurfaceConnectivities()
{

    const dynamicFvMesh& mesh(*meshPtr);

    // Mesh and conmnectivities ^^^^^^^^^^^^^^^^^
    const faceList& faces = mesh.faces();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    //-------------------------------------------

    // Temporary Vectors ^^^^^^^^^^^^^^^^^^^^^^^^
    std::vector<int> vecTmpInt;
    //-------------------------------------------

    ca_nPatches = patches.size();

    // Face-mapping vectors ^^^^^^^^^^^^^^^^^^^^^
    vecPatchFaceToFaceMap_localToLocal.clear();
    vecPatchFaceToFaceMap_localToLocal.resize(ca_nPatches);

    int faceToPointTypeSize=4;
    forAll(patches, ipatch)
    {
        vecPatchFaceToFaceMap_localToLocal[ipatch].resize(faceToPointTypeSize);
        for (int itype=0; itype<faceToPointTypeSize; itype++)
        {
            vecPatchFaceToFaceMap_localToLocal[ipatch][itype].clear();
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
            vecPatchFaceToFaceMap_localToLocal[ipatch][nPointsInFace-1]
                .push_back(iface);
        }
    }
    //-------------------------------------------

    // PointToPoint mapping Vector ^^^^^^^^^^^^^^
    vecPatchPointToPointMap_localToGlobal.clear();

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
        vecPatchPointToPointMap_localToGlobal.push_back(vecTmpInt);
    }

    // FaceToPoint connectivity vectors ^^^^^^^^^
    vecPatchFaceToPointConn_localToLocal.clear();
    vecPatchFaceToPointConn_localToLocal.resize(ca_nPatches);
    for (int ipatch=0; ipatch<ca_nPatches; ipatch++)
    {
        vecPatchFaceToPointConn_localToLocal[ipatch].resize(faceToPointTypeSize);
        for (int itype=0; itype<faceToPointTypeSize; itype++)
        {
            vecPatchFaceToPointConn_localToLocal[ipatch][itype].clear();
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
                                    vecPatchPointToPointMap_localToGlobal[ipatch].begin(),
                                    vecPatchPointToPointMap_localToGlobal[ipatch].end(),
                                    pointID
                                  );
                if (indexPtr == vecPatchPointToPointMap_localToGlobal[ipatch].end())
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
                                    vecPatchPointToPointMap_localToGlobal[ipatch].begin(),
                                    indexPtr
                                );

                vecTmpInt.push_back(indexVal);
            }
            vecPatchFaceToPointConn_localToLocal[ipatch][nPointsInFace-1].push_back(vecTmpInt);
        }
    }
    //-------------------------------------------

    // Create faceToFace mapping arrays ^^^^^^^^^
    patchFaceToFaceMap_types = new int[ca_nPatches];
    patchFaceToFaceMap_map   = new int*[ca_nPatches];
    patchFaceToFaceMap_size  = new int*[ca_nPatches];
    patchFaceToFaceMap_localToLocal = new int**[ca_nPatches];

    // face local to local
    forAll(patches, ipatch)
    {
        int totalnFaceTypes = vecPatchFaceToFaceMap_localToLocal[ipatch].size();
        int typeCount = 0;
        for (int itype=0; itype<totalnFaceTypes; itype++)
        {
            if (vecPatchFaceToFaceMap_localToLocal[ipatch][itype].size()>0)
            {
                typeCount++;
            }
        }
        patchFaceToFaceMap_types[ipatch] = typeCount;

        patchFaceToFaceMap_map[ipatch]  = new int[typeCount];
        patchFaceToFaceMap_size[ipatch] = new int[typeCount];
        patchFaceToFaceMap_localToLocal[ipatch] = new int*[typeCount];

        typeCount = 0;
        for (int itype=0; itype<totalnFaceTypes; itype++)
        {
            int nfaces = vecPatchFaceToFaceMap_localToLocal[ipatch][itype].size();
            if (nfaces>0)
            {
                patchFaceToFaceMap_map[ipatch][typeCount]  = itype+1;
                patchFaceToFaceMap_size[ipatch][typeCount] = nfaces;

                patchFaceToFaceMap_localToLocal[ipatch][typeCount] = new int[nfaces];
            
                for(int iface=0; iface<nfaces; iface++)
                {
                    patchFaceToFaceMap_localToLocal[ipatch][typeCount][iface] = 
                        vecPatchFaceToFaceMap_localToLocal[ipatch][itype][iface];
                }
                typeCount++;
            }
        }
    }
    //-------------------------------------------

    // Create pointToPoint mapping arrays ^^^^^^^
    patchPointToPointMap_size  = new int[ca_nPatches];
    patchPointToPointMap_localToGlobal = new int*[ca_nPatches];

    forAll(patches, ipatch)
    {
        int npoints = vecPatchPointToPointMap_localToGlobal[ipatch].size();
        patchPointToPointMap_size[ipatch] = npoints;
    
        patchPointToPointMap_localToGlobal[ipatch] = new int[npoints];
        
        for(int ipoint=0; ipoint< npoints; ipoint++)
        {
            patchPointToPointMap_localToGlobal[ipatch][ipoint] =
                vecPatchPointToPointMap_localToGlobal[ipatch][ipoint];
        }
    }
    //-------------------------------------------

    //  Create faceToPoint connectivity arrays ^^    
    ca_patchFaceToPointConn_types = new int[ca_nPatches];
    ca_patchFaceToPointConn_map   = new int*[ca_nPatches];
    ca_patchFaceToPointConn_size  = new int*[ca_nPatches];
    ca_patchFaceToPointConn_localToLocal = new int**[ca_nPatches];

    forAll(patches, ipatch)
    {
        int totalnFaceTypes = vecPatchFaceToPointConn_localToLocal[ipatch].size();
        int typeCount = 0;
        for (int itype=0; itype<totalnFaceTypes; itype++)
        {
            int nfaces = vecPatchFaceToPointConn_localToLocal[ipatch][itype].size();
            if (nfaces>0)
            {
                typeCount++;
            }
        }
        ca_patchFaceToPointConn_types[ipatch] = typeCount;

        ca_patchFaceToPointConn_map[ipatch]   = new int[typeCount];
        ca_patchFaceToPointConn_size[ipatch]  = new int[typeCount];
        ca_patchFaceToPointConn_localToLocal[ipatch] = new int*[typeCount];

        typeCount = 0;
        for (int itype=0; itype<totalnFaceTypes; itype++)
        {
            int nfaces = vecPatchFaceToPointConn_localToLocal[ipatch][itype].size();
            if (nfaces>0)
            {
                int npoints = itype+1;
                ca_patchFaceToPointConn_map[ipatch][itype]  = npoints;
                ca_patchFaceToPointConn_size[ipatch][itype] = nfaces;
        
                int nTypeConn = npoints * nfaces;
                
                ca_patchFaceToPointConn_localToLocal[ipatch][typeCount] = new int[nTypeConn];
                
                for(int iface=0; iface<nfaces; iface++)
                {
                    for(int ipoint=0; ipoint<npoints; ipoint++)
                    {
                        int index = ipoint+iface*npoints;
                        
                        ca_patchFaceToPointConn_localToLocal[ipatch][typeCount][index] =
                            vecPatchFaceToPointConn_localToLocal[ipatch][itype]
                                                                 [iface][ipoint];
                    }
                }

                typeCount++;
            }
        }
    }
    //-------------------------------------------

    return 0;
}


int comFoam::deleteSurfaceData()
{
    //  Create patch general data  arrays ^^^^^^^
    
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


    // Create faceToFace mapping arrays ^^^^^^^^^
    for(int ipatch=0; ipatch<ca_nPatches; ipatch++)
    {
        int ntypes = patchFaceToFaceMap_types[ipatch];
        
        for(int itype=0; itype<ntypes; itype++)
        {
            if (patchFaceToFaceMap_localToLocal[ipatch][itype] != NULL)
            {
                delete [] patchFaceToFaceMap_localToLocal[ipatch][itype];
                patchFaceToFaceMap_localToLocal[ipatch][itype] = NULL;
            }
        }

        if (patchFaceToFaceMap_localToLocal[ipatch] != NULL)
        {
            delete [] patchFaceToFaceMap_localToLocal[ipatch];
            patchFaceToFaceMap_localToLocal[ipatch] = NULL;
        }

        if (patchFaceToFaceMap_map[ipatch] != NULL)
        {
            delete [] patchFaceToFaceMap_map[ipatch];
            patchFaceToFaceMap_map[ipatch] = NULL;
        }

        if (patchFaceToFaceMap_size[ipatch] != NULL)
        {
            delete [] patchFaceToFaceMap_size[ipatch];
            patchFaceToFaceMap_size[ipatch] = NULL;
        }
    }

    if (patchFaceToFaceMap_localToLocal != NULL)
    {
        delete [] patchFaceToFaceMap_localToLocal;
        patchFaceToFaceMap_localToLocal = NULL;
    }

    if (patchFaceToFaceMap_map != NULL)
    {
        delete [] patchFaceToFaceMap_map;
        patchFaceToFaceMap_map = NULL;
    }

    if (patchFaceToFaceMap_size != NULL)
    {
        delete [] patchFaceToFaceMap_size;
        patchFaceToFaceMap_size = NULL;
    }

    if (patchFaceToFaceMap_types != NULL)
    {
        delete [] patchFaceToFaceMap_types;
        patchFaceToFaceMap_types = NULL;
    }
    //-------------------------------------------



    // Create pointToPoint mapping arrays ^^^^^^^
    for(int ipatch=0; ipatch<ca_nPatches; ipatch++)
    {
        if (patchPointToPointMap_localToGlobal[ipatch] != NULL)
        {
            delete [] patchPointToPointMap_localToGlobal[ipatch];
            patchPointToPointMap_localToGlobal[ipatch] = NULL;
        }
    }


    if (patchPointToPointMap_localToGlobal != NULL)
    {
        delete [] patchPointToPointMap_localToGlobal;
        patchPointToPointMap_localToGlobal = NULL;
    }

    if (patchPointToPointMap_size != NULL)
    {
        delete [] patchPointToPointMap_size;
        patchPointToPointMap_size = NULL;
    }
    //-------------------------------------------


    //  Create faceToPoint connectivity arrays ^^    
    
    for(int ipatch=0; ipatch<ca_nPatches; ipatch++)
    {
        int ntypes = ca_patchFaceToPointConn_types[ipatch];
        
        for(int itype=0; itype<ntypes; itype++)
        {
            if (ca_patchFaceToPointConn_localToLocal[ipatch][itype] != NULL)
            {
                delete [] ca_patchFaceToPointConn_localToLocal[ipatch][itype];
                ca_patchFaceToPointConn_localToLocal[ipatch][itype] = NULL;
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

    if (ca_patchFaceToPointConn_localToLocal != NULL)
    {
        delete [] ca_patchFaceToPointConn_localToLocal;
        ca_patchFaceToPointConn_localToLocal = NULL;
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

    vecPatchFaceToFaceMap_localToLocal.clear();
    vecPatchPointToPointMap_localToGlobal.clear();
    vecPatchFaceToPointConn_localToLocal.clear();

    return 0;
}




int comFoam::createSurfaceData()
{

    const dynamicFvMesh& mesh(*meshPtr);
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    //  Create patch general data  arrays ^^^^^^^
    ca_patchName = new std::string[ca_nPatches];
    ca_patchType = new std::string[ca_nPatches];
    ca_patchInGroup = new wordList[ca_nPatches];
    ca_patchStart = new int[ca_nPatches];
    ca_patchSize = new int[ca_nPatches];

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
 
    ca_patchPoints = new double*[ca_nPatches];
    ca_patchVel = new double*[ca_nPatches];
    ca_patchRho = new double*[ca_nPatches];
    ca_patchP = new double*[ca_nPatches];
    ca_patchT = new double*[ca_nPatches];
 
    forAll(patches, ipatch)
    {
        int npoints = patchPointToPointMap_size[ipatch];
        int nTotal = npoints * nComponents;
        
        ca_patchPoints[ipatch] = new double[nTotal];
        
        int nfaces = ca_patchSize[ipatch];
        nTotal = nfaces * nComponents;

        ca_patchVel[ipatch] = new double[nTotal];
        ca_patchRho[ipatch] = new double[nfaces];
        ca_patchP[ipatch] = new double[nfaces];
        ca_patchT[ipatch] = new double[nfaces];
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
        int npoints = patchPointToPointMap_size[ipatch];

        int localIndex = 0;
        for(int ipoint=0; ipoint<npoints; ipoint++)
        {
            int globalPointID = patchPointToPointMap_localToGlobal[ipatch][ipoint];

            for(int jcomp=0; jcomp<nComponents; jcomp++)
            {
                ca_patchPoints[ipatch][localIndex]
                    = points[globalPointID][jcomp];

                localIndex++;
            }
        }
    }

    // Cell-centered data ^^^^^^^^^^^^^^^^^^^^^^^
    const volScalarField& p(*pPtr);
    const volVectorField& U(*UPtr);
    const volScalarField& T(*TPtr);
    const volScalarField& rho(*rhoPtr);

    forAll(patches, ipatch)
    {
        if (ca_patchType[ipatch] == string("empty")) continue;

        int ntypes = patchFaceToFaceMap_types[ipatch];
        
        int localIndex = 0;
        for(int itype=0; itype<ntypes; itype++)
        {
            int nfaces = patchFaceToFaceMap_size[ipatch][itype];
            
            for(int iface=0; iface<nfaces; iface++)
            {
                int localFaceID = patchFaceToFaceMap_localToLocal[ipatch][itype][iface];
                
                for(int jcomp=0; jcomp<nComponents; jcomp++)
                {
                    int localComp = jcomp + localIndex*nComponents;
                
                    ca_patchVel[ipatch][localComp] = 
                        U.boundaryField()[ipatch][localFaceID].component(jcomp);
                }

                ca_patchP[ipatch][localIndex] = p.boundaryField()[ipatch][localFaceID];
                ca_patchT[ipatch][localIndex] = T.boundaryField()[ipatch][localFaceID];
                ca_patchRho[ipatch][localIndex] = rho.boundaryField()[ipatch][localFaceID];

                localIndex++;
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

    // Genral data registered with window
    std::string dataName = surfName+string(".patchName");
    COM_new_dataitem( dataName, 'w', COM_STRING, 1, "");
    COM_set_size(     dataName, 0, ca_nPatches);
    COM_set_array(    dataName, 0, ca_patchName);

    dataName = surfName+string(".patchType");
    COM_new_dataitem( dataName, 'w', COM_STRING, 1, "");
    COM_set_size(     dataName, 0, ca_nPatches);
    COM_set_array(    dataName, 0, ca_patchType);
    
    dataName = surfName+string(".patchInGroup");
    COM_new_dataitem( dataName, 'w', COM_STRING, 1, "");
    COM_set_size(     dataName, 0, ca_nPatches);
    COM_set_array(    dataName, 0, ca_patchInGroup);

    dataName = surfName+string(".patchInGroup");
    COM_new_dataitem( dataName, 'w', COM_STRING, 1, "");
    COM_set_size(     dataName, 0, ca_nPatches);
    COM_set_array(    dataName, 0, ca_patchInGroup);

    dataName = surfName+string(".patchStart");
    COM_new_dataitem( dataName, 'w', COM_INT, 1, "");
    COM_set_size(     dataName, 0, ca_nPatches);
    COM_set_array(    dataName, 0, ca_patchStart);

    dataName = surfName+string(".patchSize");
    COM_new_dataitem( dataName, 'w', COM_INT, 1, "");
    COM_set_size(     dataName, 0, ca_nPatches);
    COM_set_array(    dataName, 0, ca_patchSize);

    // Element data registered with window
    dataName = surfName+string(".vel");
    COM_new_dataitem( dataName, 'e', COM_DOUBLE, nComponents, "m/s");

    dataName = surfName+string(".pres");
    COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "Pa");

    dataName = surfName+string(".temp");
    COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "K");

    dataName = surfName+string(".rho");
    COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "kg/m^3");


    for(int ipatch=0; ipatch<ca_nPatches; ipatch++)
    {
        int paneID = ipatch+1;
        
        // points
        int npoints = patchPointToPointMap_size[ipatch];
        
        dataName = surfName+string(".nc");
        COM_set_size( dataName, paneID, npoints);
        COM_set_array(dataName, paneID, ca_patchPoints[ipatch], nComponents);

        // connectivity
        int ntypes = patchFaceToFaceMap_types[ipatch];
        
        for(int itype=0; itype<ntypes; itype++)
        {
            int typeID = patchFaceToFaceMap_map[ipatch][itype];
            int nfaces = patchFaceToFaceMap_size[ipatch][itype];

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
                           ca_patchFaceToPointConn_localToLocal[ipatch][itype],
                           typeID
                         );
        }
            
        // Field variables
        dataName = surfName+string(".vel");
        COM_set_array(    dataName, paneID, ca_patchVel[ipatch], nComponents);

        dataName = surfName+string(".pres");
        COM_set_array(    dataName, paneID, ca_patchP[ipatch], 1);

        dataName = surfName+string(".temp");
        COM_set_array(    dataName, paneID, ca_patchT[ipatch], 1);

        dataName = surfName+string(".rho");
        COM_set_array(    dataName, paneID, ca_patchRho[ipatch], 1);
    }

    COM_window_init_done(surfName); 

    return 0;
}

