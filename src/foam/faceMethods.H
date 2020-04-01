int comFoam::createFaceConnectivities()
{

    const dynamicFvMesh& mesh(*meshPtr);

    // Mesh and conmnectivities ^^^^^^^^^^^^^^^^^
    const faceList& faces = mesh.faces();
    //-------------------------------------------

    // Temporary Vectors ^^^^^^^^^^^^^^^^^^^^^^^^
    std::vector<int> vecTmpInt;
    std::vector< std::vector<int> > vecFaceToFaceMap;
    std::vector< std::vector< std::vector<int> >> vecFaceToPointConn;
    //-------------------------------------------

    ca_nFaces  = mesh.nFaces();

    // Face-mapping vectors ^^^^^^^^^^^^^^^^^^^^^
    vecFaceToFaceMap.clear();
    vecFaceToFaceMap.resize(faceToPointTypeSize);
    for (int itype=0; itype<faceToPointTypeSize; itype++)
    {
        vecFaceToFaceMap[itype].clear();
    }

    forAll(faces, iface)
    {
        const labelList& pointsList = faces[iface];
                    
        int nPointsInFace = pointsList.size();
        vecFaceToFaceMap[nPointsInFace-1]
            .push_back(iface);
    }
    //-------------------------------------------

    // FaceToPoint connectivity vectors ^^^^^^^^^
    vecFaceToPointConn.clear();
    vecFaceToPointConn.resize(faceToPointTypeSize);
    for (int itype=0; itype<faceToPointTypeSize; itype++)
    {
        vecFaceToPointConn[itype].clear();
    }

    forAll(faces, iface)
    {
        const labelList& pointsList = faces[iface];
        int nPointsInFace = pointsList.size();
            
        vecTmpInt.clear();
        forAll(pointsList, ipoint)
        {
            const int& pointID = pointsList[ipoint];

            vecTmpInt.push_back(pointID);
        }
        vecFaceToPointConn[nPointsInFace-1].push_back(vecTmpInt);
    }
    //-------------------------------------------

    // Create faceToFace mapping arrays ^^^^^^^^^
    ca_faceToFaceMap = new int[ca_nFaces];
    
    int sortedFaceIndex = 0;
    int typeCount = 0;
    int totalnFaceTypes = vecFaceToFaceMap.size();
    for (int itype=0; itype<totalnFaceTypes; itype++)
    {
        int nfaces = vecFaceToFaceMap[itype].size();
        if (nfaces>0)
        {
            for(int iface=0; iface<nfaces; iface++)
            {
                ca_faceToFaceMap[sortedFaceIndex] =
                    vecFaceToFaceMap[itype][iface];

                sortedFaceIndex++;
            }
            typeCount++;
        }
    }
    //-------------------------------------------

    //  Create faceToPoint connectivity arrays ^^    
    ca_faceToPointConn_types = typeCount;
    ca_faceToPointConn_map   = new int[typeCount];
    ca_faceToPointConn_size  = new int[typeCount];
    ca_faceToPointConn = new int*[typeCount];

    typeCount = 0;
    for (int itype=0; itype<totalnFaceTypes; itype++)
    {
        int nfaces = vecFaceToPointConn[itype].size();
        if (nfaces>0)
        {
            int npoints = itype+1;
            ca_faceToPointConn_map[typeCount]  = npoints;
            ca_faceToPointConn_size[typeCount] = nfaces;
    
            int nTypeConn = npoints * nfaces;
            ca_faceToPointConn[typeCount] = new int[nTypeConn];
            
            for(int iface=0; iface<nfaces; iface++)
            {
                for(int ipoint=0; ipoint<npoints; ipoint++)
                {
                    int index = ipoint+iface*npoints;
                    
                    ca_faceToPointConn[typeCount][index] =
                        vecFaceToPointConn[itype][iface][ipoint];
                }
            }

            typeCount++;
        }
    }
    //-------------------------------------------

    return 0;
}


int comFoam::deleteFaceData()
{
    //  faceToPoint connectivity arrays ^^^^^^^^^
    int ntypes = ca_faceToPointConn_types;
    for(int itype=0; itype<ntypes; itype++)
    {
        if (ca_faceToPointConn[itype] != NULL)
        {
            delete [] ca_faceToPointConn[itype];
            ca_faceToPointConn[itype] = NULL;
        }
    }

    if (ca_faceToPointConn != NULL)
    {
        delete [] ca_faceToPointConn;
        ca_faceToPointConn = NULL;
    }

    if (ca_faceToPointConn_map != NULL)
    {
        delete [] ca_faceToPointConn_map;
        ca_faceToPointConn_map = NULL;
    }

    if (ca_faceToPointConn_size != NULL)
    {
        delete [] ca_faceToPointConn_size;
        ca_faceToPointConn_size = NULL;
    }
    //-------------------------------------------

    if (ca_faceOwner != NULL)
    {
        delete [] ca_faceOwner;
        ca_faceOwner = NULL;
    }

    if (ca_faceNeighb != NULL)
    {
        delete [] ca_faceNeighb;
        ca_faceNeighb = NULL;
    }

    // Connectivity-map
    if (ca_faceToFaceMap != NULL)
    {
        delete [] ca_faceToFaceMap;
        ca_faceToFaceMap = NULL;
    }

    return 0;
}


int comFoam::createFaceData()
{
    ca_faceOwner = new int[ca_nFaces];
    ca_faceNeighb = new int[ca_nFaces];

    return 0;
}

int comFoam::updateFaceData()
{
    const dynamicFvMesh& mesh(*meshPtr);
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    int nPatches = patches.size();
    int nNeighbFaces=ca_nFaces;
    if (nPatches>0)
    {
        nNeighbFaces = patches[0].start()-1;
    }


    const labelList&  faceOwner  = mesh.faceOwner();
    const labelList&  faceNeighb = mesh.faceNeighbour();    

    int faceIndex = 0;
    for (int itype=0; itype<ca_faceToPointConn_types; itype++)
    {
        int nfaces = ca_faceToPointConn_size[itype];
        for(int iface=0; iface<nfaces; iface++)
        {
            int faceID = ca_faceToFaceMap[faceIndex];

            ca_faceOwner[faceIndex] = faceOwner[faceID];
            
            if (faceIndex>nNeighbFaces)
            {
                ca_faceNeighb[faceIndex] = -1;
            }
            else
            {                
                ca_faceNeighb[faceIndex] = faceNeighb[faceID];
            }

            faceIndex++;
        }
    }

    return 0;
}


int comFoam::registerFaceData(const char *name)
{
    Foam::Info << "rocFoam.registerFaceData: "
               << "Registering flow data with name "
               << name
               << endl;

    std::string volName = name+string("VOL");
    std::string dataName = string("");


    dataName = volName+string(".nFaces");
    COM_new_dataitem( dataName, 'w', COM_INT, 1, "");
    COM_set_size(     dataName, 0, 1);
    COM_set_array(    dataName, 0, &(ca_nFaces));
    Foam::Info << dataName << " registered." << endl;

    dataName = volName+string(".faceToPointConn_types");
    COM_new_dataitem( dataName, 'w', COM_INT, 1, "");
    COM_set_size(     dataName, 0, 1);
    COM_set_array(    dataName, 0, &(ca_faceToPointConn_types));
    Foam::Info << dataName << " registered." << endl;

    dataName = volName+string(".faceToPointConn_map");
    COM_new_dataitem(dataName, 'w', COM_INT, 1, "");
    COM_set_size( dataName, 0, ca_faceToPointConn_types);
    COM_set_array(dataName, 0, ca_faceToPointConn_map);
    Foam::Info << dataName << " registered." << endl;

    dataName = volName+string(".faceToPointConn_size");
    COM_new_dataitem(dataName, 'w', COM_INT, 1, "");
    COM_set_size( dataName, 0, ca_faceToPointConn_types);
    COM_set_array(dataName, 0, ca_faceToPointConn_size);
    Foam::Info << dataName << " registered." << endl;

    int paneID = 2; // Use this paneID for face connectivity

    // Connectivity ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    int ntypes = ca_faceToPointConn_types;
    for(int itype=0; itype<ntypes; itype++)
    {
        int typeID = ca_faceToPointConn_map[itype];
        int nfaces = ca_faceToPointConn_size[itype];

        if (typeID == 3)
        { // Triangle
            dataName = volName+string(".:t3");
        }
        else if (typeID == 4)
        { // Quad
            dataName = volName+string(".:q4");
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
                       ca_faceToPointConn[itype],
                       typeID
                     );
        Foam::Info << dataName << " registered." << endl;
    }
    
    dataName = volName+string(".faceToFaceMap");
    COM_new_dataitem(dataName, 'e', COM_INT, 1, "");
    COM_set_array(dataName, paneID, ca_faceToFaceMap, 1);
    Foam::Info << dataName << " registered." << endl;
    // ------------------------------------------    

    // Field variables
    dataName = volName+string(".owner");
    COM_new_dataitem( dataName, 'e', COM_INT, 1, "");
    COM_set_array(dataName, paneID, ca_faceOwner, 1);
    Foam::Info << dataName << " registered." << endl;

    dataName = volName+string(".neighbor");
    COM_new_dataitem( dataName, 'e', COM_INT, 1, "");
    COM_set_array(dataName, paneID, ca_faceNeighb, 1);
    Foam::Info << dataName << " registered." << endl;

    COM_window_init_done(volName);

    return 0;
}

