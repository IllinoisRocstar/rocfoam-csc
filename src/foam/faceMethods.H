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
    int totalnFaceTypes = vecFaceToFaceMap.size();
    int typeCount = 0;
    for (int itype=0; itype<totalnFaceTypes; itype++)
    {
        if (vecFaceToFaceMap[itype].size()>0)
        {
            typeCount++;
        }
    }
    faceToFaceMap_types = typeCount;
    faceToFaceMap_map  = new int[typeCount];
    faceToFaceMap_size = new int[typeCount];
    ca_faceToFaceMap = new int*[typeCount];

    typeCount = 0;
    for (int itype=0; itype<totalnFaceTypes; itype++)
    {
        int nfaces = vecFaceToFaceMap[itype].size();
        if (nfaces>0)
        {
            faceToFaceMap_map[typeCount]  = itype+1;
            faceToFaceMap_size[typeCount] = nfaces;

            ca_faceToFaceMap[typeCount] = new int[nfaces];
        
            for(int iface=0; iface<nfaces; iface++)
            {
                ca_faceToFaceMap[typeCount][iface] = 
                    vecFaceToFaceMap[itype][iface];
            }
            typeCount++;
        }
    }
    //-------------------------------------------

    //  Create faceToPoint connectivity arrays ^^    
    totalnFaceTypes = vecFaceToFaceMap.size();
    typeCount = 0;
    for (int itype=0; itype<totalnFaceTypes; itype++)
    {
        int nfaces = vecFaceToPointConn[itype].size();
        if (nfaces>0)
        {
            typeCount++;
        }
    }
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
            ca_faceToPointConn_map[itype] = npoints;
            ca_faceToPointConn_size[itype] = nfaces;
    
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
    // faceToFace mapping arrays ^^^^^^^^^^^^^^^^
    int ntypes = faceToFaceMap_types;
    for(int itype=0; itype<ntypes; itype++)
    {
        if (ca_faceToFaceMap[itype] != NULL)
        {
            delete [] ca_faceToFaceMap[itype];
            ca_faceToFaceMap[itype] = NULL;
        }
    }

    if (ca_faceToFaceMap != NULL)
    {
        delete [] ca_faceToFaceMap;
        ca_faceToFaceMap = NULL;
    }

    if (faceToFaceMap_map != NULL)
    {
        delete [] faceToFaceMap_map;
        faceToFaceMap_map = NULL;
    }

    if (faceToFaceMap_size != NULL)
    {
        delete [] faceToFaceMap_size;
        faceToFaceMap_size = NULL;
    }
    //-------------------------------------------

    //  faceToPoint connectivity arrays ^^^^^^^^^
    ntypes = ca_faceToPointConn_types;
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
    const labelList&  faceOwner  = mesh.faceOwner();
    const labelList&  faceNeighb = mesh.faceNeighbour();    

    int index = 0;
    for (int itype=0; itype<faceToFaceMap_types; itype++)
    {
        int nfaces = faceToFaceMap_size[itype];
        for(int iface=0; iface<nfaces; iface++)
        {
            int faceID = ca_faceToFaceMap[itype][iface];

            ca_faceOwner[index] = faceOwner[faceID];
            ca_faceNeighb[index] = faceNeighb[faceID];

            index++;
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
    int ntypes = faceToFaceMap_types;
    for(int itype=0; itype<ntypes; itype++)
    {
        int typeID = faceToFaceMap_map[itype];
        int nfaces = faceToFaceMap_size[itype];

        dataName = volName+string(".faceToFaceMap_")+std::to_string(itype);
        COM_new_dataitem(dataName, 'p', COM_INT, 1, "");
        COM_set_size( dataName, paneID, nfaces);
        COM_set_array(dataName, paneID, ca_faceToFaceMap[itype]);
        Foam::Info << dataName << " registered." << endl;

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

