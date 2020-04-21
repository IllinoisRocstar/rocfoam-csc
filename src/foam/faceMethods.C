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

    ca_nFaces = new int(mesh.nFaces());

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
    ca_faceToFaceMap = new int[*ca_nFaces];
    
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
    ca_faceToPointConn_types = new int(typeCount);
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

int comFoam::createFaceData()
{
    ca_faceOwner = new int[*ca_nFaces];
    ca_faceNeighb = new int[*ca_nFaces];

    return 0;
}

int comFoam::updateFaceData()
{
    const dynamicFvMesh& mesh(*meshPtr);
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    int nPatches = patches.size();
    int nNeighbFaces = *ca_nFaces;
    if (nPatches>0)
    {
        nNeighbFaces = patches[0].start()-1;
    }


    const labelList&  faceOwner  = mesh.faceOwner();
    const labelList&  faceNeighb = mesh.faceNeighbour();    

    int ntypes = *ca_faceToPointConn_types;
    int faceIndex = 0;
    for (int itype=0; itype<ntypes; itype++)
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

    std::string volName = name+std::string("VOL");
    std::string dataName = std::string("");


    dataName = volName+std::string(".nFaces");
    COM_new_dataitem( dataName, 'w', COM_INT, 1, "");
    COM_set_size(     dataName, 0, 1);
    COM_set_array(    dataName, 0, ca_nFaces);
    Foam::Info << dataName << " registered." << endl;

    dataName = volName+std::string(".faceToPointConn_types");
    COM_new_dataitem( dataName, 'w', COM_INT, 1, "");
    COM_set_size(     dataName, 0, 1);
    COM_set_array(    dataName, 0, ca_faceToPointConn_types);
    Foam::Info << dataName << " registered." << endl;

    int ntypes = *ca_faceToPointConn_types;
    dataName = volName+std::string(".faceToPointConn_map");
    COM_new_dataitem(dataName, 'w', COM_INT, 1, "");
    COM_set_size( dataName, 0, ntypes);
    COM_set_array(dataName, 0, ca_faceToPointConn_map);
    Foam::Info << dataName << " registered." << endl;

    dataName = volName+std::string(".faceToPointConn_size");
    COM_new_dataitem(dataName, 'w', COM_INT, 1, "");
    COM_set_size( dataName, 0, ntypes);
    COM_set_array(dataName, 0, ca_faceToPointConn_size);
    Foam::Info << dataName << " registered." << endl;

    int paneID = 2; // Use this paneID for face connectivity

    // Connectivity ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    for(int itype=0; itype<ntypes; itype++)
    {
        int typeID = ca_faceToPointConn_map[itype];
        int nfaces = ca_faceToPointConn_size[itype];

        if (typeID == 3)
        { // Triangle
            dataName = volName+std::string(".:t3");
        }
        else if (typeID == 4)
        { // Quad
            dataName = volName+std::string(".:q4");
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
    
    dataName = volName+std::string(".faceToFaceMap");
    COM_new_dataitem(dataName, 'e', COM_INT, 1, "");
    COM_set_array(dataName, paneID, ca_faceToFaceMap, 1);
    Foam::Info << dataName << " registered." << endl;
    // ------------------------------------------    

    // Field variables
    dataName = volName+std::string(".owner");
    COM_new_dataitem( dataName, 'e', COM_INT, 1, "");
    COM_set_array(dataName, paneID, ca_faceOwner, 1);
    Foam::Info << dataName << " registered." << endl;

    dataName = volName+std::string(".neighbor");
    COM_new_dataitem( dataName, 'e', COM_INT, 1, "");
    COM_set_array(dataName, paneID, ca_faceNeighb, 1);
    Foam::Info << dataName << " registered." << endl;

    COM_window_init_done(volName);

    return 0;
}

int comFoam::reconstCaFaceData(const char *name)
{

    std::string volName = name+std::string("VOL");
/*    Info << "rocFoam.main: Retreiving face data form window "*/
/*         << volName << "."*/
/*         << endl;*/

    std::string regNames;
    int numDataItems=0;
    
    COM_get_dataitems(volName.c_str(), &numDataItems, regNames);
    //Info << "  numDataItems = " << numDataItems << endl;

    std::vector<std::string> dataItemNames;
    dataItemNames.clear();
    std::istringstream Istr(regNames);
    for (int i=0; i<numDataItems; ++i)
    {
        std::string nameTmp;
        Istr >> nameTmp;
        dataItemNames.push_back(nameTmp);
        //Info << "  DataItem[" << i << "] = " << nameTmp << endl;
    }
    Info << endl;

    
    // Face data ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    std::string dataName = std::string("nFaces");
    nameExists(dataItemNames, dataName);
    std::string regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), 0, &ca_nFaces);
    Info << "  " << dataName.c_str() << " = " << *ca_nFaces << endl;

    dataName = std::string("faceToPointConn_types");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), 0, &ca_faceToPointConn_types);
    Info << "  " << dataName.c_str() << " = " << *ca_faceToPointConn_types << endl;

    dataName = std::string("faceToPointConn_map");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    int nComp;

    COM_get_array(regName.c_str(), 0, &ca_faceToPointConn_map);
    COM_get_size(regName.c_str(), 0, &nComp);
    for(int icomp=0; icomp<nComp; icomp++)
    {
        Info << "  " << dataName.c_str() << "[" << icomp << "] = "
             << ca_faceToPointConn_map[icomp] << endl;
    }

    dataName = std::string("faceToPointConn_size");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;

    COM_get_array(regName.c_str(), 0, &ca_faceToPointConn_size);
    COM_get_size(regName.c_str(), 0, &nComp);
    for(int icomp=0; icomp<nComp; icomp++)
    {
        Info << "  " << dataName.c_str() << "[" << icomp << "] = "
             << ca_faceToPointConn_size[icomp] << endl;
    }
    //-------------------------------------------

    // Primary allocation ^^^^^^^^^^^^^^^^^^^^^^^
    ca_faceToPointConn = new int*[*ca_faceToPointConn_types];
    //-------------------------------------------
    
    //  List of panes in this window ^^^^^^^^^^^^
    int ipane = 1;
    int paneID = 2;

    Info << "  Pane[" << ipane << "], paneID = " << paneID
         << " ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;

    // connectivity stuff ^^^^^^^^^^^^^^^^^^^^^^^
    int nConn;
    int numElem;
    std::string connNames;
    COM_get_connectivities(volName.c_str(), paneID, &nConn, connNames);
    std::istringstream connISS(connNames);
    for (int icon=0; icon<nConn; ++icon)
    {
        std::string connName;
        connISS >> connName;
        //connNames.push_back(connName);

        dataName = volName+std::string(".")+connName;
        //nameExists(dataItemNames, dataName);
        COM_get_array(dataName.c_str(), paneID, &ca_faceToPointConn[icon], &nComp);
        COM_get_size(dataName.c_str(), paneID, &numElem);

        Info << "    Connectivity[" << icon << "] = " << connName
             << ", elements = " << numElem
             << ", components =" << nComp << endl;
/*        for(int icell=0; icell<numCells; icell++)*/
/*        {*/
/*            Info << "Cell " << icell << " velocity = ";*/
/*            for(int icomp=0; icomp<nComp; icomp++)*/
/*            {*/
/*                Info << *(cellVel+icell*nComp+icomp) << " ";*/
/*            }*/
/*            Info << endl;*/
/*        }*/
    }
    //---------------------------------------

    // Mapping data ^^^^^^^^^^^^^^^^^^^^^^^^^
    dataName = std::string("faceToFaceMap");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    
    COM_get_array(regName.c_str(), paneID, &ca_faceToFaceMap, &nComp);
    COM_get_size(regName.c_str(), paneID, &numElem);
    Info << "    " << dataName.c_str() << " elements = " << numElem
         << ", components = " << nComp << endl;
/*    for(int icell=0; icell<numCells; icell++)*/
/*    {*/
/*        Info << "Cell " << icell << " velocity = ";*/
/*        for(int icomp=0; icomp<nComp; icomp++)*/
/*        {*/
/*            Info << *(cellVel+icell*nComp+icomp) << " ";*/
/*        }*/
/*        Info << endl;*/
/*    }*/

    dataName = std::string("owner");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;

    COM_get_array(regName.c_str(), paneID, &ca_faceOwner, &nComp);
    COM_get_size(regName.c_str(), paneID, &numElem);
    Info << "    " << dataName.c_str() << " elements = " << numElem
         << ", components = " << nComp << endl;
/*    for(int iface=0; iface<numElem; iface++)*/
/*    {*/
/*        Info << "Face " << iface << " owner = ";*/
/*        for(int icomp=0; icomp<nComp; icomp++)*/
/*        {*/
/*            Info << *(Owner+iface*nComp+icomp) << " ";*/
/*        }*/
/*        Info << endl;*/
/*    }*/

    dataName = std::string("neighbor");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;

    COM_get_array(regName.c_str(), paneID, &ca_faceNeighb, &nComp);
    COM_get_size(regName.c_str(), paneID, &numElem);
    Info << "    " << dataName.c_str() << " elements = " << numElem
         << ", components = " << nComp << endl;
/*    for(int iface=0; iface<numElem; iface++)*/
/*    {*/
/*        Info << "Face " << iface << " neighb = ";*/
/*        for(int icomp=0; icomp<nComp; icomp++)*/
/*        {*/
/*            Info << *(ca_faceNeighb+iface*nComp+icomp) << " ";*/
/*        }*/
/*        Info << endl;*/
/*    }*/

    Info << "  --------------------------------------------------"
         << endl;

    Info << "----------------------------------------------------"
         << endl;

    return 0;
}

int comFoam::deleteFaceData()
{
    //  faceToPoint connectivity arrays ^^^^^^^^^
    int ntypes = *ca_faceToPointConn_types;

    if (ca_faceToPointConn != NULL)
    {
        for(int itype=0; itype<ntypes; itype++)
        {
            if (ca_faceToPointConn[itype] != NULL)
            {
                delete [] ca_faceToPointConn[itype];
                ca_faceToPointConn[itype] = NULL;
            }
        }
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

    if (ca_nFaces != NULL)
    {
        delete[] ca_nFaces;
        ca_nFaces = NULL;
    }

    if (ca_faceToPointConn_types != NULL)
    {
        delete[] ca_faceToPointConn_types;
        ca_faceToPointConn_types = NULL;
    }

    return 0;
}


