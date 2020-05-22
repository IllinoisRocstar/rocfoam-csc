int comFoam::createFaceConnectivities()
{

    const dynamicFvMesh& mesh(*meshPtr);

    // Mesh and conmnectivities ^^^^^^^^^^^^^^^^^
    const faceList& faces = mesh.faces();
    //-------------------------------------------

    // Temporary STLs ^^^^^^^^^^^^^^^^^^^^^^^^^^^
    std::vector<int> vecTmpInt;
    std::map<int, std::vector<int>> mapFaceToFaceMap;
    std::map<int, std::vector< std::vector<int> >> mapFaceToPointConn;
    //-------------------------------------------

    ca_nFaces = new int(mesh.nFaces());

    // Face-mapping vectors ^^^^^^^^^^^^^^^^^^^^^
    forAll(faces, iface)
    {
        const labelList& pointsList = faces[iface];
        int nPointsInFace = pointsList.size();
        
        mapFaceToFaceMap[nPointsInFace].push_back(iface);
    }
    //-------------------------------------------

    // FaceToPoint connectivity vectors ^^^^^^^^^
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
        mapFaceToPointConn[nPointsInFace].push_back(vecTmpInt);
    }
    //-------------------------------------------

    // Create faceToFace mapping arrays ^^^^^^^^^
    ca_faceToFaceMap = new int[*ca_nFaces];
    ca_faceToFaceMap_inverse = new int[*ca_nFaces];
    
    int sortedFaceIndex = 0;
    for (auto it=mapFaceToFaceMap.begin(); it!=mapFaceToFaceMap.end(); it++)
    {
        const auto& vecFaces = it->second;
        int nfaces = vecFaces.size();
        for(int iface=0; iface<nfaces; iface++)
        {
            ca_faceToFaceMap[sortedFaceIndex] = vecFaces[iface];
            ca_faceToFaceMap_inverse[vecFaces[iface]] = sortedFaceIndex;
            
            sortedFaceIndex++;
        }
    }
    //-------------------------------------------

    //  Create faceToPoint connectivity arrays ^^    
    int nTypes = mapFaceToFaceMap.size();
    ca_faceToPointConn_types = new int(nTypes);
    ca_faceToPointConn_map   = new int[nTypes];
    ca_faceToPointConn_size  = new int[nTypes];
    ca_faceToPointConn = new int*[nTypes];

    for (auto it=mapFaceToPointConn.begin(); it!=mapFaceToPointConn.end(); it++)
    {
        const auto& nPoints = it->first;
        const auto& vecFaceToPointConn = it->second;
        int nFaces = vecFaceToPointConn.size();
        int itype = std::distance(mapFaceToPointConn.begin(), it);

        ca_faceToPointConn_map[itype]  = nPoints;
        ca_faceToPointConn_size[itype] = nFaces;

        int nTypeConn = nPoints * nFaces;
        ca_faceToPointConn[itype] = new int[nTypeConn];
        
        for(int iface=0; iface<nFaces; iface++)
        {
            for(int ipoint=0; ipoint<nPoints; ipoint++)
            {
                int index = ipoint+iface*nPoints;
                
                ca_faceToPointConn[itype][index] =
                    vecFaceToPointConn[iface][ipoint];
            }
        }
    }
    //-------------------------------------------

    return 0;
}

int comFoam::createFaceData()
{
    ca_faceOwner = new int[*ca_nFaces];
    ca_faceNeighb = new int[*ca_nFaces];

    if (phiPtr != NULL)
        ca_Phi = new double[*ca_nFaces];

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
    const surfaceScalarField& phi(*phiPtr);    

    int ntypes = *ca_faceToPointConn_types;
    int faceIndex = 0;
    for (int itype=0; itype<ntypes; itype++)
    {
        int nfaces = ca_faceToPointConn_size[itype];
        for(int iface=0; iface<nfaces; iface++)
        {
            int faceID = ca_faceToFaceMap[faceIndex];

            ca_faceOwner[faceIndex] = faceOwner[faceID];

            if (faceID>nNeighbFaces)
            {
                ca_faceNeighb[faceIndex] = -1;
            }
            else
            {                
                ca_faceNeighb[faceIndex] = faceNeighb[faceID];
            }

            if (ca_Phi != NULL)
                ca_Phi[faceIndex] = phi[faceID];

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


    int paneID = Pstream::myProcNo()+1 + Pstream::nProcs();
    // Use this paneID for face connectivity

    dataName = volName+std::string(".nFaces");
    COM_new_dataitem( dataName, 'p', COM_INT, 1, "");
    COM_set_size(     dataName, paneID, 1);
    COM_set_array(    dataName, paneID, ca_nFaces);
    Info << "  " << dataName.c_str() << " registered." << endl;

    dataName = volName+std::string(".faceToPointConn_types");
    COM_new_dataitem( dataName, 'p', COM_INT, 1, "");
    COM_set_size(     dataName, paneID, 1);
    COM_set_array(    dataName, paneID, ca_faceToPointConn_types);
    Info << "  " << dataName.c_str() << " registered." << endl;

    int ntypes = *ca_faceToPointConn_types;
    dataName = volName+std::string(".faceToPointConn_map");
    COM_new_dataitem(dataName, 'p', COM_INT, 1, "");
    COM_set_size( dataName, paneID, ntypes);
    COM_set_array(dataName, paneID, ca_faceToPointConn_map);
    Info << "  " << dataName.c_str() << " registered." << endl;

    dataName = volName+std::string(".faceToPointConn_size");
    COM_new_dataitem(dataName, 'p', COM_INT, 1, "");
    COM_set_size( dataName, paneID, ntypes);
    COM_set_array(dataName, paneID, ca_faceToPointConn_size);
    Info << "  " << dataName.c_str() << " registered." << endl;

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
        Info << "  " << dataName.c_str() << " registered." << endl;
    }
    
    dataName = volName+std::string(".faceToFaceMap");
    COM_new_dataitem(dataName, 'e', COM_INT, 1, "");
    COM_set_array(dataName, paneID, ca_faceToFaceMap, 1);
    Info << "  " << dataName.c_str() << " registered." << endl;

    dataName = volName+std::string(".faceToFaceMap_inverse");
    COM_new_dataitem(dataName, 'e', COM_INT, 1, "");
    COM_set_array(dataName, paneID, ca_faceToFaceMap_inverse, 1);
    Info << "  " << dataName.c_str() << " registered." << endl;
    // ------------------------------------------    

    // Field variables
    dataName = volName+std::string(".owner");
    COM_new_dataitem( dataName, 'e', COM_INT, 1, "");
    COM_set_array(dataName, paneID, ca_faceOwner, 1);
    Info << "  " << dataName.c_str() << " registered." << endl;

    dataName = volName+std::string(".neighbor");
    COM_new_dataitem( dataName, 'e', COM_INT, 1, "");
    COM_set_array(dataName, paneID, ca_faceNeighb, 1);
    Info << "  " << dataName.c_str() << " registered." << endl;

    if (ca_Phi != NULL)
    {
        dataName = volName+std::string(".phi");
        COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "kg/s");
        COM_set_array(    dataName, paneID, ca_Phi, 1);
        Info << "  " << dataName.c_str() << " registered." << endl;
    }

    COM_window_init_done(volName);

    return 0;
}

int comFoam::reconstCaFaceData(const char *name)
{
    std::string volName = name+std::string("VOL");
/*    std::cout << "rocFoam.main: Retreiving face data form window "*/
/*         << volName << "."*/
/*         << std::endl;*/

    std::string regNames;
    int numDataItems=0;
    
    COM_get_dataitems(volName.c_str(), &numDataItems, regNames);
    //std::cout << "  numDataItems = " << numDataItems << std::endl;

    std::vector<std::string> dataItemNames;
    dataItemNames.clear();
    std::istringstream Istr(regNames);
    for (int i=0; i<numDataItems; ++i)
    {
        std::string nameTmp;
        Istr >> nameTmp;
        dataItemNames.push_back(nameTmp);
        //std::cout << "  DataItem[" << i << "] = " << nameTmp << std::endl;
    }
    std::cout << std::endl;

    
    //  List of panes in this window ^^^^^^^^^^^^
    int paneID = Pstream::myProcNo()+1 + Pstream::nProcs();
    // Use this paneID for face connectivity

    std::cout << "  procID = " << Pstream::myProcNo()
         << ", paneID = " << paneID
         << " ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << std::endl;

    // Face data ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    std::string dataName = std::string("nFaces");
    nameExists(dataItemNames, dataName);
    std::string regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), paneID, &ca_nFaces);
    std::cout << "  " << dataName.c_str() << " = " << *ca_nFaces << std::endl;

    dataName = std::string("faceToPointConn_types");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), paneID, &ca_faceToPointConn_types);
    std::cout << "  " << dataName.c_str() << " = " << *ca_faceToPointConn_types << std::endl;

    dataName = std::string("faceToPointConn_map");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    int nComp;

    COM_get_array(regName.c_str(), paneID, &ca_faceToPointConn_map);
    COM_get_size(regName.c_str(), paneID, &nComp);
    for(int icomp=0; icomp<nComp; icomp++)
    {
        std::cout << "  " << dataName.c_str() << "[" << icomp << "] = "
             << ca_faceToPointConn_map[icomp] << std::endl;
    }

    dataName = std::string("faceToPointConn_size");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;

    COM_get_array(regName.c_str(), paneID, &ca_faceToPointConn_size);
    COM_get_size(regName.c_str(), paneID, &nComp);
    for(int icomp=0; icomp<nComp; icomp++)
    {
        std::cout << "  " << dataName.c_str() << "[" << icomp << "] = "
             << ca_faceToPointConn_size[icomp] << std::endl;
    }
    //-------------------------------------------

    // Primary allocation ^^^^^^^^^^^^^^^^^^^^^^^
    ca_faceToPointConn = new int*[*ca_faceToPointConn_types];
    //-------------------------------------------
    
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

        std::cout << "    Connectivity[" << icon << "] = " << connName
             << ", elements = " << numElem
             << ", components =" << nComp << std::endl;
/*        for(int icell=0; icell<numCells; icell++)*/
/*        {*/
/*            std::cout << "Cell " << icell << " velocity = ";*/
/*            for(int icomp=0; icomp<nComp; icomp++)*/
/*            {*/
/*                std::cout << *(cellVel+icell*nComp+icomp) << " ";*/
/*            }*/
/*            std::cout << std::endl;*/
/*        }*/
    }
    //---------------------------------------

    // Mapping data ^^^^^^^^^^^^^^^^^^^^^^^^^
    dataName = std::string("faceToFaceMap");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    
    COM_get_array(regName.c_str(), paneID, &ca_faceToFaceMap, &nComp);
    COM_get_size(regName.c_str(), paneID, &numElem);
    std::cout << "    " << dataName.c_str() << " elements = " << numElem
         << ", components = " << nComp << std::endl;
/*    for(int icell=0; icell<numCells; icell++)*/
/*    {*/
/*        std::cout << "Cell " << icell << " velocity = ";*/
/*        for(int icomp=0; icomp<nComp; icomp++)*/
/*        {*/
/*            std::cout << *(cellVel+icell*nComp+icomp) << " ";*/
/*        }*/
/*        std::cout << std::endl;*/
/*    }*/


    dataName = std::string("faceToFaceMap_inverse");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;

    COM_get_array(regName.c_str(), paneID, &ca_faceToFaceMap_inverse, &nComp);
    COM_get_size(regName.c_str(), paneID, &numElem);
    std::cout << "    " << dataName.c_str() << " elements = " << numElem
         << ", components = " << nComp << std::endl;



    dataName = std::string("owner");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;

    COM_get_array(regName.c_str(), paneID, &ca_faceOwner, &nComp);
    COM_get_size(regName.c_str(), paneID, &numElem);
    std::cout << "    " << dataName.c_str() << " elements = " << numElem
         << ", components = " << nComp << std::endl;
/*    for(int iface=0; iface<numElem; iface++)*/
/*    {*/
/*        std::cout << "Face " << iface << " owner = ";*/
/*        for(int icomp=0; icomp<nComp; icomp++)*/
/*        {*/
/*            std::cout << *(Owner+iface*nComp+icomp) << " ";*/
/*        }*/
/*        std::cout << std::endl;*/
/*    }*/

    dataName = std::string("neighbor");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;

    COM_get_array(regName.c_str(), paneID, &ca_faceNeighb, &nComp);
    COM_get_size(regName.c_str(), paneID, &numElem);
    std::cout << "    " << dataName.c_str() << " elements = " << numElem
         << ", components = " << nComp << std::endl;
/*    for(int iface=0; iface<numElem; iface++)*/
/*    {*/
/*        std::cout << "Face " << iface << " neighb = ";*/
/*        for(int icomp=0; icomp<nComp; icomp++)*/
/*        {*/
/*            std::cout << *(ca_faceNeighb+iface*nComp+icomp) << " ";*/
/*        }*/
/*        std::cout << std::endl;*/
/*    }*/

    dataName = std::string("phi");
    if (nameExists(dataItemNames, dataName))
    {
        regName = volName+std::string(".")+dataName;

        COM_get_array(regName.c_str(), paneID, &ca_Phi, &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);
        std::cout << "    " << dataName.c_str() << " elements = " << numElem
             << ", components = " << nComp << std::endl;
    }

    std::cout << "  --------------------------------------------------"
         << std::endl;

    std::cout << "----------------------------------------------------"
         << std::endl;

    return 0;
}

int comFoam::deleteFaceData()
{
    //  faceToPoint connectivity arrays ^^^^^^^^^
    
    if (ca_faceToPointConn_types != NULL)
    {
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

        delete[] ca_faceToPointConn_types;
        ca_faceToPointConn_types = NULL;
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

    if (ca_Phi != NULL)
    {
        delete[] ca_Phi;
        ca_Phi = NULL;
    }

    // Connectivity-map
    if (ca_faceToFaceMap != NULL)
    {
        delete [] ca_faceToFaceMap;
        ca_faceToFaceMap = NULL;
    }

    // Connectivity-map
    if (ca_faceToFaceMap_inverse != NULL)
    {
        delete [] ca_faceToFaceMap_inverse;
        ca_faceToFaceMap_inverse = NULL;
    }

    if (ca_nFaces != NULL)
    {
        delete[] ca_nFaces;
        ca_nFaces = NULL;
    }

    return 0;
}


