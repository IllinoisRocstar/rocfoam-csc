#include "comFoam.H"

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
                    vecFaceToPointConn[iface][ipoint] + 1; //CGNS starts from ID 1
            }
        }
    }
    //-------------------------------------------

    return 0;
}

int comFoam::createFaceData()
{
    ca_faceOwner = new int[*ca_nFaces]{0};
    ca_faceNeighb = new int[*ca_nFaces]{0};

    if (phiPtr != nullptr)
        ca_Phi = new double[*ca_nFaces]{0};

    autoPtr<surfaceVectorField>& rhoUf(rhoUfPtr);
    if (rhoUf.valid())
    {
        int nTotal = *ca_nFaces * nComponents;
        ca_RhoUf = new double[nTotal]{0};
    }

    return 0;
}

int comFoam::updateFaceData_outgoing()
{
    const dynamicFvMesh& mesh(*meshPtr);
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const labelList&  faceOwner  = mesh.faceOwner();
    const labelList&  faceNeighb = mesh.faceNeighbour();
    const surfaceScalarField& phi(*phiPtr);
    autoPtr<surfaceVectorField>& rhoUf(rhoUfPtr);

    /*
    int nInternalFaces{0};
    {
        nInternalFaces = returnReduce
                    (
                        mesh.faceNeighbour().size(),
                        sumOp<label>()
                    );
    }
    */
    int nPatches = patches.size();
    int nInternalFaces = *ca_nFaces;
    if (nPatches>0)
    {
        nInternalFaces = patches[0].start();
    }

    int ntypes    = *ca_faceToPointConn_types;
    int faceIndex = 0;
    for (int itype=0; itype<ntypes; itype++)
    {
        int nfaces = ca_faceToPointConn_size[itype];
        for(int iface=0; iface<nfaces; iface++)
        {
            int faceID = ca_faceToFaceMap[faceIndex];

            ca_faceOwner[faceIndex] = faceOwner[faceID];

            if (faceID<nInternalFaces)
            {
                ca_faceNeighb[faceIndex] = faceNeighb[faceID];                

                if (ca_Phi != nullptr)
                    ca_Phi[faceIndex] = phi[faceID];

                if (ca_RhoUf != nullptr)
                {
                    for (int i=0; i<nComponents; i++)
                    {
                        int localComp = i + faceID*nComponents;
                        ca_RhoUf[localComp] = rhoUf()[faceID][i];
                    }
                }
            }

            faceIndex++;
        }
    }
    return 0;
}

int comFoam::registerFaceData(const char *name)
{
    std::string volName = name+std::string("VOL");
    std::stringstream output{};
    output << "rocFoam.registerFaceData: "
            << "Registering flow data with name "
            << volName;
    verbose_message(output.str(), true);

    // Use this paneID for face connectivity
    int paneID = Pstream::myProcNo()+1;
    output = std::stringstream{};
    output << "procID = " << Pstream::myProcNo()
         << ", paneID = " << paneID
         << " ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^";
    verbose_message(output.str(), true);

    std::string dataName = volName+std::string(".nFaces");
    COM_new_dataitem( dataName, 'p', COM_INT, 1, "");
    COM_set_size(     dataName, paneID, 1);
    COM_set_array(    dataName, paneID, ca_nFaces);
    output = std::stringstream{};
    output << "  " << dataName.c_str() << " registered.";
    verbose_message(output.str(), true);

    dataName = volName+std::string(".faceToPointConn_types");
    COM_new_dataitem( dataName, 'p', COM_INT, 1, "");
    COM_set_size(     dataName, paneID, 1);
    COM_set_array(    dataName, paneID, ca_faceToPointConn_types);
    output = std::stringstream{};
    output << "  " << dataName.c_str() << " registered.";
    verbose_message(output.str(), true);

    int ntypes = *ca_faceToPointConn_types;
    dataName = volName+std::string(".faceToPointConn_map");
    COM_new_dataitem(dataName, 'p', COM_INT, 1, "");
    COM_set_size( dataName, paneID, ntypes);
    COM_set_array(dataName, paneID, ca_faceToPointConn_map);
    output = std::stringstream{};
    output << "  " << dataName.c_str() << " registered.";
    verbose_message(output.str(), true);

    dataName = volName+std::string(".faceToPointConn_size");
    COM_new_dataitem(dataName, 'p', COM_INT, 1, "");
    COM_set_size( dataName, paneID, ntypes);
    COM_set_array(dataName, paneID, ca_faceToPointConn_size);
    output = std::stringstream{};
    output << "  " << dataName.c_str() << " registered.";
    verbose_message(output.str(), true);
    //-----------------------------------------------------

    // Connectivity ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    for (int itype=0; itype<ntypes; itype++)
    {
        int typeID = ca_faceToPointConn_map[itype];
        int nfaces = ca_faceToPointConn_size[itype];

        if (typeID == 3)
        { // Triangle
            dataName = volName+std::string(".facet3");
        }
        else if (typeID == 4)
        { // Quad
            dataName = volName+std::string(".faceq4");
        }
        else
        { // Type not identified
            FatalErrorInFunction
                   << "=================== ERROR ===================" << endl
                   << " Face typeID " << typeID << " with size = "
                   << nfaces << " not identified!"
                   << nl << exit(FatalError);
        }

        COM_new_dataitem(dataName, 'p', COM_INT, typeID, "");
        COM_set_size( dataName, paneID, nfaces);
        COM_set_array(dataName,
                      paneID,
                      ca_faceToPointConn[itype],
                      typeID
                     );
        output = std::stringstream{};
        output << "  " << dataName.c_str() << " registered.";
        verbose_message(output.str(), true);
    }

    dataName = volName+std::string(".faceToFaceMap");
    COM_new_dataitem(dataName, 'p', COM_INT, 1, "");
    COM_set_size( dataName, paneID, *ca_nFaces);
    COM_set_array(dataName, paneID, ca_faceToFaceMap, 1);
    output = std::stringstream{};
    output << "  " << dataName.c_str() << " registered.";
    verbose_message(output.str(), true);

    dataName = volName+std::string(".faceToFaceMap_inverse");
    COM_new_dataitem(dataName, 'p', COM_INT, 1, "");
    COM_set_size( dataName, paneID, *ca_nFaces);
    COM_set_array(dataName, paneID, ca_faceToFaceMap_inverse, 1);
    output = std::stringstream{};
    output << "  " << dataName.c_str() << " registered.";
    verbose_message(output.str(), true);
    // ------------------------------------------    

    // Field variables
    dataName = volName+std::string(".owner");
    COM_new_dataitem( dataName, 'p', COM_INT, 1, "");
    COM_set_size( dataName, paneID, *ca_nFaces);
    COM_set_array(dataName, paneID, ca_faceOwner, 1);
    output = std::stringstream{};
    output << "  " << dataName.c_str() << " registered.";
    verbose_message(output.str(), true);

    dataName = volName+std::string(".neighbor");
    COM_new_dataitem( dataName, 'p', COM_INT, 1, "");
    COM_set_size( dataName, paneID, *ca_nFaces);
    COM_set_array(dataName, paneID, ca_faceNeighb, 1);
    output = std::stringstream{};
    output << "  " << dataName.c_str() << " registered.";
    verbose_message(output.str(), true);

    if (ca_Phi != nullptr)
    {
        dataName = volName+std::string(".phi");
        COM_new_dataitem( dataName, 'p', COM_DOUBLE, 1, "kg/s");
        COM_set_size( dataName, paneID, *ca_nFaces);
        COM_set_array(    dataName, paneID, ca_Phi, 1);
        output = std::stringstream{};
        output << "  " << dataName.c_str() << " registered.";
        verbose_message(output.str(), true);
    }

    if (ca_RhoUf != nullptr)
    {
        dataName = volName+std::string(".rhoUf");
        COM_new_dataitem( dataName, 'p', COM_DOUBLE, nComponents, "kg/m^2*s");
        COM_set_size( dataName, paneID, *ca_nFaces);
        COM_set_array(    dataName, paneID, ca_RhoUf, nComponents);
        output = std::stringstream{};
        output << "  " << dataName.c_str() << " registered.";
        verbose_message(output.str(), true);
    }

    COM_window_init_done(volName);

    return 0;
}

int comFoam::reconstFaceData(const char *name)
{
    std::string volName = name+std::string("VOL");

    std::stringstream output{};
    output << "rocFoam.reconstFaceData, procID = "
              << ca_myRank
              << ", Retreiving surface data form window "
              << volName << ".";
    verbose_message(output.str(), true);

    std::string regNames;
    int numDataItems=0;
    COM_get_dataitems(volName.c_str(), &numDataItems, regNames);

    std::vector<std::string> dataItemNames;
    dataItemNames.clear();
    std::istringstream Istr(regNames);
    for (int i=0; i<numDataItems; ++i)
    {
        std::string nameTmp;
        Istr >> nameTmp;

        if (nameTmp == "nFaces" ||
            nameTmp == "faceToPointConn_types" ||
            nameTmp == "faceToPointConn_map" ||
            nameTmp == "faceToPointConn_size" ||
            nameTmp == "facet3" ||
            nameTmp == "faceq4" ||
            nameTmp == "faceToFaceMap" ||
            nameTmp == "faceToFaceMap_inverse" ||
            nameTmp == "owner" ||
            nameTmp == "neighbor" ||
            nameTmp == "phi" ||
            nameTmp == "rhoUf"
            )
        {
            dataItemNames.push_back(nameTmp);

            output = std::stringstream{};;
            output << "  DataItem[" << i << "] = " << nameTmp;
            verbose_message(output.str(), true);
        }
    }

    output = std::stringstream{};;
    output << "  Number of items = " << dataItemNames.size()
            << std::endl;
    verbose_message(output.str(), true);

    //  List of panes in this window ^^^^^^^^^^^^
    int paneID = ca_myRank+1;
    // Use this paneID for face connectivity

    output = std::stringstream{};
    output << "  procID = " << ca_myRank
            << ", paneID = " << paneID
            << " ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^";
    verbose_message(output.str(), true);

    // Face data ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    std::string dataName = std::string("nFaces");
    if (nameExists(dataItemNames, dataName))
    {
        std::string regName = volName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_nFaces);

        output = std::stringstream{};
        output << "  " << dataName.c_str() << " = " << *ca_nFaces;
        verbose_message(output.str(), true);
    }

    dataName = std::string("faceToPointConn_types");
    if (nameExists(dataItemNames, dataName))
    {
        std::string regName = volName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_faceToPointConn_types);

        output = std::stringstream{};;
        output << "  " << dataName.c_str()
                  << " = " << *ca_faceToPointConn_types;
        verbose_message(output.str(), true);
    }

    dataName = std::string("faceToPointConn_map");
    if (nameExists(dataItemNames, dataName))
    {
        std::string regName = volName+std::string(".")+dataName;
        int nComp;
        COM_get_array(regName.c_str(), paneID, &ca_faceToPointConn_map);
        COM_get_size(regName.c_str(), paneID, &nComp);
        for(int icomp=0; icomp<nComp; icomp++)
        {
            output = std::stringstream{};;
            output << "  " << dataName.c_str() << "[" << icomp << "] = "
                    << ca_faceToPointConn_map[icomp];
            verbose_message(output.str(), true);
        }
    }

    dataName = std::string("faceToPointConn_size");
    if (nameExists(dataItemNames, dataName))
    {
        std::string regName = volName+std::string(".")+dataName;
        int nComp;
        COM_get_array(regName.c_str(), paneID, &ca_faceToPointConn_size);
        COM_get_size(regName.c_str(), paneID, &nComp);
        for(int icomp=0; icomp<nComp; icomp++)
        {
            output = std::stringstream{};
            output << "  " << dataName.c_str() << "[" << icomp << "] = "
                    << ca_faceToPointConn_size[icomp];
            verbose_message(output.str(), true);
        }
    }
    //-------------------------------------------

    // Primary allocation ^^^^^^^^^^^^^^^^^^^^^^^
    ca_faceToPointConn = new int*[*ca_faceToPointConn_types];
    //-------------------------------------------
    
    // connectivity stuff ^^^^^^^^^^^^^^^^^^^^^^^
    int connID{0};
    int numElem{0};

    for (int i=0; i<*ca_faceToPointConn_types; i++)
    {
        int typeID = ca_faceToPointConn_map[i];
        if (typeID == 3)
        {
            dataName = std::string("facet3");
        }
        else if (typeID == 4)
        {
            dataName = std::string("faceq4");
        }
        else
        {
            FatalErrorInFunction
                   << "Warning: Unforeseen face connectivity."
                   << nl << exit(FatalError);
        }

        if (nameExists(dataItemNames, dataName))
        {
            std::string regName = volName+std::string(".")+dataName;
            int nComp;
            COM_get_array(regName.c_str(), paneID, &ca_faceToPointConn[connID], &nComp);
            COM_get_size(regName.c_str(), paneID, &numElem);

            output = std::stringstream{};
            output << "    Connectivity[" << connID << "] = " << dataName
                    << ", elements = " << numElem
                    << ", components =" << nComp;
            verbose_message(output.str(), true);

            connID++;
        }
    }

    if (connID != *ca_faceToPointConn_types)
    {
        FatalErrorInFunction
                  << "Error: connID != ca_faceToPointConn_types: "
                  << connID << " vs " << *ca_faceToPointConn_types
                  << nl << exit(FatalError);
    }
    //---------------------------------------

    // Mapping data ^^^^^^^^^^^^^^^^^^^^^^^^^
    dataName = std::string("faceToFaceMap");
    if (nameExists(dataItemNames, dataName))
    {
        std::string regName = volName+std::string(".")+dataName;
        int nComp;
        COM_get_array(regName.c_str(), paneID, &ca_faceToFaceMap, &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);

        output = std::stringstream{};;
        output << "    " << dataName.c_str() << " elements = " << numElem
                << ", components = " << nComp;
        verbose_message(output.str(), true);
    }

    dataName = std::string("faceToFaceMap_inverse");
    if (nameExists(dataItemNames, dataName))
    {
        std::string regName = volName+std::string(".")+dataName;
        int nComp;
        COM_get_array(regName.c_str(), paneID, &ca_faceToFaceMap_inverse, &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);

        output = std::stringstream{};
        output << "    " << dataName.c_str() << " elements = " << numElem
                << ", components = " << nComp;
        verbose_message(output.str(), true);
    }

    dataName = std::string("owner");
    if (nameExists(dataItemNames, dataName))
    {
        std::string regName = volName+std::string(".")+dataName;
        int nComp;
        COM_get_array(regName.c_str(), paneID, &ca_faceOwner, &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);

        output = std::stringstream{};
        output << "    " << dataName.c_str() << " elements = " << numElem
                << ", components = " << nComp;
        verbose_message(output.str(), true);
    }

    dataName = std::string("neighbor");
    if (nameExists(dataItemNames, dataName))
    {
        std::string regName = volName+std::string(".")+dataName;
        int nComp;
        COM_get_array(regName.c_str(), paneID, &ca_faceNeighb, &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);

        output = std::stringstream{};
        output << "    " << dataName.c_str() << " elements = " << numElem
                << ", components = " << nComp;
        verbose_message(output.str(), true);
    }

    dataName = std::string("phi");
    if (nameExists(dataItemNames, dataName))
    {
        std::string regName = volName+std::string(".")+dataName;
        int nComp;
        COM_get_array(regName.c_str(), paneID, &ca_Phi, &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);

        output = std::stringstream{};
        output << "    " << dataName.c_str() << " elements = " << numElem
                << ", components = " << nComp << std::endl;
        verbose_message(output.str(), true);
    }

    dataName = std::string("rhoUf");
    if (nameExists(dataItemNames, dataName))
    {
        std::string regName = volName+std::string(".")+dataName;
        int nComp;
        COM_get_array(regName.c_str(), paneID, &ca_RhoUf, &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);

        output = std::stringstream{};
        output << "    " << dataName.c_str() << " elements = " << numElem
                << ", components = " << nComp;
        verbose_message(output.str(), true);
    }

    output = std::stringstream{};
    output << "  --------------------------------------------------"
            << std::endl
            << "----------------------------------------------------";
    verbose_message(output.str(), true);
 
    return 0;
}

int comFoam::deleteFaceData()
{
    //  faceToPoint connectivity arrays ^^^^^^^^^
    
    if (ca_faceToPointConn_types != nullptr)
    {
        if (ca_faceToPointConn != nullptr)
        {
            int ntypes = *ca_faceToPointConn_types;
            for(int itype=0; itype<ntypes; itype++)
            {
                if (ca_faceToPointConn[itype] != nullptr)
                {
                    delete [] ca_faceToPointConn[itype];
                    ca_faceToPointConn[itype] = nullptr;
                }
            }
            delete [] ca_faceToPointConn;
            ca_faceToPointConn = nullptr;
        }

        delete[] ca_faceToPointConn_types;
        ca_faceToPointConn_types = nullptr;
    }

    if (ca_faceToPointConn_map != nullptr)
    {
        delete [] ca_faceToPointConn_map;
        ca_faceToPointConn_map = nullptr;
    }

    if (ca_faceToPointConn_size != nullptr)
    {
        delete [] ca_faceToPointConn_size;
        ca_faceToPointConn_size = nullptr;
    }
    //-------------------------------------------

    if (ca_faceOwner != nullptr)
    {
        delete [] ca_faceOwner;
        ca_faceOwner = nullptr;
    }

    if (ca_faceNeighb != nullptr)
    {
        delete [] ca_faceNeighb;
        ca_faceNeighb = nullptr;
    }

    if (ca_Phi != nullptr)
    {
        delete[] ca_Phi;
        ca_Phi = nullptr;
    }

    if (ca_RhoUf != nullptr)
    {
        delete[] ca_RhoUf;
        ca_RhoUf = nullptr;
    }

    // Connectivity-map
    if (ca_faceToFaceMap != nullptr)
    {
        delete [] ca_faceToFaceMap;
        ca_faceToFaceMap = nullptr;
    }

    // Connectivity-map
    if (ca_faceToFaceMap_inverse != nullptr)
    {
        delete [] ca_faceToFaceMap_inverse;
        ca_faceToFaceMap_inverse = nullptr;
    }

    if (ca_nFaces != nullptr)
    {
        delete ca_nFaces;
        ca_nFaces = nullptr;
    }

    return 0;
}


