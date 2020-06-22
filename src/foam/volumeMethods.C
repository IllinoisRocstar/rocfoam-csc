int comFoam::createVolumeConnectivities()
{
    const dynamicFvMesh& mesh(*meshPtr);

    // Mesh and conmnectivities ^^^^^^^^^^^^^^^^^
    const labelListList& cellPoints = mesh.cellPoints();
    const cellShapeList& cellShapes = mesh.cellShapes();
    //-------------------------------------------

    ca_nPoints = new int(mesh.nPoints());
    ca_nCells  = new int(mesh.nCells());

    // Temporary STLs ^^^^^^^^^^^^^^^^^^^^^^^^^^^
    std::map<int, std::vector<int> > mapCellToCellMap;
    //-------------------------------------------

    // cellToPoint connectivities ^^^^^^^^^^^^^^^
    forAll(cellPoints, icell)
    {
        const labelList& pointsList = cellPoints[icell];
        int nPointsInCell = pointsList.size();
        mapCellToCellMap[nPointsInCell].push_back(icell);
    }
    //-------------------------------------------

    // Create cellToCell mapping ^^^^^^^^^^^^^^^^
    ca_cellToCellMap = new int[*ca_nCells];
    ca_cellToCellMap_inverse = new int[*ca_nCells];

    int sortedCellIndex = 0;
    for(auto it: mapCellToCellMap)
    {
        const auto& vecCells = it.second;
        int nCells = vecCells.size();
        for(int icell=0; icell<nCells; icell++)
        {
            ca_cellToCellMap[sortedCellIndex] = vecCells[icell];
            ca_cellToCellMap_inverse[vecCells[icell]] = sortedCellIndex;

            sortedCellIndex++;
        }
    }
    //-------------------------------------------

    if (sortedCellIndex != *ca_nCells)
    {
        Info << "========== WARNNING ===============" << endl
             << "     sortedCellIndex != ca_nCells " << endl
             << "    " << sortedCellIndex << "!=" << *ca_nCells << endl;
        return -1;
    }

    // cetToPoint Connectivity ^^^^^^^^^^^^^^^^^^
    int nTypes = mapCellToCellMap.size();
    ca_cellToPointConn_types = new int(nTypes);
    ca_cellToPointConn_map   = new int[nTypes];
    ca_cellToPointConn_size  = new int[nTypes];
    ca_cellToPointConn = new int*[nTypes];

    sortedCellIndex = 0;
    for(auto it = mapCellToCellMap.begin(); it != mapCellToCellMap.end(); it++)
    {
        const auto& nPoints = it->first;
        const auto& vecCells = it->second;
        int nCells = vecCells.size();
        int itype = std::distance(mapCellToCellMap.begin(), it);

        ca_cellToPointConn_map[itype] = nPoints;
        ca_cellToPointConn_size[itype] = nCells;

        int nTypeConn = nCells * nPoints;
        ca_cellToPointConn[itype] = new int[nTypeConn];

        for (int icell=0; icell<nCells; icell++)
        {
            int cellID = ca_cellToCellMap[sortedCellIndex];
            const cellShape& cellShape_ = cellShapes[cellID];

            for (int ipoint=0; ipoint<nPoints; ipoint++)
            {
                int index = icell*nPoints + ipoint;

                ca_cellToPointConn[itype][index] = cellShape_[ipoint];
            }
            sortedCellIndex++;
        }
    }
    //-------------------------------------------

    return 0;
}

int comFoam::createVolumeData()
{
    int nTotal = *ca_nPoints * nComponents;
    ca_Points = new double[nTotal]{0};
    ca_Disp   = new double[nTotal]{0};

    // Field-data
    nTotal = *ca_nCells * nComponents;
    ca_Vel = new double[nTotal]{0};
    ca_P   = new double[*ca_nCells]{0};
    
    if (TPtr != nullptr)
        ca_T = new double[*ca_nCells]{0};

    if (rhoPtr != nullptr)
        ca_Rho = new double[*ca_nCells]{0};

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
            ca_AlphaT = new double[*ca_nCells]{0};
            ca_Epsilon = new double[*ca_nCells]{0};
            ca_K = new double[*ca_nCells]{0};
            ca_NuT = new double[*ca_nCells]{0};
        }
    }

    return 0;
}

int comFoam::updateVolumeData_outgoing()
{
    // Point data ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    const dynamicFvMesh& mesh(*meshPtr);
    const pointField&    points = mesh.points();

    /*
    //if (mesh.moving())
        pointVectorField& pointDisplacement = const_cast<pointVectorField&>
        (
            mesh().objectRegistry::lookupObject<pointVectorField>
            (
                "pointDisplacement"
            )
        );
    */

    forAll(points, ipoint)
    {
        for(int jcomp=0; jcomp<nComponents; jcomp++)
        {
            ca_Points[ipoint*nComponents+jcomp]
                = points[ipoint][jcomp];

            /*
            if (mesh.moving())
                ca_Disp[ipoint*nComponents+jcomp]
                = pointDisplacement[ipoint][jcomp];
            */
        }
    }
    
    
    // Cell-centered data ^^^^^^^^^^^^^^^^^^^^^^^
    const volScalarField& p(*pPtr);
    const volVectorField& U(*UPtr);
    const volScalarField& T(*TPtr);
    const volScalarField& rho(*rhoPtr);
    compressible::turbulenceModel &turbulence(*turbulencePtr);

    int cellIndex = 0;
    for(int itype=0; itype<*ca_cellToPointConn_types; itype++)
    {
        int ncells = ca_cellToPointConn_size[itype];
        for(int icell=0; icell<ncells; icell++)
        {
            int cellID = ca_cellToCellMap[cellIndex];

            for(int jcomp=0; jcomp<nComponents; jcomp++)
            {
                int localComp = jcomp + cellIndex*nComponents;
            
                ca_Vel[localComp] = U[cellID][jcomp];
            }

            ca_P[cellIndex] = p[cellID];
            
            if (ca_T != nullptr)
                ca_T[cellIndex] = T[cellID];
            
            if (ca_Rho != nullptr)
                ca_Rho[cellIndex] = rho[cellID];

            // Turbulence data ^^^^^^^^^^^^^^^^^^
            if (ca_AlphaT != nullptr)
            {
                const tmp<volScalarField>& alphat = turbulence.alphat();
                ca_AlphaT[cellIndex] = alphat()[cellID];
            }

            if (ca_Epsilon != nullptr)
            {
                const tmp<volScalarField>& epsilon = turbulence.epsilon();
                ca_Epsilon[cellIndex] = epsilon()[cellID];
            }

            if (ca_K != nullptr)
            {
                const tmp<volScalarField>& k = turbulence.k();
                ca_K[cellIndex] = k()[cellID];
            }

            if (ca_NuT != nullptr)
            {
                const tmp<volScalarField>& nut = turbulence.nut();
                ca_NuT[cellIndex] = nut()[cellID];
            }
            //-----------------------------------

            cellIndex++;
        }
    }

    return 0;
}

int comFoam::registerVolumeData(const char *name)
{
    std::string volName = name+std::string("VOL");
    
    Info << endl
         << "rocFoam.registerVolumeData: "
         << "Registering flow data with name "
         << volName
         << endl;

    // grid and field data
    int paneID = Pstream::myProcNo()+1;// Use this paneID for volume connectivity

    Info << "procID = " << Pstream::myProcNo()
         << ", paneID = " << paneID
         << " ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;



    std::string dataName = volName+std::string(".nPoints");
    COM_new_dataitem( dataName, 'p', COM_INT, 1, "");
    COM_set_size(     dataName, paneID, 1);
    COM_set_array(    dataName, paneID, ca_nPoints);
    
    Info << "  " << dataName.c_str() << " registered." << endl;


std::cout << "HERE0 " << "ca_nPoints = " << ca_nPoints << std::endl;

std::cout << "HERE1 " << "*ca_nPoints = " << *ca_nPoints << std::endl;



    dataName = volName+std::string(".nCells");
    COM_new_dataitem( dataName, 'p', COM_INT, 1, "");
    COM_set_size(     dataName, paneID, 1);
    COM_set_array(    dataName, paneID, ca_nCells);
    Info << "  " << dataName.c_str() << " registered." << endl;



    dataName = volName+std::string(".cellToPointConn_types");
    COM_new_dataitem( dataName, 'p', COM_INT, 1, "");
    COM_set_size(     dataName, paneID, 1);
    COM_set_array(    dataName, paneID, ca_cellToPointConn_types);
    Info << "  " << dataName.c_str() << " registered." << endl;

    int ntypes = *ca_cellToPointConn_types;
    dataName = volName+std::string(".cellToPointConn_map");
    COM_new_dataitem(dataName, 'p', COM_INT, 1, "");
    COM_set_size( dataName, paneID, ntypes);
    COM_set_array(dataName, paneID, ca_cellToPointConn_map);
    Info << "  " << dataName.c_str() << " registered." << endl;

    dataName = volName+std::string(".cellToPointConn_size");
    COM_new_dataitem(dataName, 'p', COM_INT, 1, "");
    COM_set_size( dataName, paneID, ntypes);
    COM_set_array(dataName, paneID, ca_cellToPointConn_size);
    Info << "  " << dataName.c_str() << " registered." << endl;



    // points
    dataName = volName+std::string(".nc");
    COM_set_size( dataName, paneID, *ca_nPoints);
    COM_set_array(dataName, paneID, ca_Points, nComponents);
    Info << "  " << dataName.c_str() << " registered." << endl;

std::cout << "HERE0 " << "nPoints = " << *ca_nPoints << std::endl;

std::cout << "HERE1 " << "ntypes = " << ntypes << std::endl;


    // connectivity ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    for(int itype=0; itype<ntypes; itype++)
    {
        int typeID = ca_cellToPointConn_map[itype];
        int typeSize = ca_cellToPointConn_size[itype];

std::cout << "HERE2 " << "typeID = " << typeID << std::endl;
std::cout << "HERE3 " << "typeSize = " << typeSize << std::endl;

        if (typeID == 4)
        { // Tet
            dataName = volName+std::string(".:T4");
        }
        //else if (typeID == 5)
        //{ Type?
        //}
        else if (typeID == 6)
        { // Prism
            dataName = volName+std::string(".:P6");
        }
        //else if (typeID == 7)
        //{ // Type?
        //}
        else if (typeID == 8)
        { //Hex
            dataName = volName+std::string(".:H8");
        }
        else
        { // Type not identified

            Info << "=================== WARNING ==================="
                 << " Cell typeID " << typeID << " with size = "
                 << typeSize << " not identified!"
                 << endl;
            exit(-1);
        }

        COM_set_size( dataName, paneID, typeSize);
        COM_set_array(dataName, paneID,
                      ca_cellToPointConn[itype],
                      typeID
                     );
        Info << "  " << dataName.c_str() << " registered." << endl;
    }

    // Connectivity mapping stuff
    dataName = volName+std::string(".cellToCellMap");
    COM_new_dataitem( dataName, 'e', COM_INT, 1, "");
    COM_set_array(    dataName, paneID, ca_cellToCellMap, 1);
    Info << "  " << dataName.c_str() << " registered." << endl;

    dataName = volName+std::string(".cellToCellMap_inverse");
    COM_new_dataitem( dataName, 'e', COM_INT, 1, "");
    COM_set_array(    dataName, paneID, ca_cellToCellMap_inverse, 1);
    Info << "  " << dataName.c_str() << " registered." << endl;
    // ------------------------------------------

    // Element data registered with window
    dataName = volName+std::string(".vel");
    COM_new_dataitem( dataName, 'e', COM_DOUBLE, nComponents, "m/s");
    COM_set_array(    dataName, paneID, ca_Vel, nComponents);    
    Info << "  " << dataName.c_str() << " registered." << endl;

    dataName = volName+std::string(".pres");
    COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "Pa");
    COM_set_array(    dataName, paneID, ca_P, 1);
    Info << "  " << dataName.c_str() << " registered." << endl;

    if (ca_T != nullptr)
    {
        dataName = volName+std::string(".temp");
        COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "K");
        COM_set_array(    dataName, paneID, ca_T, 1);
        Info << "  " << dataName.c_str() << " registered." << endl;
    }

    if (ca_Rho != nullptr)
    {
        dataName = volName+std::string(".rho");
        COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "kg/m^3");
        COM_set_array(    dataName, paneID, ca_Rho, 1);
        Info << "  " << dataName.c_str() << " registered." << endl;
    }

    if (ca_Disp != nullptr)
    {
        dataName = volName+std::string(".disp");
        COM_new_dataitem( dataName, 'n', COM_DOUBLE, nComponents, "m");
        COM_set_array(    dataName, paneID, ca_Disp, nComponents);
        Info << "  " << dataName.c_str() << " registered." << endl;
    }

    // Turbulence data ^^^^^^^^^^^^^^^^^^^^^^^^^^
    if (ca_AlphaT != nullptr)
    {
        dataName = volName+std::string(".alphaT");
        COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "kg/m/s");
        COM_set_array(    dataName, paneID, ca_AlphaT, 1);
        Info << "  " << dataName.c_str() << " registered." << endl;
    }

    if (ca_Epsilon != nullptr)
    {
        dataName = volName+std::string(".epsilon");
        COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "m^2/s^3");
        COM_set_array(    dataName, paneID, ca_Epsilon, 1);
        Info << "  " << dataName.c_str() << " registered." << endl;
    }

    if (ca_K != nullptr)
    {
        dataName = volName+std::string(".k");
        COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "m^2/s^2");
        COM_set_array(    dataName, paneID, ca_K, 1);
        Info << "  " << dataName.c_str() << " registered." << endl;
    }

    if (ca_NuT != nullptr)
    {
        dataName = volName+std::string(".nuT");
        COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "m^2/s");
        COM_set_array(    dataName, paneID, ca_NuT, 1);
        Info << "  " << dataName.c_str() << " registered." << endl;
    }
    //-------------------------------------------

    COM_window_init_done(volName); 

    return 0;
}


int comFoam::reconstVolumeData(const char *name)
{
    std::string volName = name+std::string("VOL");

    std::cout << "rocFoam.reconstCaVolumeData, procID = "
              << Pstream::myProcNo()
              << ", Retreiving surface data form window "
              << volName << "."
              << std::endl;

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

        if (nameTmp == "nPoints" ||
            nameTmp == "nCells" ||
            nameTmp == "cellToPointConn_types" ||
            nameTmp == "cellToPointConn_map" ||
            nameTmp == "cellToPointConn_size" ||
            nameTmp == "nc" ||
            nameTmp == "cellToCellMap" ||
            nameTmp == "cellToCellMap_inverse" ||
            nameTmp == "vel" ||
            nameTmp == "pres" ||
            nameTmp == "temp" ||
            nameTmp == "rho" ||
            nameTmp == "alphaT" ||
            nameTmp == "epsilon" ||
            nameTmp == "k" ||
            nameTmp == "nuT"
            )
        {
            dataItemNames.push_back(nameTmp);
            std::cout << "  DataItem[" << i << "] = " << nameTmp << std::endl;
        }
    }
    std::cout << "  Number of items = " << dataItemNames.size()
              << std::endl << std::endl;

    // Volume data ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    int paneID = Pstream::myProcNo()+1;// Use this paneID for volume connectivity
    std::cout << "  procID = " << Pstream::myProcNo()
         << ", paneID = " << paneID
         << " ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << std::endl;

    std::string dataName = std::string("nPoints");
    nameExists(dataItemNames, dataName);
    std::string regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), paneID, &ca_nPoints);
    std::cout << "  " << dataName.c_str() << " = " << *ca_nPoints << std::endl;
    
    dataName = std::string("nCells");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), paneID, &ca_nCells);
    std::cout << "  " << dataName.c_str() << " = " << *ca_nCells << std::endl;

    dataName = std::string("cellToPointConn_types");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), paneID, &ca_cellToPointConn_types);
    std::cout << "  " << dataName.c_str() << " = " << *ca_cellToPointConn_types << std::endl;

    dataName = std::string("cellToPointConn_map");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    int nComp;
    COM_get_array(regName.c_str(), paneID, &ca_cellToPointConn_map);
    COM_get_size(regName.c_str(), paneID, &nComp);
    for(int icomp=0; icomp<nComp; icomp++)
    {
        std::cout << "  " << dataName.c_str() << "[" << icomp << "] = "
             << ca_cellToPointConn_map[icomp] << std::endl;
    }

    dataName = std::string("cellToPointConn_size");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), paneID, &ca_cellToPointConn_size);
    COM_get_size(regName.c_str(), paneID, &nComp);
    for(int icomp=0; icomp<nComp; icomp++)
    {
        std::cout << "  " << dataName.c_str() << "[" << icomp << "] = "
             << ca_cellToPointConn_size[icomp] << std::endl;
    }
    //-------------------------------------------

    // Primary allocation ^^^^^^^^^^^^^^^^^^^^^^^
    ca_cellToPointConn = new int*[*ca_cellToPointConn_types];
    //-------------------------------------------

    // Point and connectivity stuff ^^^^^^^^^
    dataName = std::string("nc");
    //nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    int nPoints;
    COM_get_array(regName.c_str(), paneID, &ca_Points, &nComp);
    COM_get_size(regName.c_str(), paneID, &nPoints);
    std::cout << "  " << dataName.c_str() << " nPoints = " << nPoints
              << ", components = " << nComp << std::endl;

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
        COM_get_array(dataName.c_str(), paneID, &ca_cellToPointConn[icon], &nComp);
        COM_get_size(dataName.c_str(), paneID, &numElem);

        std::cout << "    Connectivity[" << icon << "] = " << connName
                  << ", elements = " << numElem
                  << ", components =" << nComp << std::endl;
    }
    //---------------------------------------

    // Mapping data ^^^^^^^^^^^^^^^^^^^^^^^^^
    dataName = std::string("cellToCellMap");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), paneID, &ca_cellToCellMap, &nComp);
    COM_get_size(regName.c_str(), paneID, &numElem);
    std::cout << "    " << dataName.c_str() << " elements = " << numElem
         << ", components = " << nComp << std::endl;

    dataName = std::string("cellToCellMap_inverse");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), paneID, &ca_cellToCellMap_inverse, &nComp);
    COM_get_size(regName.c_str(), paneID, &numElem);
    std::cout << "    " << dataName.c_str() << " elements = " << numElem
              << ", components = " << nComp << std::endl;
    //---------------------------------------

    // Field data ^^^^^^^^^^^^^^^^^^^^^^^^^^^
    dataName = std::string("vel");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), paneID, &ca_Vel, &nComp);
    COM_get_size(regName.c_str(), paneID, &numElem);
    std::cout << "    " << dataName.c_str() << " elements = " << numElem
              << ", components = " << nComp << std::endl;
    
    dataName = std::string("pres");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), paneID, &ca_P, &nComp);
    COM_get_size(regName.c_str(), paneID, &numElem);
    std::cout << "    " << dataName.c_str() << " elements = " << numElem
              << ", components = " << nComp << std::endl;

    dataName = std::string("temp");
    if (nameExists(dataItemNames, dataName))
    {
        regName = volName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_T, &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);
        std::cout << "    " << dataName.c_str() << " elements = " << numElem
                  << ", components = " << nComp << std::endl;
    }

    dataName = std::string("rho");
    if (nameExists(dataItemNames, dataName))
    {
        regName = volName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_Rho, &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);
        std::cout << "    " << dataName.c_str() << " elements = " << numElem
                  << ", components = " << nComp << std::endl;
    }

    // Turbulence data ^^^^^^^^^^^^^^^^^^^^^^^^^^
    dataName = std::string("alphaT");
    if (nameExists(dataItemNames, dataName))
    {
        regName = volName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_AlphaT, &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);
        std::cout << "    " << dataName.c_str() << " elements = " << numElem
                  << ", components = " << nComp << std::endl;
    }

    dataName = std::string("epsilon");
    if (nameExists(dataItemNames, dataName))
    {
        regName = volName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_Epsilon, &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);
        std::cout << "    " << dataName.c_str() << " elements = " << numElem
                  << ", components = " << nComp << std::endl;
    }

    dataName = std::string("k");
    if (nameExists(dataItemNames, dataName))
    {
        regName = volName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_K, &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);
        std::cout << "    " << dataName.c_str() << " elements = " << numElem
                  << ", components = " << nComp << std::endl;
    }    

    dataName = std::string("nuT");
    if (nameExists(dataItemNames, dataName))
    {
        regName = volName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_NuT, &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);
        std::cout << "    " << dataName.c_str() << " elements = " << numElem
                  << ", components = " << nComp << std::endl;
    }
    //-------------------------------------------
    
    std::cout << "  --------------------------------------------------"
         << std::endl;

    return 0;
}

int comFoam::deleteVolumeData()
{
    if (ca_cellToPointConn_map != nullptr)
    {
        delete [] ca_cellToPointConn_map;
        ca_cellToPointConn_map = nullptr;
    }
    
    if (ca_cellToPointConn_size != nullptr)
    {
        delete [] ca_cellToPointConn_size;
        ca_cellToPointConn_size = nullptr;
    }

    if (ca_cellToCellMap != nullptr)
    {
        delete [] ca_cellToCellMap;
        ca_cellToCellMap = nullptr;
    }

    if (ca_cellToCellMap_inverse != nullptr)
    {
        delete [] ca_cellToCellMap_inverse;
        ca_cellToCellMap_inverse = nullptr;
    }
   
    if (ca_cellToPointConn_types != nullptr)
    {
        int ntype = *ca_cellToPointConn_types;
        if (ca_cellToPointConn != nullptr)
        {
            for (int itype=0; itype<ntype; itype++)
            {
                if (ca_cellToPointConn[itype] != nullptr)
                {
                    delete [] ca_cellToPointConn[itype];
                    ca_cellToPointConn[itype] = nullptr;
                }
            }
            delete [] ca_cellToPointConn;
            ca_cellToPointConn = nullptr;
        }

        delete[] ca_cellToPointConn_types;
        ca_cellToPointConn_types = nullptr;
    }

    if (ca_Points != nullptr)
    {
        delete [] ca_Points;
        ca_Points = nullptr;
    }

    if (ca_Disp != nullptr)
    {
        delete [] ca_Disp;
        ca_Disp = nullptr;
    }

    if (ca_Vel != nullptr)
    {
        delete [] ca_Vel;
        ca_Vel = nullptr;
    }

    if (ca_P != nullptr)
    {
        delete [] ca_P;
        ca_P = nullptr;
    }

    if (ca_T != nullptr)
    {
        delete [] ca_T;
        ca_T = nullptr;
    }

    if (ca_Rho != nullptr)
    {
        delete [] ca_Rho;
        ca_Rho = nullptr;
    }

    if (ca_nPoints!= nullptr)
    {
        delete ca_nPoints;
        ca_nPoints = nullptr;
    }

    if (ca_nCells != nullptr)
    {
        delete ca_nCells;
        ca_nCells = nullptr;
    }

    // Turbulence data ^^^^^^^^^^^^^^^^^^^^^^^^^^
    if (ca_AlphaT != nullptr)
    {
        delete [] ca_AlphaT;
        ca_AlphaT = nullptr;
    }

    if (ca_Epsilon != nullptr)
    {
        delete [] ca_Epsilon;
        ca_Epsilon = nullptr;
    }

    if (ca_K != nullptr)
    {
        delete [] ca_K;
        ca_K = nullptr;
    }

    if (ca_NuT != nullptr)
    {
        delete [] ca_NuT;
        ca_NuT = nullptr;
    }
    //-------------------------------------------
    
    return 0;
}


