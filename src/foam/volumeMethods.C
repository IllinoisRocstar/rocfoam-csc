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
    for(auto it=mapCellToCellMap.begin(); it!=mapCellToCellMap.end(); it++)
    {
        const auto& vecCells = it->second;
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
    for(auto it=mapCellToCellMap.begin(); it!=mapCellToCellMap.end(); it++)
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
    ca_Points = new double[nTotal];

    // Field-data
    nTotal = *ca_nCells * nComponents;
    ca_Vel = new double[nTotal];

    ca_P   = new double[*ca_nCells];
    
    if (TPtr != NULL)
        ca_T   = new double[*ca_nCells];

    if (rhoPtr != NULL)
        ca_Rho = new double[*ca_nCells];

    return 0;
}

int comFoam::updateVolumeData()
{
    // Point data ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    const dynamicFvMesh& mesh(*meshPtr);
    const pointField&    points = mesh.points();

    forAll(points, ipoint)
    {
        for(int jcomp=0; jcomp<nComponents; jcomp++)
        {
            ca_Points[ipoint*nComponents+jcomp]
                = points[ipoint][jcomp];
        }
    }
    
    // Cell-centered data ^^^^^^^^^^^^^^^^^^^^^^^
    const volScalarField& p(*pPtr);
    const volVectorField& U(*UPtr);
    const volScalarField& T(*TPtr);
    const volScalarField& rho(*rhoPtr);

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
            
            if (ca_T != NULL)
                ca_T[cellIndex] = T[cellID];
            
            if (ca_Rho != NULL)
                ca_Rho[cellIndex] = rho[cellID];

            cellIndex++;
        }
    }

    return 0;
}

int comFoam::registerVolumeData(const char *name)
{
    Info << "rocFoam.registerVolumeData: "
         << "Registering flow data with name "
         << name
         << endl;

    std::string volName = name+std::string("VOL");

    // Genral data registered with window
    // Solver data    
    std::string dataName = volName+std::string(".time");
    COM_new_dataitem( dataName, 'w', COM_DOUBLE, 1, "");
    COM_set_size(     dataName, 0, 1);
    COM_set_array(    dataName, 0, ca_time);
    Info << dataName << " registered." << endl;

    dataName = volName+std::string(".timeName");
    COM_new_dataitem( dataName, 'w', COM_CHAR, 1, "");
    COM_set_size(     dataName, 0, genCharSize);
    COM_set_array(    dataName, 0, ca_timeName);
    Info << dataName << " registered." << endl;

    dataName = volName+std::string(".deltaT");
    COM_new_dataitem( dataName, 'w', COM_DOUBLE, 1, "");
    COM_set_size(     dataName, 0, 1);
    COM_set_array(    dataName, 0, ca_deltaT);
    Info << dataName << " registered." << endl;

    dataName = volName+std::string(".deltaT0");
    COM_new_dataitem( dataName, 'w', COM_DOUBLE, 1, "");
    COM_set_size(     dataName, 0, 1);
    COM_set_array(    dataName, 0, ca_deltaT0);
    Info << dataName << " registered." << endl;

    dataName = volName+std::string(".timeIndex");
    COM_new_dataitem( dataName, 'w', COM_INT, 1, "");
    COM_set_size(     dataName, 0, 1);
    COM_set_array(    dataName, 0, ca_timeIndex);
    Info << dataName << " registered." << endl;

    dataName = volName+std::string(".runStat");
    COM_new_dataitem( dataName, 'w', COM_INT, 1, "");
    COM_set_size(     dataName, 0, 1);
    COM_set_array(    dataName, 0, ca_runStat);
    Info << dataName << " registered." << endl;

    // grid and field data
    int paneID = Pstream::myProcNo()+1;// Use this paneID for volume connectivity

    Info << "procID = " << Pstream::myProcNo()
         << ", paneID = " << paneID
         << " ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;

    dataName = volName+std::string(".nPoints");
    COM_new_dataitem( dataName, 'p', COM_INT, 1, "");
    COM_set_size(     dataName, paneID, 1);
    COM_set_array(    dataName, paneID, ca_nPoints);
    Info << "  " << dataName.c_str() << " registered." << endl;

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

    // connectivity ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    for(int itype=0; itype<ntypes; itype++)
    {
        int typeID = ca_cellToPointConn_map[itype];
        int typeSize = ca_cellToPointConn_size[itype];

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
            return -1;
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

    if (ca_T != NULL)
    {
        dataName = volName+std::string(".temp");
        COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "K");
        COM_set_array(    dataName, paneID, ca_T, 1);
        Info << "  " << dataName.c_str() << " registered." << endl;
    }

    if (ca_Rho != NULL)
    {
        dataName = volName+std::string(".rho");
        COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "kg/m^3");
        COM_set_array(    dataName, paneID, ca_Rho, 1);
        Info << "  " << dataName.c_str() << " registered." << endl;
    }

    COM_window_init_done(volName); 

    return 0;
}


int comFoam::reconstCaVolumeData(const char *name)
{
    std::string volName = name+std::string("VOL");

    std::cout << "rocFoam.reconstCaVolumeData, proID = "
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

        std::string subName = nameTmp.substr(0,4);
        if (subName != "file" && subName != "nFil")
        {
            dataItemNames.push_back(nameTmp);
            std::cout << "  DataItem[" << i << "] = " << nameTmp << std::endl;
        }
    }
    std::cout << std::endl;

    // Flow stat data ^^^^^^^^^^^^^^^^^^^^^^^^^^^
    std::string dataName = std::string("time");
    nameExists(dataItemNames, dataName);
    std::string regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), 0, &ca_time);
    std::cout << "  " << dataName.c_str() << " = " << *ca_time << std::endl;

    dataName = std::string("timeIndex");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), 0, &ca_timeIndex);
    std::cout << "  " << dataName.c_str() << " = " << *ca_timeIndex << std::endl;

    dataName = std::string("deltaT");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), 0, &ca_deltaT);
    std::cout << "  " << dataName.c_str() << " = " << *ca_deltaT << std::endl;

    dataName = std::string("deltaT0");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), 0, &ca_deltaT0);
    std::cout << "  " << dataName.c_str() << " = " << *ca_deltaT0 << std::endl;

    int nComp;
    dataName = std::string("timeName");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), 0, &ca_timeName);
    COM_get_size(regName.c_str(), 0, &nComp);
    std::cout << "  " << dataName.c_str() << " = " << ca_timeName << std::endl;

    dataName = std::string("runStat");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), 0, &ca_runStat);
    std::cout << "  " << dataName.c_str() << " = " << *ca_runStat << std::endl;
    //-------------------------------------------
    
    // Volume data ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    int paneID = Pstream::myProcNo()+1;// Use this paneID for volume connectivity
    std::cout << "  procID = " << Pstream::myProcNo()
         << ", paneID = " << paneID
         << " ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << std::endl;

    dataName = std::string("nPoints");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
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
/*    std::cout << "    " << dataName.c_str() << " points = " << *ca_nPoints*/
/*         << ", components = " << nComp << std::endl;*/
/*    for(int ipoint=0; ipoint<*ca_nPoints; ipoint++)*/
/*    {*/
/*        std::cout << "Node " << ipoint << " ca_Points = ";*/
/*        for(int icomp=0; icomp<nComp; icomp++)*/
/*        {*/
/*          std::cout << *(ca_Points+ipoint*nComp+icomp) << " ";*/
/*        }*/
/*        std::cout << std::endl;*/
/*    }*/

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
    dataName = std::string("cellToCellMap");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    
    COM_get_array(regName.c_str(), paneID, &ca_cellToCellMap, &nComp);
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
/*    for(int icell=0; icell<numCells; icell++)*/
/*    {*/
/*        std::cout << "Cell " << icell << " velocity = ";*/
/*        for(int icomp=0; icomp<nComp; icomp++)*/
/*        {*/
/*            std::cout << *(cellVel+icell*nComp+icomp) << " ";*/
/*        }*/
/*        std::cout << std::endl;*/
/*    }*/
    
    dataName = std::string("pres");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;

    COM_get_array(regName.c_str(), paneID, &ca_P, &nComp);
    COM_get_size(regName.c_str(), paneID, &numElem);
    std::cout << "    " << dataName.c_str() << " elements = " << numElem
         << ", components = " << nComp << std::endl;
/*    for(int icell=0; icell<numCells; icell++)*/
/*    {*/
/*        std::cout << "Cell " << icell << " pressure = ";*/
/*        for(int icomp=0; icomp<nComp; icomp++)*/
/*        {*/
/*            std::cout << *(cellPres+icell*nComp+icomp) << " ";*/
/*        }*/
/*        std::cout << std::endl;*/
/*    }*/

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

    std::cout << "  --------------------------------------------------"
         << std::endl;

    return 0;
}



int comFoam::deleteVolumeData()
{
    if (ca_cellToPointConn_map != NULL)
    {
        delete [] ca_cellToPointConn_map;
        ca_cellToPointConn_map = NULL;
    }
    
    if (ca_cellToPointConn_size != NULL)
    {
        delete [] ca_cellToPointConn_size;
        ca_cellToPointConn_size = NULL;
    }

    if (ca_cellToCellMap != NULL)
    {
        delete [] ca_cellToCellMap;
        ca_cellToCellMap = NULL;
    }

    if (ca_cellToCellMap_inverse != NULL)
    {
        delete [] ca_cellToCellMap_inverse;
        ca_cellToCellMap_inverse = NULL;
    }
   
    if (ca_cellToPointConn_types != NULL)
    {
        int ntype = *ca_cellToPointConn_types;
        if (ca_cellToPointConn != NULL)
        {
            for (int itype=0; itype<ntype; itype++)
            {
                if (ca_cellToPointConn[itype] != NULL)
                {
                    delete [] ca_cellToPointConn[itype];
                    ca_cellToPointConn[itype] = NULL;
                }
            }
            delete [] ca_cellToPointConn;
            ca_cellToPointConn = NULL;
        }

        delete[] ca_cellToPointConn_types;
        ca_cellToPointConn_types = NULL;
    }

    
    if (ca_Points != NULL)
    {
        delete[] ca_Points;
        ca_Points = NULL;
    }

    if (ca_Vel != NULL)
    {
        delete[] ca_Vel;
        ca_Vel = NULL;
    }

    if (ca_P != NULL)
    {
        delete[] ca_P;
        ca_P = NULL;
    }

    if (ca_T != NULL)
    {
        delete[] ca_T;
        ca_T = NULL;
    }

    if (ca_Rho != NULL)
    {
        delete[] ca_Rho;
        ca_Rho = NULL;
    }

    if (ca_nPoints!= NULL)
    {
        delete[] ca_nPoints;
        ca_nPoints = NULL;
    }

    if (ca_nCells != NULL)
    {
        delete[] ca_nCells;
        ca_nCells = NULL;
    }

    return 0;
}


