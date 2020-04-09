int comFoam::createVolumeConnectivities()
{

    const dynamicFvMesh& mesh(*meshPtr);

    // Mesh and conmnectivities ^^^^^^^^^^^^^^^^^
    const labelListList& cellPoints = mesh.cellPoints();
    const cellShapeList& cellShapes = mesh.cellShapes();
    //-------------------------------------------

    ca_nPoints = new int(mesh.nPoints());
    ca_nCells  = new int(mesh.nCells());
    // Temporary Vectors ^^^^^^^^^^^^^^^^^^^^^^^^
    std::vector<int> vecTmpInt;

    std::vector< std::vector<int> > vecCellToCellMap;
    //-------------------------------------------

    // cellToPoint connectivities ^^^^^^^^^^^^^^^
    vecCellToCellMap.clear();
    vecCellToCellMap.resize(cellToPointTypeSize);
    for (int i=0; i<cellToPointTypeSize; i++)
    {
        vecCellToCellMap[i].clear();
    }

    forAll(cellPoints, icell)
    {
        const labelList& pointsList = cellPoints[icell];
        int nPointsInCell = pointsList.size();
        vecCellToCellMap[nPointsInCell-1].push_back(icell);
    }
    //-------------------------------------------

    // Create cellToCell mapping ^^^^^^^^^^^^^^^^
    ca_cellToCellMap = new int[*ca_nCells];

    int sortedCellIndex = 0;
    int typeCount = 0;
    int totalnCellTypes = vecCellToCellMap.size();
    for(int itype=0; itype<totalnCellTypes; itype++)
    {
        int nCells = vecCellToCellMap[itype].size();
        if (nCells > 0)
        {
            for(int icell=0; icell<nCells; icell++)
            {
                ca_cellToCellMap[sortedCellIndex] =
                    vecCellToCellMap[itype][icell];

                sortedCellIndex++;
            }
            typeCount++;
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
    ca_cellToPointConn_types = new int(typeCount);
    ca_cellToPointConn_map   = new int[typeCount];
    ca_cellToPointConn_size  = new int[typeCount];
    ca_cellToPointConn = new int*[typeCount];

    typeCount = 0;
    sortedCellIndex = 0;
    for(int itype=0; itype<totalnCellTypes; itype++)
    {
        int nCells = vecCellToCellMap[itype].size();
        if (nCells > 0)
        {
            int nPoints = itype+1;
            ca_cellToPointConn_map[typeCount] = nPoints;
            ca_cellToPointConn_size[typeCount] = nCells;

            int nTypeConn = nCells * nPoints;
            ca_cellToPointConn[typeCount] = new int[nTypeConn];

            for (int icell=0; icell<nCells; icell++)
            {
                int cellID = ca_cellToCellMap[sortedCellIndex];
                const cellShape& cellShape_ = cellShapes[cellID];

                for (int ipoint=0; ipoint<nPoints; ipoint++)
                {
                    int index = icell*nPoints + ipoint;

                    ca_cellToPointConn[typeCount][index] = cellShape_[ipoint];
                }
                sortedCellIndex++;
            }
            typeCount++;
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
    ca_T   = new double[*ca_nCells];
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
            ca_T[cellIndex] = T[cellID];
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

    std::string volName = name+string("VOL");

    // Genral data registered with window
    // Solver data    
    std::string dataName = volName+string(".time");
    COM_new_dataitem( dataName, 'w', COM_DOUBLE, 1, "");
    COM_set_size(     dataName, 0, 1);
    COM_set_array(    dataName, 0, ca_time);
    Info << dataName << " registered." << endl;

    dataName = volName+string(".deltaT");
    COM_new_dataitem( dataName, 'w', COM_DOUBLE, 1, "");
    COM_set_size(     dataName, 0, 1);
    COM_set_array(    dataName, 0, ca_deltaT );
    Info << dataName << " registered." << endl;

    dataName = volName+string(".runStat");
    COM_new_dataitem( dataName, 'w', COM_INT, 1, "");
    COM_set_size(     dataName, 0, 1);
    COM_set_array(    dataName, 0, ca_runStat);
    Info << dataName << " registered." << endl;

    // grid and field data
    dataName = volName+string(".nPoints");
    COM_new_dataitem( dataName, 'w', COM_INT, 1, "");
    COM_set_size(     dataName, 0, 1);
    COM_set_array(    dataName, 0, ca_nPoints);
    Info << dataName << " registered." << endl;

    dataName = volName+string(".nCells");
    COM_new_dataitem( dataName, 'w', COM_INT, 1, "");
    COM_set_size(     dataName, 0, 1);
    COM_set_array(    dataName, 0, ca_nCells);
    Info << dataName << " registered." << endl;

    dataName = volName+string(".cellToPointConn_types");
    COM_new_dataitem( dataName, 'w', COM_INT, 1, "");
    COM_set_size(     dataName, 0, 1);
    COM_set_array(    dataName, 0, ca_cellToPointConn_types);
    Info << dataName << " registered." << endl;

    int ntypes = *ca_cellToPointConn_types;
    dataName = volName+string(".cellToPointConn_map");
    COM_new_dataitem(dataName, 'w', COM_INT, 1, "");
    COM_set_size( dataName, 0, ntypes);
    COM_set_array(dataName, 0, ca_cellToPointConn_map);
    Info << dataName << " registered." << endl;

    dataName = volName+string(".cellToPointConn_size");
    COM_new_dataitem(dataName, 'w', COM_INT, 1, "");
    COM_set_size( dataName, 0, ntypes);
    COM_set_array(dataName, 0, ca_cellToPointConn_size);
    Info << dataName << " registered." << endl;

    int paneID = 1; // Use this paneID for volume connectivity
    // points
    dataName = volName+string(".nc");
    COM_set_size( dataName, paneID, *ca_nPoints);
    COM_set_array(dataName, paneID, ca_Points, nComponents);
    Info << dataName << " registered." << endl;

    // connectivity ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    for(int itype=0; itype<ntypes; itype++)
    {
        int typeID = ca_cellToPointConn_map[itype];
        int typeSize = ca_cellToPointConn_size[itype];

        if (typeID == 4)
        { // Tet
            dataName = volName+string(".:T4");
        }
        //else if (typeID == 5)
        //{ Type?
        //}
        else if (typeID == 6)
        { // Prism
            dataName = volName+string(".:P4");
        }
        //else if (typeID == 7)
        //{ // Type?
        //}
        else if (typeID == 8)
        { //Hex
            dataName = volName+string(".:H8");
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
        COM_set_array( dataName, paneID,
                       ca_cellToPointConn[itype],
                       typeID
                     );
        Info << dataName << " registered." << endl;
    }
        
    // Connectivity mapping stuff
    dataName = volName+string(".cellToCellMap");
    COM_new_dataitem( dataName, 'e', COM_INT, 1, "");
    COM_set_array(    dataName, paneID, ca_cellToCellMap, 1);
    Info << dataName << " registered." << endl;
    // ------------------------------------------

    // Element data registered with window
    dataName = volName+string(".vel");
    COM_new_dataitem( dataName, 'e', COM_DOUBLE, nComponents, "m/s");
    COM_set_array(    dataName, paneID, ca_Vel, nComponents);    
    Info << dataName << " registered." << endl;

    dataName = volName+string(".pres");
    COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "Pa");
    COM_set_array(    dataName, paneID, ca_P, 1);
    Info << dataName << " registered." << endl;

    dataName = volName+string(".temp");
    COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "K");
    COM_set_array(    dataName, paneID, ca_T, 1);
    Info << dataName << " registered." << endl;

    dataName = volName+string(".rho");
    COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "kg/m^3");
    COM_set_array(    dataName, paneID, ca_Rho, 1);
    Info << dataName << " registered." << endl;

    COM_window_init_done(volName); 

    return 0;
}


int comFoam::reconstCaVolumeData(const char *name)
{

    std::string volName = name+string("VOL");
    Info << "rocFoam.main: Retreiving volume data form window "
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

    // Flow stat data ^^^^^^^^^^^^^^^^^^^^^^^^^^^
    std::string dataName = string("time");
    nameExists(dataItemNames, dataName);
    std::string regName = volName+string(".")+dataName;
    COM_get_array(regName.c_str(), 0, &ca_time);
    Info << "  " << dataName.c_str() << " = " << *ca_time << endl;

    dataName = string("deltaT");
    nameExists(dataItemNames, dataName);
    regName = volName+string(".")+dataName;
    COM_get_array(regName.c_str(), 0, &ca_deltaT);
    std::cout << "  " << dataName.c_str() << " = " << *ca_deltaT << std::endl;

    dataName = string("runStat");
    nameExists(dataItemNames, dataName);
    regName = volName+string(".")+dataName;
    COM_get_array(regName.c_str(), 0, &ca_runStat);
    Info << "  " << dataName.c_str() << " = " << *ca_runStat << endl;
    //-------------------------------------------
    
    // Volume data ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    dataName = string("nPoints");
    nameExists(dataItemNames, dataName);
    regName = volName+string(".")+dataName;
    COM_get_array(regName.c_str(), 0, &ca_nPoints);
    Info << "  " << dataName.c_str() << " = " << *ca_nPoints << endl;
    
    dataName = string("nCells");
    nameExists(dataItemNames, dataName);
    regName = volName+string(".")+dataName;
    COM_get_array(regName.c_str(), 0, &ca_nCells);
    Info << "  " << dataName.c_str() << " = " << *ca_nCells << endl;

    dataName = string("cellToPointConn_types");
    nameExists(dataItemNames, dataName);
    regName = volName+string(".")+dataName;
    COM_get_array(regName.c_str(), 0, &ca_cellToPointConn_types);
    Info << "  " << dataName.c_str() << " = " << *ca_cellToPointConn_types << endl;

    dataName = string("cellToPointConn_map");
    nameExists(dataItemNames, dataName);
    regName = volName+string(".")+dataName;
    int nComp;
    COM_get_array(regName.c_str(), 0, &ca_cellToPointConn_map);
    COM_get_size(regName.c_str(), 0, &nComp);
    for(int icomp=0; icomp<nComp; icomp++)
    {
        Info << "  " << dataName.c_str() << "[" << icomp << "] = "
             << ca_cellToPointConn_map[icomp] << endl;
    }

    dataName = string("cellToPointConn_size");
    nameExists(dataItemNames, dataName);
    regName = volName+string(".")+dataName;
    COM_get_array(regName.c_str(), 0, &ca_cellToPointConn_size);
    COM_get_size(regName.c_str(), 0, &nComp);
    for(int icomp=0; icomp<nComp; icomp++)
    {
        Info << "  " << dataName.c_str() << "[" << icomp << "] = "
             << ca_cellToPointConn_size[icomp] << endl;
    }
    //-------------------------------------------

    // Primary allocation ^^^^^^^^^^^^^^^^^^^^^^^
    ca_cellToPointConn = new int*[*ca_cellToPointConn_types];
    //-------------------------------------------

    int paneID = 1;
    int ipane = 0;

    Info << "  Pane[" << ipane << "], paneID = " << paneID
         << " ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;

    // Point and connectivity stuff ^^^^^^^^^
    dataName = string("nc");
    //nameExists(dataItemNames, dataName);
    regName = volName+string(".")+dataName;
    int nPoints;
    
    COM_get_array(regName.c_str(), paneID, &ca_Points, &nComp);
    COM_get_size(regName.c_str(), paneID, &nPoints);
/*    Info << "    " << dataName.c_str() << " points = " << *ca_nPoints*/
/*         << ", components = " << nComp << endl;*/
/*    for(int ipoint=0; ipoint<*ca_nPoints; ipoint++)*/
/*    {*/
/*        Info << "Node " << ipoint << " ca_Points = ";*/
/*        for(int icomp=0; icomp<nComp; icomp++)*/
/*        {*/
/*          Info << *(ca_Points+ipoint*nComp+icomp) << " ";*/
/*        }*/
/*        Info << endl;*/
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

        dataName = volName+string(".")+connName;
        //nameExists(dataItemNames, dataName);
        COM_get_array(dataName.c_str(), paneID, &ca_cellToPointConn[icon], &nComp);
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
    dataName = string("cellToCellMap");
    nameExists(dataItemNames, dataName);
    regName = volName+string(".")+dataName;
    
    COM_get_array(regName.c_str(), paneID, &ca_cellToCellMap, &nComp);
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

    // Field data ^^^^^^^^^^^^^^^^^^^^^^^^^^^
    dataName = string("vel");
    nameExists(dataItemNames, dataName);
    regName = volName+string(".")+dataName;

    COM_get_array(regName.c_str(), paneID, &ca_Vel, &nComp);
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
    
    dataName = string("pres");
    nameExists(dataItemNames, dataName);
    regName = volName+string(".")+dataName;

    COM_get_array(regName.c_str(), paneID, &ca_P, &nComp);
    COM_get_size(regName.c_str(), paneID, &numElem);
    Info << "    " << dataName.c_str() << " elements = " << numElem
         << ", components = " << nComp << endl;
/*    for(int icell=0; icell<numCells; icell++)*/
/*    {*/
/*        Info << "Cell " << icell << " pressure = ";*/
/*        for(int icomp=0; icomp<nComp; icomp++)*/
/*        {*/
/*            Info << *(cellPres+icell*nComp+icomp) << " ";*/
/*        }*/
/*        Info << endl;*/
/*    }*/

    dataName = string("temp");
    nameExists(dataItemNames, dataName);
    regName = volName+string(".")+dataName;

    COM_get_array(regName.c_str(), paneID, &ca_T, &nComp);
    COM_get_size(regName.c_str(), paneID, &numElem);
    Info << "    " << dataName.c_str() << " elements = " << numElem
         << ", components = " << nComp << endl;
/*    for(int icell=0; icell<numCells; icell++)*/
/*    {*/
/*        Info << "Cell " << icell << " temperature = ";*/
/*        for(int icomp=0; icomp<nComp; icomp++)*/
/*        {*/
/*            Info << *(cellTemp+icell*nComp+icomp) << " ";*/
/*        }*/
/*        Info << endl;*/
/*    }*/

    dataName = string("rho");
    nameExists(dataItemNames, dataName);
    regName = volName+string(".")+dataName;

    COM_get_array(regName.c_str(), paneID, &ca_Rho, &nComp);
    COM_get_size(regName.c_str(), paneID, &numElem);
    Info << "    " << dataName.c_str() << " elements = " << numElem
         << ", components = " << nComp << endl;
/*    for(int icell=0; icell<numCells; icell++)*/
/*    {*/
/*        Info << "Cell " << icell << " density = ";*/
/*        for(int icomp=0; icomp<nComp; icomp++)*/
/*        {*/
/*            Info << *(cellRho+icell*nComp+icomp) << " ";*/
/*        }*/
/*        Info << endl;*/
/*    }*/

    Info << "  --------------------------------------------------"
         << endl;

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
   
    if (ca_cellToPointConn != NULL)
    {
        for (int itype=0; itype<*ca_cellToPointConn_types; itype++)
        {
            if (ca_cellToPointConn[itype] != NULL)
            {
                delete [] ca_cellToPointConn[itype];
            }
        }
        ca_cellToPointConn = NULL;
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

    if (ca_cellToPointConn_types != NULL)
    {
        delete[] ca_cellToPointConn_types;
        ca_cellToPointConn_types = NULL;
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


