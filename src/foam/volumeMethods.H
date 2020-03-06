int comFoam::createVolumeConnectivities()
{

    const dynamicFvMesh& mesh(*meshPtr);

    // Mesh and conmnectivities ^^^^^^^^^^^^^^^^^
    const labelListList& cellPoints = mesh.cellPoints();
    const cellShapeList& cellShapes = mesh.cellShapes();
    //-------------------------------------------

    ca_nPoints = mesh.nPoints();
    ca_nCells  = mesh.nCells();

    // Temporary Vectors ^^^^^^^^^^^^^^^^^^^^^^^^
    std::vector<int> vecTmpInt;

    std::vector< std::vector<int> > vecCellToCellMap_localToGlobal;
    //-------------------------------------------

    // cellToPoint connectivities ^^^^^^^^^^^^^^^
    vecCellToCellMap_localToGlobal.clear();
    vecCellToCellMap_localToGlobal.resize(cellToPointTypeSize);
    for (int i=0; i<cellToPointTypeSize; i++)
    {
        vecCellToCellMap_localToGlobal[i].clear();
    }

    forAll(cellPoints, icell)
    {
        const labelList& pointsList = cellPoints[icell];
        int nPointsInCell = pointsList.size();
        vecCellToCellMap_localToGlobal[nPointsInCell-1].push_back(icell);
    }
    //-------------------------------------------

    // Create local-to-global cell map according to types
    int totalnCellTypes = vecCellToCellMap_localToGlobal.size();

    for(int itype=0; itype<totalnCellTypes; itype++)
    {
        if (vecCellToCellMap_localToGlobal[itype].size()>0)
        {
            ca_cellToPointConn_types++;
        }
    }

    ca_cellToPointConn_map  = new int[ca_cellToPointConn_types];
    ca_cellToPointConn_size = new int[ca_cellToPointConn_types];
    ca_cellToCellMap_localToGlobal = new int*[ca_cellToPointConn_types];
    
    int checkSum = 0;
    ca_cellToPointConn_types = 0;
    for(int itype=0; itype<totalnCellTypes; itype++)
    {
        if (vecCellToCellMap_localToGlobal[itype].size()>0)
        {
            ca_cellToPointConn_map[ca_cellToPointConn_types] = itype+1;
            
            int typeSize = vecCellToCellMap_localToGlobal[itype].size();
            ca_cellToPointConn_size[ca_cellToPointConn_types] = typeSize;

            ca_cellToCellMap_localToGlobal[ca_cellToPointConn_types] = new int[typeSize];

            for(int icell=0; icell<typeSize; icell++)
            {
                ca_cellToCellMap_localToGlobal[ca_cellToPointConn_types][icell] =
                    vecCellToCellMap_localToGlobal[itype][icell];
            
                checkSum++;
            }

            ca_cellToPointConn_types++;
        }
    }

    if (checkSum != ca_nCells)
    {
        Foam::Info << "========== WARNNING ===============" << endl
                   << "     checkSum != ca_nCells " << endl
                   << "    " << checkSum << "!=" << ca_nCells << endl;
        return -1;
    }

    ca_celltoPointConn = new int*[ca_cellToPointConn_types];
    for(int itype=0; itype<ca_cellToPointConn_types; itype++)
    {
        int nPoints_ = ca_cellToPointConn_map[itype];
        int nCells_  = ca_cellToPointConn_size[itype];
        
        int nTypeConn = nCells_ * nPoints_;
        ca_celltoPointConn[itype] = new int[nTypeConn];

        for (int icell=0; icell<nCells_; icell++)
        {
            int cellID = ca_cellToCellMap_localToGlobal[itype][icell];
            const cellShape& cellShape_ = cellShapes[cellID];

            for (int ipoint=0; ipoint<nPoints_; ipoint++)
            {
                int index = icell*nPoints_ + ipoint;
                ca_celltoPointConn[itype][index] = cellShape_[ipoint];
            }
        }
    }
    
    //-------------------------------------------

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
   
    if (ca_cellToCellMap_localToGlobal != NULL)
    {
        for (int itype=0; itype<ca_cellToPointConn_types; itype++)
        {
            if (ca_cellToCellMap_localToGlobal[itype] != NULL)
            {
                delete [] ca_cellToCellMap_localToGlobal[itype];
            }
        }
        delete [] ca_cellToCellMap_localToGlobal;
        ca_cellToCellMap_localToGlobal = NULL;
    }

    if (ca_celltoPointConn != NULL)
    {
        for (int itype=0; itype<ca_cellToPointConn_types; itype++)
        {
            if (ca_celltoPointConn[itype] != NULL)
            {
                delete [] ca_celltoPointConn[itype];
            }
        }
        ca_celltoPointConn = NULL;
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

    return 0;
}


int comFoam::createVolumeData()
{
    int nTotal = ca_nPoints * nComponents;
    ca_Points = new double[nTotal];

    nTotal = ca_nCells * nComponents;
    ca_Vel = new double[nTotal];

    ca_P = new double[ca_nCells];
    ca_T = new double[ca_nCells];
    ca_Rho = new double[ca_nCells];
    
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
    for(int itype=0; itype<ca_cellToPointConn_types; itype++)
    {
        int ncells = ca_cellToPointConn_size[itype];
        for(int icell=0; icell<ncells; icell++)
        {
            int cellID = ca_cellToCellMap_localToGlobal[itype][icell];

            for(int jcomp=0; jcomp<nComponents; jcomp++)
            {
                int localComp = jcomp + cellIndex*nComponents;
            
                ca_Vel[localComp] = U[cellID][jcomp];
            }

            ca_P[cellIndex] = p[cellIndex];
            ca_T[cellIndex] = T[cellIndex];
            ca_Rho[cellIndex] = rho[cellIndex];
            
            cellIndex++;
        }
    }

    return 0;
}

int comFoam::registerVolumeData(const char *name)
{
    Foam::Info << "rocFoam.registerVolumeData: "
               << "Registering flow data with name "
               << name
               << endl;

    std::string volName = name+string("VOL");

    // Genral data registered with window
    std::string dataName = volName+string(".winNProc");
    COM_new_dataitem( dataName, 'w', COM_INT, 1, "");
    COM_set_size(     dataName, 0, 1);
    COM_set_array(    dataName, 0, &(winNProc));

    dataName = volName+string(".winTime");
    COM_new_dataitem( dataName, 'w', COM_DOUBLE, 1, "");
    COM_set_size(     dataName, 0, 1);
    COM_set_array(    dataName, 0, &(winTime));

    dataName = volName+string(".winDeltaT");
    COM_new_dataitem( dataName, 'w', COM_DOUBLE, 1, "");
    COM_set_size(     dataName, 0, 1);
    COM_set_array(    dataName, 0, &(winDeltaT) );

    dataName = volName+string(".winRun");
    COM_new_dataitem( dataName, 'w', COM_INT, 1, "");
    COM_set_size(     dataName, 0, 1);
    COM_set_array(    dataName, 0, &(winRun));

    // points
    int paneID = 1; // Use this ID for surface pane 
                    // in volume window
    dataName = volName+string(".nc");
    COM_set_size( dataName, paneID, ca_nPoints);
    COM_set_array(dataName, paneID, ca_Points, nComponents);

    // connectivity
    for(int itype=0; itype<ca_cellToPointConn_types; itype++)
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

            Foam::Info << "=================== WARNING ==================="
                       << " Cell typeID " << typeID << " with size = "
                       << typeSize << " not identified!"
                       << endl;
            return -1;
        }

        COM_set_size( dataName, paneID, typeSize);
        COM_set_array( dataName, paneID,
                       ca_celltoPointConn[itype],
                       typeID
                     );
    }
        
    // Element data registered with window
    dataName = volName+string(".vel");
    COM_new_dataitem( dataName, 'e', COM_DOUBLE, nComponents, "m/s");
    COM_set_array(    dataName, paneID, ca_Vel, nComponents);    

    dataName = volName+string(".pres");
    COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "Pa");
    COM_set_array(    dataName, paneID, ca_P, 1);

    dataName = volName+string(".temp");
    COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "K");
    COM_set_array(    dataName, paneID, ca_T, 1);

    dataName = volName+string(".rho");
    COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "kg/m^3");
    COM_set_array(    dataName, paneID, ca_Rho, 1);

    COM_window_init_done(volName); 

    return 0;
}

