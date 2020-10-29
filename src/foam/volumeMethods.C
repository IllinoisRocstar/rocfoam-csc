#include "comFoam.H"

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
        std::stringstream output{};
        output << "========== WARNNING ===============" << endl
             << "     sortedCellIndex != ca_nCells " << endl
             << "    " << sortedCellIndex << "!=" << *ca_nCells;
        message(output.str(), true);

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

                ca_cellToPointConn[itype][index] = cellShape_[ipoint] + 1; //CGNS starts from ID 1;
            }
            sortedCellIndex++;
        }
    }
    //-------------------------------------------

    return 0;
}

int comFoam::createVolumeData()
{
    int npoints = *ca_nPoints;
    int nTotal  = npoints * nComponents;
    ca_Points = new double[nTotal]{0};

    std::string dynamicSolverType = ca_dynamicSolverType;
    if
    (
        *ca_isDynamicFvMesh == 1 &&
        (
            dynamicSolverType == "displacementLaplacian" ||
            dynamicSolverType == "solidBodyDisplacementLaplacian"
        )
    )
    {
        if (ca_Disp == nullptr)
            ca_Disp = new double[nTotal]{0};

        if (pointUpdated == nullptr)
            pointUpdated = new bool[npoints]{false};

        /*
        dynamicFvMesh& mesh(*meshPtr);
        if (pointDisplacementNewPtr == nullptr)
        {
            pointDisplacementNewPtr = new pointVectorField
            (
                IOobject
                (
                    "pointDisplacementNew",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                pointMesh::New(mesh)
            );
        }
        */
    }
    
    // Field-data
    nTotal = *ca_nCells * nComponents;
    ca_Vel = new double[nTotal]{0};
    ca_P   = new double[*ca_nCells]{0};
    
    if (TPtr != nullptr)
        ca_T = new double[*ca_nCells]{0};

    if (rhoPtr != nullptr)
        ca_Rho = new double[*ca_nCells]{0};

    const volVectorField& U(*UPtr);

#ifdef HAVE_OFE20
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

     std::string simulationType = ITstream
                                  (
                                      turbProperties.lookup("simulationType")
                                  ).toString();
#elif defined(HAVE_OF7)
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

#elif defined(HAVE_OF8)
    IOdictionary turbProperties
    ( 
        IOobject
        (
            momentumTransportModel::typeName,
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    word simulationType = turbProperties.lookup("simulationType");
#endif


    if (simulationType == "RAS")
    {
        const dictionary& subDict = turbProperties.subDict("RAS");

#ifdef HAVE_OFE20
        word RASModel = ITstream (subDict.lookup("RASModel")).toString();
#elif defined(HAVE_OF7)
        word RASModel = subDict.lookup("RASModel");
#elif defined(HAVE_OF8)
        word RASModel = subDict.lookup("model");
#endif

        if (RASModel == "kEpsilon")
        {
            ca_AlphaT = new double[*ca_nCells]{0};
            ca_K = new double[*ca_nCells]{0};
            ca_Epsilon = new double[*ca_nCells]{0};
            ca_NuT = new double[*ca_nCells]{0};
        }
        else if (RASModel == "kOmegaSST")
        {
            ca_AlphaT = new double[*ca_nCells]{0};
            ca_K = new double[*ca_nCells]{0};
            ca_Omega = new double[*ca_nCells]{0};
            ca_NuT = new double[*ca_nCells]{0};
        }
        else
        {
            FatalErrorInFunction
                << "Error: turbulence model not recongnized by the CSC Module."
                << nl << exit(FatalError);
        }
    }

    return 0;
}

int comFoam::updateVolumeData_outgoing()
{
    // Point data ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    const dynamicFvMesh& mesh(*meshPtr);
    const pointField&    points = mesh.points();
    //if (*isDynamicFvMesh == 1
    {
        forAll(points, ipoint)
        {
            for(int jcomp=0; jcomp<nComponents; jcomp++)
            {
                ca_Points[ipoint*nComponents+jcomp]
                    = points[ipoint][jcomp];
            }
        }
    }

    if (ca_Disp != nullptr)
    {
#if HAVE_OFE20
                const pointVectorField& pointDisplacement_ = mesh.lookupObject<displacementMotionSolver>
                        (
                            "dynamicMeshDict"
                        ).pointDisplacement();

#elif defined(HAVE_OF7) || defined(HAVE_OF8)

                const motionSolver& motion_ =
                    refCast<const dynamicMotionSolverFvMesh>(mesh).motion();

                const pointVectorField& pointDisplacement_ =
                        refCast<const displacementMotionSolver>(motion_).pointDisplacement();
#endif

        if (pointDisplacement_.size())
        {
            forAll(points, ipoint)
            {
                for(int jcomp=0; jcomp<nComponents; jcomp++)
                {
                    ca_Disp[ipoint*nComponents+jcomp]
                        = pointDisplacement_[ipoint][jcomp];
                }
            }
        }

        if (pointDisplacementNewPtr != nullptr)
        {
            pointVectorField &pointDisplacementNew(*pointDisplacementNewPtr);
            forAll(pointDisplacement_, ipoint)
            {
                pointDisplacementNew[ipoint] = pointDisplacement_[ipoint];
            }
        }
    }
    
    // Cell-centered data ^^^^^^^^^^^^^^^^^^^^^^^
    const volScalarField& p(*pPtr);
    const volVectorField& U(*UPtr);
    const volScalarField& T(*TPtr);
    const volScalarField& rho(*rhoPtr);

#ifdef HAVE_OFE20
    const compressible::turbulenceModel& turbulence(*turbulencePtr);
#elif defined(HAVE_OF7)
    const compressible::turbulenceModel& turbulence(*turbulencePtr);
#elif defined(HAVE_OF8)
    const compressible::momentumTransportModel& turbulence(*turbulencePtr);
    //const fluidThermophysicalTransportModel& thermoTransModel(*thermophysicalTransportPtr);
#endif

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
#ifdef HAVE_OFE20
                const tmp<volScalarField>& alphat = turbulence.alphat();
#elif defined(HAVE_OF7)
                const tmp<volScalarField>& alphat = turbulence.alphat();
#elif defined(HAVE_OF8)
                const tmp<volScalarField>& alphat = refCast<const volScalarField>
                    (
                        mesh.objectRegistry::lookupObject<volScalarField>
                        (
                            "alphat"
                        )
                    );
#endif

                ca_AlphaT[cellIndex] = alphat()[cellID];
            }

            if (ca_K != nullptr)
            {
                const tmp<volScalarField>& k = turbulence.k();
                ca_K[cellIndex] = k()[cellID];
            }

            if (ca_Epsilon != nullptr)
            {
                const tmp<volScalarField>& epsilon = turbulence.epsilon();
                ca_Epsilon[cellIndex] = epsilon()[cellID];
            }

            if (ca_Omega != nullptr)
            {
                const tmp<volScalarField>& omega = 
                    mesh.objectRegistry::lookupObject<volScalarField>
                    (
                        "omega"
                    );
                ca_Omega[cellIndex] = omega()[cellID];
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
    
    std::stringstream output{};
    output << "rocFoam.registerVolumeData: "
         << "Registering flow data with name "
         << volName;
    verbose_message(output.str(), true);

    // grid and field data
    int paneID = Pstream::myProcNo()+1;// Use this paneID for volume connectivity

    output = std::stringstream{};
    output << "procID = " << Pstream::myProcNo()
         << ", paneID = " << paneID
         << " ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^";
    verbose_message(output.str(), true);

    std::string dataName = volName+std::string(".nPoints");
    COM_new_dataitem( dataName, 'p', COM_INT, 1, "");
    COM_set_size(     dataName, paneID, 1);
    COM_set_array(    dataName, paneID, ca_nPoints);
    output = std::stringstream{};
    output << "  " << dataName.c_str() << " registered.";
    verbose_message(output.str(), true);

    dataName = volName+std::string(".nCells");
    COM_new_dataitem( dataName, 'p', COM_INT, 1, "");
    COM_set_size(     dataName, paneID, 1);
    COM_set_array(    dataName, paneID, ca_nCells);
    output = std::stringstream{};
    output << "  " << dataName.c_str() << " registered.";
    verbose_message(output.str(), true);

    dataName = volName+std::string(".cellToPointConn_types");
    COM_new_dataitem( dataName, 'p', COM_INT, 1, "");
    COM_set_size(     dataName, paneID, 1);
    COM_set_array(    dataName, paneID, ca_cellToPointConn_types);
    output = std::stringstream{};
    output << "  " << dataName.c_str() << " registered.";
    verbose_message(output.str(), true);

    int ntypes = *ca_cellToPointConn_types;
    dataName = volName+std::string(".cellToPointConn_map");
    COM_new_dataitem(dataName, 'p', COM_INT, 1, "");
    COM_set_size( dataName, paneID, ntypes);
    COM_set_array(dataName, paneID, ca_cellToPointConn_map);
    output = std::stringstream{};
    output << "  " << dataName.c_str() << " registered.";
    verbose_message(output.str(), true);

    dataName = volName+std::string(".cellToPointConn_size");
    COM_new_dataitem(dataName, 'p', COM_INT, 1, "");
    COM_set_size( dataName, paneID, ntypes);
    COM_set_array(dataName, paneID, ca_cellToPointConn_size);
    output = std::stringstream{};
    output << "  " << dataName.c_str() << " registered.";
    verbose_message(output.str(), true);

    // points
    dataName = volName+std::string(".nc");
    COM_set_size( dataName, paneID, *ca_nPoints);
    COM_set_array(dataName, paneID, ca_Points, nComponents);
    output = std::stringstream{};
    output << "  " << dataName.c_str() << " registered.";
    verbose_message(output.str(), true);

    // connectivity ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    for (int itype=0; itype<ntypes; itype++)
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
            FatalErrorInFunction
                << "=================== ERROR ===================" << endl
                << " Cell typeID " << typeID << " with size = "
                << typeSize << " not identified!"
                << nl << exit(FatalError);
        }

        COM_set_size( dataName, paneID, typeSize);
        COM_set_array(dataName,
                      paneID,
                      ca_cellToPointConn[itype],
                      typeID
                     );

        output = std::stringstream{};
        output << "  " << dataName.c_str() << " registered.";
        verbose_message(output.str(), true);
    }

    // Connectivity mapping stuff
    dataName = volName+std::string(".cellToCellMap");
    COM_new_dataitem( dataName, 'e', COM_INT, 1, "");
    COM_set_array(    dataName, paneID, ca_cellToCellMap, 1);
    output = std::stringstream{};
    output << "  " << dataName.c_str() << " registered.";
    verbose_message(output.str(), true);

    dataName = volName+std::string(".cellToCellMap_inverse");
    COM_new_dataitem( dataName, 'e', COM_INT, 1, "");
    COM_set_array(    dataName, paneID, ca_cellToCellMap_inverse, 1);
    output = std::stringstream{};
    output << "  " << dataName.c_str() << " registered.";
    verbose_message(output.str(), true);
    // ------------------------------------------

    // Element data registered with window
    dataName = volName+std::string(".vel");
    COM_new_dataitem( dataName, 'e', COM_DOUBLE, nComponents, "m/s");
    COM_set_array(    dataName, paneID, ca_Vel, nComponents);    
    output = std::stringstream{};
    output << "  " << dataName.c_str() << " registered.";
    verbose_message(output.str(), true);

    dataName = volName+std::string(".pres");
    COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "Pa");
    COM_set_array(    dataName, paneID, ca_P, 1);
    output = std::stringstream{};
    output << "  " << dataName.c_str() << " registered.";
    verbose_message(output.str(), true);

    if (ca_T != nullptr)
    {
        dataName = volName+std::string(".temp");
        COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "K");
        COM_set_array(    dataName, paneID, ca_T, 1);
        output = std::stringstream{};
        output << "  " << dataName.c_str() << " registered.";
        verbose_message(output.str(), true);
    }

    if (ca_Rho != nullptr)
    {
        dataName = volName+std::string(".rho");
        COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "kg/m^3");
        COM_set_array(    dataName, paneID, ca_Rho, 1);
        output = std::stringstream{};
        output << "  " << dataName.c_str() << " registered.";
        verbose_message(output.str(), true);
    }

    if (ca_Disp != nullptr)
    {
        dataName = volName+std::string(".disp");
        COM_new_dataitem( dataName, 'n', COM_DOUBLE, nComponents, "m");
        COM_set_array(    dataName, paneID, ca_Disp, nComponents);
        output = std::stringstream{};
        output << "  " << dataName.c_str() << " registered.";
        verbose_message(output.str(), true);
    }

    // Turbulence data ^^^^^^^^^^^^^^^^^^^^^^^^^^
    if (ca_AlphaT != nullptr)
    {
        dataName = volName+std::string(".alphaT");
        COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "kg/m/s");
        COM_set_array(    dataName, paneID, ca_AlphaT, 1);
        output = std::stringstream{};
        output << "  " << dataName.c_str() << " registered.";
        verbose_message(output.str(), true);
    }

    if (ca_K != nullptr)
    {
        dataName = volName+std::string(".k");
        COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "m^2/s^2");
        COM_set_array(    dataName, paneID, ca_K, 1);
        output = std::stringstream{};
        output << "  " << dataName.c_str() << " registered.";
        verbose_message(output.str(), true);
    }

    if (ca_Epsilon != nullptr)
    {
        dataName = volName+std::string(".epsilon");
        COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "m^2/s^3");
        COM_set_array(    dataName, paneID, ca_Epsilon, 1);
        output = std::stringstream{};
        output << "  " << dataName.c_str() << " registered.";
        verbose_message(output.str(), true);
    }

    if (ca_Omega != nullptr)
    {
        dataName = volName+std::string(".omega");
        COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "1/s");
        COM_set_array(    dataName, paneID, ca_Omega, 1);
        output = std::stringstream{};
        output << "  " << dataName.c_str() << " registered.";
        verbose_message(output.str(), true);
    }

    if (ca_NuT != nullptr)
    {
        dataName = volName+std::string(".nuT");
        COM_new_dataitem( dataName, 'e', COM_DOUBLE, 1, "m^2/s");
        COM_set_array(    dataName, paneID, ca_NuT, 1);
        output = std::stringstream{};
        output << "  " << dataName.c_str() << " registered.";
        verbose_message(output.str(), true);
    }
    //-------------------------------------------

    COM_window_init_done(volName); 

    return 0;
}


int comFoam::reconstVolumeData(const char *name)
{
    std::string volName = name+std::string("VOL");

    std::stringstream output{};
    output << "rocFoam.reconstCaVolumeData, procID = "
           << ca_myRank
           << ", Retreiving volume data form window "
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
            nameTmp == "disp" ||
            nameTmp == "alphaT" ||
            nameTmp == "epsilon" ||
            nameTmp == "omega" ||
            nameTmp == "k" ||
            nameTmp == "nuT"
            )
        {
            dataItemNames.push_back(nameTmp);

            output = std::stringstream{};
            output << "  DataItem[" << i << "] = " << nameTmp;
            verbose_message(output.str(), true);
        }
    }
    output = std::stringstream{};
    output << "  Number of items = " << dataItemNames.size()
              << std::endl;
    verbose_message(output.str(), true);

    // Volume data ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    int paneID = ca_myRank+1;// Use this paneID for volume connectivity
    output = std::stringstream{};
    output << "  procID = " << ca_myRank
         << ", paneID = " << paneID
         << " ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^";
    verbose_message(output.str(), true);

    std::string dataName = std::string("nPoints");
    nameExists(dataItemNames, dataName);
    std::string regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), paneID, &ca_nPoints);
    output = std::stringstream{};
    output << "  " << dataName.c_str() << " = " << *ca_nPoints;
    verbose_message(output.str(), true);

    dataName = std::string("nCells");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), paneID, &ca_nCells);
    output = std::stringstream{};
    output << "  " << dataName.c_str() << " = " << *ca_nCells;
    verbose_message(output.str(), true);

    dataName = std::string("cellToPointConn_types");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), paneID, &ca_cellToPointConn_types);
    output = std::stringstream{};
    output << "  " << dataName.c_str() << " = " << *ca_cellToPointConn_types;
    verbose_message(output.str(), true);

    dataName = std::string("cellToPointConn_map");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    int nComp;
    COM_get_array(regName.c_str(), paneID, &ca_cellToPointConn_map);
    COM_get_size(regName.c_str(), paneID, &nComp);
    for(int icomp=0; icomp<nComp; icomp++)
    {
        output = std::stringstream{};
        output << "  " << dataName.c_str() << "[" << icomp << "] = "
             << ca_cellToPointConn_map[icomp];
        verbose_message(output.str(), true);
    }

    dataName = std::string("cellToPointConn_size");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), paneID, &ca_cellToPointConn_size);
    COM_get_size(regName.c_str(), paneID, &nComp);
    for(int icomp=0; icomp<nComp; icomp++)
    {
        output = std::stringstream{};
        output << "  " << dataName.c_str() << "[" << icomp << "] = "
             << ca_cellToPointConn_size[icomp];
        verbose_message(output.str(), true);
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

    output = std::stringstream{};
    output << "  " << dataName.c_str() << " nPoints = " << nPoints
              << ", components = " << nComp;
    verbose_message(output.str(), true);

    int nConn;
    int numElem;
    std::string connNames;
    COM_get_connectivities(volName.c_str(), paneID, &nConn, connNames);
    std::istringstream connISS(connNames);
    for (int icon=0; icon<nConn; ++icon)
    {
        std::string connName;
        connISS >> connName;

        dataName = volName+std::string(".")+connName;
        COM_get_array(dataName.c_str(), paneID, &ca_cellToPointConn[icon], &nComp);
        COM_get_size(dataName.c_str(), paneID, &numElem);

        output = std::stringstream{};
        output << "    Connectivity[" << icon << "] = " << connName
                  << ", elements = " << numElem
                  << ", components =" << nComp;
        verbose_message(output.str(), true);
    }
    //---------------------------------------

    // Mapping data ^^^^^^^^^^^^^^^^^^^^^^^^^
    dataName = std::string("cellToCellMap");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), paneID, &ca_cellToCellMap, &nComp);
    COM_get_size(regName.c_str(), paneID, &numElem);
    output = std::stringstream{};
    output << "    " << dataName.c_str() << " elements = " << numElem
         << ", components = " << nComp;
    verbose_message(output.str(), true);

    dataName = std::string("cellToCellMap_inverse");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), paneID, &ca_cellToCellMap_inverse, &nComp);
    COM_get_size(regName.c_str(), paneID, &numElem);
    output = std::stringstream{};
    output << "    " << dataName.c_str() << " elements = " << numElem
              << ", components = " << nComp;
    verbose_message(output.str(), true);
    //---------------------------------------

    // Field data ^^^^^^^^^^^^^^^^^^^^^^^^^^^
    dataName = std::string("vel");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), paneID, &ca_Vel, &nComp);
    COM_get_size(regName.c_str(), paneID, &numElem);
    output = std::stringstream{};
    output << "    " << dataName.c_str() << " elements = " << numElem
              << ", components = " << nComp;
    verbose_message(output.str(), true);
    
    dataName = std::string("pres");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), paneID, &ca_P, &nComp);
    COM_get_size(regName.c_str(), paneID, &numElem);
    output = std::stringstream{};
    output << "    " << dataName.c_str() << " elements = " << numElem
              << ", components = " << nComp;
    verbose_message(output.str(), true);


    dataName = std::string("temp");
    if (nameExists(dataItemNames, dataName))
    {
        regName = volName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_T, &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);
        output = std::stringstream{};
        output << "    " << dataName.c_str() << " elements = " << numElem
                  << ", components = " << nComp;
        verbose_message(output.str(), true);
    }

    dataName = std::string("rho");
    if (nameExists(dataItemNames, dataName))
    {
        regName = volName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_Rho, &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);
        output = std::stringstream{};
        output << "    " << dataName.c_str() << " elements = " << numElem
                  << ", components = " << nComp;
        verbose_message(output.str(), true);
    }

    dataName = std::string("disp");
    if (nameExists(dataItemNames, dataName))
    {
        regName = volName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_Disp, &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);
        output = std::stringstream{};
        output << "    " << dataName.c_str() << " elements = " << numElem
                  << ", components = " << nComp;
        verbose_message(output.str(), true);

        if (pointUpdated == nullptr)
            pointUpdated = new bool[numElem]{false};
    }

    // Turbulence data ^^^^^^^^^^^^^^^^^^^^^^^^^^
    dataName = std::string("alphaT");
    if (nameExists(dataItemNames, dataName))
    {
        regName = volName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_AlphaT, &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);
        output = std::stringstream{};
        output << "    " << dataName.c_str() << " elements = " << numElem
                  << ", components = " << nComp;
        verbose_message(output.str(), true);
    }

    dataName = std::string("k");
    if (nameExists(dataItemNames, dataName))
    {
        regName = volName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_K, &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);
        output = std::stringstream{};
        output << "    " << dataName.c_str() << " elements = " << numElem
                  << ", components = " << nComp;
        verbose_message(output.str(), true);
    }    

    dataName = std::string("epsilon");
    if (nameExists(dataItemNames, dataName))
    {
        regName = volName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_Epsilon, &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);
        output = std::stringstream{};
        output << "    " << dataName.c_str() << " elements = " << numElem
                  << ", components = " << nComp;
        verbose_message(output.str(), true);
    }

    dataName = std::string("omega");
    if (nameExists(dataItemNames, dataName))
    {
        regName = volName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_Omega, &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);
        output = std::stringstream{};
        output << "    " << dataName.c_str() << " elements = " << numElem
                  << ", components = " << nComp;
        verbose_message(output.str(), true);
    }

    dataName = std::string("nuT");
    if (nameExists(dataItemNames, dataName))
    {
        regName = volName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_NuT, &nComp);
        COM_get_size(regName.c_str(), paneID, &numElem);
        output = std::stringstream{};
        output << "    " << dataName.c_str() << " elements = " << numElem
                  << ", components = " << nComp;
        verbose_message(output.str(), true);
    }
    //-------------------------------------------

    output = std::stringstream{};
    output << "  --------------------------------------------------";
    verbose_message(output.str(), true);

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

    if (ca_Disp != nullptr)
    {
        delete [] ca_Disp;
        ca_Disp = nullptr;
    }

    if (pointUpdated != nullptr)
    {
        delete [] pointUpdated;
        pointUpdated = nullptr;
    }

    if (pointDisplacementNewPtr != nullptr)
    {
        delete pointDisplacementNewPtr;
        pointDisplacementNewPtr = nullptr;
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

    if (ca_K != nullptr)
    {
        delete [] ca_K;
        ca_K = nullptr;
    }

    if (ca_Epsilon != nullptr)
    {
        delete [] ca_Epsilon;
        ca_Epsilon = nullptr;
    }

    if (ca_Omega != nullptr)
    {
        delete [] ca_Omega;
        ca_Omega = nullptr;
    }

    if (ca_NuT != nullptr)
    {
        delete [] ca_NuT;
        ca_NuT = nullptr;
    }
    //-------------------------------------------
    
    return 0;
}


