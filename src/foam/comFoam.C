#include "comFoam.H"
#include "cellShape.H"

comFoam::comFoam()
    : nPoints(0),
      nFaces(0),
      nCells(0),
      nPatches(0),
      winVolName(""),
      winSrfName(""),
      solverType(""),
      winComm(NULL),
      winNProc(0),
      winRank(0),
      winTime(0.0),
      winDeltaT(0.0),
      winRun(1)
{};

comFoam::comFoam(int *pargc, void **pargv, const char *name)
    : nPoints(0),
      nFaces(0),
      nCells(0),
      nPatches(0),
      winVolName(""),
      winSrfName(""),
      solverType(""),
      winComm(NULL),
      winNProc(0),
      winRank(0),
      winTime(0.0),
      winDeltaT(0.0),
      winRun(1)
{
    flowInit(pargc, pargv, name);
}

//^^^ DEFINITION OF COM-RELATED MTHODS ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
int comFoam::flowInit(int *pargc, void **pargv, const char *name)
{
    MPI_Comm tmpComm = COM_get_default_communicator();  

    int tmpRank;
    MPI_Comm_rank(tmpComm, &tmpRank);
    
    if (tmpRank == 0)
    {
        std::cout << "rocFoam.flowInit: Initializing flow solver."
                  << std::endl;
    }
    //  OpenFOAM initializer ^^^^^^^^^^^^^^^^^^^^
    comFoam *comFoamPtr = NULL;

    std::string volName = name+string("VOL");
    std::string objectName = volName+string(".object");

    COM_get_object(objectName.c_str(), 0, &comFoamPtr);

    int argc = *pargc;
    char** argv = reinterpret_cast<char**>(pargv);

    comFoamPtr->initialize(argc, argv);

    //  Other initializations ^^^^^^^^^^^^^^^^^^^
    // extractData can be called here, or in
    // rocFoam driver with comExtractData

    comFoamPtr->createVolumeConnectivities();
    comFoamPtr->createVolumeData();
    comFoamPtr->updateVolumeData();
    comFoamPtr->registerData(name);

    return 0;
}


int comFoam::flowLoop(const char *name)
{
    Foam::Info << "rocFoam.flowLoop: Iterating flow solver." << Foam::endl;

    //  Call the flow iterator ^^^^^^^^^^^^^^^^^^
    comFoam *comFoamPtr = NULL;
    std::string volName = name+string("VOL");
    std::string objectName = volName+string(".object");
    COM_get_object(objectName.c_str(), 0, &comFoamPtr);

    comFoamPtr->loop();
    
    return 0;
}

int comFoam::flowStep(const char *name)
{

    Foam::Info << "rocFoam.flowStep: Stepping flow solver." << Foam::endl;

    //  Call the flow iterator ^^^^^^^^^^^^^^^^^^
    comFoam *comFoamPtr = NULL;
    std::string volName = name+string("VOL");
    std::string objectName = volName+string(".object");
    COM_get_object(objectName.c_str(), 0, &comFoamPtr);

    comFoamPtr->step();
    
    return 0;
}

//int comFoam::flowExtractData(const char *name)
//{
//    Foam::Info << "rocFoam.extractData: Extracting flow data." << Foam::endl;
//    //  Call the flow iterator ^^^^^^^^^^^^^^^^^^
//    comFoam *comFoamPtr = NULL;
//    std::string volName = name+string("VOL");
//    std::string objectName = volName+string(".object");
//    COM_get_object(objectName.c_str(), 0, &comFoamPtr);
//    comFoamPtr->extractData();
//    
//    return 0;
//}

int comFoam::flowRegisterData(const char *name)
{

    Foam::Info << "rocFoam.registerData: Extracting flow data." << Foam::endl;

    //  Call the flow iterator ^^^^^^^^^^^^^^^^^^
    comFoam *comFoamPtr = NULL;
    std::string volName = name+string("VOL");
    std::string objectName = volName+string(".object");
    COM_get_object(objectName.c_str(), 0, &comFoamPtr);

    comFoamPtr->registerData(name);
    
    return 0;
}

//^^^^^ REGISTER FUNCTIONS ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
int comFoam::registerFunctions(const char *name)
{
    //  Anouncing default communicator  ^^^^^^^^^^^^^^^^^^^
    MPI_Comm tmpComm = COM_get_default_communicator();  

    int tmpRank;
    MPI_Comm_rank(tmpComm, &tmpRank);
    
    if (tmpRank == 0)
    {
        std::cout << "rocFoam.flowRegister: "
                  << "Registering flow functions with name "
                  << name
                  << std::endl;
    }
    
    //  Register module with COM ^^^^^^^^^^^^^^^^^^^^^^^^^^
    comFoam *comFoamPtr = NULL;

    //std::string name="ROCFOAM";
    std::string objectName = name+string(".object");
    COM_get_object(objectName.c_str(), 0, &comFoamPtr);

    /// Register functions ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    std::vector<COM_Type> types(13,COM_INT);

    types[0] = COM_RAWDATA;
    types[2] = COM_VOID;

    std::string functionName = name+string(".flowInit");
    COM_set_member_function
    (
        functionName.c_str(),
        reinterpret_cast<Member_func_ptr>(&comFoam::flowInit),
        objectName.c_str(),
        "biii",
        &types[0]
    );

    functionName = name+string(".flowLoop");
    COM_set_member_function
    (
        functionName.c_str(),
        reinterpret_cast<Member_func_ptr>(&comFoam::flowLoop),
        objectName.c_str(),
        "bi",
        &types[0]
    );

    functionName = name+string(".flowStep");
    COM_set_member_function
    (
        functionName.c_str(),
        reinterpret_cast<Member_func_ptr>(&comFoam::flowStep),
        objectName.c_str(),
        "bi",
        &types[0]
    );

//    functionName = name+string(".flowExtractData");
//    COM_set_member_function
//    (
//        functionName.c_str(),
//        reinterpret_cast<Member_func_ptr>(&comFoam::flowExtractData),
//        objectName.c_str(),
//        "bi",
//        &types[0]
//    );

//    functionName = name+string(".flowRegisterDate");
//    COM_set_member_function
//    (
//        functionName.c_str(),
//        reinterpret_cast<Member_func_ptr>(&comFoam::flowExtractData),
//        objectName.c_str(),
//        "bi",
//        &types[0]
//    );


    //COM_set_member_function
    //(
    //    (name + string(".flowFin")).c_str(),
    //    reinterpret_cast<Member_func_ptr>(&rhoCentral::flowFin),
    //    objectName.c_str(), "b", &types[0]
    //);

//    //  Registering data of this module to COM ^^^^^^^^^^^^
//    std::string dataName = name+string(".winNProc");
//    COM_new_dataitem( dataName.c_str(), 'w', COM_INT, 1, "");
//    COM_set_size(     dataName.c_str(), 0, 1);
//    COM_set_array(    dataName.c_str(), 0, &(comFoamPtr->winNProc));

//    dataName = name+string(".winTime");
//    COM_new_dataitem( dataName.c_str(), 'w', COM_DOUBLE, 1, "");
//    COM_set_size(     dataName.c_str(), 0, 1);
//    COM_set_array(    dataName.c_str(), 0, &(comFoamPtr->winTime));

//    dataName = name+string(".winDeltaT");
//    COM_new_dataitem( dataName.c_str(), 'w', COM_DOUBLE, 1, "");
//    COM_set_size(     dataName.c_str(), 0, 1);
//    COM_set_array(    dataName.c_str(), 0, &(comFoamPtr->winDeltaT) );

//    dataName = name+string(".winRun");
//    COM_new_dataitem( dataName.c_str(), 'w', COM_INT, 1, "");
//    COM_set_size(     dataName.c_str(), 0, 1);
//    COM_set_array(    dataName.c_str(), 0, &(comFoamPtr->winRun));

    COM_window_init_done(name); 

    return 0;
}
//---------------------------------------------------------

//===================================================================
int comFoam::createVolumeConnectivities()
{

    const dynamicFvMesh& mesh(*meshPtr);
    const Foam::Time&    runTime(*runTimePtr);


    // Mesh and conmnectivities ^^^^^^^^^^^^^^^^^
    const pointField& points     = mesh.points();
    const cellList&   cells      = mesh.cells();
    const faceList&   faces      = mesh.faces();
    const labelList&  faceOwner  = mesh.faceOwner();
    const labelList&  faceNeighb = mesh.faceNeighbour();    

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const labelListList& cellPoints = mesh.cellPoints();
    const cellShapeList& cellShapes = mesh.cellShapes();
    //-------------------------------------------

    // Temporary Vectors ^^^^^^^^^^^^^^^^^^^^^^^^
    std::vector< std::vector<double> > vecVecTmpDouble;
    std::vector<double> vecTmpDouble;
    
    std::vector<int> vecTmpInt;
    std::vector< std::vector<int> > vecVecTmpInt;
    //-------------------------------------------

    nPoints = mesh.nPoints();
    nFaces  = mesh.nFaces();
    nCells  = mesh.nCells();

    // cellToPoint connectivities ^^^^^^^^^^^^^^^
    vecCellToPointConn.clear();
    //vecCellToPointConnVTK.clear();
    vecCellTypeLocalToGlobalMap.clear();
    
    vecTmpInt.clear();
    vecVecTmpInt.clear();
    
    //vecTmpInt.push_back(0);
    //vecVecTmpInt.push_back(vecTmpInt);
    // predefine the number of faces per cell
    int cellPointTypeSize=8;
    vecCellToPointConn.resize(cellPointTypeSize);
    //vecCellToPointConnVTK.resize(cellPointTypeSize);
    vecCellTypeLocalToGlobalMap.resize(cellPointTypeSize);
    for (int i=0; i<cellPointTypeSize; i++)
    {
        vecCellToPointConn[i].clear();
        //vecCellToPointConnVTK[i].clear();
        vecCellTypeLocalToGlobalMap[i].clear();
    }

    forAll(cellPoints, icell)
    {
        const labelList& pointsList = cellPoints[icell];
        int nPointsInCell = pointsList.size();

        vecTmpInt.clear();
        forAll(pointsList, ipoint)
        {
            vecTmpInt.push_back(pointsList[ipoint]);
        }
        
        vecCellToPointConn[nPointsInCell-1].push_back(vecTmpInt);
        vecCellTypeLocalToGlobalMap[nPointsInCell-1].push_back(icell);
    }
    
    // Connectivity for VTK format
    /*
    forAll(cellShapes, icell)
    {
        const cellShape& cellShape = cellShapes[icell];
        int nPointsInCell = cellShape.size();

        vecTmpInt.clear();
        forAll(cellShape, ipoint)
        {
            vecTmpInt.push_back(cellShape[ipoint]);
        }
        
        vecCellToPointConnVTK[nPointsInCell-1].push_back(vecTmpInt);        
    }
    */
    //-------------------------------------------


    // Create local-to-global cell map according to types
    int totalnCellTypes = vecCellToPointConn.size();

    for(int itype=0; itype<totalnCellTypes; itype++)
    {
        if (vecCellToPointConn[itype].size()>0)
        {
            nCellTypes++;
        }
    }

    nCellTypeMap  = new int[nCellTypes];
    nCellTypeSize = new int[nCellTypes];
    cellTypeLocalToGlobalMap = new int*[nCellTypes];
    
    int checkSum = 0;
    nCellTypes = 0;
    for(int itype=0; itype<totalnCellTypes; itype++)
    {
        if (vecCellToPointConn[itype].size()>0)
        {
            nCellTypeMap[nCellTypes] = itype+1;
            
            int typeSize = vecCellToPointConn[itype].size();
            nCellTypeSize[nCellTypes] = typeSize;

            cellTypeLocalToGlobalMap[nCellTypes] = new int[typeSize];

            for(int icell=0; icell<typeSize; icell++)
            {
                cellTypeLocalToGlobalMap[nCellTypes][icell] =
                    vecCellTypeLocalToGlobalMap[itype][icell];
            
                checkSum++;
            }

            nCellTypes++;
        }
    }

    if (checkSum != nCells)
    {
        Foam::Info << "========== WARNNING ==============="
                   << "     checkSum != nCells " << endl;
        return -1;
    }

    comArrCelltoPointConn = new int*[nCellTypes];
    for(int itype=0; itype<nCellTypes; itype++)
    {
        int nPoints_ = nCellTypeMap[itype];
        int nCells_  = nCellTypeSize[itype];
        
        int nTypeConn = nCells_ * nPoints_;
        comArrCelltoPointConn[itype] = new int[nTypeConn];

        for (int icell=0; icell<nCells_; icell++)
        {
            int cellID = cellTypeLocalToGlobalMap[itype][icell];
            const cellShape& cellShape_ = cellShapes[cellID];

            for (int ipoint=0; ipoint<nPoints_; ipoint++)
            {
                int index = icell*nPoints_ + ipoint;
                comArrCelltoPointConn[itype][index] = cellShape_[ipoint];
            }
        }
    }
    
    //-------------------------------------------


    // Save faceToCell and faceToPoint conn ^^^^^
    vecOwners.clear();
    vecNeighbs.clear();
    vecFaceToPointConn.clear();
    forAll(faces, i)
    {
        const labelList& pointsList = faces[i];
        
        vecOwners.push_back(faceOwner[i]);
        vecNeighbs.push_back(faceNeighb[i]);

        vecTmpInt.clear();
        forAll(pointsList, j)
        {
            vecTmpInt.push_back(pointsList[j]);
        }
        
        vecFaceToPointConn.push_back(vecTmpInt);
    }
    //-------------------------------------------


    // Save cellToFace connectivity ^^^^^^^^^^^^
    // Not sure if it is needed!
    /*
    vecCellToFaceConn.clear();
    vecTmpInt.clear();
    vecVecTmpInt.clear();
    // predefine the number of faces per cell
    int cellFaceTypeSize=6;
    vecCellToFaceConn.resize(cellFaceTypeSize)
    for (int i=0; i<cellTypeSize; i++)
    {
        vecCellToFaceConn[i].clear();
    }

    forAll(cells, icell)
    {
        const labelList& faceList = cells[icell];
        int nFacesInCell = faceList.size();
        
        vecTmpInt.clear();
        forAll(faceList, iface)
        {
            vecTmpInt.push_back(faceList[iface]);
        }
        vecCellToFaceConn[nFacesInCell-1].push_back(vecTmpInt);
    }
    */
    //-------------------------------------------
    return 0;
}

int comFoam::deleteVolumeData()
{

    if (nCellTypeMap != NULL)
    {
        delete [] nCellTypeMap;
        nCellTypeMap = NULL;
    }
    
    if (nCellTypeSize != NULL)
    {
        delete [] nCellTypeSize;
        nCellTypeSize = NULL;
    }
   
    if (cellTypeLocalToGlobalMap != NULL)
    {
        for (int itype=0; itype<nCellTypes; itype++)
        {
            if (cellTypeLocalToGlobalMap[itype] != NULL)
            {
                delete [] cellTypeLocalToGlobalMap[itype];
            }
        }
        delete [] cellTypeLocalToGlobalMap;
        cellTypeLocalToGlobalMap = NULL;
    }

    if (comArrCelltoPointConn != NULL)
    {
        for (int itype=0; itype<nCellTypes; itype++)
        {
            if (comArrCelltoPointConn[itype] != NULL)
            {
                delete [] comArrCelltoPointConn[itype];
            }
        }
    }
    comArrCelltoPointConn = NULL;

    if (comArrPoints != NULL)
    {
        delete[] comArrPoints;
        comArrPoints = NULL;
    }

    if (comArrVel != NULL)
    {
        delete[] comArrVel;
        comArrVel = NULL;
    }

    if (comArrP != NULL)
    {
        delete[] comArrP;
        comArrP = NULL;
    }

    if (comArrT != NULL)
    {
        delete[] comArrT;
        comArrT = NULL;
    }

    if (comArrRho != NULL)
    {
        delete[] comArrRho;
        comArrRho = NULL;
    }
}


int comFoam::createVolumeData()
{
    int nTotal = nPoints * nComponents;
    comArrPoints = new double[nTotal];

    nTotal = nCells * nComponents;
    comArrVel = new double[nTotal];

    comArrP = new double[nCells];
    comArrT = new double[nCells];
    comArrRho = new double[nCells];
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
            comArrPoints[ipoint*nComponents+jcomp]
                = points[ipoint][jcomp];
        }
    }
    
    // Cell-centered data ^^^^^^^^^^^^^^^^^^^^^^^
    const volScalarField& p(*pPtr);
    const volVectorField& U(*UPtr);
    const volScalarField& T(*TPtr);
    const volScalarField& rho(*rhoPtr);

    for(int itype=0; itype<nCellTypes; itype++)
    {
        int size = nCellTypeSize[itype];
        for(int icell=0; icell<size; icell++)
        {
            int cellID = cellTypeLocalToGlobalMap[itype][icell];

            for(int jcomp=0; jcomp<nComponents; jcomp++)
            {
                int index = itype*size*nComponents +
                            icell*nComponents +
                            jcomp;

                comArrVel[index] = U[cellID][jcomp];
            }

            int index = itype*size + icell;
            comArrP[index] = p[cellID];
            comArrT[index] = T[cellID];
            comArrRho[index] = rho[cellID];
        }
    }

    return 0;
}


int comFoam::createSurfaceConnectivities()
{

    const dynamicFvMesh& mesh(*meshPtr);
    const Foam::Time&     runTime(*runTimePtr);
    const volScalarField& p(*pPtr);
    const volVectorField& U(*UPtr);
    const volScalarField& T(*TPtr);
    const volScalarField& rho(*rhoPtr);

    // Mesh and conmnectivities ^^^^^^^^^^^^^^^^^
    const pointField& points     = mesh.points();
    const cellList&   cells      = mesh.cells();
    const faceList&   faces      = mesh.faces();
    const labelList&  faceOwner  = mesh.faceOwner();
    const labelList&  faceNeighb = mesh.faceNeighbour();    

    const polyBoundaryMesh& patches    = mesh.boundaryMesh();
    const labelListList&    cellPoints = mesh.cellPoints();
    const cellShapeList&    cellShapes = mesh.cellShapes();
    //-------------------------------------------

    // Temporary Vectors ^^^^^^^^^^^^^^^^^^^^^^^^
    std::vector< std::vector<double> > vecVecTmpDouble;
    std::vector<double> vecTmpDouble;
    
    std::vector<int> vecTmpInt;
    std::vector< std::vector<int> > vecVecTmpInt;
    //-------------------------------------------

    nPatches = patches.size();
    
    
    // Patch Vectors ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    patchGlobalPointIndex.clear();
    patchLocalConnectivity.clear();

    vecPatchVel.clear();
    vecPatchRho.clear();
    vecPatchP.clear();
    vecPatchT.clear();

    vecPatchName.clear();
    vecPatchType.clear();
    vecPatchInGroup.clear();

    vecPatchStart.clear();
    vecPatchNFaces.clear();
    vecPatchNPoints.clear();
    //-------------------------------------------
    
    // Save Patch faceToNode connectivity  ^^^^^^
    // Local and Global
    forAll(patches, ipatch)
    {
        const polyPatch& patch = patches[ipatch];

        const word& patchName = patch.name();
        const word& patchType = patch.type();
        const wordList& patchInGroup = patch.inGroups();

        const label& patchStart = patch.start();
        const int& patchSize = patch.size();

        vecPatchName   .push_back(patchName);
        vecPatchType   .push_back(patchType);
        vecPatchInGroup.push_back(patchInGroup);
        vecPatchStart  .push_back(patchStart);
        vecPatchNFaces .push_back(patchSize);

        std::vector<int> vecTmpInt_patchGlobalPointIndex;
        std::vector<std::vector<int>> vecTmpInt_patchLocalConnectivity;

        for(int iface=0; iface<patchSize; iface++)
        {
            const label& face = patchStart + iface;
            const labelList& pointsList = faces[face];

            vecTmpInt.clear();
            forAll(pointsList, ipoint)
            {
                const int& point = pointsList[ipoint];
                std::vector<int>::iterator index = std::find
                                  (
                                    vecTmpInt_patchGlobalPointIndex.begin(),
                                    vecTmpInt_patchGlobalPointIndex.end(),
                                    point
                                  );

                if  ( index == vecTmpInt_patchGlobalPointIndex.end() )
                {
                    vecTmpInt_patchGlobalPointIndex.push_back(point);

                    index = vecTmpInt_patchGlobalPointIndex.end() - 1;
                }

                int indexVal = std::distance
                               (
                                vecTmpInt_patchGlobalPointIndex.begin(),
                                index
                               );

                vecTmpInt.push_back( indexVal );
            }
            vecTmpInt_patchLocalConnectivity.push_back(vecTmpInt);
        }

        patchGlobalPointIndex.push_back(vecTmpInt_patchGlobalPointIndex);
        patchLocalConnectivity.push_back(vecTmpInt_patchLocalConnectivity);

        vecPatchNPoints.push_back(vecTmpInt_patchGlobalPointIndex.size()) ;
    }
    //-------------------------------------------

    //  Save Flow Quantities on Patches ^^^^^^^^^
    vecPatchVel.clear();
    vecPatchRho.clear();
    vecPatchP.clear();
    vecPatchT.clear();
    for(int ipatch=0; ipatch<nPatches; ipatch++)
    {
        const int& nFaces = vecPatchNFaces[ipatch];

        //  Patch Face Velocity
        vecVecTmpDouble.clear();
        if (vecPatchType[ipatch] != string("empty"))
        {
            for(int iface=0; iface<nFaces; iface++)
            {
                vecTmpDouble.clear();
                for(int k=0; k<nComponents; k++)
                {
                    vecTmpDouble.push_back
                    ( 
                        U.boundaryField()[ipatch][iface].component(k)
                    );
                }
                vecVecTmpDouble.push_back(vecTmpDouble);
            }
        }
        vecPatchVel.push_back(vecVecTmpDouble);

        //  Patch Face Density
        vecTmpDouble.clear();
        if (vecPatchType[ipatch] != string("empty"))
        {
            for(int iface=0; iface<nFaces; iface++)
            {
                vecTmpDouble.push_back( rho.boundaryField()[ipatch][iface] );
            }
        }
        vecPatchRho.push_back(vecTmpDouble);

        //  Patch Face Pressure
        vecTmpDouble.clear();
        if (vecPatchType[ipatch] != string("empty"))
        {
            for(int iface=0; iface<nFaces; iface++)
            {
                vecTmpDouble.push_back( p.boundaryField()[ipatch][iface] );
            }
        }
        vecPatchP.push_back(vecTmpDouble);

        //  Patch Face Temperature
        vecTmpDouble.clear();
        if (vecPatchType[ipatch] != string("empty"))
        {
            for(int iface=0; iface<nFaces; iface++)
            {
                vecTmpDouble.push_back( T.boundaryField()[ipatch][iface] );
            }
        }
        vecPatchT.push_back(vecTmpDouble);
    }
    //-------------------------------------------
}


int comFoam::registerData(const char *name)
{
    Foam::Info << "rocFoam.registerData: "
               << "Registering flow data with name "
               << name
               << endl;

    std::string volName = name+string("VOL");

    std::string dataName = volName+string(".winNProc");
    COM_new_dataitem( dataName.c_str(), 'w', COM_INT, 1, "");
    COM_set_size(     dataName.c_str(), 0, 1);
    COM_set_array(    dataName.c_str(), 0, &(winNProc));

    dataName = volName+string(".winTime");
    COM_new_dataitem( dataName.c_str(), 'w', COM_DOUBLE, 1, "");
    COM_set_size(     dataName.c_str(), 0, 1);
    COM_set_array(    dataName.c_str(), 0, &(winTime));

    dataName = volName+string(".winDeltaT");
    COM_new_dataitem( dataName.c_str(), 'w', COM_DOUBLE, 1, "");
    COM_set_size(     dataName.c_str(), 0, 1);
    COM_set_array(    dataName.c_str(), 0, &(winDeltaT) );

    dataName = volName+string(".winRun");
    COM_new_dataitem( dataName.c_str(), 'w', COM_INT, 1, "");
    COM_set_size(     dataName.c_str(), 0, 1);
    COM_set_array(    dataName.c_str(), 0, &(winRun));

    // points
    dataName = volName+string(".nc");
    COM_set_size( dataName.c_str(), 1, nPoints);
    COM_set_array(dataName.c_str(), 1, comArrPoints, nComponents);

    // connectivity
    for(int itype=0; itype<nCellTypes; itype++)
    {
        int typeID = nCellTypeMap[itype];
        int typeSize = nCellTypeSize[itype];
        
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
            dataName = volName+string(".:q8");
        }
        else
        { // Type not identified

            Foam::Info << "=================== WARNING ==================="
                       << " Cell typeID " << typeID << " with size = "
                       << typeSize << " not identified!"
                       << endl;
            return -1;
        }

        COM_set_size( dataName.c_str(), 1, typeSize);
        COM_set_array(
                        dataName.c_str(),
                        1,
                        comArrCelltoPointConn[itype],
                        typeID
                     );
    }
        
    //dataName = volName+string(".:q8");
    //COM_set_size( dataName.c_str(), 1, nCells);
    //COM_set_array(dataName.c_str(), 1, comFoamPtr->comCONN8, 8);

    dataName = volName+string(".vel");
    COM_new_dataitem( dataName.c_str(), 'e', COM_DOUBLE, nComponents, "m/s");
    COM_set_array(    dataName.c_str(), 1, comArrVel, nComponents);    

    dataName = volName+string(".pres");
    COM_new_dataitem( dataName.c_str(), 'e', COM_DOUBLE, 1, "Pa");
    COM_set_array(    dataName.c_str(), 1, comArrP, 1);

    dataName = volName+string(".temp");
    COM_new_dataitem( dataName.c_str(), 'e', COM_DOUBLE, 1, "K");
    COM_set_array(    dataName.c_str(), 1, comArrT, 1);

    dataName = volName+string(".rho");
    COM_new_dataitem( dataName.c_str(), 'e', COM_DOUBLE, 1, "kg/m^3");
    COM_set_array(    dataName.c_str(), 1, comArrRho, 1);

    COM_window_init_done(volName); 

    return 0;
}

