#include "comFoam.H"
#include "cellShape.H"

comFoam::comFoam()
{
    initSet();
};

comFoam::comFoam(int *pargc, void **pargv, const char *name)
{
    initSet();
    flowInit(pargc, pargv, name);
}

int comFoam::initSet()
{
    // Variables ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ca_nPoints = NULL;  //single value
    ca_nCells = NULL;   //single value
    ca_nFaces = NULL;   //single value
    ca_nPatches = NULL; //single value
    //-------------------------------------------

    // COM Volume Arrays^^^^^^^^^^^^^^^^^^^^^^^^^^
    // Mapping
    ca_cellToCellMap = NULL;

    // Connectivity
    ca_cellToPointConn_types = NULL; //single value
    ca_cellToPointConn_map = NULL;
    ca_cellToPointConn_size = NULL;
    ca_cellToPointConn = NULL;

    // Field Data
    ca_Points = NULL;
    ca_Vel = NULL;
    ca_P = NULL;
    ca_T = NULL;
    ca_Rho = NULL;
    //-------------------------------------------
    
    // COM Face Arrays^^^^^^^^^^^^^^^^^^^^^^^^^^^
    // Mapping
    ca_faceToFaceMap = NULL;

    // Connectivity
    ca_faceToPointConn_types = NULL; //single value
    ca_faceToPointConn_map = NULL;
    ca_faceToPointConn_size = NULL;
    ca_faceToPointConn = NULL;

    // Field data
    ca_faceOwner = NULL;
    ca_faceNeighb = NULL;
    //-------------------------------------------

    // COM Patch Arrays^^^^^^^^^^^^^^^^^^^^^^^^^^
    // General data

    patchNameStr = NULL; //single value for the last
    patchTypeStr = NULL; //single value for the last

    ca_patchName = NULL; //single value for the last
    ca_patchType = NULL; //single value for the last
    ca_patchInGroup = NULL; //single value for the last
    ca_patchStart = NULL; //single value for the last
    ca_patchSize = NULL; //single value for the last

    // PointToPoint Mapping
    ca_patchPointToPointMap_size = NULL; //single value for the last
    ca_patchPointToPointMap = NULL;

    // FaceToFace Mapping
    ca_patchFaceToFaceMap = NULL;

    // FaceToPoint Mapping
    ca_patchFaceToPointConn_types = NULL; //single value for the last
    ca_patchFaceToPointConn_map = NULL;
    ca_patchFaceToPointConn_size = NULL;
    ca_patchFaceToPointConn = NULL;

    // Field data
    ca_patchPoints = NULL;
    ca_patchVel = NULL;
    ca_patchP = NULL;
    ca_patchT = NULL;
    ca_patchRho = NULL;

    // File data
    ca_nFiles = NULL;
    ca_fileSize = NULL;
    ca_fileName = NULL;
    ca_filePath = NULL;
    ca_fileContent = NULL;

    //-------------------------------------------

    //  Window data ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    //std::string winVolName; /// Tracks *this* volume window name.
    //std::string winSurfName; /// Tracks *this* volume window name.

    solverType = "";
    winComm = NULL;
    
    ca_nProc = 1;
    ca_myRank = 0;

    // registered data set during the simulation
    ca_runStat = NULL;
    ca_time = NULL;
    ca_deltaT = NULL;

    return 0;
}



//^^^ DEFINITION OF COM-RELATED MTHODS ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
int comFoam::flowInit(int *pargc, void **pargv, const char *name)
{

    MPI_Comm tmpComm = COM_get_default_communicator();  
    int tmpRank;
    MPI_Comm_rank(tmpComm, &tmpRank);
    
    if (tmpRank == 0)
    {
        std::cout << "rocFoam.flowInit: Initializing flow solver with name "
                  << name << std::endl;
    }
    //  OpenFOAM initializer ^^^^^^^^^^^^^^^^^^^^
    comFoam *comFoamPtr = NULL;

    std::string volName = name+std::string("VOL");
    std::string surfName = name+std::string("SURF");
    std::string objectName = volName+std::string(".object");
    COM_get_object(objectName.c_str(), 0, &comFoamPtr);

    //char** argv = reinterpret_cast<char**>(pargv);
    int argc = *pargc+2;
    char** argv;
    argv = new char*[*pargc+2];
    for (int i=0; i<*pargc; i++)
    {
        std::string strTmp = reinterpret_cast<char*>(pargv[i]);
        argv[i] = new char[strTmp.length()+1];
        std::strcpy(argv[i], strTmp.c_str());
    }

    std::string strTmp = "-case";
    argv[*pargc] = new char[strTmp.length()+1];
    std::strcpy(argv[*pargc], strTmp.c_str());

    strTmp = "./";
    argv[*pargc+1] = new char[strTmp.length()+1];
    std::strcpy(argv[*pargc+1], strTmp.c_str());
    
    comFoamPtr->initialize(argc, argv);

    //  Other initializations ^^^^^^^^^^^^^^^^^^^
    // extractData can be called here, or in
    // rocFoam driver with comExtractData

    /*
    comFoamPtr->createVolumeConnectivities();
    comFoamPtr->createVolumeData();
    comFoamPtr->updateVolumeData();
    comFoamPtr->registerVolumeData(name);

    comFoamPtr->createFaceConnectivities();
    comFoamPtr->createFaceData();
    comFoamPtr->updateFaceData();
    comFoamPtr->registerFaceData(name);

    comFoamPtr->createSurfaceConnectivities();
    comFoamPtr->createSurfaceData();
    comFoamPtr->updateSurfaceData();
    comFoamPtr->registerSurfaceData(name);
    */

    createVolumeConnectivities();
    createVolumeData();
    updateVolumeData();
    registerVolumeData(name);

    createFaceConnectivities();
    createFaceData();
    updateFaceData();
    registerFaceData(name);

    createSurfaceConnectivities();
    createSurfaceData();
    updateSurfaceData();
    registerSurfaceData(name);

    readInitFiles(strTmp);
    registerInitFiles(name);

    return 0;
}



int comFoam::reconstCaData(int *pargc, void **pargv, const char *name)

{
    MPI_Comm tmpComm = COM_get_default_communicator();
    int tmpRank;
    MPI_Comm_rank(tmpComm, &tmpRank);
    
    if (tmpRank == 0)
    {
        std::cout << "rocFoam.reconstCaData: Initializing CA "
                  << " reconstructions for window "
                  << name << std::endl;
    }

    //  OpenFOAM initializer ^^^^^^^^^^^^^^^^^^^^
    comFoam *comFoamPtr = NULL;
    std::string volName = name+std::string("VOL");
    std::string surfName = name+std::string("SURF");
    std::string dataName = volName+std::string(".object");
    COM_get_object(dataName.c_str(), 0, &comFoamPtr);

    reconstCaVolumeData(name);
    reconstCaFaceData(name);
    reconstCaSurfaceData(name);
    std::string tmpDir = "./fluidTmp";
    deleteInitFiles(tmpDir);
    reconstCaInitFiles(name, tmpDir);

    //int argc = *pargc;
    //char** argv = reinterpret_cast<char**>(pargv);
    int argc = *pargc+2;
    char** argv;
    argv = new char*[*pargc+2];
        
    std::string strTmp = "-case";
    argv[*pargc] = new char[strTmp.length()+1];
    std::strcpy(argv[*pargc], strTmp.c_str());

    strTmp = tmpDir;
    argv[*pargc+1] = new char[strTmp.length()+1];
    std::strcpy(argv[*pargc+1], strTmp.c_str());

    comFoamPtr->initialize(argc, argv, true);

//    comFoamPtr->reconstDynamicFvMesh();
//    comFoamPtr->createFields_COM();

    //deleteInitFiles(tmpDir);

    return 0;
}

int comFoam::flowLoop(const char *name)
{
    Foam::Info << "rocFoam.flowLoop: Iterating flow solver." << Foam::endl;

    //  Call the flow iterator ^^^^^^^^^^^^^^^^^^
    comFoam *comFoamPtr = NULL;
    std::string volName = name+std::string("VOL");
    std::string objectName = volName+std::string(".object");
    COM_get_object(objectName.c_str(), 0, &comFoamPtr);

    comFoamPtr->loop();
    
    return 0;
}

int comFoam::flowStep(const char *name)
{

    Foam::Info << "rocFoam.flowStep: Stepping flow solver." << Foam::endl;

    //  Call the flow iterator ^^^^^^^^^^^^^^^^^^
    comFoam *comFoamPtr = NULL;
    std::string volName = name+std::string("VOL");
    std::string objectName = volName+std::string(".object");
    COM_get_object(objectName.c_str(), 0, &comFoamPtr);

    comFoamPtr->step();
    
    return 0;
}

//int comFoam::flowExtractData(const char *name)
//{
//    Foam::Info << "rocFoam.extractData: Extracting flow data." << Foam::endl;
//    //  Call the flow iterator ^^^^^^^^^^^^^^^^^^
//    comFoam *comFoamPtr = NULL;
//    std::string volName = name+std::string("VOL");
//    std::string objectName = volName+std::string(".object");
//    COM_get_object(objectName.c_str(), 0, &comFoamPtr);
//    comFoamPtr->extractData();
//    
//    return 0;
//}

//int comFoam::flowRegisterVolumeData(const char *name)
//{

//    Foam::Info << "rocFoam.registerVolumeData: Extracting flow data." << Foam::endl;
//    //  Call the flow iterator ^^^^^^^^^^^^^^^^^^
//    comFoam *comFoamPtr = NULL;
//    std::string volName = name+std::string("VOL");
//    std::string objectName = volName+std::string(".object");
//    COM_get_object(objectName.c_str(), 0, &comFoamPtr);
//    comFoamPtr->registerVolumeData(name);
//    
//    return 0;
//}

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
    
    std::string volName = name+std::string("VOL");
    std::string surfName = name+std::string("SURF");

    //  Register module with COM ^^^^^^^^^^^^^^^^^^^^^^^^^^
    comFoam *comFoamPtr = NULL;

    //std::string name="ROCFOAM";
    std::string objectName = volName+std::string(".object");
    COM_get_object(objectName.c_str(), 0, &comFoamPtr);

    /// Register functions ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    /*
    std::vector<COM_Type> types(13,COM_INT);
    types[0] = COM_RAWDATA;
    types[2] = COM_VOID;
    */

    std::vector<COM_Type> types(13,COM_VOID);
    types[0] = COM_RAWDATA;
    types[1] = COM_INT;

    std::string functionName = volName+std::string(".flowInit");
    COM_set_member_function
    (
        functionName.c_str(),
        reinterpret_cast<Member_func_ptr>(&comFoam::flowInit),
        objectName.c_str(),
        "biii",
        &types[0]
    );

    functionName = volName+std::string(".flowLoop");
    COM_set_member_function
    (
        functionName.c_str(),
        reinterpret_cast<Member_func_ptr>(&comFoam::flowLoop),
        objectName.c_str(),
        "bi",
        &types[0]
    );

    functionName = volName+std::string(".flowStep");
    COM_set_member_function
    (
        functionName.c_str(),
        reinterpret_cast<Member_func_ptr>(&comFoam::flowStep),
        objectName.c_str(),
        "bi",
        &types[0]
    );

    functionName = volName+std::string(".flowReconstCaData");
    COM_set_member_function
    (
        functionName.c_str(),
        reinterpret_cast<Member_func_ptr>(&comFoam::reconstCaData),
        objectName.c_str(),
        "biii",
        &types[0]
    );

//    functionName = name+std::string(".flowRegisterDate");
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
    //    (name + std::string(".flowFin")).c_str(),
    //    reinterpret_cast<Member_func_ptr>(&rhoCentral::flowFin),
    //    objectName.c_str(), "b", &types[0]
    //);

//    //  Registering data of this module to COM ^^^^^^^^^^^^
//    std::string dataName = name+std::string(".winNProc");
//    COM_new_dataitem( dataName.c_str(), 'w', COM_INT, 1, "");
//    COM_set_size(     dataName.c_str(), 0, 1);
//    COM_set_array(    dataName.c_str(), 0, &(comFoamPtr->winNProc));

//    dataName = name+std::string(".winTime");
//    COM_new_dataitem( dataName.c_str(), 'w', COM_DOUBLE, 1, "");
//    COM_set_size(     dataName.c_str(), 0, 1);
//    COM_set_array(    dataName.c_str(), 0, &(comFoamPtr->winTime));

//    dataName = name+std::string(".winDeltaT");
//    COM_new_dataitem( dataName.c_str(), 'w', COM_DOUBLE, 1, "");
//    COM_set_size(     dataName.c_str(), 0, 1);
//    COM_set_array(    dataName.c_str(), 0, &(comFoamPtr->winDeltaT) );

//    dataName = name+std::string(".winRun");
//    COM_new_dataitem( dataName.c_str(), 'w', COM_INT, 1, "");
//    COM_set_size(     dataName.c_str(), 0, 1);
//    COM_set_array(    dataName.c_str(), 0, &(comFoamPtr->winRun));

    COM_window_init_done(volName); 

    return 0;
}
//---------------------------------------------------------

//===================================================================

#include "volumeMethods.C"
#include "faceMethods.C"
#include "surfaceMethods.C"
#include "reconstMethods.C"
#include "tempFiles.C"

comFoam::~comFoam()
{   
    deleteVolumeData();
    deleteFaceData();
    deleteSurfaceData();
    deleteFilesData();

    if (ca_runStat != NULL){
        delete [] ca_runStat;
        ca_runStat = NULL;
    }
    
    if (ca_time != NULL){
        delete [] ca_time;
        ca_time = NULL;
    }
    
    if (ca_deltaT != NULL){
        delete [] ca_deltaT;
        ca_deltaT = NULL;
    }
}

