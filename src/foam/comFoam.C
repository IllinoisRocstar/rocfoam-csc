#include "comFoam.H"
#include "cellShape.H"

comFoam::comFoam()
{
    //initSet();
};

comFoam::comFoam(int *pargc, void **pargv, const char *name)
{
    //initSet();
    flowInit(pargc, pargv, name);
}

/*int comFoam::initSet()
{
    // Variables ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ca_nPoints = nullptr;  //single value
    ca_nCells = nullptr;   //single value
    ca_nFaces = nullptr;   //single value
    ca_nPatches = nullptr; //single value
    //-------------------------------------------

    // COM Volume Arrays^^^^^^^^^^^^^^^^^^^^^^^^^^
    // Mapping
    ca_cellToCellMap = nullptr;
    ca_cellToCellMap_inverse = nullptr;

    // Connectivity
    ca_cellToPointConn_types = nullptr; //single value
    ca_cellToPointConn_map = nullptr;
    ca_cellToPointConn_size = nullptr;
    ca_cellToPointConn = nullptr;

    // Field Data
    ca_Points = nullptr;
    ca_Vel = nullptr;
    ca_P = nullptr;
    ca_T = nullptr;
    ca_Rho = nullptr;
    ca_patchSf = nullptr;
    //-------------------------------------------
    
    // COM Face Arrays^^^^^^^^^^^^^^^^^^^^^^^^^^^
    // Mapping
    ca_faceToFaceMap = nullptr;
    ca_faceToFaceMap_inverse = nullptr;

    // Connectivity
    ca_faceToPointConn_types = nullptr; //single value
    ca_faceToPointConn_map = nullptr;
    ca_faceToPointConn_size = nullptr;
    ca_faceToPointConn = nullptr;

    // Field data
    ca_faceOwner = nullptr;
    ca_faceNeighb = nullptr;
    //-------------------------------------------

    // COM Patch Arrays^^^^^^^^^^^^^^^^^^^^^^^^^^
    // General data

    patchNameStr = nullptr; //single value for the last
    patchTypeStr = nullptr; //single value for the last

    ca_patchName = nullptr; //single value for the last
    ca_patchType = nullptr; //single value for the last
    ca_patchInGroup = nullptr; //single value for the last
    ca_patchStart = nullptr; //single value for the last
    ca_patchSize = nullptr; //single value for the last

    // PointToPoint Mapping
    ca_patchPointToPointMap_size = nullptr; //single value for the last
    ca_patchPointToPointMap = nullptr;

    // FaceToFace Mapping
    ca_patchFaceToFaceMap = nullptr;
    ca_patchFaceToFaceMap_inverse = nullptr;

    // FaceToPoint Mapping
    ca_patchFaceToPointConn_types = nullptr; //single value for the last
    ca_patchFaceToPointConn_map = nullptr;
    ca_patchFaceToPointConn_size = nullptr;
    ca_patchFaceToPointConn = nullptr;

    // Field data
    ca_patchPoints = nullptr;
    ca_patchVel = nullptr;
    ca_patchP = nullptr;
    ca_patchT = nullptr;
    ca_patchRho = nullptr;

    // File data
    ca_nFiles = nullptr;
    ca_fileSize = nullptr;
    ca_fileName = nullptr;
    ca_filePath = nullptr;
    ca_fileContent = nullptr;

    //-------------------------------------------

    //  Window data ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    //std::string winVolName; /// Tracks *this* volume window name.
    //std::string winSurfName; /// Tracks *this* volume window name.

    solverType = "";
    winComm = nullptr;
    
    ca_nProc = 1;
    ca_myRank = 0;

    // registered data set during the simulation
    ca_runStat = nullptr;
    ca_time = nullptr;
    ca_deltaT = nullptr;

    ca_timeIndex = nullptr;
    ca_timeName = nullptr;
    ca_deltaT0 = nullptr;

    return 0;
} */



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
    comFoam *comFoamPtr = nullptr;

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
    
    initialize(argc, argv);
    //  Other initializations ^^^^^^^^^^^^^^^^^^^
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

    //std::string tmpDir = strTmp+tmpFluidDir;
    deleteInitFiles(tmpFluidDir);
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
                  << "reconstructions for window "
                  << name << std::endl;
    }

    //  OpenFOAM initializer ^^^^^^^^^^^^^^^^^^^^
    comFoam *comFoamPtr = nullptr;
    std::string volName = name+std::string("VOL");
    std::string surfName = name+std::string("SURF");
    std::string dataName = volName+std::string(".object");
    COM_get_object(dataName.c_str(), 0, &comFoamPtr);

    deleteInitFiles(tmpFluidDir);
    for(int iproc=0; iproc<ca_nProc; iproc++)
    {
        if(iproc==ca_myRank)
        {
            reconstCaVolumeData(name);
            reconstCaFaceData(name);
            reconstCaSurfaceData(name);
            reconstCaInitFiles(name, tmpFluidDir);
        }
        MPI_Barrier(winComm);
    }

    //int argc = *pargc;
    //char** argv = reinterpret_cast<char**>(pargv);
//    int argc = *pargc+2;
//    char** argv;
//    argv = new char*[*pargc+2];
//        
//    std::string strTmp = "-case";
//    argv[*pargc] = new char[strTmp.length()+1];
//    std::strcpy(argv[*pargc], strTmp.c_str());

//    strTmp = tmpFluidDir;
//    argv[*pargc+1] = new char[strTmp.length()+1];
//    std::strcpy(argv[*pargc+1], strTmp.c_str());

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

    strTmp = tmpFluidDir;
    argv[*pargc+1] = new char[strTmp.length()+1];
    std::strcpy(argv[*pargc+1], strTmp.c_str());


    comFoamPtr->initialize(argc, argv);

    //deleteInitFiles(tmpFluidDir);

    return 0;
}

int comFoam::flowLoop(const char *name)
{
    Foam::Info << "rocFoam.flowLoop: Iterating flow solver." << Foam::endl;

    //  Call the flow iterator ^^^^^^^^^^^^^^^^^^
    comFoam *comFoamPtr = nullptr;
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
    comFoam *comFoamPtr = nullptr;
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
//    comFoam *comFoamPtr = nullptr;
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
//    comFoam *comFoamPtr = nullptr;
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
    comFoam *comFoamPtr = nullptr;

    //std::string name="ROCFOAM";
    std::string objectName = volName+std::string(".object");
    COM_get_object(objectName.c_str(), 0, &comFoamPtr);

    /// Register functions ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    /*
    std::vector<COM_Type> types(13,COM_INT);
    types[0] = COM_RAWDATA;
    types[2] = COM_VOID;
    */

    COM_Type types[13]={COM_VOID};
    types[0] = COM_RAWDATA;
    types[1] = COM_INT;

    std::string functionName = volName+std::string(".flowInit");
    COM_set_member_function
    (
        functionName.c_str(),
        reinterpret_cast<Member_func_ptr>(&comFoam::flowInit),
        objectName.c_str(),
        "biii",
        types
    );

    functionName = volName+std::string(".flowLoop");
    COM_set_member_function
    (
        functionName.c_str(),
        reinterpret_cast<Member_func_ptr>(&comFoam::flowLoop),
        objectName.c_str(),
        "bi",
        types
    );

    functionName = volName+std::string(".flowStep");
    COM_set_member_function
    (
        functionName.c_str(),
        reinterpret_cast<Member_func_ptr>(&comFoam::flowStep),
        objectName.c_str(),
        "bi",
        types
    );

    functionName = volName+std::string(".flowReconstCaData");
    COM_set_member_function
    (
        functionName.c_str(),
        reinterpret_cast<Member_func_ptr>(&comFoam::reconstCaData),
        objectName.c_str(),
        "biii",
        types
    );

    COM_window_init_done(volName);

    return 0;
}
//---------------------------------------------------------

//===================================================================

#include "volumeMethods.C"
#include "faceMethods.C"
#include "surfaceMethods.C"
#include "reconstMethods.C"
#include "initFiles.C"

comFoam::~comFoam()
{   
    deleteVolumeData();
    deleteFaceData();
    deleteSurfaceData();
    deleteFilesData();

    if (ca_runStat != nullptr){delete ca_runStat; ca_runStat = nullptr;}
    if (ca_time != nullptr){delete ca_time; ca_time = nullptr;}
    if (ca_deltaT != nullptr){delete ca_deltaT; ca_deltaT = nullptr;}
    if (ca_deltaT0 != nullptr){delete ca_deltaT0; ca_deltaT0 = nullptr;}
    if (ca_timeIndex != nullptr){delete ca_timeIndex; ca_timeIndex = nullptr;}
    
    if (ca_timeName != nullptr)
    {
        delete [] ca_timeName;
        ca_timeName = nullptr;
    }
}

