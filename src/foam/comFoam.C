#include "comFoam.H"
#include "cellShape.H"

comFoam::comFoam()
    : ca_nPoints(0),
      ca_nFaces(0),
      ca_nCells(0),
      ca_nPatches(0),
      winVolName(""),
      winSurfName(""),
      solverType(""),
      winComm(NULL),
      winNProc(0),
      winRank(0),
      winTime(0.0),
      winDeltaT(0.0),
      winRun(1)
{};

comFoam::comFoam(int *pargc, void **pargv, const char *name)
    : ca_nPoints(0),
      ca_nFaces(0),
      ca_nCells(0),
      ca_nPatches(0),
      winVolName(""),
      winSurfName(""),
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
    comFoamPtr->registerVolumeData(name);

    /*
    comFoamPtr->createFaceConnectivities();
    comFoamPtr->createFaceData();
    comFoamPtr->updateFaceData();
    comFoamPtr->registerSurfaceData(name);
    */

    comFoamPtr->createSurfaceConnectivities();
    comFoamPtr->createSurfaceData();
    comFoamPtr->updateSurfaceData();
    comFoamPtr->registerSurfaceData(name);

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

//int comFoam::flowRegisterVolumeData(const char *name)
//{

//    Foam::Info << "rocFoam.registerVolumeData: Extracting flow data." << Foam::endl;
//    //  Call the flow iterator ^^^^^^^^^^^^^^^^^^
//    comFoam *comFoamPtr = NULL;
//    std::string volName = name+string("VOL");
//    std::string objectName = volName+string(".object");
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

#   include "volumeMethods.H"
#   include "faceMethods.H"
#   include "surfaceMethods.H"
