#include "comFoam.H"
#include "cellShape.H"

comFoam::comFoam()
{}

// Collective calls ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
int comFoam::createCSCdata()
{
    createStatusData();
    createVolumeConnectivities();
    createVolumeData();
    createFaceConnectivities();
    createFaceData();
    createSurfaceConnectivities();
    createSurfaceData();
    deleteTempFiles(tmpFluidDir);
    std::string strTmp = "./";
    readFilesData(strTmp);

    return 0;
}

int comFoam::updateCSCdata()
{
    updateStatusData();

    if (ca_nCells != nullptr)
        updateVolumeData_outgoing();
    if (ca_nFaces != nullptr)
        updateFaceData_outgoing();
    if (ca_nPatches != nullptr)
        updateSurfaceData_outgoing();

    return 0;
}

int comFoam::registerCSCdata(const char *name)
{

    registerStatusData(name);
    registerFilesData(name);
    registerVolumeData(name);
    registerFaceData(name);
    registerSurfaceData(name);

    return 0;
}

int comFoam::deleteCSCdata()
{
    deleteFilesData();
    deleteSurfaceData();
    deleteFaceData();
    deleteVolumeData();
    deleteStatusData();

    return 0;
}

int comFoam::reconstCSCdata(const char *name)
{
    deleteTempFiles(tmpFluidDir);
    for (int iproc=0; iproc<ca_nProc; iproc++)
    {
        MPI_Barrier(winComm);
        if(iproc!=ca_myRank)
            continue;

        reconstStatusData(name);
        reconstVolumeData(name);
        reconstFaceData(name);
        reconstSurfaceData(name);
        reconstFilesData(name);
    }
    MPI_Barrier(winComm);

    return 0;
}
//-----------------------------------------------

//^^^ DEFINITION OF COM-RELATED MTHODS ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
int comFoam::flowInit(int *pargc, void **pargv, const char *name)
{
    if (ca_myRank == 0)
    {
        std::cout << "rocFoam.flowInit: Initializing flow solver with name "
                  << name << std::endl;
    }
    //  OpenFOAM initializer ^^^^^^^^^^^^^^^^^^^^
    comFoam *comFoamPtr = nullptr;

    std::string winName = name;
    std::string objectName = winName+std::string(".object");
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
    
    initFOAM(argc, argv);
    //  Other initializations ^^^^^^^^^^^^^^^^^^^
    createCSCdata();
    updateCSCdata();
    registerCSCdata(name);

    if (argv != nullptr)
    {
        for (int i=0; i<argc; i++)
        {
            if (argv[i] != nullptr)
            {
                delete [] argv[i];
                argv[i] = nullptr;
            }
        }
        delete [] argv;
        argv = nullptr;
    }

    return 0;
}


int comFoam::restartInit(int *pargc, void **pargv, const char *name)

{
    if (ca_myRank == 0)
    {
        std::cout << "rocFoam.restartInit: Initializing CSC "
                  << "reconstructions for window "
                  << name << std::endl;
    }

    reconstCSCdata(name);

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

    //comFoamPtr->initFOAM(argc, argv);
    initFOAM(argc, argv);

    //deleteInitFiles(tmpFluidDir);

    if (argv != nullptr)
    {
        for (int i=0; i<argc; i++)
        {
            if (argv[i] != nullptr)
            {
                delete [] argv[i];
                argv[i] = nullptr;
            }
        }
        delete [] argv;
        argv = nullptr;
    }

    return 0;
}

int comFoam::flowLoop()
{
    Foam::Info << "rocFoam.flowLoop: Iterating flow solver."
               << Foam::endl;
    loop();
    return 0;
}

int comFoam::flowStep()
{
    Foam::Info << "rocFoam.flowStep: Stepping flow solver."
               << Foam::endl;

    updateSurfaceData_incoming();
    step();
    updateCSCdata();

    return 0;
}

// Rocstar Agent methods ^^^^^^^^^^^^^^^^^^^^^^^^
void comFoam::initialize
(
    const double *initTime,
    const MPI_Comm* flowComm,
    const int* manInitHandle,
    const char* surfName_,
    const char* volName_,
    const int* obtainHandle
)
{
    std::string volName = std::string(volName_);
    std::string surfName = std::string(surfName_);

    
    /*
    std::string volTmp = "_vol";
    std::string surfTmp = "_srf";
    size_t volStart = volName.find(volTmp, 0);
    size_t surfStart = surfName.find(surfTmp, 0);


    
    if (volStart == std::string::npos ||
        surfStart == std::string::npos)
    {
        std::cout << "WARNING: The name of input windows are not"
                  << " consistent with what required names."
                  << std::endl;
        exit(-1);    
    }
    
    if (volStart != surfStart)
    {
        std::cout << "WARNING: temp volume and surface windows"
                  << " do not follow naming rule."
                  << std::endl;
        exit(-1);
    }
        
    std::string subType = winName; //volName.substr(0, volStart);
    char* name = const_cast<char*>(subType.c_str());
    
    if (ca_myRank == 0)
    {
        std::cout << "rocFoam.initialize: Initializing "
                  << "reconstructions of windows for "
                  << name << std::endl;

        std::cout << "ManInitHandle is "
                  << std::string((*manInitHandle < 0) ? ("not set") : ("set"));
        std::cout << "ObtainHandle is "
                  << std::string((*obtainHandle < 0) ? ("not set") : ("set"));
    }
    */

    //loadInternal(name);

    std::string newVolName = winName+"VOL";
    std::string newSurfName = winName+"SURF";

    copyWindow(volName.c_str(),
               newVolName.c_str());
    copyWindow(surfName.c_str(),
               newSurfName.c_str());

    reconstCSCdata(winName.c_str());

    int argc{3};
    if (ca_nProc>1) argc++;

    char** argv = new char*[argc];
    
    std::string strTmp = "-rocFoam";
    argv[0] = new char[strTmp.length()+1];
    std::strcpy(argv[0], strTmp.c_str());

    strTmp = "-case";
    argv[1] = new char[strTmp.length()+1];
    std::strcpy(argv[1], strTmp.c_str());

    strTmp = tmpFluidDir;
    argv[2] = new char[strTmp.length()+1];
    std::strcpy(argv[2], strTmp.c_str());

    if (ca_nProc>1)
    {
        strTmp = "-parallel";
        argv[3] = new char[strTmp.length()+1];
        std::strcpy(argv[3], strTmp.c_str());
    }

    initFOAM(argc, argv);
    
    if (*initTime != *ca_time)
    {
        std::cout << "WARNING: initTime!=ca_time, "
                  << "initTime = " << *initTime
                  << ", ca_time = " << *ca_time
                  << std::endl;
    }

    if (argv != nullptr)
    {
        for (int i=0; i<argc; i++)
        {
            if (argv[i] != nullptr)
            {
                delete [] argv[i];
                argv[i] = nullptr;
            }
        }
        delete [] argv;
        argv = nullptr;
    }

    //if (*manInitHandle > 0)
    //    COM_call_function(*manInitHandle,
    //                      newSurfName.c_str(),
    //                      newVolName.c_str());
}

void comFoam::update_solution
(
    double* currentTime,
    double* timeStep,
    int* bcHandle,
    int* gmHandle
)
{
    Info << "rocFoam.flowStepRocStar: Stepping flow solver."
         << endl;
         
    Info << "  bcHandle is "
         << std::string((*bcHandle < 0) ? ("not set") : ("set"))
         << endl;

    Info << "  gmHandle is "
         << std::string((*gmHandle < 0) ? ("not set") : ("set"))
         << endl;

    if (*currentTime != *ca_time)
        Info << "  Flow solver time and the input time"
             << " are not the same " << *ca_time 
             << " vs " << *currentTime << endl;

    if (*gmHandle>=0)
    {
        double alpha{1};
        COM_call_function(*gmHandle, &alpha);
    }
        
    updateSurfaceData_incoming();
    step(timeStep);
    updateCSCdata();
}

void comFoam::finalize()
{
    Info << "rocFoam.finalize: "
         << "Finalizing flow solver."
         << endl;

    finalizeFoam();
}
//-----------------------------------------------



//^^^^^ REGISTER FUNCTIONS ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
int comFoam::registerFunctions(const char *name)
{
    if (ca_myRank == 0)
    {
        std::cout << "rocFoam.flowRegister: "
                  << "Registering flow functions with name "
                  << name
                  << std::endl;
    }
    
    std::string winName = name;

    //  Register module with COM ^^^^^^^^^^^^^^^^^^^^^^^^^^
    comFoam *comFoamPtr = nullptr;

    //std::string name="ROCFOAM";
    std::string objectName = winName+std::string(".object");
    COM_get_object(objectName.c_str(), 0, &comFoamPtr);

    /// Register functions ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    // Stand-alone driver functions
    COM_Type types[13]={COM_VOID};
    types[0] = COM_RAWDATA;
    types[1] = COM_INT;

    std::string functionName = winName+std::string(".flowInit");
    COM_set_member_function
    (
        functionName.c_str(),
        reinterpret_cast<Member_func_ptr>(&comFoam::flowInit),
        objectName.c_str(),
        "biii",
        types
    );

    functionName = winName+std::string(".flowLoop");
    COM_set_member_function
    (
        functionName.c_str(),
        reinterpret_cast<Member_func_ptr>(&comFoam::flowLoop),
        objectName.c_str(),
        "b",
        types
    );

    functionName = winName+std::string(".flowStep");
    COM_set_member_function
    (
        functionName.c_str(),
        reinterpret_cast<Member_func_ptr>(&comFoam::flowStep),
        objectName.c_str(),
        "b",
        types
    );


    //types[0] = COM_RAWDATA;
    //types[1] = COM_VOID;
    functionName = winName+std::string(".flowRestartInit");
    COM_set_member_function
    (
        functionName.c_str(),
        reinterpret_cast<Member_func_ptr>(&comFoam::restartInit),
        objectName.c_str(),
        "biii",
        types
    );


    // RocStar Agent driver functions ^^^^^^^^^^^^^^^^^^^^^
    COM_Type init_types[]
    {
        COM_RAWDATA,      // G
        COM_DOUBLE,       // initialTime
        COM_MPI_COMM,     // communicator
        COM_INT,          // manInitHandle,
        COM_STRING,       // win_surf
        COM_STRING,       // win_vol
        COM_INT           // obtainHandle
    };
    functionName = winName+std::string(".initialize");
    COM_set_member_function
    (
        functionName.c_str(),
        reinterpret_cast<Member_func_ptr>(&comFoam::initialize),
        objectName.c_str(),
        "biiiiii",
        init_types
    );

    COM_Type update_types[]
    {
        COM_RAWDATA,    // G
        COM_DOUBLE,     // currentTime
        COM_DOUBLE,     // initTime
        COM_INT,        // handle1
        COM_INT         // handle2
    };
    functionName = winName+std::string(".update_solution");
    COM_set_member_function
    (
        functionName.c_str(),
        reinterpret_cast<Member_func_ptr>(&comFoam::update_solution),
        objectName.c_str(),
        "biiii",
        update_types
    );

    COM_Type fin_types[]
    {
        COM_RAWDATA // G
    };
    functionName = winName+std::string(".finalize");
    COM_set_member_function
    (
        functionName.c_str(),
        reinterpret_cast<Member_func_ptr>(&comFoam::finalize),
        objectName.c_str(),
        "b",
        fin_types
    );
    //-----------------------------------------------------

    COM_window_init_done(winName);

    return 0;
}

//---------------------------------------------------------

//===================================================================

#include "statusMethods.C"
#include "volumeMethods.C"
#include "faceMethods.C"
#include "surfaceMethods.C"
#include "reconstMethods.C"
#include "filesMethods.C"

int comFoam::finalizeFoam()
{
    deleteCSCdata();
    
    return 0;
}

comFoam::~comFoam()
{
    finalizeFoam();
}


