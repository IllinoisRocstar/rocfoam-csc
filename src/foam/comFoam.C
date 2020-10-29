#include "comFoam.H"

comFoam::comFoam()
{}

// Collective calls ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
int comFoam::createCSCdata()
{
    createStatusData();
    updateStatusData();
    createVolumeConnectivities();
    createVolumeData();
    createFaceConnectivities();
    createFaceData();
    createSurfaceConnectivities();
    createSurfaceData();

    std::string strTmp = "./";
    readFilesData(strTmp);
    deleteTempFiles(tmpFluidDir);
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
    MPI_Barrier(winComm);
    //for (int iproc=0; iproc<ca_nProc; iproc++)
    //{
    //    MPI_Barrier(winComm);
    //    if(iproc!=ca_myRank)
    //        continue;

        reconstStatusData(name);
        reconstVolumeData(name);
        reconstFaceData(name);
        reconstSurfaceData(name);
        reconstFilesData(name);
    //}
    //MPI_Barrier(winComm);

    return 0;
}
//-----------------------------------------------

//^^^ DEFINITION OF COM-RELATED MTHODS ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
int comFoam::flowInit(int *pargc, void **pargv, const char *name)
{
    Foam::Info << solverType.c_str() << ": Initializing flow solver." << endl;
    //  OpenFOAM initializer ^^^^^^^^^^^^^^^^^^^^
    comFoam *comFoamPtr = nullptr;

    std::string winName_ = name;
    std::string objectName = winName_+std::string(".object");
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
    Info << solverType.c_str() << ": Initializing CSC restart." << endl;

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
    Info << solverType.c_str() << ": Looping flow solver." << endl;

    loop();
    return 0;
}

int comFoam::flowStep()
{
    Info << solverType.c_str() << ": Stepping flow solver." << endl;

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
    Info << solverType.c_str() << ": Initializing Rocstar restart." << endl;

    std::string volName = std::string(volName_);
    std::string surfName = std::string(surfName_);

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
        WarningInFunction
            << "WARNING: initTime!=ca_time, "
            << "initTime = " << *initTime
            << ", ca_time = " << *ca_time
            << endl;
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

    if (*manInitHandle > 0)
        COM_call_function(*manInitHandle,
                          newSurfName.c_str(),
                          newVolName.c_str());
}

void comFoam::update_solution
(
    double* currentTime,
    double* timeStep,
    int* bcHandle,
    int* gmHandle
)
{
    Info << solverType.c_str() << ": Stepping flow solver." << endl;
         
    std::stringstream output{};
    output << "  bcHandle is "
         << std::string((*bcHandle < 0) ? ("not set") : ("set"))
         << std::endl
         << "  gmHandle is "
         << std::string((*gmHandle < 0) ? ("not set") : ("set"));
    verbose_message(output.str());


    if (*currentTime != *ca_time)
    {
        std::ostringstream doubleToOs;
        doubleToOs << std::scientific 
                   << std::setprecision(IODigits);
        doubleToOs << std::abs(*ca_time - *currentTime);

        output = std::stringstream{};
        output << "  Flow solver time and the input time"
             << " are not the same " << *ca_time 
             << " vs " << *currentTime
             << ", diff = " << doubleToOs.str();
        verbose_message(output.str());
    }

    step(timeStep, gmHandle);
    updateCSCdata();
}

void comFoam::finalize()
{
    Info << solverType.c_str() << ": Finalizing flow solver." << endl;

    finalizeFoam();
}
//-----------------------------------------------

void comFoam::message(std::string message, bool parallel)
{
    if (!parallel)
    {
        if (ca_myRank == 0)
        {
            std::cout << message << std::endl;
        }
    }
    else
    {
        std::cout << message << std::endl;
    }
}

void comFoam::verbose_message(std::string message, bool parallel)
{
#ifdef VERBOSE
    if (!parallel)
    {
        if (ca_myRank == 0)
        {
            std::cout << message << std::endl;
        }
    }
    else
    {
        std::cout << message << std::endl;
    }
#endif
}

//^^^^^ REGISTER FUNCTIONS ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
int comFoam::registerFunctions(const char *name)
{
    Info << solverType.c_str() << ": Registering solver methods." << endl;

    std::string winName_ = name;

    //  Register module with COM ^^^^^^^^^^^^^^^^^^^^^^^^^^
    comFoam *comFoamPtr = nullptr;

    //std::string name="ROCFOAM";
    std::string objectName = winName_+std::string(".object");
    COM_get_object(objectName.c_str(), 0, &comFoamPtr);

    /// Register functions ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    // Stand-alone driver functions
    COM_Type types[13]={COM_VOID};
    types[0] = COM_RAWDATA;
    types[1] = COM_INT;

    std::string functionName = winName_+std::string(".flowInit");
    COM_set_member_function
    (
        functionName.c_str(),
        reinterpret_cast<Member_func_ptr>(&comFoam::flowInit),
        objectName.c_str(),
        "biii",
        types
    );

    functionName = winName_+std::string(".flowLoop");
    COM_set_member_function
    (
        functionName.c_str(),
        reinterpret_cast<Member_func_ptr>(&comFoam::flowLoop),
        objectName.c_str(),
        "b",
        types
    );

    functionName = winName_+std::string(".flowStep");
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
    functionName = winName_+std::string(".flowRestartInit");
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
    functionName = winName_+std::string(".initialize");
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
    functionName = winName_+std::string(".update_solution");
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
    functionName = winName_+std::string(".finalize");
    COM_set_member_function
    (
        functionName.c_str(),
        reinterpret_cast<Member_func_ptr>(&comFoam::finalize),
        objectName.c_str(),
        "b",
        fin_types
    );
    //-----------------------------------------------------

    COM_window_init_done(winName_);

    return 0;
}

//---------------------------------------------------------

//===================================================================

int comFoam::finalizeFoam()
{
    deleteCSCdata();
    
    return 0;
}

comFoam::~comFoam()
{
    finalizeFoam();
}


