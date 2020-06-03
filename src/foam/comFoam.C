#include "comFoam.H"
#include "cellShape.H"

comFoam::comFoam()
{}

/*
comFoam::comFoam(int *pargc, void **pargv, const char *name)
{
    flowInit(pargc, pargv, name);
}
*/


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
    if (ca_myRank == 0)
    {
        std::cout << "rocFoam.reconstCaData: Initializing CA "
                  << "reconstructions for window "
                  << name << std::endl;
    }

    //  OpenFOAM initializer ^^^^^^^^^^^^^^^^^^^^
    /*
    comFoam *comFoamPtr = nullptr;
    std::string volName = name+std::string("VOL");
    std::string surfName = name+std::string("SURF");
    std::string dataName = volName+std::string(".object");
    COM_get_object(dataName.c_str(), 0, &comFoamPtr);
    */

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


    //comFoamPtr->initialize(argc, argv);
    initialize(argc, argv);

    //deleteInitFiles(tmpFluidDir);

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
    step();
    return 0;
}

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
    
    std::string volName = name+std::string("VOL");
    std::string surfName = name+std::string("SURF");

    //  Register module with COM ^^^^^^^^^^^^^^^^^^^^^^^^^^
    comFoam *comFoamPtr = nullptr;

    //std::string name="ROCFOAM";
    std::string objectName = volName+std::string(".object");
    COM_get_object(objectName.c_str(), 0, &comFoamPtr);

    /// Register functions ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    // Stand-alone driver functions
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


    // RocStar Agent driver functions ^^^^^^^^^^^^^^^^^^^^^
    COM_Type init_types[]
    {
        COM_RAWDATA,          // G
        COM_DOUBLE_PRECISION, // initialTime
        COM_MPI_COMM,         // communicator
        COM_INTEGER,          // manInitHandle,
        COM_STRING,           // win_surf
        COM_STRING,           // win_vol
        COM_INTEGER           // obtainHandle
    };
    functionName = volName+std::string(".initialize");
    COM_set_member_function
    (
        functionName.c_str(),
        reinterpret_cast<Member_func_ptr>(&comFoam::flowInitRocStar),
        objectName.c_str(),
        "biiiiii",
        init_types
    );

    COM_Type update_types[]
    {
        COM_RAWDATA,          // G
        COM_DOUBLE_PRECISION, // currentTime
        COM_DOUBLE_PRECISION, // initTime
        COM_INTEGER           // handle1
    };
    functionName = volName+std::string(".update_solution");
    COM_set_member_function
    (
        functionName.c_str(),
        reinterpret_cast<Member_func_ptr>(&comFoam::flowUpdateRocStar),
        objectName.c_str(),
        "biii",
        update_types
    );

    COM_Type fin_types[]
    {
        COM_RAWDATA // G
    };
    functionName = volName+std::string(".finalize");
    COM_set_member_function
    (
        functionName.c_str(),
        reinterpret_cast<Member_func_ptr>(&comFoam::flowFinRocStar),
        objectName.c_str(),
        "b",
        fin_types
    );
    //-----------------------------------------------------

    COM_window_init_done(volName);

    return 0;
}

void comFoam::flowInitRocStar
(
    const double& initTime,
    const MPI_Comm& flowComm,
    const int& manInitHandle,
    const std::string& volName,
    const std::string& surfName,
    const int& obtainHandle
)
{
    std::string volTmp = "VOL";
    std::string surfTmp = "SURF";
    size_t volStart = volName.find(volTmp);
    size_t surfStart = surfName.find(surfTmp);
    
    if (volStart != surfStart)
    {
        std::cout << "WARNING: Volume and surface windows"
                  << " do not follow naming rule."
                  << std::endl;
        exit(-1);
    }
    
    std::string subType = volName.substr(0, volStart);
    char* name = const_cast<char*>(subType.c_str());
    
    /* if (initTime == 0)
    {
        if (ca_myRank == 0)
        {
            std::cout << "rocFoam.flowInitRocStar: Initializing flow solver with name "
                      << name << std::endl;
        }
        //  OpenFOAM initializer ^^^^^^^^^^^^^^^^^^^^
        int argc{3};
        if (ca_nProc>1)
            argc++;

        char** argv = char*[argc];
        argv[0] = reinterpret_cast<char*>("rocFoam");
        argv[1] = reinterpret_cast<char*>("-case");
        argv[2] = reinterpret_cast<char*>("./");
        
        if (ca_nProc>1)
            argv[3] = reinterpret_cast<char*>("-parallel");
        
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
    }
    else
    { */
        if (ca_myRank == 0)
        {
            std::cout << "rocFoam.flowInitRocStar: Initializing CA "
                      << "reconstructions for window "
                      << name << std::endl;
        }

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

        int argc{3};
        if (ca_nProc>1)
            argc++;

        char** argv = new char*[argc];
        
        //argv[0] = "rocFoam";
        //argv[1] = "-case";
        //argv[2] = tmpFluidDir; //.c_str();
        //if (ca_nProc>1)
        //    argv[3] = "-parallel";

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

        initialize(argc, argv);
        
        if (initTime != *ca_time)
        {
            std::cout << "WARNING: initTime!=ca_time, "
                      << "initTime = " << initTime
                      << ", ca_time = " << *ca_time
                      << std::endl;
        }
        //deleteInitFiles(tmpFluidDir);
    //}
}

void comFoam::flowUpdateRocStar
(
    double& currentTime,
    double& timeStep,
    int& handles
)
{

    Foam::Info << "rocFoam.flowStepRocStar: Stepping flow solver."
               << Foam::endl;
    step(&timeStep);
}


void comFoam::flowFinRocStar()
{
    if (ca_myRank == 0)
    {
        std::cout << "rocFoam.flowFinRocStar: "
                  << "Registering flow functions with name "
                  << std::endl;
    }

    finalize();
}




//---------------------------------------------------------

//===================================================================

#include "volumeMethods.C"
#include "faceMethods.C"
#include "surfaceMethods.C"
#include "reconstMethods.C"
#include "initFiles.C"

int comFoam::finalize()
{
    deleteVolumeData();
    deleteFaceData();
    deleteSurfaceData();
    deleteFilesData();

    if (ca_runStat != nullptr) {delete ca_runStat; ca_runStat = nullptr;}
    if (ca_time != nullptr) {delete ca_time; ca_time = nullptr;}
    if (ca_deltaT != nullptr) {delete ca_deltaT; ca_deltaT = nullptr;}
    if (ca_deltaT0 != nullptr) {delete ca_deltaT0; ca_deltaT0 = nullptr;}
    if (ca_timeIndex != nullptr) {delete ca_timeIndex; ca_timeIndex = nullptr;}
    
    if (ca_timeName != nullptr)
    {
        delete [] ca_timeName;
        ca_timeName = nullptr;
    }
    
    return 0;
}



comFoam::~comFoam()
{
    finalize();
}


