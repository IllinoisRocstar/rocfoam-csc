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

    if (false)
    {
        for (int count=0; count<3; count++)
        {
            std::string volName;
            if (count == 0)
            {
                volName = name;
            }
            else if (count ==1)
            {
                volName = name+std::string("VOL");
            }
            else if (count ==2)
            {
                volName = name+std::string("SURF");
            }

            std::string regNames;
            int numDataItems=0;
            COM_get_dataitems(volName.c_str(), &numDataItems, regNames);
            //std::cout << "  numDataItems = " << numDataItems << std::endl;

            std::vector<std::string> dataItemNames;
            dataItemNames.clear();
            std::istringstream Istr(regNames);
            for (int i=0; i<numDataItems; ++i)
            {
                std::string nameTmp;
                Istr >> nameTmp;
                dataItemNames.push_back(nameTmp);
            }

            dataItemNames.push_back("nc");
            numDataItems++;


            int nPanes;
            int* paneList;
            COM_get_panes(volName.c_str(), &nPanes, &paneList);
            std::cout << "  Number of Panes = " << nPanes << std::endl;

            for (int ipane=0; ipane<=nPanes; ipane++)
            {
                int paneID;
                if ( ipane == 0 )
                {
                    paneID = 0;
                }
                else
                {
                    paneID = paneList[ipane-1];
                }
                
                std::cout << "WindoName = " << volName
                          << ", Pane[" << ipane << "], paneID = " << paneID
                          << " ^^^^^^^^^^^^^^^^^^^^^^^" << std::endl;

                for (int i=0; i<numDataItems; i++)
                {
                
                    std::string nameTmp = dataItemNames[i];

                    char loc[50];
                    COM_Type type;
                    int ncomp;
                    std::string unit;
                    int size{0};


                    std::string regName = volName+std::string(".")+nameTmp;
                    if (ipane == 0)
                    {
                        COM_get_dataitem(regName.c_str(), loc, &type, &ncomp, &unit);
                    }
                    else
                    {
                        COM_get_dataitem(regName.c_str(), loc, &type, &ncomp, &unit);
                    }

                    if 
                    (
                        (*loc == 'w' && ipane == 0) ||
                        (*loc == 'p' && ipane != 0)
                    )
                    {
                        COM_get_size(regName.c_str(), paneID, &size);
                        if (size>0)
                        {
                            std::cout << "  DataItem[" << i << "] = " << nameTmp
                                      << ", loc = " << *loc
                                      << ", type = " << type
                                      << ", paneID = " << paneID
                                      << ", size = " << size
                                      << ", ncomp = " << ncomp
                                      << ", unit = " << unit
                                      << std::endl;
                        }
                    }
                }
                
                if (ipane != 0)
                {
                    std::string  dataName = std::string("nc");
                    std::string regName = volName+std::string(".")+dataName;
                    int nPoints;
                    int nComp;
                    COM_get_array(regName.c_str(), paneID, &ca_Points, &nComp);
                    COM_get_size(regName.c_str(), paneID, &nPoints);
                    std::cout << "  " << dataName.c_str() << " nPoints = " << nPoints
                              << ", components = " << nComp << std::endl;                

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

                        dataName = volName+std::string(".")+connName;
                        //nameExists(dataItemNames, dataName);
                        COM_get_array(dataName.c_str(), paneID, &ca_cellToPointConn[icon], &nComp);
                        COM_get_size(dataName.c_str(), paneID, &numElem);

                        std::cout << "    Connectivity[" << icon << "] = " << connName
                                  << ", elements = " << numElem
                                  << ", components =" << nComp << std::endl;
                    }
                }
            }
        }
    }

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
        if(iproc==ca_myRank)
        {
            reconstStatusData(name);
            reconstVolumeData(name);
            reconstFaceData(name);
            reconstSurfaceData(name);
            reconstFilesData(name);
        }
        MPI_Barrier(winComm);
    }

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
    step();
    updateCSCdata();
    return 0;
}

// Rocstar Agent methods ^^^^^^^^^^^^^^^^^^^^^^^^
void comFoam::initialize
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
    
    if (ca_myRank == 0)
    {
        std::cout << "rocFoam.initialize: Initializing "
                  << "reconstructions of windows for "
                  << name << std::endl;
    }

    reconstCSCdata(name);

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
    
    if (initTime != *ca_time)
    {
        std::cout << "WARNING: initTime!=ca_time, "
                  << "initTime = " << initTime
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
}

void comFoam::update_solution
(
    double& currentTime,
    double& timeStep,
    int& handles
)
{

    Foam::Info << "rocFoam.flowStepRocStar: Stepping flow solver."
               << Foam::endl;
    step(&timeStep);
    updateCSCdata();
}

void comFoam::finalize()
{
    if (ca_myRank == 0)
    {
        std::cout << "rocFoam.finalize: "
                  << "Finalizing flow solver."
                  << std::endl;
    }

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
        COM_RAWDATA,          // G
        COM_DOUBLE_PRECISION, // initialTime
        COM_MPI_COMM,         // communicator
        COM_INTEGER,          // manInitHandle,
        COM_STRING,           // win_surf
        COM_STRING,           // win_vol
        COM_INTEGER           // obtainHandle
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
        COM_RAWDATA,          // G
        COM_DOUBLE_PRECISION, // currentTime
        COM_DOUBLE_PRECISION, // initTime
        COM_INTEGER           // handle1
    };
    functionName = winName+std::string(".update_solution");
    COM_set_member_function
    (
        functionName.c_str(),
        reinterpret_cast<Member_func_ptr>(&comFoam::update_solution),
        objectName.c_str(),
        "biii",
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


