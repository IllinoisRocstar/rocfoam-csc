#include "comFoam.H"

comFoam::comFoam()
    : winNameVol(""),
      solverType(""),
      winComm(NULL),
      winNProc(0),
      winRank(0),
      winTime(0.0),
      winDeltaT(0.0),
      winRun(1)
{};

comFoam::comFoam(int *pargc, void **pargv, int *verbIn)
    : winNameVol(""),
      solverType(""),
      winComm(NULL),
      winNProc(0),
      winRank(0),
      winTime(0.0),
      winDeltaT(0.0),
      winRun(1)
{
    flowInit(pargc, pargv, verbIn);
}

//^^^ DEFINITION OF COM-RELATED MTHODS ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
int comFoam::flowInit(int *pargc, void **pargv, int *verbIn)
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
    std::string name="ROCFOAM";
    std::string globalName(name+".global");
    COM_get_object(globalName.c_str(), 0, &comFoamPtr);

    int argc = *pargc;
    char** argv = reinterpret_cast<char**>(pargv);

    comFoamPtr->initialize(argc, argv);
    
    //  Other initializations ^^^^^^^^^^^^^^^^^^^

    return 0;
}


int comFoam::flowLoop()
{

    Foam::Info << "rocFoam.flowLoop: Iterating flow solver." << Foam::endl;

    //  Call the flow iterator ^^^^^^^^^^^^^^^^^^
    comFoam *comFoamPtr = NULL;

    std::string name="ROCFOAM";
    std::string globalName(name+".global");
    COM_get_object(globalName.c_str(), 0, &comFoamPtr);

    comFoamPtr->loop();
    
    return 0;
}

int comFoam::flowStep()
{

    Foam::Info << "rocFoam.flowStep: Stepping flow solver." << Foam::endl;

    //  Call the flow iterator ^^^^^^^^^^^^^^^^^^
    comFoam *comFoamPtr = NULL;

    std::string name="ROCFOAM";
    std::string globalName(name+".global");
    COM_get_object(globalName.c_str(), 0, &comFoamPtr);

    comFoamPtr->step();
    
    return 0;
}

//^^^^^ REGISTER FUNCTIONS ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
int comFoam::flowRegister()
{
    //  Anouncing default communicator  ^^^^^^^^^^^^^^^^^^^
    MPI_Comm tmpComm = COM_get_default_communicator();  

    int tmpRank;
    MPI_Comm_rank(tmpComm, &tmpRank);
    
    if (tmpRank == 0)
    {
        std::cout << "rocFoam.flowRegister: Registering flow functions."
                  << std::endl << std::endl;
    }
    
    //  Register module with COM ^^^^^^^^^^^^^^^^^^^^^^^^^^
    comFoam *comFoamPtr = NULL;

    std::string name="ROCFOAM";
    std::string globalName(name+".global");
    COM_get_object(globalName.c_str(), 0, &comFoamPtr);

    /// Register functions ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    std::vector<COM_Type> types(13,COM_INT);

    types[0] = COM_RAWDATA;
    types[2] = COM_VOID;

    COM_set_member_function
    (
        (name + string(".flowInit")).c_str(),
        reinterpret_cast<Member_func_ptr>(&comFoam::flowInit),
        globalName.c_str(), "biii", &types[0]
    );


    COM_set_member_function
    (
        (name + string(".flowLoop")).c_str(),
        reinterpret_cast<Member_func_ptr>(&comFoam::flowLoop),
        globalName.c_str(), "b", &types[0]
    );

    COM_set_member_function
    (
        (name + string(".flowStep")).c_str(),
        reinterpret_cast<Member_func_ptr>(&comFoam::flowStep),
        globalName.c_str(), "b", &types[0]
    );

    //COM_set_member_function
    //(
    //    (name + string(".flowFin")).c_str(),
    //    reinterpret_cast<Member_func_ptr>(&rhoCentral::flowFin),
    //    globalName.c_str(), "b", &types[0]
    //);

    //  Registering data of this module to COM ^^^^^^^^^^^^
    std::string dataName="";

    dataName = name+string(".winNProc");
    COM_new_dataitem( dataName.c_str(), 'w', COM_INT, 1, "");
    COM_set_size(     dataName.c_str(), 0, 1);
    COM_set_array(    dataName.c_str(), 0, &(comFoamPtr->winNProc));

    dataName = name+string(".winTime");
    COM_new_dataitem( dataName.c_str(), 'w', COM_DOUBLE, 1, "");
    COM_set_size(     dataName.c_str(), 0, 1);
    COM_set_array(    dataName.c_str(), 0, &(comFoamPtr->winTime));

    dataName = name+string(".winDeltaT");
    COM_new_dataitem( dataName.c_str(), 'w', COM_DOUBLE, 1, "");
    COM_set_size(     dataName.c_str(), 0, 1);
    COM_set_array(    dataName.c_str(), 0, &(comFoamPtr->winDeltaT) );

    dataName = name+string(".winRun");
    COM_new_dataitem( dataName.c_str(), 'w', COM_INT, 1, "");
    COM_set_size(     dataName.c_str(), 0, 1);
    COM_set_array(    dataName.c_str(), 0, &(comFoamPtr->winRun));


 
    
    COM_window_init_done(name); 

    return 0;
}
//---------------------------------------------------------

//===================================================================


