#include "comFoam.H"

comFoam::comFoam()
    : solverType(""),
      winName(""),
      winComm(NULL),
      winNProc(0),
      winRank(0)
{};

comFoam::comFoam(int *pargc, void **pargv, int *verbIn)
    : solverType(""),
      winName(""),
      winComm(NULL),
      winNProc(0),
      winRank(0)
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

    //COM_set_member_function
    //(
    //    (name + string(".flowFin")).c_str(),
    //    reinterpret_cast<Member_func_ptr>(&rhoCentral::flowFin),
    //    globalName.c_str(), "b", &types[0]
    //);

    //  Registering nproc for this module to COM ^^^^^^^^^^
    COM_new_dataitem( (name+string(".winNProc")).c_str(), 'w', COM_INT, 1, "");
    COM_set_size(     (name+string(".winNProc")).c_str(), 0, 1);
    COM_set_array(    (name+string(".winNProc")).c_str(), 0, &(comFoamPtr->winNProc));

    COM_window_init_done(name); 

    return 0;
}
//---------------------------------------------------------

//===================================================================


