#include "comFoam.H"
#include "PstreamGlobals.H"

comFoamModule::comFoamModule(){};

comFoamModule::comFoamModule(int *pargc, void **pargv, int *verbIn)
{
    flowInit(pargc, pargv, verbIn);
}


int comFoamModule::flowInit(int *pargc, void **pargv, int *verbIn)
{
    int argc = *pargc;
    char** argv = (char**)(pargv);

    Foam::Info << "RFModule.flowInit: Initializing flow solver." << Foam::endl;

    //  OpenFOAM initializer ^^^^^^^^^^^^^^^^^^^^
    comFoamModule *comFoamPtr = NULL;

    std::string name="CFModule";
    std::string globalName(name+".global");
    COM_get_object(globalName.c_str(), 0, &comFoamPtr);

    comFoamPtr->initialize(argc, argv);
    
    //  Other initializations ^^^^^^^^^^^^^^^^^^^

    return 0;
}

int comFoamModule::flowLoop()
{

    Foam::Info << "RFModule.flowLoop: flow interation." << Foam::endl;

    //  Call the flow iterator ^^^^^^^^^^^^^^^^^^
    comFoamModule *comFoamPtr = NULL;

    std::string name="CFModule";
    std::string globalName(name+".global");
    COM_get_object(globalName.c_str(), 0, &comFoamPtr);

    comFoamPtr->loop();
    
    return 0;
}

//^^^^^ LOAD MODULES ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

// C/C++ bindings to load rocFoam
extern "C" void comfoam_load_module(const char *name)
{
  comFoamModule::Load(name);
}

void comFoamModule::Load(const char *name)
{
    Foam::Info << "RFModule.Load: Loading comFoamModule with name "
               << name << "." << Foam::endl;

    //  Anouncing default communicator  ^^^^^^^^^^^^^^^^^^^
    MPI_Comm tmpComm;
    tmpComm = COM_get_default_communicator();  

    int tmpRank, tmpNProc;
    MPI_Comm_rank(tmpComm, &tmpRank);
    MPI_Comm_size(tmpComm, &tmpNProc);
    
    Foam::Info << "RFModoule.Load: Rank #" << tmpRank
               << " on communicator " << tmpComm
               << " with " << tmpNProc << " processes."
               << Foam::endl;

    Foam::Info << "RFModule.Load: Rank #" << tmpRank
               << " Loading FsiFoamModule with name " 
               << name << Foam::endl;

    Foam::Info << Foam::endl;


    //  Register module with COM ^^^^^^^^^^^^^^^^^^^^^^^^^^
    comFoamModule *comFoamPtr = new comFoamModule();

    COM_new_window(name, MPI_COMM_NULL);
    //COM_new_window(name, tmpComm);

    comFoamPtr->winName = name;

    //MPI_Comm_dup(tmpComm, &(comFoamPtr->winComm));
    comFoamPtr->winComm = tmpComm;
    
    //Foam::PstreamGlobals::MPI_comFoam_to_openFoam = comFoamPtr->winComm;
    Foam::PstreamGlobals::MPI_COMM_FOAM = comFoamPtr->winComm;
    
    MPI_Comm_rank(comFoamPtr->winComm, &(comFoamPtr->winRank));
    MPI_Comm_size(comFoamPtr->winComm, &(comFoamPtr->winNProc));

    std::string globalName = name + string(".global");

    COM_new_dataitem(globalName.c_str(), 'w', COM_VOID, 1, "");

    COM_set_object(globalName.c_str(), 0, comFoamPtr);


    /// Register functions ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    std::vector<COM_Type> types(13,COM_INT);

    types[0] = COM_RAWDATA;
    types[2] = COM_VOID;

    COM_set_member_function
    (
        (name + string(".flowInit")).c_str(),
        (Member_func_ptr)(&comFoamModule::flowInit),
        globalName.c_str(), "biii", &types[0]
    );


    COM_set_member_function
    (
        (name + string(".flowLoop")).c_str(),
        (Member_func_ptr)(&comFoamModule::flowLoop),
        globalName.c_str(), "b", &types[0]
    );

    //COM_set_member_function
    //(
    //    (name + string(".flowFin")).c_str(),
    //    (Member_func_ptr)(&comFoamModule::flowFin),
    //    globalName.c_str(), "b", &types[0]
    //);

    //  Registering nproc for this module to COM ^^^^^^^^^^
    COM_new_dataitem( (name+string(".winNProc")).c_str(), 'w', COM_INT, 1, "");
    COM_set_size(     (name+string(".winNProc")).c_str(), 0, 1);
    COM_set_array(    (name+string(".winNProc")).c_str(), 0, &(comFoamPtr->winNProc));

    COM_window_init_done(name); 

    return;
}
//-------------------------------------------------------------------


//^^^^^ UNLOAD MODULES ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

// C/C++ bindings to unload rocFoam
extern "C" void comfoam_unload_module(const char *name)
{
  comFoamModule::Unload(name);
}

void comFoamModule::Unload(const std::string &name)
{
    std::cout << "RFModule.Unload: Unloading comFoamModule with name "
              << name << "." << std::endl;

    comFoamModule *comFoamPtr = NULL;
    std::string globalName(name+".global");

    COM_get_object(globalName.c_str(), 0, &comFoamPtr);

    delete comFoamPtr;

    COM_delete_window(std::string(name));
}
//-------------------------------------------------------------------





