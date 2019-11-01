#include "comFoam.H"

comFoamModule::comFoamModule(){};

comFoamModule::comFoamModule(int *pargc, void **pargv, int *verbIn)
{
    flowInit(pargc, pargv, verbIn);
}


int comFoamModule::flowInit(int *pargc, void **pargv, int *verbIn)
{
    int argc = *pargc;
    char** argv = (char**)(pargv);
    verbosity = *verbIn;

    std::cout << "RFModule.flowInit: Initializing flow solver." << std::endl;
    

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

    std::cout << "RFModule.flowLoop: flow interation." << std::endl;

    //  Call the flow iterator ^^^^^^^^^^^^^^^^^^
    comFoamModule *comFoamPtr = NULL;

    std::string name="CFModule";
    std::string globalName(name+".global");
    COM_get_object(globalName.c_str(), 0, &comFoamPtr);

    comFoamPtr->loop();
    
    return 0;
}

//^^^^^ FINALIZE MODULES ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
void comFoamModule::Unload(const std::string &name)
{
    std::cout << "RFModule.Unload: Unloading comFoamModule with name "
              << name << "." << std::endl;
              
    comFoamModule *comFoamPtr = NULL;
    std::string globalName(name+".global");

    COM_get_object(globalName.c_str(), 0, &comFoamPtr);
    
    comFoamPtr->finalize();

    COM_assertion_msg(comFoamPtr->validate_object()==0, "Invalid object");
    delete comFoamPtr;

    COM_delete_window(std::string(name));
}


void comFoamModule::Load(const char *name)
{
    std::cout << "RFModule.Load: Loading comFoamModule with name "
              << name << "." << std::endl;

    //  Register module with COM ^^^^^^^^^^^
    comFoamModule *comFoamPtr = new comFoamModule();

    COM_new_window(name, MPI_COMM_NULL);

    comFoamPtr->windowName = name;

    std::string globalName = name + string(".global");

    COM_new_dataitem(globalName.c_str(), 'w', COM_VOID, 1, "");

    COM_set_object(globalName.c_str(), 0, comFoamPtr);


    /// Register functions
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

    COM_window_init_done(name); 

    return;
}



/// @brief C/C++ bindings to load IcoFoamModule
extern "C" void comfoam_load_module(const char *name)
{

std::cout << "HERE 1001" << std::endl;

  comFoamModule::Load(name);
  
std::cout << "HERE 1002" << std::endl;
  
}

/// @brief C/C++ bindings to unload IcoFoamModule
extern "C" void comfoam_unload_module(const char *name)
{
  comFoamModule::Unload(name);
}




