#include "comFoam.H"

using namespace COM;

int rocFoamModule::initializeAll(int *pargc, void **pargv, int *verbIn)
{
    int argc = *pargc;
    char** argv = (char**)(pargv);
    verbosity = *verbIn;

    std::cout << "OF Module verbosity = " << verbosity << std::endl;

    //  OpenFOAM initializer ^^^^^^^^^^^^^^^
    initialize(argc, argv);
    
    //  Other initializations ^^^^^^^^^^^^^^

    return 0;
}


void rocFoamModule::Load(const char *name)
{
    std::cout << "RFModule.Load: Loading rocFoamModule with name "
              << name << "." << std::endl;

    //  Register module with COM ^^^^^^^^^^^
    rocFoamModule *rocFoamPtr = new rocFoamModule();

    COM_new_window(name, MPI_COMM_NULL);
    rocFoamPtr->my_window_name = name;

    std::string global_name(name+".global");

    COM_new_dataitem(global_name.c_str(), 'w', COM_VOID, 1, "");

    COM_set_object(global_name.c_str(), 0, rocFoamPtr);


    /// Register functions
    std::vector<COM_Type> types(13,COM_INT);

    types[0] = COM_RAWDATA;
    types[2] = COM_VOID;
    
    COM_set_member_function
    (
        (name+".initializeAll").c_str(),
        (Member_func_ptr)(&rocFoamModule::initializeAll),
        global_name.c_str(), "biii", &types[0]
    );

    COM_set_member_function
    (
        (name+".loop").c_str(),
        (Member_func_ptr)(&rocFoamModule::loop),
        global_name.c_str(), "b", &types[0]
    );

    COM_set_member_function
    (
        (name+".finalize").c_str(),
        (Member_func_ptr)(&rocFoamModule::finalize),
        global_name.c_str(), "b", &types[0]
    );
    

    COM_window_init_done(name); 

    return;
}


void rocFoamModule::Unload(const std::string &name)
{
    std::cout << "RFModule.Unload: Unloading rocFoamModule with name "
              << name << "." << std::endl;
              
    rocFoamModule *rocFoamPtr = NULL;
    std::string global_name(name+".global");

    COM_get_object(global_name.c_str(), 0, &rocFoamPtr);

    COM_assertion_msg(rocFoamPtr->validate_object()==0, "Invalid object");

    delete rocFoamPtr;

    COM_delete_window(std::string(name));
}



/// @brief C/C++ bindings to load IcoFoamModule
extern "C" void rocFoamModule_load_module(const char *name)
{
  rocFoamModule::Load(name);
}

/// @brief C/C++ bindings to unload IcoFoamModule
extern "C" void rocFoamModule_unload_module(const char *name)
{
  rocFoamModule::Unload(name);
}
