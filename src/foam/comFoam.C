#include "comFoam.H"

comFoamModule::comFoamModule()
    : rocFoam(),
      solverType(""),
      winName(""),
      winComm(NULL),
      winNProc(0),
      winRank(0)
{};

comFoamModule::comFoamModule(int *pargc, void **pargv, int *verbIn)
    : rocFoam(),
      solverType(""),
      winName(""),
      winComm(NULL),
      winNProc(0),
      winRank(0)
{
    flowInit(pargc, pargv, verbIn);
}


//^^^ DEFINITION OF COM-RELATED MTHODS ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
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
/*
extern "C" void comfoam_load_module(const char *name)
{
  comFoamModule::Load(name);
}


void comFoamModule::Load(const char *name)
{
    return;
}
*/

//-------------------------------------------------------------------


//^^^^^ UNLOAD MODULES ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

// C/C++ bindings to unload rocFoam
/*
extern "C" void comfoam_unload_module(const char *name)
{
  comFoamModule::Unload(name);
}

void comFoamModule::Unload(const std::string &name)
{
    return;
}
*/
//-------------------------------------------------------------------





