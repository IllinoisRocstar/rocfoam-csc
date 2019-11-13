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
    int argc = *pargc;
    char** argv = reinterpret_cast<char**>(pargv);

    Foam::Info << "RFModule.flowInit: Initializing flow solver." << Foam::endl;

    //  OpenFOAM initializer ^^^^^^^^^^^^^^^^^^^^
    comFoam *comFoamPtr = NULL;

    std::string name="ROCFOAM";
    std::string globalName(name+".global");
    COM_get_object(globalName.c_str(), 0, &comFoamPtr);

    comFoamPtr->initialize(argc, argv);
    
    //  Other initializations ^^^^^^^^^^^^^^^^^^^

    return 0;
}


int comFoam::flowLoop()
{

    Foam::Info << "RFModule.flowLoop: flow interation." << Foam::endl;

    //  Call the flow iterator ^^^^^^^^^^^^^^^^^^
    comFoam *comFoamPtr = NULL;

    std::string name="ROCFOAM";
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
  comFoam::Load(name);
}


void comFoam::Load(const char *name)
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
  comFoam::Unload(name);
}

void comFoam::Unload(const std::string &name)
{
    return;
}
*/
//-------------------------------------------------------------------





