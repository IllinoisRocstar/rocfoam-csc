#include "comFoam.H"
#include "PstreamGlobals.H"

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





