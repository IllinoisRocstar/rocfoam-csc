#include "comRhoPimple.H"
#include "comRhoCentral.H"

#ifndef COMLOADUNLOAD_
#define COMLOADUNLOAD_

using namespace COM;

COM_EXTERN_MODULE(comfoam);

//^^^^^ (UN)LOAD METHOD ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
// C/C++ bindings to load rocFoam
extern "C" void comfoam_load_module(const char *name, const char *solverType)
{
    if (solverType == "rhoCentral")
    {
        comRhoCentralModule::Load(name);
    }
    else if (solverType == "rhoPimple")
    {
        comRhoPimpleModule::Load(name);
    }
}

// C/C++ bindings to unload rocFoam
extern "C" void comfoam_unload_module(const char *name, const char *solverType)
{
    if (solverType == "rhoCentral")
    {
        comRhoCentralModule::Unload(name);
    }
    else if (solverType == "rhoPimple")
    {
        comRhoPimpleModule::Unload(name);
    }
}//===================================================================


#endif
