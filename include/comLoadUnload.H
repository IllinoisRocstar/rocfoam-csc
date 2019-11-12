#include "rocRhoPimple.H"
#include "rocRhoCentral.H"

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
        rhoCentral::Load(name);
    }
    else if (solverType == "rhoPimple")
    {
        rhoPimple::Load(name);
    }
}

// C/C++ bindings to unload rocFoam
extern "C" void comfoam_unload_module(const char *name, const char *solverType)
{
    if (solverType == "rhoCentral")
    {
        rhoCentral::Unload(name);
    }
    else if (solverType == "rhoPimple")
    {
        rhoPimple::Unload(name);
    }
}//===================================================================


#endif
