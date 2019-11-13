#include "rocRhoPimple.H"
#include "rocRhoCentral.H"

#ifndef COMLOADUNLOAD_
#define COMLOADUNLOAD_

using namespace COM;

COM_EXTERN_MODULE(comfoam);

//^^^^^ (UN)LOAD METHOD ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
// C/C++ bindings to load rocFoam
extern "C" void comfoam_load_module(const char *name, char *solverType)
{
    if (solverType == "rocRhoCentral")
    {
        rhoCentral::load(name);
    }
    else if (solverType == "rocRhoPimple")
    {
        rhoPimple::load(name);
    }
    else
    {
        std::cout << "comFoam.load: Undefined solverType " << solverType << std::endl;
        throw -4;
    }
}

// C/C++ bindings to unload rocFoam
extern "C" void comfoam_unload_module(const char *name, const char *solverType)
{
    if (solverType == "rocRhoCentral")
    {
        rhoCentral::unload(name);
    }
    else if (solverType == "rocRhoPimple")
    {
        rhoPimple::unload(name);
    }
    else
    {
        std::cout << "comFoam.load: Undefined solverType " << solverType << std::endl;
        throw -4;
    }
}//===================================================================


#endif
