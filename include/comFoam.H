#ifndef COMFOAM_
#define COMFOAM_

#include "com.h"
#include "rocRhoCentral.H"


using namespace COM;

class comFoamModule : public rocRhoCentral, public COM_Object
{
public:
    std::string winName; /// Tracks *this* window name.

    MPI_Comm winComm; /// Tracks *this* window communicator.
    
    int winNProc;
    int winRank;

    //  Object Functions ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    static void Load(const char *name);
    static void Unload(const std::string &name);

    //  Constructor/Deconstructor ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    comFoamModule();
    comFoamModule(int *pargc, void **pargv, int *verbIn);

    ~comFoamModule(){};
private:

protected:
    int flowInit(int *pargc, void **pargv, int *verbIn);
    int flowStat();
    int flowLoop();
    int flowFin();
};

#endif
