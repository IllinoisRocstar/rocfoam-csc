#ifndef COMFOAM_
#define COMFOAM_

//#include "com.h"

#include "rocFoam.H"

class comFoamModule : public rocFoam
{

public:

    std::string solverType;
    std::string winName; /// Tracks *this* window name.
    MPI_Comm winComm; /// Tracks *this* window communicator.
    
    int winNProc;
    int winRank;

    //  Object Functions ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    //static void Load(const char *name);
    //static void Unload(const std::string &name);

    //void Load(const char *name);
    //void Unload(const char *name);

    //  Constructor/Deconstructor ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    comFoamModule();
    comFoamModule(int *pargc, void **pargv, int *verbIn);

    ~comFoamModule() {}; 
private:

protected:
    int flowInit(int *pargc, void **pargv, int *verbIn);
    int flowLoop();

    //virtual int flowFin() { return 0; }
    //virtual int flowStat() { return 0; }
};

#endif
