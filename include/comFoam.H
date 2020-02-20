#ifndef COMFOAM_
#define COMFOAM_

#include "rocFoam.H"
#include "com.h"

class comFoam : public COM_Object, public rocFoam
{

public:

    //  Object Functions ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    //static void Load(const char *name);
    //static void Unload(const std::string &name);

    //void Load(const char *name);
    //void Unload(const char *name);

    //  Constructor/Deconstructor ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    comFoam();
    comFoam(int *pargc, void **pargv, int *verbIn);

    ~comFoam() {}; 

protected:

    //  Window data ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    std::string winNameVol; /// Tracks *this* volume window name.
    std::string winNameSurf; /// Tracks *this* volume window name.

    std::string solverType;
    MPI_Comm winComm; /// Tracks *this* window communicator.
    
    int winNProc;
    int winRank;

    scalar winTime;
    scalar winDeltaT;
    int winRun;
    
    //  Nodal Coordinates
    
    double *comXYZ = NULL; // [x1,y1,z1, x2,y2,z2, ... ,xN,yN,zN]
    int nComponents = 3;
    int nNodes = 0;

    int *comFaces = NULL;
    int *comOwner = NULL;
    int *comNeighb = NULL;
    int nFaces = 0;


protected:
    int flowInit(int *pargc, void **pargv, int *verbIn);
    int flowLoop();
    int flowStep();
    
    int flowRegister();

    //virtual int flowFin() { return 0; }
    //virtual int flowStat() { return 0; }
};

#endif
