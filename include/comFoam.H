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
    comFoam(int *pargc, void **pargv, const char *name);

    ~comFoam() {}; 

protected:

    int nPoints;
    int nFaces;
    int nCells;
    int nPatches;

    // Field Vectors ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    std::vector< std::vector<double> > vecPoints;
    std::vector< std::vector<int> > vecFaceToPointConn;
    std::vector<int> vecOwners;
    std::vector<int> vecNeighbs;
    
    std::vector< std::vector<double> > vecFieldVel;
    std::vector<double> vecFieldRho;
    std::vector<double> vecFieldP;
    std::vector<double> vecFieldT;
    //-------------------------------------------

    // Patch Vectors ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    std::vector< std::vector<int> > patchGlobalPointIndex;
    std::vector< std::vector< std::vector<int> >> patchLocalConnectivity;

    std::vector< std::vector< std::vector<double> >> vecPatchVel;
    std::vector< std::vector<double> > vecPatchRho;
    std::vector< std::vector<double> > vecPatchP;
    std::vector< std::vector<double> > vecPatchT;

    std::vector<std::string> vecPatchName;
    std::vector<std::string> vecPatchType;
    // Note: I need to change this
    std::vector< wordList > vecPatchInGroup;

    std::vector<int> vecPatchStart;
    std::vector<int> vecPatchNFaces;
    std::vector<int> vecPatchNPoints;
    //-------------------------------------------

    //  Window data ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    std::string winVolName; /// Tracks *this* volume window name.
    std::string winSrfName; /// Tracks *this* volume window name.

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


protected:
    int flowInit(int *pargc, void **pargv, const char *name);
    int flowLoop(const char *name);
    int flowStep(const char *name);

    int flowExtractData(const char *name);
    
    
    //virtual int flowFin() { return 0; }
    //virtual int flowStat() { return 0; }

    int registerFunctions(const char *name);
    int extractData();
};

#endif
