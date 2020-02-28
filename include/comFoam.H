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

    ~comFoam()
    {   
        deleteVolumeData();
    }

protected:

    int nPoints;
    int nFaces;
    int nCells;
    int nPatches;

    // Volume Connectivity Vectors ^^^^^^^^^^^^^^
    std::vector< std::vector< std::vector<int> >> vecCellToFaceConn;
    std::vector< std::vector< std::vector<int> >> vecCellToPointConn;
    //std::vector<std::vector< std::vector<int> >> vecCellToPointConnVTK;
    std::vector< std::vector<int> > vecCellTypeLocalToGlobalMap;
    std::vector< std::vector<int> > vecFaceToPointConn;
    std::vector<int> vecOwners;
    std::vector<int> vecNeighbs;
    
    //-------------------------------------------

    // COM VOLUME ARRAYS^^^^^^^^^^^^^^^^^^^^^^^^^
    // Connectivity ^^^^^^^^^^^^^^^^^^^
    int   nCellTypes = 0;
    int   nComponents = 3;
    int*  nCellTypeMap = NULL;
    int*  nCellTypeSize = NULL;
    
    int** cellTypeLocalToGlobalMap = NULL;
    int** comArrCelltoPointConn = NULL;

    // Field Data ^^^^^^^^^^^^^^^^^^^^^
    double* comArrPoints = NULL;
    double* comArrVel = NULL;
    double* comArrP = NULL;
    double* comArrT = NULL;
    double* comArrRho = NULL;
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
    

    

    int *comFaces = NULL;
    int *comOwner = NULL;
    int *comNeighb = NULL;


protected:
    int flowInit(int *pargc, void **pargv, const char *name);
    int flowLoop(const char *name);
    int flowStep(const char *name);

    int flowExtractData(const char *name);
    int flowRegisterData(const char *name);
    
    //virtual int flowFin() { return 0; }
    //virtual int flowStat() { return 0; }

    int registerFunctions(const char *name);

    int createVolumeConnectivities();
    int createVolumeData();
    int deleteVolumeData();
    int updateVolumeData();
    int registerData(const char *name);


    int createSurfaceConnectivities();


};

#endif
