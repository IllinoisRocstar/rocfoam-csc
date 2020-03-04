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
        deleteSurfaceData();
    }

private:

    // Variables ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    int ca_nPoints;
    int ca_nFaces;
    int ca_nCells;
    int ca_nPatches;

    // Volume Connectivity Vectors ^^^^^^^^^^^^^^
    std::vector< std::vector<int> > vecCellToCellMap_localToGlobal;
    std::vector<int> vecOwners;
    std::vector<int> vecNeighbs;
    //-------------------------------------------

    // COM Volume Arrays^^^^^^^^^^^^^^^^^^^^^^^^^^
    // Connectivity ^^^^^^^^^^^^^^^^^^^
    const int nComponents = 3;
    int   ca_cellToPointConn_types = 0;
    int*  ca_cellToPointConn_map = NULL;
    int*  ca_cellToPointConn_size = NULL;

    int** ca_cellToCellMap_localToGlobal = NULL;
    int** ca_celltoPointConn = NULL;

    int* ca_faces = NULL;
    int* ca_owner = NULL;
    int* ca_neighb = NULL;

    // Field Data ^^^^^^^^^^^^^^^^^^^^^
    double* ca_Points = NULL;
    double* ca_Vel = NULL;
    double* ca_P = NULL;
    double* ca_T = NULL;
    double* ca_Rho = NULL;
    //-------------------------------------------
    
    // Patch Connectivity Vectors ^^^^^^^^^^^^^^
    std::vector< std::vector< std::vector<int> >> vecPatchFaceToFaceMap_localToLocal;
    std::vector< std::vector<int> >  vecPatchPointToPointMap_localToGlobal;
    std::vector< std::vector< std::vector< std::vector<int> >>>
        vecPatchFaceToPointConn_localToLocal;
    //-------------------------------------------

    // COM Patch Arrays^^^^^^^^^^^^^^^^^^^^^^^^^^
    // Connectivity ^^^^^^^^^^^^^^^^^^^
    int* patchFaceToFaceMap_types = NULL;
    int** patchFaceToFaceMap_map = NULL;
    int** patchFaceToFaceMap_size = NULL;
    int*** patchFaceToFaceMap_localToLocal = NULL;

    int* patchPointToPointMap_size = NULL;
    int** patchPointToPointMap_localToGlobal = NULL;

    int* ca_patchFaceToPointConn_types = NULL;
    int** ca_patchFaceToPointConn_map = NULL;
    int** ca_patchFaceToPointConn_size = NULL;
    int*** ca_patchFaceToPointConn_localToLocal = NULL;

    // General data ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    std::string* ca_patchName = NULL;
    std::string* ca_patchType = NULL;
    wordList* ca_patchInGroup = NULL;
    int* ca_patchStart = NULL;
    int* ca_patchSize = NULL;

    // Field data ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    double** ca_patchPoints = NULL;
    double** ca_patchVel = NULL;
    double** ca_patchRho = NULL;
    double** ca_patchP = NULL;
    double** ca_patchT = NULL;
    //-------------------------------------------

    // Methods ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    int createVolumeConnectivities();
    int createVolumeData();
    int deleteVolumeData();
    int updateVolumeData();
    int registerVolumeData(const char *name);

    int createSurfaceConnectivities();
    int createSurfaceData();
    int deleteSurfaceData();
    int updateSurfaceData();
    int registerSurfaceData(const char *name);    
    
    
protected:
    //  Window data ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    std::string winVolName; /// Tracks *this* volume window name.
    std::string winSurfName; /// Tracks *this* volume window name.

    std::string solverType;
    MPI_Comm winComm; /// Tracks *this* window communicator.
    
    int winNProc;
    int winRank;

    scalar winTime;
    scalar winDeltaT;
    int winRun;
    
    // Methods ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    int registerFunctions(const char *name);

    int flowInit(int *pargc, void **pargv, const char *name);
    int flowLoop(const char *name);
    int flowStep(const char *name);

};

#endif
