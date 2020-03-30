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
        deleteFaceData();
        deleteSurfaceData();
    }

private:
    const int nComponents = 3;
    const int cellToPointTypeSize=8;
    const int faceToPointTypeSize=4;

    // Methods ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    int createVolumeConnectivities();
    int createVolumeData();
    int deleteVolumeData();
    int updateVolumeData();
    int registerVolumeData(const char *name);

    int createFaceConnectivities();
    int createFaceData();
    int deleteFaceData();
    int updateFaceData();
    int registerFaceData(const char *name);

    int createSurfaceConnectivities();
    int createSurfaceData();
    int deleteSurfaceData();
    int updateSurfaceData();
    int registerSurfaceData(const char *name);    

protected:

    // Variables ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    int ca_nPoints;
    int ca_nFaces;
    int ca_nCells;
    int ca_nPatches;
    //-------------------------------------------

    // COM Volume Arrays^^^^^^^^^^^^^^^^^^^^^^^^^^
    // Connectivity
    int   ca_cellToPointConn_types = 0;
    int*  ca_cellToPointConn_map = NULL;
    int*  ca_cellToPointConn_size = NULL;

    // Mapping
    int*  ca_cellToCellMap = NULL;
    int**  ca_celltoPointConn = NULL;

    // Field Data
    double* ca_Points = NULL;
    double* ca_Vel = NULL;
    double* ca_P = NULL;
    double* ca_T = NULL;
    double* ca_Rho = NULL;
    //-------------------------------------------
    
    // COM Face Arrays^^^^^^^^^^^^^^^^^^^^^^^^^^^
    // Mapping
    int*  ca_faceToFaceMap = NULL;

    // Connectivity
    int   ca_faceToPointConn_types = 0;
    int*  ca_faceToPointConn_map = NULL;
    int*  ca_faceToPointConn_size = NULL;
    int** ca_faceToPointConn = NULL;

    // Field data
    int* ca_faceOwner = NULL;
    int* ca_faceNeighb = NULL;
    //-------------------------------------------

    // COM Patch Arrays^^^^^^^^^^^^^^^^^^^^^^^^^^
    // General data
    std::string* ca_patchName = NULL;
    std::string* ca_patchType = NULL;
    wordList*    ca_patchInGroup = NULL;
    int*         ca_patchStart = NULL;
    int*         ca_patchSize = NULL;

    // PointToPoint Mapping
    int*  ca_patchPointToPointMap_size = NULL;
    int** ca_patchPointToPointMap = NULL;

    // FaceToFace Mapping
    int** ca_patchFaceToFaceMap = NULL;

    // FaceToPoint Mapping
    int*   ca_patchFaceToPointConn_types = NULL;
    int**  ca_patchFaceToPointConn_map = NULL;
    int**  ca_patchFaceToPointConn_size = NULL;
    int*** ca_patchFaceToPointConn = NULL;

    // Field data
    double** ca_patchPoints = NULL;
    double** ca_patchVel = NULL;
    double** ca_patchRho = NULL;
    double** ca_patchP = NULL;
    double** ca_patchT = NULL;
    //-------------------------------------------

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
