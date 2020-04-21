#ifndef COMFOAM_
#define COMFOAM_

#include "rocFoam.H"
#include "com.h"
#include "fileContainer.H"
#include <filesystem.hpp>

class comFoam : public COM_Object, public rocFoam
{

public:
    //  Constructor/Deconstructor ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    comFoam();
    comFoam(int *pargc, void **pargv, const char *name);
    ~comFoam();
private:
    int initSet();

protected:
    const int nComponents = 3;
private:
    const int cellToPointTypeSize=8;
    const int faceToPointTypeSize=4;

    // Methods ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    int createVolumeConnectivities();
    int createVolumeData();
    int deleteVolumeData();
    int updateVolumeData();
    int registerVolumeData(const char *name);
    int reconstCaVolumeData(const char *name);

    int createFaceConnectivities();
    int createFaceData();
    int deleteFaceData();
    int updateFaceData();
    int registerFaceData(const char *name);
    int reconstCaFaceData(const char *name);

    int createSurfaceConnectivities();
    int createSurfaceData();
    int deleteSurfaceData();
    int updateSurfaceData();
    int registerSurfaceData(const char *name);
    int reconstCaSurfaceData(const char *name);

    // tempFiles-specific
    int readInitFiles(const std::string& rootAddr);
    int createInitFiles(const std::string& rootAddr);
    int deleteInitFiles(const std::string& addr);
    int deleteFilesData();
    int registerInitFiles(const char *name);
    int reconstCaInitFiles(const char *name, const std::string& rootAddr);
    
    
    int readRecursive(const std::string& rootAddr,
                      std::string fullAddr,
                      std::vector<fileContainer>& vecFile,
                      int& fileCount);


    bool nameExists(const std::vector<std::string>& dataItemNames,
                          const std::string& dataName);
    int reconstMesh();

protected:
    // Variables ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    int* ca_nPoints;  //single value
    int* ca_nCells;   //single value
    int* ca_nFaces;   //single value
    int* ca_nPatches; //single value
    //-------------------------------------------

    // COM Volume Arrays^^^^^^^^^^^^^^^^^^^^^^^^^^
    // Mapping
    int*  ca_cellToCellMap;

    // Connectivity
    int*  ca_cellToPointConn_types; //single value
    int*  ca_cellToPointConn_map;
    int*  ca_cellToPointConn_size;
    int** ca_cellToPointConn;

    // Field Data
    double* ca_Points;
    double* ca_Vel;
    double* ca_P;
    double* ca_T;
    double* ca_Rho;
    //-------------------------------------------
    
    // COM Face Arrays^^^^^^^^^^^^^^^^^^^^^^^^^^^
    // Mapping
    int*  ca_faceToFaceMap;

    // Connectivity
    int*  ca_faceToPointConn_types; //single value
    int*  ca_faceToPointConn_map;
    int*  ca_faceToPointConn_size;
    int** ca_faceToPointConn;

    // Field data
    int* ca_faceOwner;
    int* ca_faceNeighb;
    //-------------------------------------------

    // COM Patch Arrays^^^^^^^^^^^^^^^^^^^^^^^^^^
    // General data

    std::string* patchNameStr; //single value for the last
    std::string* patchTypeStr; //single value for the last

    char** ca_patchName; //single value for the last
    char** ca_patchType; //single value for the last
    wordList**    ca_patchInGroup; //single value for the last
    int**         ca_patchStart; //single value for the last
    int**         ca_patchSize; //single value for the last

    // PointToPoint Mapping
    int** ca_patchPointToPointMap_size; //single value for the last
    int** ca_patchPointToPointMap;

    // FaceToFace Mapping
    int** ca_patchFaceToFaceMap;

    // FaceToPoint Mapping
    int**  ca_patchFaceToPointConn_types; //single value for the last
    int**  ca_patchFaceToPointConn_map;
    int**  ca_patchFaceToPointConn_size;
    int*** ca_patchFaceToPointConn;

    // Field data
    double** ca_patchPoints;
    double** ca_patchVel;
    double** ca_patchP;
    double** ca_patchT;
    double** ca_patchRho;
    //-------------------------------------------
    
    // Files
    int*   ca_nFiles;
    int*   ca_fileSize;
    char** ca_fileName;
    char** ca_filePath;
    char** ca_fileContent;

    //-------------------------------------------

protected:
    //  Window data ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    //std::string winVolName; /// Tracks *this* volume window name.
    //std::string winSurfName; /// Tracks *this* volume window name.

    std::string solverType;
    MPI_Comm winComm; /// Tracks *this* window communicator.
    
    int ca_nProc;
    int ca_myRank;

    // registered data set during the simulation
    int* ca_runStat;
    double* ca_time;
    double* ca_deltaT;
    
    
    // Methods ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    
    // Solver
    int flowInit(int *pargc, void **pargv, const char *name);
    int flowLoop(const char *name);
    int flowStep(const char *name);

    // COM-specific
    int registerFunctions(const char *name);
    int reconstCaData(int* pargc, void** pargv, const char* name);
    int reconstDynamicFvMesh();

public: 
    static void copyWindow(const char *name1, const char *name2);
};

#endif

