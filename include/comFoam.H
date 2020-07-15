#ifndef COMFOAM_
#define COMFOAM_

#include "rocFoam.H"

#include "displacementMotionSolver.H"
#include "dynamicMotionSolverFvMesh.H"
#include "surfaceDisplacementPointPatchVectorFieldFSI.H"

#include "com.h"
#include "fileContainer.H"
#include <boost/filesystem.hpp>
#include <iomanip>

class comFoam : public COM_Object, public rocFoam
{
private:
    const int faceToPointTypeSize{4};
    const int panePadding{100};
    const int IODigits{10};
    const std::string tmpFluidDir{"./fluidTmp"};
    const std::string movingWallName = "flap";
    
protected:
    const int nComponents{3};
    const int genCharSize{20};

    pointVectorField* pointDisplacementNewPtr{nullptr};

public:
    //  Constructor/Deconstructor ^^^^^^^^^^^^^^^^^^^^^^^^^
    comFoam();
    //comFoam(int *pargc, void **pargv, const char *name);
    ~comFoam();
    
    virtual int finalizeFoam();

private:
    // CSC collective methods
    int createCSCdata();
    int updateCSCdata();
    int registerCSCdata(const char *name);
    int deleteCSCdata();
    int reconstCSCdata(const char *name);
    

    // CSC status methods
    int createStatusData();
    int updateStatusData();
    int registerStatusData(const char *name);
    int reconstStatusData(const char *name);
    int deleteStatusData();

    // CSC volume methods
    int createVolumeConnectivities();
    int createVolumeData();
    int deleteVolumeData();
    int updateVolumeData_outgoing();
    int registerVolumeData(const char *name);
    int reconstVolumeData(const char *name);

    // CSC face methods
    int createFaceConnectivities();
    int createFaceData();
    int deleteFaceData();
    int updateFaceData_outgoing();
    int registerFaceData(const char *name);
    int reconstFaceData(const char *name);

    // CSC surface methods
    int createSurfaceConnectivities();
    int createSurfaceData();
    int deleteSurfaceData();
    int updateSurfaceData_incoming();
    int updateSurfaceData_outgoing();
    int registerSurfaceData(const char *name);
    int reconstSurfaceData(const char *name);

    // CSC files methods
    int readFilesData(const std::string& rootAddr);
    int createFilesData();
    int deleteFilesData();
    int registerFilesData(const char *name);
    int reconstFilesData(const char *name);

    int deleteTempFiles(const std::string& addr);

    int createSysConstFiles
    (
        const std::string& rootAddr,
        const std::vector<fileContainer>& vecFile
    );
    int createFieldFiles
    (
        const std::string& rootAddr,
        const std::vector<fileContainer>& vecFile
    );
    int createBoundaryFile
    (
        const std::string& rootAddr,
        const std::vector<fileContainer>& vecFile
    );
    int createConnectivityFiles
    (
        const std::string& rootAddr,
        const std::vector<fileContainer>& vecFile
    );
    int createUniformTimeFile(const std::string& rootAddr);
    int createPointsFile(const std::string& rootAddr);
    int createOwnerFile(const std::string& rootAddr);
    int createNeighborFile(const std::string& rootAddr);
    int createFacesFile(const std::string& rootAddr);
    int deleteInitFiles(const std::string& addr);

    
    std::string createBaseFile
    (
        std::string name,
        std::string loc,
        std::string type,
        std::string dim
    );

    bool fileShouldBeRead
    (
        const std::string& locaParAddr,
        const std::string& localAddr,
        const std::string& fileName
    );
    
    int readRecursive(const std::string& rootAddr,
                      std::string fullAddr,
                      std::vector<fileContainer>& vecFile,
                      int& fileCount);


    bool nameExists(const std::vector<std::string>& dataItemNames,
                          const std::string& dataName);
    int reconstMesh();


    // Variables ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    int* ca_nPoints{nullptr};  //single value
    int* ca_nCells{nullptr};   //single value
    int* ca_nFaces{nullptr};   //single value
    int* ca_nPatches{nullptr}; //single value
    //-------------------------------------------

    // COM Volume Arrays^^^^^^^^^^^^^^^^^^^^^^^^^^
    // Mapping
    int*  ca_cellToCellMap{nullptr};
    int*  ca_cellToCellMap_inverse{nullptr};

    // Connectivity
    int*  ca_cellToPointConn_types{nullptr}; //single value
    int*  ca_cellToPointConn_map{nullptr};
    int*  ca_cellToPointConn_size{nullptr};
    int** ca_cellToPointConn{nullptr};

    // Field Data
    double* ca_Points{nullptr};
    double* ca_Disp{nullptr};

    double* ca_Vel{nullptr};
    double* ca_P{nullptr};
    double* ca_T{nullptr};
    double* ca_Rho{nullptr};
    double* ca_Phi{nullptr};
    double* ca_RhoUf{nullptr};

    // Turbulence data
    double* ca_AlphaT{nullptr};
    double* ca_Epsilon{nullptr};
    double* ca_K{nullptr};
    double* ca_NuT{nullptr};
    //-------------------------------------------
    
    // COM Face Arrays^^^^^^^^^^^^^^^^^^^^^^^^^^^
    // Mapping
    int*  ca_faceToFaceMap{nullptr};
    int*  ca_faceToFaceMap_inverse{nullptr};

    // Connectivity
    int*  ca_faceToPointConn_types{nullptr}; //single value
    int*  ca_faceToPointConn_map{nullptr};
    int*  ca_faceToPointConn_size{nullptr};
    int** ca_faceToPointConn{nullptr};

    // Field data
    int* ca_faceOwner{nullptr};
    int* ca_faceNeighb{nullptr};
    //-------------------------------------------

    // COM Patch Arrays^^^^^^^^^^^^^^^^^^^^^^^^^^
    // General data

    std::string* patchNameStr{nullptr}; //single value for the last
    std::string* patchTypeStr{nullptr}; //single value for the last

    char** ca_patchName{nullptr}; //single value for the last
    char** ca_patchType{nullptr}; //single value for the last
    wordList** ca_patchInGroup{nullptr}; //single value for the last
    int**      ca_patchStart{nullptr}; //single value for the last
    int**      ca_patchSize{nullptr}; //single value for the last

    // PointToPoint Mapping
    int** ca_patchPointToPointMap_size{nullptr}; //single value for the last
    int** ca_patchPointToPointMap{nullptr};

    // FaceToFace Mapping
    int** ca_patchFaceToFaceMap{nullptr};
    int** ca_patchFaceToFaceMap_inverse{nullptr};

    // FaceToPoint Mapping
    int**  ca_patchFaceToPointConn_types{nullptr}; //single value for the last
    int**  ca_patchFaceToPointConn_map{nullptr};
    int**  ca_patchFaceToPointConn_size{nullptr};
    int*** ca_patchFaceToPointConn{nullptr};

    // Field data
    double** ca_patchPoints{nullptr};
    double** ca_patchVel{nullptr};
    double** ca_patchP{nullptr};
    double** ca_patchT{nullptr};
    double** ca_patchRho{nullptr};
    double** ca_patchPhi{nullptr};
    double** ca_patchRhoUf{nullptr};

    // Turbulence data
    double** ca_patchAlphaT{nullptr};
    double** ca_patchEpsilon{nullptr};
    double** ca_patchK{nullptr};
    double** ca_patchNuT{nullptr};
    
    // Agent data
    int**    ca_bcflag{nullptr}; //single value for the last
    double** ca_patchNf{nullptr};
    double** ca_patchSf{nullptr};
    double** ca_patchTrac{nullptr};
    double** ca_patchDisp{nullptr};     //Displacement
    double** ca_patchMassFlux{nullptr}; //Mass flux (scalar)
    double** ca_patchMomentum{nullptr}; //Momentum flux (vector)
    double** ca_patchFlameT{nullptr};   //Flame Temperature (scalar)
    //-------------------------------------------
    
    // Files
    int*   ca_nFiles{nullptr};
    int*   ca_fileSize{nullptr};
    char** ca_fileName{nullptr};
    char** ca_filePath{nullptr};
    char** ca_fileContent{nullptr};

    //-------------------------------------------

    // Status data
    int*    ca_runStat{nullptr};
    int*    ca_timeIndex{nullptr};
    double* ca_time{nullptr};
    char*   ca_timeName{nullptr};
    double* ca_deltaT{nullptr};
    double* ca_deltaT0{nullptr};
    
    int* ca_movingMesh{nullptr};
    char* ca_dynamicFvMesh{nullptr};
    

protected:

    //  Window data ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    //std::string winVolName; /// Tracks *this* volume window name.
    //std::string winSurfName; /// Tracks *this* volume window name.
    std::string solverType{""};
    MPI_Comm winComm{};
    
    int ca_nProc{0};
    int ca_myRank{0};
    std::string winName{""};
    
    // Methods ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    // rocFoam driver methods
    int flowInit(int *pargc, void **pargv, const char *name);
    int flowLoop();
    int flowStep();

    // COM-specific
    int registerFunctions(const char *name);
    int restartInit(int* pargc, void** pargv, const char* name);

    virtual int loadInternal(const char* name) {return 0;}
    //-------------------------------------------

public:

    // RocStar-specific methods ^^^^^^^^^^^^^^^^^
    void initialize
    (
        const double *initTime,
        const MPI_Comm* flowComm,
        const int* manInitHandle,
        const char* volName,
        const char* surfName,
        const int* obtainHandle
    );

    void update_solution
    (
        double* currentTime,
        double* timeStep,
        int* bcHandle,
        int* gmHandle
    );

    void finalize();
    //-------------------------------------------

    static void copyWindow(const string& name1, const string& name2);
    int findGlobalIndex(int* arr, const int& size,  const int& elem);
    static std::string removeTrailZero(std::string in);

    size_t findLoc
    (
        const std::string& fullAddr,
        const std::string& content,
        const std::string& exp,
        size_t start = 0
    );

    size_t findLoc
    (
        const std::string& fullAddr,
        const std::string& content,
        const std::string& exp,
        size_t start,
        size_t end
    );

    size_t findOnlyLoc
    (
        const std::string& fullAddr,
        const std::string& content,
        const std::string& exp,
        size_t start=0,
        size_t end=std::string::npos
    );

};

#endif

