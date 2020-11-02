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
    const std::string movingWallTypeName = "rocstarFSI";

protected:
    const int nComponents{3};
    const size_t genCharSize{50};

    pointVectorField* pointDisplacementNewPtr{nullptr};

    // control timeStepping when called by Rocstar
    bool modifiedDeltaT{false};
    double unmodifiedDeltaTvalue{0};
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

    //CSC zones methods
    int createZonesData();
    int deleteZonesData();
    int registerZonesData(const char *name);
    int reconstZonesData(const char *name);

protected:
    int updateSurfaceData_incoming(const int& count = 1);

private:
    int updateSurfaceData_outgoing();
    int registerSurfaceData(const char *name);
    int reconstSurfaceData(const char *name);

    void readMeshDict(const dictionary& currentDict, bool& foundSolver);

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

    int createCellZonesFile(const std::string& rootAddr);
    int createFaceZonesFile(const std::string& rootAddr);
    int createPointZonesFile(const std::string& rootAddr);
    
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
    bool*   pointUpdated{nullptr}; //Keep track of updated FSI points

    double* ca_Vel{nullptr};
    double* ca_P{nullptr};
    double* ca_T{nullptr};
    double* ca_Rho{nullptr};
    double* ca_Phi{nullptr};
    double* ca_RhoUf{nullptr};

    // Turbulence data
    double* ca_AlphaT{nullptr};
    double* ca_K{nullptr};
    double* ca_Epsilon{nullptr};
    double* ca_Omega{nullptr};
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
    int* ca_maxNameLength{nullptr};
    int* ca_maxTypeLength{nullptr};

    std::string* patchNameStr{nullptr};
    std::string* patchTypeStr{nullptr};

    char* ca_patchName{nullptr};
    char* ca_patchType{nullptr};
    wordList** ca_patchInGroup{nullptr}; //single value for the last


    int* ca_patchStart{nullptr};
    int* ca_patchSize{nullptr};

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
    double** ca_patchK{nullptr};
    double** ca_patchEpsilon{nullptr};
    double** ca_patchOmega{nullptr};
    double** ca_patchNuT{nullptr};
    
    // Agent data
    int**    ca_bcflag{nullptr}; //single value for the last
    double** ca_patchNf{nullptr};
    double** ca_patchSf{nullptr};
    double** ca_patchTrac{nullptr};

    double** ca_patchDisp{nullptr};      //Displacement
    double** patchDispOld{nullptr};      //Displacement
    //bool**   patchPointUpdated{nullptr}; //Keep track of updated pathcPoints

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
    
    int* ca_isDynamicFvMesh{nullptr};
    char* ca_dynamicFvMeshType{nullptr};
    char* ca_dynamicSolverType{nullptr};
    //-------------------------------------------
    
    // cellZones
    std::string* cellZonesTypeStr{nullptr};
    std::string* cellZonesNameStr{nullptr};

    int* ca_nCellZones{nullptr};
    int* ca_cellZonesCount{nullptr};
    int** ca_cellZonesList{nullptr};

    int* ca_cellZonesTypeMaxLength{nullptr};
    int* ca_cellZonesNameMaxLength{nullptr};
    char* ca_cellZonesType{nullptr};
    char* ca_cellZonesName{nullptr};
    //-------------------------------------------

    // faceZones
    std::string* faceZonesTypeStr{nullptr};
    std::string* faceZonesNameStr{nullptr};

    int* ca_nFaceZones{nullptr};
    int* ca_faceZonesCount{nullptr};
    int** ca_faceZonesList{nullptr};

    int* ca_faceZonesTypeMaxLength{nullptr};
    int* ca_faceZonesNameMaxLength{nullptr};
    char* ca_faceZonesType{nullptr};
    char* ca_faceZonesName{nullptr};
    //-------------------------------------------

    // pointZones
    std::string* pointZonesTypeStr{nullptr};
    std::string* pointZonesNameStr{nullptr};

    int* ca_nPointZones{nullptr};
    int* ca_pointZonesCount{nullptr};
    int** ca_pointZonesList{nullptr};

    int* ca_pointZonesTypeMaxLength{nullptr};
    int* ca_pointZonesNameMaxLength{nullptr};
    char* ca_pointZonesType{nullptr};
    char* ca_pointZonesName{nullptr};
    //-------------------------------------------

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

    void verbose_message(std::string message, bool parrallel=false);
    void message(std::string message, bool parrallel=false);
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
    void compareWarningExit(
            const int& val1,
            const int& val2,
            const string& name1,
            const string& name2);

    size_t findChar
    (
        const std::string& fullAddr,
        const std::string& content,
        const std::string& exp,
        size_t start = 0,
        size_t end = std::string::npos
    );

    size_t findWord
    (
        const std::string& fullAddr,
        const std::string& content,
        const std::string& exp,
        size_t start = 0,
        size_t end = std::string::npos
    );

    size_t findWordOnly
    (
        const std::string& fullAddr,
        const std::string& content,
        const std::string& exp,
        size_t start = 0,
        size_t end = std::string::npos
    );

};




#endif

