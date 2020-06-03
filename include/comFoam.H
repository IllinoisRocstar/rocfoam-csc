#ifndef COMFOAM_
#define COMFOAM_

#include "rocFoam.H"
#include "com.h"
#include "fileContainer.H"
#include <filesystem.hpp>
#include <iomanip>

class comFoam : public COM_Object, public rocFoam
{
private:
    const int faceToPointTypeSize{4};
    const int panePadding{100};
    const int IODigits{10};
    const std::string tmpFluidDir{"./fluidTmp"};
protected:
    const int nComponents{3};
    const int genCharSize{15};

public:
    //  Constructor/Deconstructor ^^^^^^^^^^^^^^^^^^^^^^^^^
    comFoam();
    //comFoam(int *pargc, void **pargv, const char *name);
    ~comFoam();
    
protected:
    virtual int finalize();

    // Methods ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
protected:
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

private:
    // tempFiles-specific
    int readInitFiles(const std::string& rootAddr);
    int createInitFiles(const std::string& rootAddr);
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
    int deleteFilesData();
    int registerInitFiles(const char *name);
    int reconstCaInitFiles(const char *name, const std::string& rootAddr);
    
    std::string createBaseScalar(std::string name);

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

protected:
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
    double* ca_Vel{nullptr};
    double* ca_P{nullptr};
    double* ca_T{nullptr};
    double* ca_Rho{nullptr};
    double* ca_Phi{nullptr};
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
    wordList**    ca_patchInGroup{nullptr}; //single value for the last
    int**         ca_patchStart{nullptr}; //single value for the last
    int**         ca_patchSize{nullptr}; //single value for the last

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
    
    double** ca_patchSf{nullptr};
    //-------------------------------------------
    
    // Files
    int*   ca_nFiles{nullptr};
    int*   ca_fileSize{nullptr};
    char** ca_fileName{nullptr};
    char** ca_filePath{nullptr};
    char** ca_fileContent{nullptr};

    //-------------------------------------------

    //  Window data ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    //std::string winVolName; /// Tracks *this* volume window name.
    //std::string winSurfName; /// Tracks *this* volume window name.

    std::string solverType{""};
    MPI_Comm winComm{};
    
    int ca_nProc{0};
    int ca_myRank{0};

    // registered data set during the simulation
    int*    ca_runStat{nullptr};
    int*    ca_timeIndex{nullptr};
    double* ca_time{nullptr};
    char*   ca_timeName{nullptr};
    double* ca_deltaT{nullptr};
    double* ca_deltaT0{nullptr};
    
    
    // Methods ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    // rocFoam driver methods
    int flowInit(int *pargc, void **pargv, const char *name);
    int flowLoop();
    int flowStep();

    // COM-specific
    int registerFunctions(const char *name);
    int reconstCaData(int* pargc, void** pargv, const char* name);
    //-------------------------------------------

public:

    // RocStar-specific methods ^^^^^^^^^^^^^^^^^
    void flowInitRocStar
    (
        const double& initTime,
        const MPI_Comm& flowComm,
        const int& manInitHandle,
        const std::string& volName,
        const std::string& surfName,
        const int& obtainHandle
    );

    void flowUpdateRocStar
    (
        double& currentTime,
        double& timeStep,
        int& handles
    );

    void flowFinRocStar();
    //-------------------------------------------

    static void copyWindow(const char *name1, const char *name2);
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
        const std::string& exp
    );

    size_t findOnlyLoc
    (
        const std::string& fullAddr,
        const std::string& content,
        const std::string& exp,
        size_t start
    );

    size_t findOnlyLoc
    (
        const std::string& fullAddr,
        const std::string& content,
        const std::string& exp,
        size_t start,
        size_t end
    );

};

#endif

