#ifndef COMFOAM_
#define COMFOAM_

#include "com.h"
//#include "com_devel.hpp"
#include "rocRhoCentral.H"


using namespace COM;

class comFoamModule : public rocRhoCentral, public COM_Object
{
public:
    std::string windowName; /// Tracks *this* window name.
    std::string globalName;


    //  Function Handlers ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    int flowInitHandle;
    int flowStatHandle;
    int flowLoopHandle;
    int flowFinHandle;

    int verbosity;    
    

    //  Status Variables ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    int numDataItems=0;
    std::vector<std::string> dataItemNames;
    int numPanes;
    int* paneList;

    double* Coord;
    int numNodes = 0;

    int numConn;
    int numElem;
    int* Conn;
    std::vector<std::string> connNames;

    char getDataItemLoc;
    COM_Type getDataItemType;
    int numElementNodes;
    std::string getDataItemUnits;

    int timeArrayLength;


    //  Object Functions ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    static void Load(const char *name);
    static void Unload(const std::string &name, comFoamModule &object);

    int comInitialize(int argc, char *argv[]);
    int comStatus();
    int comLoop();
    int comFinalize();

    int flowInit(int *pargc, void **pargv, int *verbIn);
    int flowStat();
    int flowLoop();
    int flowFin();


    //  Constructor/Deconstructor ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    comFoamModule();
    comFoamModule(int argc, char *argv[]);

    ~comFoamModule()
    {
        std::cout << "comFoam decunstructor = "
         << test_bool << " " << listOptions_ << std::endl;
    };


private:


protected:
    // ====== Data to register with IMPACT =======
    //  Surface Mesh
    /*
    int numPointsSurface;
    int numElementsSurface;
    std::vector<double> surfaceCoordinates;
    std::vector<unsigned int> surfaceConnectivity;
    std::map<int, int> surfaceNodeMap;
    std::map<int, int> interfaceToFoamNodeMap;
    std::map<int, int> foamGlobalToPatchNodeMap;
    std::vector<double> surfacePressure;
    std::vector<double> surfaceTraction;
    std::vector<double> surfaceDisplacement;
    std::vector<double> solidDisplacement;
    std::vector<double> time;
    std::vector<double> deltaTime;
    */
    
    


    // function status variables
    /*
    std::vector<int> initStatus;
    std::vector<int> runStatus;
    */

    // timing variables
    // std::vector<double> endTime;

    // ====== OpenFOAM stuff ========
    //
    //    General, init stuff
    // Foam::argList *argsPtr;
    // Foam::Time *runTimePtr;
    //
    //    Fluid data
    /*
    autoPtr<dynamicFvMesh> meshPtr;
    IOdictionary *transportPropertiesPtr;
    dimensionedScalar *nuPtr;
    dimensionedScalar *rhoFluidPtr;
    volScalarField *pPtr;
    volVectorField *UPtr;
    surfaceScalarField *phiPtr;
    label pRefCell;
    scalar pRefValue;
    */   
};

#endif
