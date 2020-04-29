#ifndef ROCFOAM_
#define ROCFOAM_

#include "dynamicFvMesh.H"
#include "fvCFD.H"
#include "fvcSmooth.H"
#include "localEulerDdtScheme.H"
#include "turbulentFluidThermoModel.H"
#include "PstreamGlobals.H"

class rocFoam
{
public:
    //  Initializations ^^^^^^^^^^^^^^
    Foam::argList *argsPtr;
    Foam::Time *runTimePtr;
    //--------------------------------

    //  Constructor/deconstructor ^^^^
    rocFoam();
    virtual ~rocFoam();
    //--------------------------------

    //  Methods status ^^^^^^^^^^^^^^^^^^^^^
    int initializeStat;
    int loopStat;
    int stepStat;
    int finalizeStat;
    
    double testStat;
    virtual double errorEvaluate(int argc, char *argv[]) { return 0.0; };
    //--------------------------------------


private:
    int initSet();

    char *solverType;
protected:
    //  Control variables ^^^^^^^^^^^^
    // static scalar constexpr small=1.0e-10;
    // static scalar constexpr great=1.0e+10;

    bool listOptions;
    bool LTS;
    bool adjustTimeStep;
    bool overrideTimeStep;

    scalar maxCo;
    scalar maxDeltaT;
    scalar overrideDeltaT;
    scalar CoNum;
    scalar meanCoNum;
    //--------------------------------

    //  Solver pointers ^^^^^^^^^^^^^^
    volScalarField *pPtr;
    const volScalarField *TPtr;
    const volScalarField *psiPtr;

    volScalarField *ePtr;
    volScalarField *rhoPtr;
    volVectorField *UPtr;
    volVectorField *rhoUPtr;
    volScalarField *rhoEPtr;
    surfaceScalarField *phiPtr;

    autoPtr<dynamicFvMesh> meshPtr;
    autoPtr<compressible::turbulenceModel> turbulencePtr;
    tmp<volScalarField> trDeltaT;
    //--------------------------------

    //  Main functions ^^^^^^^^^^^^^^^^^^^^^
    virtual int initialize() { return 0; };
    virtual int initialize(int argc, char *argv[]) { return 0; };
    virtual int loop() { return 0; };
    virtual int step() { return 0; };
    virtual int finalize();
    //--------------------------------------

    //  Main functions ^^^^^^^^^^^^^^^^^^^^^
    virtual int createFields() { return 0; };
    virtual int setRDeltaT() { return 0; };
    virtual int createControl() { return 0; };
    //--------------------------------------

    //  Execution functions ^^^^^^^^^^^^^^^^
    int createArgs(int argc, char *argv[]);

    int PostProcess(int argc, char *argv[]);
    int addRegionOption();
    int addDictOption();
    int addFunctionObjectOptions();
    int setRootCase();

    int setRootCaseLists();
    int listOptionsFunc();
    int listOutput();

    int createTime();
    int createTimeControls();
    int readTimeControls();

    int setDeltaT();
    int setDeltaT(scalar &overrideDeltaT);

    // int createMesh();
    int createDynamicFvMesh();

    int createRDeltaT();
    //--------------------------------------
};

#endif
