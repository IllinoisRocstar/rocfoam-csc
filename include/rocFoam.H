#ifndef ROCFOAM_
#define ROCFOAM_

#include "dynamicFvMesh.H"
#include "fvCFD.H"
#include "fvcSmooth.H"
#include "localEulerDdtScheme.H"

#ifdef HAVE_OF7
    #include "turbulentFluidThermoModel.H"
#elif defined(HAVE_OF8)
    #include "fluidThermoMomentumTransportModel.H"
    #include "fluidThermophysicalTransportModel.H"
#endif

#include "PstreamGlobals.H"

class rocFoam
{
public:
    //  Initializations ^^^^^^^^^^^^^^
    Foam::argList *argsPtr{nullptr};
    Foam::Time *runTimePtr{nullptr};
    //--------------------------------

    //  Constructor/deconstructor ^^^^
    rocFoam();
    virtual ~rocFoam();
    //--------------------------------

    //  Methods status ^^^^^^^^^^^^^^^
    int initializeStat{0};
    int loopStat{0};
    int stepStat{0};
    int finalizeStat{0};
    
    double testStat;
    virtual double errorEvaluate(int argc, char *argv[]) { return 0.0; };
    //--------------------------------

protected:
    //  Control variables ^^^^^^^^^^^^
    // static scalar constexpr small=1.0e-10;
    // static scalar constexpr great=1.0e+10;

    bool listOptions{false};
    bool LTS{false};
    bool adjustTimeStep{false};
    //bool overrideTimeStep{false};

    scalar maxCo{0};
    scalar maxDeltaT{0};
    //scalar overrideDeltaT{0};
    scalar CoNum{0};
    scalar meanCoNum{0};
    //--------------------------------

    //  Solver pointers ^^^^^^^^^^^^^^
    volScalarField* pPtr{nullptr};
    const volScalarField* TPtr{nullptr};
    const volScalarField* psiPtr{nullptr};

    volScalarField* ePtr{nullptr};
    volScalarField* rhoPtr{nullptr};
    volVectorField* UPtr{nullptr};
    volVectorField* rhoUPtr{nullptr};
    volScalarField* rhoEPtr{nullptr};
    surfaceScalarField* phiPtr{nullptr};
    autoPtr<surfaceVectorField> rhoUfPtr{nullptr};
    
    autoPtr<dynamicFvMesh> meshPtr{nullptr};

#ifdef HAVE_OF7
    autoPtr<compressible::turbulenceModel> turbulencePtr{nullptr};
#elif  defined(HAVE_OF8)
    autoPtr<compressible::momentumTransportModel> turbulencePtr{nullptr};
    autoPtr<fluidThermophysicalTransportModel> thermophysicalTransportPtr{nullptr};
#endif

    tmp<volScalarField> trDeltaT{nullptr};
    //--------------------------------

    //  Main functions ^^^^^^^^^^^^^^^^^^^^^
    virtual int initFOAM() { return 0; };
    virtual int initFOAM(int argc, char *argv[]) { return 0; };
    virtual int loop() { return 0; };
    //virtual int step(double* newDeltaT = nullptr) { return 0; };
    virtual int step(double* newDeltaT = nullptr,
                     int* handle = nullptr){ return 0; };
    virtual int finalizeFoam();
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
    //int setDeltaT(double* newDeltaT = nullptr);

    // int createMesh();
    int createDynamicFvMesh();

    int createRDeltaT();
    //--------------------------------------
};

#endif
