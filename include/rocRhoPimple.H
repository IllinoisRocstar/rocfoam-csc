#ifndef COMRHOPIMPLE_
#define COMRHOPIMPLE_

#include "comFoam.H"

#include "CorrectPhi.H"
#include "fluidThermo.H"
#include "fvOptions.H"
#include "pimpleControl.H"
#include "pressureControl.H"

class rhoPimple : public comFoam
{
public:
    //  Solver pointers ^^^^^^^^^^^^^^^^^^^^^
    pimpleControl *pimplePtr;
    pressureControl *pressureControlPtr;
    volScalarField *dpdtPtr;
    volScalarField *KPtr;
    fv::options *fvOptionsPtr;
    IOMRFZoneList *MRFPtr;
    fvVectorMatrix *UEqnPtr;

    autoPtr<fluidThermo> pThermoPtr;
    autoPtr<surfaceVectorField> rhoUfPtr;
    autoPtr<volScalarField> divrhoUPtr;
    tmp<fvVectorMatrix> tUEqnPtr;
    //--------------------------------------

    //  Control variables ^^^^^^^^^^^^^^^^^^
    bool correctPhi;
    bool checkMeshCourantNo;
    bool moveMeshOuterCorrectors;

    scalar cumulativeContErr;
    //--------------------------------------


    //  Object Functions ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    static void load(const char *name);
    static void unload(const std::string &name);


    //  Constructor, Deconstructor ^^^^^^^^^
    rhoPimple();
    rhoPimple(int argc, char *argv[]);

    ~rhoPimple();
    //--------------------------------------

    //  Main functions ^^^^^^^^^^^^^^^^^^^^^
    int initialize(int argc, char *argv[]);
    int loop();
    int finalize();
    //--------------------------------------

    //  Execution functions ^^^^^^^^^^^^^^^^
    int createControl();
    int createDyMControls();
    int readTimeControls();
    int initContinuityErrs();

    int createFields();
    int createMRF();
    int compressibleCreatePhi();
    int createFvOptions();

    int createFieldRefs();
    int createRhoUfIfPresent();

    int compressibleCourantNo();
    int setInitialDeltaT();

    int readDyMControls();
    int setRDeltaT();
    int correctPhi_();
    int meshCourantNo();
    int rhoEqn_();
    int UEqn_();
    int EEqn_();
    int pcEqn();
    int compressibleContinuityErrs();
    int pEqn_();
    //--------------------------------------

private:

protected:
    //int flowInit(int *pargc, void **pargv, int *verbIn);
    //int flowLoop();

    //int flowStat();
    //int flowFin();
};

#endif
