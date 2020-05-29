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
    //  Constructor, Deconstructor ^^^^^^^^^
    rhoPimple();
    rhoPimple(int argc, char *argv[]);

    ~rhoPimple();
    //--------------------------------------

    //  COM-specific methods ^^^^^^^^^^^^^^^
    static void load(const char *name);
    static void unload(const char *name);
    //--------------------------------------

private:
    //int initSet();

    //  Solver pointers ^^^^^^^^^^^^^^^^^^^^
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

    //  Solver methods  ^^^^^^^^^^^^^^^^^^^^
    int createControl();
    int createDyMControls();
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

    double errorEvaluate(int argc, char *argv[]);
    //--------------------------------------
    
public:
    //  Main functions ^^^^^^^^^^^^^^^^^^^^^
    int initialize(int argc, char *argv[]);
    int loop();
    int step();
    int finalize();
    //--------------------------------------
};

#endif
