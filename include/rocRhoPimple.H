#ifndef ROCRHOPIMPLE_
#define ROCRHOPIMPLE_

#include "rocFoam.H"

#include "CorrectPhi.H"
#include "fluidThermo.H"
#include "fvOptions.H"
#include "pimpleControl.H"
#include "pressureControl.H"

class rocRhoPimple : public rocFoam
{
private:
protected:
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

    scalar cumulativeContErr = 0;
    //--------------------------------------

    //  Constructor, Deconstructor ^^^^^^^^^
    rocRhoPimple();
    rocRhoPimple(int argc, char *argv[]);
    ~rocRhoPimple();
    //--------------------------------------

    //  Main/virtual functions ^^^^^^^^^^^^^^^^^^^^^
    int initialize(int argc, char *argv[]);
    int loop();
    int finalize();

    int createControl();
    //--------------------------------------

    //  Execution functions ^^^^^^^^^^^^^^^^
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
};

#endif
