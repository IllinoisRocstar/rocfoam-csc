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

    pimpleControl& getPimpleControl()
    {
        return *pimplePtr;
    }

private:
    pimpleControl* pimplePtr{nullptr};
    //int initSet();

    //  Solver pointers ^^^^^^^^^^^^^^^^^^^^
    pressureControl* pressureControlPtr{nullptr};
    volScalarField* dpdtPtr{nullptr};
    volScalarField* KPtr{nullptr};
    fv::options* fvOptionsPtr{nullptr};
    IOMRFZoneList* MRFPtr{nullptr};
    fvVectorMatrix* UEqnPtr{nullptr};

    autoPtr<fluidThermo> pThermoPtr{nullptr};
    autoPtr<volScalarField> divrhoUPtr{nullptr};
    tmp<fvVectorMatrix> tUEqnPtr{nullptr};

#ifdef HAVE_OF8
    dimensionedScalar* initialMassPtr{nullptr};
#endif
    //--------------------------------------

    //  Control variables ^^^^^^^^^^^^^^^^^^
    bool correctPhi{false};
    bool checkMeshCourantNo{false};
    bool moveMeshOuterCorrectors{false};

    scalar cumulativeContErr{0};

    bool modifiedDeltaT{false};
    double unmodifiedDeltaTvalue{0};
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

#ifdef HAVE_OF7
    int pcEqn();
#endif

    int compressibleContinuityErrs();
    int incompressibleContinuityErrs();
    int pEqn_();

    double errorEvaluate(int argc, char *argv[]);
    //--------------------------------------

protected:
    int loadInternal(const char* name);
    
public:
    //  Main functions ^^^^^^^^^^^^^^^^^^^^^
    int initFOAM(int argc, char *argv[]);
    int loop();
    int step(double* newDeltaT = nullptr,
             int* handle = nullptr);
    int finalizeFoam();
    //--------------------------------------
};

#endif
