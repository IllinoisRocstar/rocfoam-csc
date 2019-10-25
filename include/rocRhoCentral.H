#ifndef ROCRHOCENTRAL_
#define ROCCENTRAL_

#include "rocFoam.H"

#include "directionInterpolate.H"
#include "fixedRhoFvPatchScalarField.H"
#include "psiThermo.H"

class rocRhoCentral : public rocFoam
{
private:
protected:
public:
    //  Solver pointers ^^^^^^^^^^^^^^^^^^^^^
    surfaceScalarField *posPtr;
    surfaceScalarField *negPtr;
    surfaceScalarField *amaxSfPtr;

    autoPtr<psiThermo> pThermoPtr;
    //--------------------------------------

    //  Control variables ^^^^^^^^^^^^^^^^^^
    word fluxScheme;
    bool inviscid;
    //--------------------------------------

    //  Constructor, Deconstructor ^^^^^^^^^
    rocRhoCentral();
    rocRhoCentral(int argc, char *argv[]);
    ~rocRhoCentral() {}
    //--------------------------------------

    //  Main functions ^^^^^^^^^^^^^^^^^^^^^
    int initialize(int argc, char *argv[]);
    int loop();
    int finalize();
    //--------------------------------------

    //  Execution functions ^^^^^^^^^^^^^^^^
    int createFields();
    int createFieldRefs();

    int readFluxScheme();
    int readTimeControls();

    int centralCourantNo();
    int setRDeltaT();
    //--------------------------------------
};

#endif
