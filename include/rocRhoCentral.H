#ifndef RHOCENTRAL_
#define RHOCENTRAL_

#include "comFoam.H"

#include "directionInterpolate.H"
#include "fixedRhoFvPatchScalarField.H"
#include "psiThermo.H"

class rhoCentral : public comFoam
{
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

    //  Object Functions ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    static void load(const char *name);
    static void unload(const std::string &name);


    //  Constructor, Deconstructor ^^^^^^^^^
    rhoCentral();
    rhoCentral(int argc, char *argv[]);
    //rhoCentral(int *pargc, void **pargv, int *verbIn);

    ~rhoCentral();
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

private:

protected:
    //int flowInit(int *pargc, void **pargv, int *verbIn);
    //int flowLoop();

    //int flowStat();
    //int flowFin();
};

#endif
