#ifndef COMRHOCENTRAL_
#define COMRHOCENTRAL_

//#include "com.h"

#include "comFoam.H"
#include "rocFoam.H"

#include "directionInterpolate.H"
#include "fixedRhoFvPatchScalarField.H"
#include "psiThermo.H"

class comRhoCentralModule : public comFoamModule
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
    static void Load(const char *name);
    static void Unload(const std::string &name);


    //  Constructor, Deconstructor ^^^^^^^^^
    comRhoCentralModule();
    comRhoCentralModule(int argc, char *argv[]);
    comRhoCentralModule(int *pargc, void **pargv, int *verbIn);
    ~comRhoCentralModule();
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
    int flowInit(int *pargc, void **pargv, int *verbIn);
    int flowLoop();

    //int flowStat();
    //int flowFin();
};

#endif
