#ifndef RHOCENTRAL_
#define RHOCENTRAL_

#include "comFoam.H"

#include "directionInterpolate.H"
#include "fixedRhoFvPatchScalarField.H"
#include "psiThermo.H"

class rhoCentral : public comFoam
{
public:
    //  Constructor, Deconstructor ^^^^^^^^^
    rhoCentral();
    rhoCentral(int argc, char *argv[]);
    //rhoCentral(int *pargc, void **pargv, int *verbIn);

    ~rhoCentral();
    //--------------------------------------

    //  COM-specific methods ^^^^^^^^^^^^^^^
    static void load(const char *name);
    static void unload(const char *name);
    //--------------------------------------

private:
    //int initSet();

    //  Solver pointers ^^^^^^^^^^^^^^^^^^^^^
    surfaceScalarField *posPtr{nullptr};
    surfaceScalarField *negPtr{nullptr};
    surfaceScalarField *amaxSfPtr{nullptr};

    autoPtr<psiThermo> pThermoPtr{nullptr};
    //--------------------------------------

    //  Control variables ^^^^^^^^^^^^^^^^^^
    word fluxScheme{""};
    bool inviscid{false};
    //--------------------------------------

   
    // Solver methods ^^^^^^^^^^^^^^^^^^^^^^
    int createFields();
    int createFieldRefs();

    int readFluxScheme();

    int centralCourantNo();
    int setRDeltaT();
    double errorEvaluate(int argc, char *argv[]);
    //--------------------------------------

public:
    //  Main methods ^^^^^^^^^^^^^^^^^^^^^^^
    int initialize(int argc, char *argv[]);
    int loop();
    int step();
    int finalize();
    //--------------------------------------
};

#endif
