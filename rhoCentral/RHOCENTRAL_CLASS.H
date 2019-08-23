#ifndef RHOCENTRAL_
#define RHOCENTRAL_

#include "FOAM_CLASS.H"

#include "psiThermo.H"
#include "fixedRhoFvPatchScalarField.H"
#include "directionInterpolate.H"


class RHOCENTRAL_CLASS : public FOAM_CLASS
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
   RHOCENTRAL_CLASS();
   RHOCENTRAL_CLASS(int argc,char *argv[]); 
   ~RHOCENTRAL_CLASS() {}
   //--------------------------------------

   //  Main functions ^^^^^^^^^^^^^^^^^^^^^
   int Initialize(int argc,char *argv[]);
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
