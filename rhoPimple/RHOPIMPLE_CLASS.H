#ifndef RHOPIMPLE_
#define RHOPIMPLE_

#include "FOAM_CLASS.H"

#include "fluidThermo.H"
#include "pimpleControl.H"
#include "pressureControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"

class RHOPIMPLE_CLASS : public FOAM_CLASS
{
private:

protected:

public:


   //  Solver pointers ^^^^^^^^^^^^^^^^^^^^^
   pimpleControl *pimplePtr;
   
   autoPtr<fluidThermo> pThermoPtr;
   
   pressureControl *pressureControlPtr;
   volScalarField *dpdtPtr;
   volScalarField *KPtr;
   fv::options *fvOptionsPtr;
   autoPtr<surfaceVectorField> rhoUfPtr;
   autoPtr<volScalarField> divrhoUPtr;
   IOMRFZoneList *MRFPtr;
   //--------------------------------------
      
   //  Control variables ^^^^^^^^^^^^^^^^^^
   bool correctPhi;
   bool checkMeshCourantNo;
   bool moveMeshOuterCorrectors;

   scalar cumulativeContErr = 0;
   //--------------------------------------
   
   //  Constructor, Deconstructor ^^^^^^^^^
   RHOPIMPLE_CLASS();
   RHOPIMPLE_CLASS(int argc,char *argv[]); 
   ~RHOPIMPLE_CLASS() {}
   //--------------------------------------

   //  Main/virtual functions ^^^^^^^^^^^^^^^^^^^^^
   int Initialize(int argc,char *argv[]);
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
