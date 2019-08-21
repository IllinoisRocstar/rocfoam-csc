#ifndef ROCFOAM_
#define ROCFOAM_

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "fixedRhoFvPatchScalarField.H"
#include "directionInterpolate.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"


class RHOCENTRAL_CLASS
{
private:

protected:

public:

   //  Initializations ^^^^^^^^^^^^^^
   Foam::argList *argsPtr;
   Foam::Time *runTimePtr;
   //--------------------------------
   
   //  Solver pointers ^^^^^^^^^^^^^^
   volScalarField *pPtr;
   const volScalarField *TPtr;
   const volScalarField *psiPtr;
   //const volScalarField *muPtr;
   //tmp<volScalarField> muPtr;
   
   volScalarField *ePtr;
   volScalarField *rhoPtr;
   volVectorField *UPtr;
   volVectorField *rhoUPtr;
   volScalarField *rhoEPtr;
   surfaceScalarField *posPtr;
   surfaceScalarField *negPtr;
   surfaceScalarField *phiPtr;
   surfaceScalarField *amaxSfPtr;

   autoPtr<psiThermo> thermoPtr;
   autoPtr<dynamicFvMesh> meshPtr;   
   autoPtr<compressible::turbulenceModel> turbulencePtr;
   tmp<volScalarField> trDeltaT;
   //--------------------------------
      
   //  Control variables ^^^^^^^^^^^^
   const scalar small=1e-10;
   
   word fluxScheme;
   bool inviscid;
   bool adjustTimeStep;
   bool listOptions;
   bool LTS;
   scalar maxCo;
   scalar maxDeltaT;
   scalar CoNum;
   scalar meanCoNum;
   //--------------------------------
   
   //  Constructor, Deconstructor ........ 
   RHOCENTRAL_CLASS();
   RHOCENTRAL_CLASS(int argc,char *argv[]); 
   ~RHOCENTRAL_CLASS(){}

   //  Main functions ....................
   int Initialize(int argc,char *argv[]);
   int loop();
   int finalize();

   //  Execution functions ...............
   int PostProcess(int argc,char *argv[]);
   int addRegionOption();
   int addDictOption();
   int addFunctionObjectOptions();
   int setRootCase();

   int setRootCaseLists();
   int listOptions_Func();
   int listOutput();
   
   int createTime();
   int createMesh();
   int createDynamicFvMesh();
   int createFields();
   int createControl();
   
   int createFieldRefs();
   int createRDeltaT();
   int createTimeControls();
   int readFluxScheme();
   int readTimeControls();
   int setDeltaT();
   int setRDeltaT();
   int centralCourantNo();

};

#endif
