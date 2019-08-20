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

   // ====== OpenFOAM stuff ========
   //
   //    General, init stuff
   Foam::argList *argsPtr;
   Foam::Time *runTimePtr;
   
   // postProcess usses fvMesh while the solver uses
   // dynamicFvMesh. For the moment, use the latter.
   //Foam::fvMesh *meshPtr;
   
   //Foam::dynamicFvMesh *meshPtr;
   autoPtr<dynamicFvMesh> meshPtr;
   
   
   //  Reference variables. For the moment, using pointers insted
   /*
   volScalarField& p;
   const volScalarField& T;
   const volScalarField& psi;
   const volScalarField& mu;
   */

   volScalarField *pPtr;
   const volScalarField *TPtr;
   const volScalarField *psiPtr;
   //const volScalarField *muPtr;
   //tmp<volScalarField> muPtr;
   

   //psiThermo *thermoPtr;
   autoPtr<psiThermo> thermoPtr;
   
   volScalarField *ePtr;
   volScalarField *rhoPtr;
   volVectorField *UPtr;
   volVectorField *rhoUPtr;
   volScalarField *rhoEPtr;
   surfaceScalarField *posPtr;
   surfaceScalarField *negPtr;
   surfaceScalarField *phiPtr;
   
   surfaceScalarField *amaxSfPtr = NULL;
   
   //compressible::turbulenceModel *turbulencePtr;
   autoPtr<compressible::turbulenceModel> turbulencePtr;
   tmp<volScalarField> trDeltaT;
      
   word fluxScheme;
   
   bool inviscid;
   bool adjustTimeStep;
   scalar maxCo;
   scalar maxDeltaT;
   
   scalar CoNum=0.0;
   scalar meanCoNum=0.0;
   scalar small=1e-10;
   
   
   bool listOptions = false;
   bool LTS = false;


   RHOCENTRAL_CLASS();
   RHOCENTRAL_CLASS(int argc,char *argv[]); 

   int loop();
   
   int Initialize(int argc,char *argv[]);
   int PostProcess(int argc,char *argv[]);
   int addRegionOption();
   int addDictOption();
   int addFunctionObjectOptions();
   int setRootCase(int argc,char *argv[]);

   int setRootCaseLists(int argc,char *argv[]);
   int listOptions_Func();
   int listOutput();
   
   int createTime(int argc,char *argv[]);
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

   ~RHOCENTRAL_CLASS()
   {
      delete pPtr;
      delete TPtr;
      delete psiPtr;
      //delete muPtr;
      
      //delete thermoPtr;
      
      delete ePtr;
      delete rhoPtr;
      delete UPtr;
      delete rhoUPtr;
      delete rhoEPtr;
      delete posPtr;
      delete negPtr;
      delete phiPtr;
      
       delete amaxSfPtr;
      
      //delete turbulencePtr;
      //delete trDeltaT;
   }

};

#endif
