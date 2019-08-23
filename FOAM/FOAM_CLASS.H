#ifndef FOAM_
#define FOAM_


#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "turbulentFluidThermoModel.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

class FOAM_CLASS
{
private:

protected:

public:
   //  Initializations ^^^^^^^^^^^^^^
   Foam::argList *argsPtr;
   Foam::Time *runTimePtr;
   //--------------------------------

   //  Control variables ^^^^^^^^^^^^
   static constexpr scalar small=1e-10;
   static constexpr scalar great=1e+10;

   bool listOptions;
   bool LTS;
   bool adjustTimeStep;

   scalar maxCo;
   scalar maxDeltaT;
   scalar CoNum;
   scalar meanCoNum;
   //--------------------------------

   //  Solver pointers ^^^^^^^^^^^^^^
   volScalarField *pPtr;
   const volScalarField *TPtr;
   const volScalarField *psiPtr;
   
   volScalarField *ePtr;
   volScalarField *rhoPtr;
   volVectorField *UPtr;
   volVectorField *rhoUPtr;
   volScalarField *rhoEPtr;
   surfaceScalarField *phiPtr;

   autoPtr<dynamicFvMesh> meshPtr;   
   autoPtr<compressible::turbulenceModel> turbulencePtr;
   tmp<volScalarField> trDeltaT;
   //--------------------------------

   //  Constructor/deconstructor ^^^^
   FOAM_CLASS();
   virtual ~FOAM_CLASS();
   //--------------------------------


   //  Main functions ^^^^^^^^^^^^^^^^^^^^^
   virtual int Initialize() { return 0;};
   virtual int loop() { return 0;};
   virtual int finalize() { return 0;};
   virtual int createFields() { return 0;};
   virtual int setRDeltaT() { return 0;};
   virtual int createControl();
   //--------------------------------------

   //  Execution functions ^^^^^^^^^^^^^^^^
   int PostProcess(int argc,char *argv[]);
   int addRegionOption();
   int addDictOption();
   int addFunctionObjectOptions();
   int setRootCase();
   

   int setRootCaseLists();
   int listOptions_Func();
   int listOutput();
   

   int createTime();
   int createTimeControls();
   //int createMesh();
   int createDynamicFvMesh();

   int createRDeltaT();
   int setDeltaT();
   //--------------------------------------

};

#endif
