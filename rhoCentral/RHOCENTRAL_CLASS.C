#include "RHOCENTRAL_CLASS.H"

RHOCENTRAL_CLASS::RHOCENTRAL_CLASS() {}
RHOCENTRAL_CLASS::RHOCENTRAL_CLASS(int argc,char *argv[])
{
   Initialize(argc, argv);
}


int RHOCENTRAL_CLASS::Initialize(int argc,char *argv[])
{
   #define NO_CONTROL

   // Mohammad: Not quite sure where this line should be
   argsPtr = new Foam::argList(argc, argv);

   //  postProcess.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   PostProcess(argc, argv);
   // ---------------------------------------------------

   //  setRootCaseLists.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   setRootCaseLists(argc, argv);
   // ---------------------------------------------------

   //  createTime.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   createTime(argc, argv);
   // ---------------------------------------------------

   //  createDynamicFvMesh.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^
   createDynamicFvMesh();
   // ---------------------------------------------------

   //  createFields.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   createFields();
   
   // ---------------------------------------------------

   //  createFieldRefs.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   createFieldRefs();
   // ---------------------------------------------------

   //  createTimeControls.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   createTimeControls();
   // ---------------------------------------------------

   //  readFluxScheme.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   readFluxScheme();
   // ---------------------------------------------------

   turbulencePtr->validate();

   Foam::Info << "End of initialization of openFoamPar module." << endl;

   return 0;
}


int RHOCENTRAL_CLASS::loop()
{

    dimensionedScalar v_zero("v_zero", dimVolume/dimTime, 0.0);

    // Courant numbers used to adjust the time-step
    //scalar CoNum = 0.0;
    //scalar meanCoNum = 0.0;

    Info<< "\nStarting time loop\n" << endl;

    while (runTimePtr->run())
    {
        //  readTimeControls.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        readTimeControls();

Foam::Info << "LTS = " << LTS << endl;


        // ---------------------------------------------------

        if (!LTS)
        {
            //  setDeltaT.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            setDeltaT();
            // ---------------------------------------------
            &(*runTimePtr)++;

            // Do any mesh changes
            meshPtr->update();
        }

        // --- Directed interpolation of primitive fields onto faces

        surfaceScalarField rho_pos(interpolate(*rhoPtr, *posPtr));
        surfaceScalarField rho_neg(interpolate(*rhoPtr, *negPtr));

        surfaceVectorField rhoU_pos(interpolate(*rhoUPtr, *posPtr, UPtr->name()));
        surfaceVectorField rhoU_neg(interpolate(*rhoUPtr, *negPtr, UPtr->name()));

        volScalarField rPsi("rPsi", 1.0/ *psiPtr);
        surfaceScalarField rPsi_pos(interpolate(rPsi, *posPtr, TPtr->name()));
        surfaceScalarField rPsi_neg(interpolate(rPsi, *negPtr, TPtr->name()));

        surfaceScalarField e_pos(interpolate(*ePtr, *posPtr, TPtr->name()));
        surfaceScalarField e_neg(interpolate(*ePtr, *negPtr, TPtr->name()));

        surfaceVectorField U_pos("U_pos", rhoU_pos/rho_pos);
        surfaceVectorField U_neg("U_neg", rhoU_neg/rho_neg);

        surfaceScalarField p_pos("p_pos", rho_pos*rPsi_pos);
        surfaceScalarField p_neg("p_neg", rho_neg*rPsi_neg);

        surfaceScalarField phiv_pos("phiv_pos", U_pos & meshPtr->Sf());
        surfaceScalarField phiv_neg("phiv_neg", U_neg & meshPtr->Sf());

        // Make fluxes relative to mesh-motion
        if (meshPtr->moving())
        {
            phiv_pos -= meshPtr->phi();
            phiv_neg -= meshPtr->phi();
        }

        volScalarField c("c", sqrt(thermoPtr->Cp() / thermoPtr->Cv() * rPsi));
        surfaceScalarField cSf_pos
        (
            "cSf_pos",
            interpolate(c, *posPtr, TPtr->name()) * meshPtr->magSf()
        );

        surfaceScalarField cSf_neg
        (
            "cSf_neg",
            interpolate(c, *negPtr, TPtr->name()) * meshPtr->magSf()
        );

        surfaceScalarField ap
        (
            "ap",
            max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero)
        );
        surfaceScalarField am
        (
            "am",
            min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero)
        );

        surfaceScalarField a_pos("a_pos", ap/(ap - am));

        //surfaceScalarField amaxSf("amaxSf", max(mag(am), mag(ap)));
        if (amaxSfPtr == NULL)
        {
            amaxSfPtr = new surfaceScalarField("amaxSf", max(mag(am), mag(ap)));
        }
        //else
        //{
            *amaxSfPtr = max(mag(am), mag(ap));
        //}

        surfaceScalarField aSf("aSf", am*a_pos);

        if (fluxScheme == "Tadmor")
        {
            aSf = -0.5 * *amaxSfPtr;
            a_pos = 0.5;
        }

        surfaceScalarField a_neg("a_neg", 1.0 - a_pos);

        phiv_pos *= a_pos;
        phiv_neg *= a_neg;

        surfaceScalarField aphiv_pos("aphiv_pos", phiv_pos - aSf);
        surfaceScalarField aphiv_neg("aphiv_neg", phiv_neg + aSf);

        // Reuse amaxSf for the maximum positive and negative fluxes
        // estimated by the central scheme
        *amaxSfPtr = max(mag(aphiv_pos), mag(aphiv_neg));

        //  centralCourantNo.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        centralCourantNo();
        // ---------------------------------------------------

        if (LTS)
        {
            // setRDeltaT.H
            // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            setRDeltaT();
            // -------------------------------
            &(*runTimePtr)++;
        }

        Info<< "Time = " << runTimePtr->timeName() << nl << endl;

        *phiPtr = aphiv_pos*rho_pos + aphiv_neg*rho_neg;

        surfaceVectorField phiUp
        (
            (aphiv_pos*rhoU_pos + aphiv_neg*rhoU_neg)
          + (a_pos*p_pos + a_neg*p_neg) * meshPtr->Sf()
        );

        surfaceScalarField phiEp
        (
            "phiEp",
            aphiv_pos*(rho_pos*(e_pos + 0.5*magSqr(U_pos)) + p_pos)
          + aphiv_neg*(rho_neg*(e_neg + 0.5*magSqr(U_neg)) + p_neg)
          + aSf*p_pos - aSf*p_neg
        );

        // Make flux for pressure-work absolute
        if (meshPtr->moving())
        {
            phiEp += meshPtr->phi()*(a_pos*p_pos + a_neg*p_neg);
        }

        volScalarField muEff("muEff", turbulencePtr->muEff());
        volTensorField tauMC("tauMC", muEff*dev2(Foam::T(fvc::grad(*UPtr))));

        // --- Solve density
        solve(fvm::ddt(*rhoPtr) + fvc::div(*phiPtr));

        // --- Solve momentum
        solve(fvm::ddt(*rhoUPtr) + fvc::div(phiUp));

        //U.ref() = rhoU()/rho()
        //UPtr->ref() = rhoUPtr->ref() / rhoPtr->ref();
        UPtr->ref() = *rhoUPtr / *rhoPtr;

       
        UPtr->correctBoundaryConditions();
        rhoUPtr->boundaryFieldRef() == rhoPtr->boundaryField() * UPtr->boundaryField();

        if (!inviscid)
        {
            solve
            (
                fvm::ddt(*rhoPtr, *UPtr) - fvc::ddt(*rhoPtr, *UPtr)
              - fvm::laplacian(muEff, *UPtr)
              - fvc::div(tauMC)
            );
            *rhoUPtr = *rhoPtr * *UPtr;
        }

        // --- Solve energy
        surfaceScalarField sigmaDotU
        (
            "sigmaDotU",
            (
                fvc::interpolate(muEff) * meshPtr->magSf() * fvc::snGrad(*UPtr)
              + fvc::dotInterpolate(meshPtr->Sf(), tauMC)
            )
          & (a_pos*U_pos + a_neg*U_neg)
        );

        solve
        (
            fvm::ddt(*rhoEPtr)
          + fvc::div(phiEp)
          - fvc::div(sigmaDotU)
        );

        *ePtr = *rhoEPtr / *rhoPtr - 0.5*magSqr(*UPtr);
        ePtr->correctBoundaryConditions();
        thermoPtr->correct();
        rhoEPtr->boundaryFieldRef() ==
            rhoPtr->boundaryField()*
            (
                ePtr->boundaryField() + 0.5*magSqr(UPtr->boundaryField())
            );

        if (!inviscid)
        {
            solve
            (
                fvm::ddt(*rhoPtr, *ePtr) - fvc::ddt(*rhoPtr, *ePtr)
              - fvm::laplacian(turbulencePtr->alphaEff(), *ePtr)
            );
            thermoPtr->correct();
            *rhoEPtr = *rhoPtr * (*ePtr + 0.5*magSqr(*UPtr));
        }

        //pPtr->ref() = *rhoPtr() / *psiPtr();
        //pPtr->ref() = rhoPtr()->ref() / psiPtr()->ref();
        pPtr->ref() = *rhoPtr / *psiPtr;

        pPtr->correctBoundaryConditions();
        rhoPtr->boundaryFieldRef() == psiPtr->boundaryField() * pPtr->boundaryField();

        turbulencePtr->correct();

        runTimePtr->write();

        Info<< "ExecutionTime = " << runTimePtr->elapsedCpuTime() << " s"
            << "  ClockTime = " << runTimePtr->elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;


    return 0;
}








int RHOCENTRAL_CLASS::PostProcess(int argc, char *argv[])
{

   //  createTime.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#ifndef CREATE_TIME
   #define CREATE_TIME createTime(argc, argv);
#endif
   // ---------------------------------------------------

   //  createMesh.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#ifndef CREATE_MESH
   //#define CREATE_MESH createMesh();
   #define CREATE_MESH createDynamicFvMesh();
#endif
   // ---------------------------------------------------

   //  createFields.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#ifndef CREATE_FIELDS
   #define CREATE_FIELDS createFields();
#endif
   // ---------------------------------------------------

   //  createControl.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#ifndef CREATE_CONTROL
   #define CREATE_CONTROL createControl();
#endif
   // ---------------------------------------------------

   #define INCLUDE_FILE(X) INCLUDE_FILE2(X)
   #define INCLUDE_FILE2(X) X

   Foam::argList::addBoolOption
   (
       argList::postProcessOptionName,
       "Execute functionObjects only"
   );

Foam::Info << "FROM POSTPROCESS00" << endl;

   if (argList::postProcess(argc, argv))
   {
      Foam::timeSelector::addOptions();


std::cout << "FROM POSTPROCESS01" << endl;

      //  addRegionOption.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      addRegionOption();
      // ---------------------------------------------------

std::cout << "FROM POSTPROCESS02" << endl;

      //  addFunctionObjectOptions.H  ^^^^^^^^^^^^^^^^^^^^^^
      addFunctionObjectOptions();
      // ---------------------------------------------------

std::cout << "FROM POSTPROCESS03" << endl;

      // Set functionObject post-processing mode
      functionObject::postProcess = true;

std::cout << "FROM POSTPROCESS03" << endl;

      //  setRootCase.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      //setRootCase(argc, argv);
      // --------------------------------------------------

std::cout << "FROM POSTPROCESS04" << endl;


      if (argsPtr->optionFound("list"))
      {
         functionObjectList::list();

std::cout << "FROM POSTPROCESS05" << endl;
         
         return 0;
      }

std::cout << "FROM POSTPROCESS06" << endl;


      INCLUDE_FILE(CREATE_TIME)
      Foam::instantList timeDirs = Foam::timeSelector::select0(*runTimePtr, *argsPtr);
      INCLUDE_FILE(CREATE_MESH)

      #ifndef NO_CONTROL
         INCLUDE_FILE(CREATE_CONTROL)
      #endif

      forAll(timeDirs, timei)
      {
         runTimePtr->setTime(timeDirs[timei], timei);

         Info<< "Time = " << runTimePtr->timeName() << endl;

         FatalIOError.throwExceptions();

         try
         {
            INCLUDE_FILE(CREATE_FIELDS)

            #ifdef CREATE_FIELDS_2
               #include INCLUDE_FILE(CREATE_FIELDS_2)
            #endif

            #ifdef CREATE_FIELDS_3
               #include INCLUDE_FILE(CREATE_FIELDS_3)
            #endif


Foam::Info << "FROM POSTPROCESS" << endl;


            // Externally stored dictionary for functionObjectList
            // if not constructed from runTime
            dictionary functionsControlDict("controlDict");

            HashSet<word> selectedFields;

            // Construct functionObjectList
            autoPtr<functionObjectList> functionsPtr
            (
               functionObjectList::New
               (
                  *argsPtr,
                  *runTimePtr,
                  functionsControlDict,
                  selectedFields
               )
            );

            functionsPtr->execute();
         }
         catch (IOerror& err)
         {
            Warning<< err << endl;
         }

         // Clear the objects owned by the mesh
         meshPtr->objectRegistry::clear();

         Info<< endl;
      }

      Info<< "End\n" << endl;

      return 0;
   }


Foam::Info << "FROM POSTPROCESS000" << endl;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

   #undef INCLUDE_FILE
   #undef INCLUDE_FILE2

   #undef CREATE_MESH
   #undef CREATE_FIELDS
   #undef CREATE_CONTROL

   return 0;
}


int RHOCENTRAL_CLASS::addRegionOption()
{
   Foam::argList::addOption
   (
       "region",
       "name",
       "specify alternative mesh region"
   );
   
   return 0;
}


int RHOCENTRAL_CLASS::addDictOption()
{
   Foam::argList::addOption
   (
       "dict",
       "file",
       "read control dictionary from specified location"
   );
   return 0;
}


int RHOCENTRAL_CLASS::addFunctionObjectOptions()
{

   //  addDictOption.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   addDictOption();
   // --------------------------------------------------
   
   Foam::argList::addOption
   (
       "field",
       "name",
       "Specify the name of the field to be processed, e.g. U"
   );
   Foam::argList::addOption
   (
       "fields",
       "list",
       "Specify a list of fields to be processed, e.g. '(U T p)' - "
       "regular expressions not currently supported"
   );
   Foam::argList::addOption
   (
       "func",
       "name",
       "Specify the name of the functionObject to execute, e.g. Q"
   );
   Foam::argList::addOption
   (
       "funcs",
       "list",
       "Specify the names of the functionObjects to execute, e.g. '(Q div(U))'"
   );
   Foam::argList::addBoolOption
   (
       "list",
       "List the available configured functionObjects"
   );
   
   return 0;
}


int RHOCENTRAL_CLASS::setRootCase(int argc,char *argv[])
{

Foam::Info << "In setRootCase ..." << endl;

   Foam::argList args(argc, argv);
   if (!argsPtr->checkRootCase())
   {
       Foam::FatalError.exit();
   }
   return 0;
}



int RHOCENTRAL_CLASS::setRootCaseLists(int argc,char *argv[])
{

   //  listOptions.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   listOptions_Func();
   // --------------------------------------------------

Foam::Info << "listOptions_Func Done" << endl;

   
   //  setRootCase.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   //setRootCase(argc, argv);
   // --------------------------------------------------

Foam::Info << "setRootCase Done" << endl;


   //  listOutput.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   listOutput();

Foam::Info << "listOutput Done" << endl;   
   
   // --------------------------------------------------

   return 0;
}


int RHOCENTRAL_CLASS::listOptions_Func()
{
   argList::addBoolOption
   (
       "listSwitches",
       "List switches declared in libraries but not set in etc/controlDict"
   );
   argList::addBoolOption
   (
       "listRegisteredSwitches",
       "List switches registered for run-time modification"
   );
   argList::addBoolOption
   (
       "listUnsetSwitches",
       "List switches declared in libraries but not set in etc/controlDict"
   );

   #ifdef fvPatchField_H
   argList::addBoolOption
   (
       "listScalarBCs",
       "List scalar field boundary conditions (fvPatchField<scalar>)"
   );
   argList::addBoolOption
   (
       "listVectorBCs",
       "List vector field boundary conditions (fvPatchField<vector>)"
   );
   #endif

   #ifdef functionObject_H
   argList::addBoolOption
   (
       "listFunctionObjects",
       "List functionObjects"
   );
   #endif

   #ifdef fvOption_H
   argList::addBoolOption
   (
       "listFvOptions",
       "List fvOptions"
   );
   #endif

   #if defined(turbulentTransportModel_H) || defined(turbulentFluidThermoModel_H)
   argList::addBoolOption
   (
       "listTurbulenceModels",
       "List turbulenceModels"
   );
   #endif
   
   return 0;
}


int RHOCENTRAL_CLASS::listOutput()
{
   //bool listOptions = false ;

   if
   (
       argsPtr->optionFound("listSwitches")
   )
   {
       debug::listSwitches(argsPtr->optionFound("includeUnsetSwitches"));
       listOptions = true;
   }

   if
   (
       argsPtr->optionFound("listRegisteredSwitches")
   )
   {
       debug::listRegisteredSwitches(argsPtr->optionFound("includeUnsetSwitches"));
       listOptions = true;
   }

   #ifdef fvPatchField_H
   if (argsPtr->optionFound("listScalarBCs"))
   {
       Info<< "scalarBCs"
           << fvPatchField<scalar>::dictionaryConstructorTablePtr_->sortedToc()
           << endl;
       listOptions = true;
   }

   if (argsPtr->optionFound("listVectorBCs"))
   {
       Info<< "vectorBCs"
           << fvPatchField<vector>::dictionaryConstructorTablePtr_->sortedToc()
           << endl;
       listOptions = true;
   }
   #endif

   #ifdef functionObject_H
   if (argsPtr->optionFound("listFunctionObjects"))
   {
       Info<< "functionObjects"
           << functionObject::dictionaryConstructorTablePtr_->sortedToc()
           << endl;
       listOptions = true;
   }
   #endif

   #ifdef fvOption_H
   if (argsPtr->optionFound("listFvOptions"))
   {
       Info<< "fvOptions"
           << fv::option::dictionaryConstructorTablePtr_->sortedToc()
           << endl;
       listOptions = true;
   }
   #endif

   #ifdef turbulentTransportModel_H
   if (argsPtr->optionFound("listTurbulenceModels"))
   {
       Info<< "Turbulence models"
           << incompressible::turbulenceModel::
              dictionaryConstructorTablePtr_->sortedToc()
           << endl;

       Info<< "RAS models"
           << incompressible::RASModel::
              dictionaryConstructorTablePtr_->sortedToc()
           << endl;

       Info<< "LES models"
           << incompressible::LESModel::
              dictionaryConstructorTablePtr_->sortedToc()
           << endl;
       listOptions = true;
   }
   #elif defined(turbulentFluidThermoModel_H)
   if (argsPtr->optionFound("listTurbulenceModels"))
   {
       Info<< "Turbulence models"
           << compressible::turbulenceModel::
              dictionaryConstructorTablePtr_->sortedToc()
           << endl;

       Info<< "RAS models"
           << compressible::RASModel::
              dictionaryConstructorTablePtr_->sortedToc()
           << endl;

       Info<< "LES models"
           << compressible::LESModel::
              dictionaryConstructorTablePtr_->sortedToc()
           << endl;
       listOptions = true;
   }
   #endif

   if (listOptions)
   {
       exit(0);
   }
   
   return 0;
}





int RHOCENTRAL_CLASS::createTime(int argc,char *argv[])
{

   Foam::Info<< "Create time\n" << Foam::endl;
   //Foam::Time runTimePtr(Foam::Time::controlDictName, *argsPtr);
   runTimePtr = new Foam::Time(Foam::Time::controlDictName, *argsPtr);
   // Mohammad: Not quite sure where this line should be
   
   return 0;
}

/*
// Depreciate this for the moment. Only used in postPtocess.H
int RHOCENTRAL_CLASS::createMesh()
{
   Foam::Info
      << "Create mesh for time = "
      << runTimePtr->timeName() << Foam::nl << Foam::endl;

   Foam::fvMesh meshPtr
   (
      Foam::IOobject
      (
         Foam::fvMesh::defaultRegion,
         runTimePtr->timeName(),
         *runTimePtr,
         Foam::IOobject::MUST_READ
      )
   );

   return 0;
}
*/

int RHOCENTRAL_CLASS::createDynamicFvMesh()
{
    Info<< "Create mesh for time = "
        << runTimePtr->timeName() << nl << endl;

    /*autoPtr<dynamicFvMesh> meshPtr
    (
        dynamicFvMesh::New
        (
            IOobject
            (
                dynamicFvMesh::defaultRegion,
                runTimePtr->timeName(),
                *runTimePtr,
                IOobject::MUST_READ
            )
        )
    );*/

    meshPtr = dynamicFvMesh::New
        (
            IOobject
            (
                dynamicFvMesh::defaultRegion,
                runTimePtr->timeName(),
                *runTimePtr,
                IOobject::MUST_READ
            )
        );
   
   return 0;
}



int RHOCENTRAL_CLASS::createFields()
{
   // createRDeltaT.H
   // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   createRDeltaT();
   // -------------------------------

   Info<< "Reading thermophysical properties\n" << endl;

/*
   autoPtr<psiThermo> thermo
   (
      psiThermo::New(*meshPtr)
   );
*/

    thermoPtr = autoPtr<psiThermo>
    (
        psiThermo::New(*meshPtr)
    );


   //volScalarField& e = thermoPtr->he();
   ePtr = &thermoPtr->he();

   Info<< "Reading field U\n" << endl;

   /*volVectorField U
   (
      IOobject
      (
         "U",
         runTimePtr->timeName(),
         *meshPtr,
         IOobject::MUST_READ,
         IOobject::AUTO_WRITE
      ),
      *meshPtr
   );*/
   
   UPtr = new volVectorField
   (
      IOobject
      (
         "U",
         runTimePtr->timeName(),
         *meshPtr,
         IOobject::MUST_READ,
         IOobject::AUTO_WRITE
      ),
      *meshPtr
   );

   /*volScalarField rho
   (
      IOobject
      (
         "rho",
         runTimePtr->timeName(),
         *meshPtr,
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
      ),
      thermoPtr->rho()
   );*/
   
   rhoPtr = new volScalarField
   (
      IOobject
      (
         "rho",
         runTimePtr->timeName(),
         *meshPtr,
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
      ),
      thermoPtr->rho()
   );

   /*volVectorField rhoU
   (
      IOobject
      (
         "rhoU",
         runTimePtr->timeName(),
         *meshPtr,
         IOobject::NO_READ,
         IOobject::NO_WRITE
      ),
      rho*U
   ); */

   rhoUPtr = new volVectorField
   (
      IOobject
      (
         "rhoU",
         runTimePtr->timeName(),
         *meshPtr,
         IOobject::NO_READ,
         IOobject::NO_WRITE
      ),
      *rhoPtr * *UPtr
   );

   /*volScalarField rhoE
   (
      IOobject
      (
         "rhoE",
         runTimePtr->timeName(),
         *meshPtr,
         IOobject::NO_READ,
         IOobject::NO_WRITE
      ),
      rho*(e + 0.5*magSqr(U))
   ); */

   rhoEPtr = new volScalarField
   (
      IOobject
      (
         "rhoE",
         runTimePtr->timeName(),
         *meshPtr,
         IOobject::NO_READ,
         IOobject::NO_WRITE
      ),
      *rhoPtr * (*ePtr + 0.5*magSqr(*UPtr))
   );

   /*surfaceScalarField pos
   (
      IOobject
      (
         "pos",
         runTimePtr->timeName(),
         *meshPtr
      ),
      *meshPtr,
      dimensionedScalar(dimless, 1.0)
   );*/
   
   posPtr = new surfaceScalarField
   (
      IOobject
      (
         "pos",
         runTimePtr->timeName(),
         *meshPtr
      ),
      *meshPtr,
      dimensionedScalar(dimless, 1.0)
   );

   /*surfaceScalarField neg
   (
      IOobject
      (
         "neg",
         runTimePtr->timeName(),
         *meshPtr
      ),
      *meshPtr,
      dimensionedScalar(dimless, -1.0)
   );*/
   
   negPtr = new surfaceScalarField
   (
      IOobject
      (
         "neg",
         runTimePtr->timeName(),
         *meshPtr
      ),
      *meshPtr,
      dimensionedScalar(dimless, -1.0)
   );

   //surfaceScalarField phi("phi", fvc::flux(rhoU));
   
   phiPtr = new surfaceScalarField("phi", fvc::flux(*rhoUPtr));

   Info<< "Creating turbulence model\n" << endl;

   /*autoPtr<compressible::turbulenceModel> turbulence
   (
      compressible::turbulenceModel::New
      (
         rho,
         U,
         phi,
         thermo
      )
   ); */

/*   turbulencePtr = new compressible::turbulenceModel::New
      (
         *rhoPtr,
         *UPtr,
         *phiPtr,
         *thermoPtr
      );
*/

    turbulencePtr = autoPtr<compressible::turbulenceModel>
    (
        compressible::turbulenceModel::New
        (
            *rhoPtr,
            *UPtr,
            *phiPtr,
            *thermoPtr
        )
    );

   
    return 0;
}


int RHOCENTRAL_CLASS::createControl()
{
   #if defined(NO_CONTROL)
   #elif defined(PISO_CONTROL)
       #include "createPisoControl.H"
   #elif defined(PIMPLE_CONTROL)
       #include "createPimpleControl.H"
   #elif defined(SIMPLE_CONTROL)
       #include "createSimpleControl.H"
   #endif
   
   return 0;
}


int RHOCENTRAL_CLASS::createRDeltaT()
{
   //bool LTS = fv::localEulerDdt::enabled(*meshPtr);
   LTS = fv::localEulerDdt::enabled(*meshPtr);

   if (LTS)
   {
      Info<< "Using LTS" << endl;

      trDeltaT = tmp<volScalarField>
      (
         new volScalarField
         (
            IOobject
            (
               fv::localEulerDdt::rDeltaTName,
               runTimePtr->timeName(),
               *meshPtr,
               IOobject::READ_IF_PRESENT,
               IOobject::AUTO_WRITE
            ),
            *meshPtr,
            dimensionedScalar(dimless/dimTime, 1),
            extrapolatedCalculatedFvPatchScalarField::typeName
         )
      );
   }
   
   return 0;
}


int RHOCENTRAL_CLASS::createFieldRefs()
{
   /*
   volScalarField& p = thermo.p();
   const volScalarField& T = thermo.T();
   const volScalarField& psi = thermo.psi(); */
   const volScalarField& mu = thermoPtr->mu();
   
   pPtr   = &thermoPtr->p();
   TPtr   = &thermoPtr->T();
   psiPtr = &thermoPtr->psi();
   //muPtr  = &thermoPtr->mu();
   
   //bool inviscid(true);
   inviscid = true;
   if (max(mu.primitiveField()) > 0.0)
   //if (max( thermoPtr->mu.primitiveField() ) > 0.0)
   {
       inviscid = false;
   }

   return 0;
}


int RHOCENTRAL_CLASS::createTimeControls()
{
   //bool adjustTimeStep =
   adjustTimeStep = 
      runTimePtr->controlDict().lookupOrDefault("adjustTimeStep", false);

   //scalar maxCo =
   maxCo =
      runTimePtr->controlDict().lookupOrDefault<scalar>("maxCo", 1.0);

   //scalar maxDeltaT =
   maxDeltaT =
      runTimePtr->controlDict().lookupOrDefault<scalar>("maxDeltaT", great);
         
   return 0;
}


int RHOCENTRAL_CLASS::readFluxScheme()
{
   //word fluxScheme("Kurganov");
   word fluxScheme("Kurganov");
   if (meshPtr->schemesDict().readIfPresent("fluxScheme", fluxScheme))
   {
      if ((fluxScheme == "Tadmor") || (fluxScheme == "Kurganov"))
      {
         Info<< "fluxScheme: " << fluxScheme << endl;
      }
      else
      {
         FatalErrorInFunction
         << "fluxScheme: " << fluxScheme
         << " is not a valid choice. "
         << "Options are: Tadmor, Kurganov"
         << abort(FatalError);
      }
   }
   
   return 0;
}

int RHOCENTRAL_CLASS::readTimeControls()
{
   adjustTimeStep =
       runTimePtr->controlDict().lookupOrDefault("adjustTimeStep", false);

   maxCo =
       runTimePtr->controlDict().lookupOrDefault<scalar>("maxCo", 1.0);

   maxDeltaT =
       runTimePtr->controlDict().lookupOrDefault<scalar>("maxDeltaT", great);
   

Foam::Info << "adjustTimeStep = " << adjustTimeStep << endl;
Foam::Info << "maxCo = " << maxCo << endl;
Foam::Info << "maxDeltaT = " << maxDeltaT << endl;


   return 0;
}

int RHOCENTRAL_CLASS::setDeltaT()
{
   if (adjustTimeStep)
   {
       scalar maxDeltaTFact = maxCo/(CoNum + small);
       scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

       runTimePtr->setDeltaT
       (
           min
           (
               deltaTFact * runTimePtr->deltaTValue(),
               maxDeltaT
           )
       );

       Info<< "deltaT = " <<  runTimePtr->deltaTValue() << endl;
   }
   return 0;
}


int RHOCENTRAL_CLASS::setRDeltaT()
{
    volScalarField& rDeltaT = trDeltaT.ref();

    scalar rDeltaTSmoothingCoeff
    (
        runTimePtr->controlDict().lookupOrDefault<scalar>
        (
            "rDeltaTSmoothingCoeff",
            0.02
        )
    );

    // Set the reciprocal time-step from the local Courant number
    rDeltaT.ref() = max
    (
        1/dimensionedScalar(dimTime, maxDeltaT),
        fvc::surfaceSum(*amaxSfPtr)()()
       /((2*maxCo) * meshPtr->V())
    );

    // Update tho boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);

    Info<< "Flow time scale min/max = "
        << gMin(1/rDeltaT.primitiveField())
        << ", " << gMax(1/rDeltaT.primitiveField()) << endl;
        
        
    return 0;
}



int RHOCENTRAL_CLASS::centralCourantNo()
{
   if (meshPtr->nInternalFaces())
   {
       scalarField sumAmaxSf(fvc::surfaceSum(*amaxSfPtr)().primitiveField());

       CoNum = 0.5*gMax(sumAmaxSf/meshPtr->V().field()) * runTimePtr->deltaTValue();

       meanCoNum =
           0.5*(gSum(sumAmaxSf)/gSum(meshPtr->V().field())) * runTimePtr->deltaTValue();
   }

   Info<< "Mean and max Courant Numbers = "
       << meanCoNum << " " << CoNum << endl;

   return 0;
}




