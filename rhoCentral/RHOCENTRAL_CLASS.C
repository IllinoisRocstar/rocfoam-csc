#include "RHOCENTRAL_CLASS.H"

RHOCENTRAL_CLASS::RHOCENTRAL_CLASS() :
                argsPtr(NULL), runTimePtr(NULL), pPtr(NULL), TPtr(NULL), psiPtr(NULL),
                ePtr(NULL), rhoPtr(NULL), UPtr(NULL), rhoUPtr(NULL), rhoEPtr(NULL), 
                posPtr(NULL), negPtr(NULL), phiPtr(NULL), amaxSfPtr(NULL),
                thermoPtr(NULL), meshPtr(NULL), turbulencePtr(NULL), 
                trDeltaT(NULL),
                fluxScheme(""), inviscid(false), adjustTimeStep(false), listOptions(false),
                LTS(false), maxCo(0.0), maxDeltaT(0.0), CoNum(0.0), meanCoNum(0.0)
{}

RHOCENTRAL_CLASS::RHOCENTRAL_CLASS(int argc,char *argv[]) :
                argsPtr(NULL), runTimePtr(NULL), pPtr(NULL), TPtr(NULL), psiPtr(NULL),
                ePtr(NULL), rhoPtr(NULL), UPtr(NULL), rhoUPtr(NULL), rhoEPtr(NULL), 
                posPtr(NULL), negPtr(NULL), phiPtr(NULL), amaxSfPtr(NULL),
                thermoPtr(NULL), meshPtr(NULL), turbulencePtr(NULL), 
                trDeltaT(NULL),
                fluxScheme(""), inviscid(false), adjustTimeStep(false), listOptions(false),
                LTS(false), maxCo(0.0), maxDeltaT(0.0), CoNum(0.0), meanCoNum(0.0)
{
   Initialize(argc, argv);
}


int RHOCENTRAL_CLASS::Initialize(int argc,char *argv[]){
    #define NO_CONTROL

    // Mohammad: Not quite sure where this line should be
    argsPtr = new Foam::argList(argc, argv);

    //  postProcess.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    PostProcess(argc, argv);
    // ---------------------------------------------------

    //  setRootCaseLists.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    setRootCaseLists();
    // ---------------------------------------------------

    //  createTime.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    createTime();
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

    compressible::turbulenceModel &turbulence(*turbulencePtr);

    turbulence.validate();

    Foam::Info << "End of initialization of openFoamPar module." << endl;

    return 0;
}


int RHOCENTRAL_CLASS::loop()
{

    dynamicFvMesh &mesh(*meshPtr);
    Foam::Time &runTime(*runTimePtr);
    volScalarField &p(*pPtr);
    volVectorField &U(*UPtr);
    volVectorField &rhoU(*rhoUPtr);
    const volScalarField &T(*TPtr);
    const volScalarField &psi(*psiPtr);
    volScalarField &e(*ePtr);
    volScalarField &rho(*rhoPtr);
    volScalarField &rhoE(*rhoEPtr);
    surfaceScalarField &pos(*posPtr);
    surfaceScalarField &neg(*negPtr);
    surfaceScalarField &phi(*phiPtr);
    Foam::psiThermo &thermo(*thermoPtr);
    compressible::turbulenceModel &turbulence(*turbulencePtr);



    dimensionedScalar v_zero("v_zero", dimVolume/dimTime, 0.0);

    // Courant numbers used to adjust the time-step
    //scalar CoNum = 0.0;
    //scalar meanCoNum = 0.0;

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        //  readTimeControls.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        readTimeControls();
        // ---------------------------------------------------

        if (!LTS)
        {
            //  setDeltaT.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            setDeltaT();
            // ---------------------------------------------
            runTime++;

            // Do any mesh changes
            mesh.update();
        }

        // --- Directed interpolation of primitive fields onto faces

        surfaceScalarField rho_pos(interpolate(rho, pos));
        surfaceScalarField rho_neg(interpolate(rho, neg));

        surfaceVectorField rhoU_pos(interpolate(rhoU, pos, U.name()));
        surfaceVectorField rhoU_neg(interpolate(rhoU, neg, U.name()));

        volScalarField rPsi("rPsi", 1.0/ psi);
        surfaceScalarField rPsi_pos(interpolate(rPsi, pos, T.name()));
        surfaceScalarField rPsi_neg(interpolate(rPsi, neg, T.name()));

        surfaceScalarField e_pos(interpolate(e, pos, T.name()));
        surfaceScalarField e_neg(interpolate(e, neg, T.name()));

        surfaceVectorField U_pos("U_pos", rhoU_pos/rho_pos);
        surfaceVectorField U_neg("U_neg", rhoU_neg/rho_neg);

        surfaceScalarField p_pos("p_pos", rho_pos*rPsi_pos);
        surfaceScalarField p_neg("p_neg", rho_neg*rPsi_neg);

        surfaceScalarField phiv_pos("phiv_pos", U_pos & mesh.Sf());
        surfaceScalarField phiv_neg("phiv_neg", U_neg & mesh.Sf());

        // Make fluxes relative to mesh-motion
        if (mesh.moving())
        {
            phiv_pos -= mesh.phi();
            phiv_neg -= mesh.phi();
        }

        volScalarField c("c", sqrt(thermo.Cp() / thermo.Cv() * rPsi));
        surfaceScalarField cSf_pos
        (
            "cSf_pos",
            interpolate(c, pos, T.name()) * mesh.magSf()
        );

        surfaceScalarField cSf_neg
        (
            "cSf_neg",
            interpolate(c, neg, T.name()) * mesh.magSf()
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
        surfaceScalarField &amaxSf(*amaxSfPtr);

        amaxSf = max(mag(am), mag(ap));

        surfaceScalarField aSf("aSf", am*a_pos);

        if (fluxScheme == "Tadmor")
        {
            aSf = -0.5 * amaxSf;
            a_pos = 0.5;
        }

        surfaceScalarField a_neg("a_neg", 1.0 - a_pos);

        phiv_pos *= a_pos;
        phiv_neg *= a_neg;

        surfaceScalarField aphiv_pos("aphiv_pos", phiv_pos - aSf);
        surfaceScalarField aphiv_neg("aphiv_neg", phiv_neg + aSf);

        // Reuse amaxSf for the maximum positive and negative fluxes
        // estimated by the central scheme
        amaxSf = max(mag(aphiv_pos), mag(aphiv_neg));

        //  centralCourantNo.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        centralCourantNo();
        // ---------------------------------------------------

        if (LTS)
        {
            // setRDeltaT.H
            // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            setRDeltaT();
            // -------------------------------
            runTime++;
        }

        Info<< "Time = " << runTime.timeName() << nl << endl;

        phi = aphiv_pos*rho_pos + aphiv_neg*rho_neg;

        surfaceVectorField phiUp
        (
            (aphiv_pos*rhoU_pos + aphiv_neg*rhoU_neg)
          + (a_pos*p_pos + a_neg*p_neg) * mesh.Sf()
        );

        surfaceScalarField phiEp
        (
            "phiEp",
            aphiv_pos*(rho_pos*(e_pos + 0.5*magSqr(U_pos)) + p_pos)
          + aphiv_neg*(rho_neg*(e_neg + 0.5*magSqr(U_neg)) + p_neg)
          + aSf*p_pos - aSf*p_neg
        );

        // Make flux for pressure-work absolute
        if (mesh.moving())
        {
            phiEp += mesh.phi()*(a_pos*p_pos + a_neg*p_neg);
        }

        volScalarField muEff("muEff", turbulence.muEff());
        volTensorField tauMC("tauMC", muEff*dev2(Foam::T(fvc::grad(U))));

        // --- Solve density
        solve(fvm::ddt(rho) + fvc::div(phi));

        // --- Solve momentum
        solve(fvm::ddt(rhoU) + fvc::div(phiUp));

        U.ref() = rhoU()/rho();
       
        U.correctBoundaryConditions();
        rhoU.boundaryFieldRef() == rho.boundaryField() * U.boundaryField();

        if (!inviscid)
        {
            solve
            (
                fvm::ddt(rho, U) - fvc::ddt(rho, U)
              - fvm::laplacian(muEff, U)
              - fvc::div(tauMC)
            );
            rhoU = rho * U;
        }

        // --- Solve energy
        surfaceScalarField sigmaDotU
        (
            "sigmaDotU",
            (
                fvc::interpolate(muEff) * mesh.magSf() * fvc::snGrad(U)
              + fvc::dotInterpolate(mesh.Sf(), tauMC)
            )
          & (a_pos*U_pos + a_neg*U_neg)
        );

        solve
        (
            fvm::ddt(rhoE)
          + fvc::div(phiEp)
          - fvc::div(sigmaDotU)
        );

        e = rhoE / rho - 0.5*magSqr(U);
        e.correctBoundaryConditions();
        thermo.correct();
        rhoE.boundaryFieldRef() ==
            rho.boundaryField()*
            (
                e.boundaryField() + 0.5*magSqr(U.boundaryField())
            );

        if (!inviscid)
        {
            solve
            (
                fvm::ddt(rho, e) - fvc::ddt(rho, e)
              - fvm::laplacian(turbulence.alphaEff(), e)
            );
            thermo.correct();
            rhoE = rho * (e + 0.5*magSqr(U));
        }

        p.ref() = rho() / psi();

        p.correctBoundaryConditions();
        rho.boundaryFieldRef() == psi.boundaryField() * p.boundaryField();

        turbulence.correct();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;


    return 0;
}








int RHOCENTRAL_CLASS::PostProcess(int argc,char *argv[])
{

   //  createTime.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#ifndef CREATE_TIME
   #define CREATE_TIME createTime();
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


    Foam::argList &args(*argsPtr);
    Foam::Time &runTime(*runTimePtr);
    

   Foam::argList::addBoolOption
   (
       argList::postProcessOptionName,
       "Execute functionObjects only"
   );

Foam::Info << "FROM POSTPROCESS00" << endl;

   if (argList::postProcess(argc, argv))
   {
      Foam::timeSelector::addOptions();

      //  addRegionOption.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      addRegionOption();
      // ---------------------------------------------------

      //  addFunctionObjectOptions.H  ^^^^^^^^^^^^^^^^^^^^^^
      addFunctionObjectOptions();
      // ---------------------------------------------------

      // Set functionObject post-processing mode
      functionObject::postProcess = true;
      //  setRootCase.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      setRootCase();
      // --------------------------------------------------

      if (args.optionFound("list"))
      {
         functionObjectList::list();
         return 0;
      }

      INCLUDE_FILE(CREATE_TIME)
      Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);

      INCLUDE_FILE(CREATE_MESH)
      dynamicFvMesh &mesh(*meshPtr);

      #ifndef NO_CONTROL
         INCLUDE_FILE(CREATE_CONTROL)
      #endif

      forAll(timeDirs, timei)
      {
         runTime.setTime(timeDirs[timei], timei);

         Info<< "Time = " << runTime.timeName() << endl;

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

            // Externally stored dictionary for functionObjectList
            // if not constructed from runTime
            dictionary functionsControlDict("controlDict");

            HashSet<word> selectedFields;

            // Construct functionObjectList
            autoPtr<functionObjectList> functionsPtr
            (
               functionObjectList::New
               (
                  args,
                  runTime,
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
         mesh.objectRegistry::clear();

         Info<< endl;
      }

      Info<< "End\n" << endl;

      return 0;
   }

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


int RHOCENTRAL_CLASS::setRootCase()
{

    Foam::argList &args(*argsPtr);

   //Foam::argList args(argc, argv);
   if (!args.checkRootCase())
   {
       Foam::FatalError.exit();
   }
   return 0;
}



int RHOCENTRAL_CLASS::setRootCaseLists()
{

   //  listOptions.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   listOptions_Func();
   // --------------------------------------------------
   
   //  setRootCase.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   setRootCase();
   // --------------------------------------------------

   //  listOutput.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   listOutput();
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

    Foam::argList &args(*argsPtr);

   //bool listOptions = false ;

   if
   (
       args.optionFound("listSwitches")
   )
   {
       debug::listSwitches(args.optionFound("includeUnsetSwitches"));
       listOptions = true;
   }

   if
   (
       args.optionFound("listRegisteredSwitches")
   )
   {
       debug::listRegisteredSwitches(args.optionFound("includeUnsetSwitches"));
       listOptions = true;
   }

   #ifdef fvPatchField_H
   if (args.optionFound("listScalarBCs"))
   {
       Info<< "scalarBCs"
           << fvPatchField<scalar>::dictionaryConstructorTablePtr_->sortedToc()
           << endl;
       listOptions = true;
   }

   if (args.optionFound("listVectorBCs"))
   {
       Info<< "vectorBCs"
           << fvPatchField<vector>::dictionaryConstructorTablePtr_->sortedToc()
           << endl;
       listOptions = true;
   }
   #endif

   #ifdef functionObject_H
   if (args.optionFound("listFunctionObjects"))
   {
       Info<< "functionObjects"
           << functionObject::dictionaryConstructorTablePtr_->sortedToc()
           << endl;
       listOptions = true;
   }
   #endif

   #ifdef fvOption_H
   if (args.optionFound("listFvOptions"))
   {
       Info<< "fvOptions"
           << fv::option::dictionaryConstructorTablePtr_->sortedToc()
           << endl;
       listOptions = true;
   }
   #endif

   #ifdef turbulentTransportModel_H
   if (args.optionFound("listTurbulenceModels"))
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
   if (args.optionFound("listTurbulenceModels"))
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



int RHOCENTRAL_CLASS::createTime()
{

    Foam::argList &args(*argsPtr);


   Foam::Info<< "Create time\n" << Foam::endl;
   //Foam::Time runTimePtr(Foam::Time::controlDictName, *argsPtr);
   runTimePtr = new Foam::Time(Foam::Time::controlDictName, args);
   // Mohammad: Not quite sure where this line should be
   
   return 0;
}


int RHOCENTRAL_CLASS::createDynamicFvMesh()
{
    Foam::Time &runTime(*runTimePtr);

    Info<< "Create mesh for time = "
        << runTime.timeName() << nl << endl;

    meshPtr = dynamicFvMesh::New
        (
            IOobject
            (
                dynamicFvMesh::defaultRegion,
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
        );
   
   return 0;
}



int RHOCENTRAL_CLASS::createFields()
{

    Foam::Time &runTime(*runTimePtr);
    dynamicFvMesh &mesh(*meshPtr);


   // createRDeltaT.H
   // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   createRDeltaT();
   // -------------------------------

   Info<< "Reading thermophysical properties\n" << endl;

    thermoPtr = autoPtr<psiThermo>
    (
        psiThermo::New(mesh)
    );
    Foam::psiThermo &thermo(*thermoPtr);


   ePtr = &thermo.he();
   volScalarField &e(*ePtr);

   Info<< "Reading field U\n" << endl;

   UPtr = new volVectorField
   (
      IOobject
      (
         "U",
         runTime.timeName(),
         mesh,
         IOobject::MUST_READ,
         IOobject::AUTO_WRITE
      ),
      mesh
   );
   volVectorField &U(*UPtr);
   
   

   rhoPtr = new volScalarField
   (
      IOobject
      (
         "rho",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
      ),
      thermo.rho()
   );
   volScalarField &rho(*rhoPtr);

   rhoUPtr = new volVectorField
   (
      IOobject
      (
         "rhoU",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
      ),
      rho * U
   );
   volVectorField &rhoU(*rhoUPtr);

   rhoEPtr = new volScalarField
   (
      IOobject
      (
         "rhoE",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
      ),
      rho * (e + 0.5*magSqr(U))
   );

   posPtr = new surfaceScalarField
   (
      IOobject
      (
         "pos",
         runTime.timeName(),
         mesh
      ),
      mesh,
      dimensionedScalar(dimless, 1.0)
   );

   negPtr = new surfaceScalarField
   (
      IOobject
      (
         "neg",
         runTime.timeName(),
         mesh
      ),
      mesh,
      dimensionedScalar(dimless, -1.0)
   );

   phiPtr = new surfaceScalarField("phi", fvc::flux(rhoU));
   surfaceScalarField &phi(*phiPtr);

   Info<< "Creating turbulence model\n" << endl;

    turbulencePtr = autoPtr<compressible::turbulenceModel>
    (
        compressible::turbulenceModel::New
        (
            rho,
            U,
            phi,
            thermo
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
    dynamicFvMesh &mesh(*meshPtr);
    Foam::Time &runTime(*runTimePtr);


   LTS = fv::localEulerDdt::enabled(mesh);

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
               runTime.timeName(),
               mesh,
               IOobject::READ_IF_PRESENT,
               IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(dimless/dimTime, 1),
            extrapolatedCalculatedFvPatchScalarField::typeName
         )
      );
   }
   
   return 0;
}


int RHOCENTRAL_CLASS::createFieldRefs()
{

    Foam::psiThermo &thermo(*thermoPtr);

   /*
   volScalarField& p = thermo.p();
   const volScalarField& T = thermo.T();
   const volScalarField& psi = thermo.psi(); */
   const volScalarField& mu = thermo.mu();
   
   pPtr   = &thermo.p();
   TPtr   = &thermo.T();
   psiPtr = &thermo.psi();
   
   //bool inviscid(true);
   inviscid = true;
   if (max(mu.primitiveField()) > 0.0)
   {
       inviscid = false;
   }

   return 0;
}


int RHOCENTRAL_CLASS::createTimeControls()
{

    Foam::Time &runTime(*runTimePtr);

   //bool adjustTimeStep =
   adjustTimeStep = 
      runTime.controlDict().lookupOrDefault("adjustTimeStep", false);

   //scalar maxCo =
   maxCo =
      runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);

   //scalar maxDeltaT =
   maxDeltaT =
      runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", great);
         
   return 0;
}


int RHOCENTRAL_CLASS::readFluxScheme()
{

    dynamicFvMesh &mesh(*meshPtr);

   //word fluxScheme("Kurganov");
   word fluxScheme("Kurganov");
   if (mesh.schemesDict().readIfPresent("fluxScheme", fluxScheme))
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
    Foam::Time &runTime(*runTimePtr);

   adjustTimeStep =
       runTime.controlDict().lookupOrDefault("adjustTimeStep", false);

   maxCo =
       runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);

   maxDeltaT =
       runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", great);
   

Foam::Info << "adjustTimeStep = " << adjustTimeStep << endl;
Foam::Info << "maxCo = " << maxCo << endl;
Foam::Info << "maxDeltaT = " << maxDeltaT << endl;


   return 0;
}

int RHOCENTRAL_CLASS::setDeltaT()
{
    Foam::Time &runTime(*runTimePtr);

   if (adjustTimeStep)
   {
       scalar maxDeltaTFact = maxCo/(CoNum + small);
       scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

       runTime.setDeltaT
       (
           min
           (
               deltaTFact * runTime.deltaTValue(),
               maxDeltaT
           )
       );

       Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
   }
   return 0;
}


int RHOCENTRAL_CLASS::setRDeltaT()
{
    Foam::Time &runTime(*runTimePtr);
    dynamicFvMesh &mesh(*meshPtr);
    surfaceScalarField &amaxSf(*amaxSfPtr);
    
    

    volScalarField& rDeltaT = trDeltaT.ref();

    scalar rDeltaTSmoothingCoeff
    (
        runTime.controlDict().lookupOrDefault<scalar>
        (
            "rDeltaTSmoothingCoeff",
            0.02
        )
    );

    // Set the reciprocal time-step from the local Courant number
    rDeltaT.ref() = max
    (
        1/dimensionedScalar(dimTime, maxDeltaT),
        fvc::surfaceSum(amaxSf)()()
       /((2*maxCo) * mesh.V())
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
    Foam::Time &runTime(*runTimePtr);
    dynamicFvMesh &mesh(*meshPtr);
    surfaceScalarField &amaxSf(*amaxSfPtr);
    

   if (mesh.nInternalFaces())
   {
       scalarField sumAmaxSf(fvc::surfaceSum(amaxSf)().primitiveField());

       CoNum = 0.5*gMax(sumAmaxSf/mesh.V().field()) * runTime.deltaTValue();

       meanCoNum =
           0.5*(gSum(sumAmaxSf)/gSum(mesh.V().field())) * runTime.deltaTValue();
   }

   Info<< "Mean and max Courant Numbers = "
       << meanCoNum << " " << CoNum << endl;

   return 0;
}


int RHOCENTRAL_CLASS::finalize()
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
    
    return 0;
}

