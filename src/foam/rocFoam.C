#include "rocFoam.H"
#include <vector>

rocFoam::rocFoam()
{ initSet(); }

int rocFoam::initSet()
{
    solverType = NULL;
    argsPtr = NULL;
    runTimePtr = NULL;
    listOptions = false;
    LTS = false;
    adjustTimeStep = false;
    overrideTimeStep = false;
    maxCo = 0.0;
    maxDeltaT = 0.0;
    overrideDeltaT = small;
    CoNum = 0.0;
    meanCoNum = 0.0;
    pPtr = NULL;
    TPtr = NULL;
    psiPtr = NULL;
    ePtr = NULL;
    rhoPtr = NULL;
    UPtr = NULL;
    rhoUPtr = NULL;
    rhoEPtr = NULL;
    phiPtr = NULL;
    meshPtr = NULL;
    turbulencePtr = NULL;
    trDeltaT = NULL;
    initializeStat = -1;
    loopStat = -1;
    stepStat = -1;
    finalizeStat = -1;
    testStat = -1.0;
    
    return 0;
}


rocFoam::~rocFoam()
{
    finalize();
}

int rocFoam::finalize()
{
    // Delete thing that are allocated here
    if (runTimePtr != NULL)
    {
      delete runTimePtr;
      runTimePtr = NULL;
    }

    if (argsPtr != NULL)
    {
         // Once this is deleted, all information
         // about the openfoam-related stuff is gone.
         // Currently have no better place to delete
         // this. One option is not to delete this.

         delete argsPtr;
         argsPtr = NULL;
    }

    return 0;
}

int rocFoam::PostProcess(int argc, char *argv[])
{
    //  createTime.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^
#ifndef CREATE_TIME
#define CREATE_TIME createTime();
#endif
    // ------------------------------------------

    //  createMesh.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^
#ifndef CREATE_MESH
//#define CREATE_MESH createMesh();
#define CREATE_MESH createDynamicFvMesh();
#endif
    // ------------------------------------------

    //  createFields.H  ^^^^^^^^^^^^^^^^^^^^^^^^^
#ifndef CREATE_FIELDS
#define CREATE_FIELDS createFields();
#endif
    // ------------------------------------------

    //  createControl.H  ^^^^^^^^^^^^^^^^^^^^^^^^
#ifndef CREATE_CONTROL
#define CREATE_CONTROL createControl();
#endif
    // ------------------------------------------

#define INCLUDE_FILE(X) INCLUDE_FILE2(X)
#define INCLUDE_FILE2(X) X

    Foam::argList &args(*argsPtr);
    Foam::Time &runTime(*runTimePtr);

    Foam::argList::addBoolOption(argList::postProcessOptionName,
                                 "Execute functionObjects only");

    if (argList::postProcess(argc, argv))
    {
        Foam::timeSelector::addOptions();

        //  addRegionOption.H  ^^^^^^^^^^^^^^^^^^
        addRegionOption();
        // --------------------------------------

        //  addFunctionObjectOptions.H  ^^^^^^^^^
        addFunctionObjectOptions();
        // --------------------------------------

        // Set functionObject post-processing mode
        functionObject::postProcess = true;
        //  setRootCase.H  ^^^^^^^^^^^^^^^^^^^^^^
        setRootCase();
        // --------------------------------------

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

            Info << "Time = " << runTime.timeName() << endl;

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
            catch (IOerror &err)
            {
                Warning << err << endl;
            }

            // Clear the objects owned by the mesh
            mesh.objectRegistry::clear();

            Info << endl;
        }

        Info << "End\n" << endl;

        return 0;
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    // //

#undef INCLUDE_FILE
#undef INCLUDE_FILE2

#undef CREATE_MESH
#undef CREATE_FIELDS
#undef CREATE_CONTROL

    return 0;
}

int rocFoam::addRegionOption()
{
    Foam::argList::addOption
    (
        "region",
        "name",
        "specify alternative mesh region"
    );

    return 0;
}

int rocFoam::addDictOption()
{
    Foam::argList::addOption
    (
        "dict",
        "file",
        "read control dictionary from specified location"
    );
    
    return 0;
}

int rocFoam::addFunctionObjectOptions()
{
    //  addDictOption.H  ^^^^^^^^^^^^^^^^^^^^^^^^
    addDictOption();
    // ------------------------------------------

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
        "list", "List the available configured functionObjects"
    );

    return 0;
}

int rocFoam::createArgs(int argc, char *argv[])
{
    argsPtr = new Foam::argList(argc, argv);
    return 0;
}


int rocFoam::setRootCase()
{
    Foam::argList &args(*argsPtr);

    // Foam::argList args(argc, argv);
    if (!args.checkRootCase())
    {
        Foam::FatalError.exit();
    }
    return 0;
}

int rocFoam::setRootCaseLists()
{

    //  listOptions.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^
    listOptionsFunc();
    // ------------------------------------------

    //  setRootCase.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^
    setRootCase();
    // ------------------------------------------

    //  listOutput.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^
    listOutput();
    // ------------------------------------------

    return 0;
}

int rocFoam::listOptionsFunc()
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

int rocFoam::listOutput()
{
    Foam::argList &args(*argsPtr);

    // bool listOptions = false ;
    
    if (args.optionFound("listSwitches"))
    {
        debug::listSwitches(args.optionFound("includeUnsetSwitches"));
        listOptions = true;
    }

    if (args.optionFound("listRegisteredSwitches"))
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
        Info << "functionObjects"
             << functionObject::dictionaryConstructorTablePtr_->sortedToc()
             << endl;
        listOptions = true;
    }
#endif




#ifdef fvOption_H
    if (args.optionFound("listFvOptions"))
    {
        Info << "fvOptions"
             << fv::option::dictionaryConstructorTablePtr_->sortedToc() << endl;
        listOptions = true;
    }
#endif

#ifdef turbulentTransportModel_H
    if (args.optionFound("listTurbulenceModels"))
    {
        Info << "Turbulence models"
             << incompressible::turbulenceModel::dictionaryConstructorTablePtr_
                    ->sortedToc()
             << endl;

        Info << "RAS models"
             << incompressible::RASModel::dictionaryConstructorTablePtr_
                    ->sortedToc()
             << endl;

        Info << "LES models"
             << incompressible::LESModel::dictionaryConstructorTablePtr_
                    ->sortedToc()
             << endl;
        listOptions = true;
    }

#elif defined(turbulentFluidThermoModel_H)
    if (args.optionFound("listTurbulenceModels"))
    {
        Info << "Turbulence models"
             << compressible::turbulenceModel::dictionaryConstructorTablePtr_
                    ->sortedToc()
             << endl;

        Info << "RAS models"
             << compressible::RASModel::dictionaryConstructorTablePtr_
                    ->sortedToc()
             << endl;

        Info << "LES models"
             << compressible::LESModel::dictionaryConstructorTablePtr_
                    ->sortedToc()
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

int rocFoam::createTime()
{
    Foam::argList &args(*argsPtr);

    Foam::Info << "Create time\n" << Foam::endl;
    // Foam::Time runTimePtr(Foam::Time::controlDictName, *argsPtr);
    runTimePtr = new Foam::Time(Foam::Time::controlDictName, args);
    
    return 0;
}

// ^^^ The two folloing functions are the same ^^^
int rocFoam::createTimeControls()
{
    Foam::Time &runTime(*runTimePtr);

    adjustTimeStep =
        runTime.controlDict().lookupOrDefault("adjustTimeStep", false);

    maxCo = runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);

    maxDeltaT =
        runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", great);

    return 0;
}

int rocFoam::readTimeControls()
{
    Foam::Time &runTime(*runTimePtr);

    adjustTimeStep =
        runTime.controlDict().lookupOrDefault("adjustTimeStep", false);

    maxCo = runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);

    maxDeltaT =
        runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", great);

    return 0;
}
//-----------------------------------------------

int rocFoam::setDeltaT()
{
    Foam::Time &runTime(*runTimePtr);

    if (adjustTimeStep && !overrideTimeStep)
    {
        scalar maxDeltaTFact = maxCo / (CoNum + small);
        scalar deltaTFact = min
        (
            min(maxDeltaTFact, 1.0 + 0.1 * maxDeltaTFact),
            1.2
        );

        runTime.setDeltaT(min(deltaTFact * runTime.deltaTValue(), maxDeltaT));

        Info << "deltaT = " << runTime.deltaTValue() << endl;
    }
    /* else if (overrideTimeStep && !adjustTimeStep)
    {
        FatalErrorInFunction
            << "AdjustTimeStep should be set in control dictionary"
            << Foam::endl;    
    } */
    else if (overrideTimeStep)
    {
        scalar maxDeltaTFact = maxCo / (CoNum + small);
        scalar deltaTFact = min
        (
            min(maxDeltaTFact, 1.0 + 0.1 * maxDeltaTFact),
            1.2
        );

        scalar mainDeltaT = min(deltaTFact * runTime.deltaTValue(), maxDeltaT);
        if (overrideDeltaT > mainDeltaT)
        {
            Info << "Warning: deltaT = " << mainDeltaT 
                 << ", overrideDeltaT = " << overrideDeltaT << endl;
        }

        runTime.setDeltaT(overrideDeltaT);

        Info << "deltaT = " << runTime.deltaTValue() << endl;
    }
    
    return 0;
}

/*
int rocFoam::setDeltaT(scalar &overrideDeltaT)
{
    Foam::Time &runTime(*runTimePtr);

    if (!adjustTimeStep)
    {
        FatalErrorInFunction
            << "AdjustTimeStep should be set in control dictionary"
            << Foam::endl;
    }
    else
    {
        scalar maxDeltaTFact = maxCo / (CoNum + small);
        scalar deltaTFact = min
        (
            min(maxDeltaTFact, 1.0 + 0.1 * maxDeltaTFact),
            1.2
        );

        scalar mainDeltaT = min(deltaTFact * runTime.deltaTValue(), maxDeltaT);
        if (overrideDeltaT > mainDeltaT)
        {
            Info << "Warning: deltaT = " << mainDeltaT 
                 << ", overrideDeltaT = " << overrideDeltaT << endl;
        }

        runTime.setDeltaT(overrideDeltaT);

        Info << "deltaT = " << runTime.deltaTValue() << endl;
    }
    return 0;
}
*/

int rocFoam::createRDeltaT()
{
    dynamicFvMesh &mesh(*meshPtr);
    Foam::Time &runTime(*runTimePtr);

    LTS = fv::localEulerDdt::enabled(mesh);
    if (LTS)
    {
        Info << "Using LTS" << endl;

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
                dimensionedScalar(dimless / dimTime, 1),
                extrapolatedCalculatedFvPatchScalarField::typeName
            )
        );
    }

    return 0;
}

int rocFoam::createDynamicFvMesh()
{
    Foam::Time &runTime(*runTimePtr);

    Info << "Create mesh for time = " << runTime.timeName() << nl << endl;

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


