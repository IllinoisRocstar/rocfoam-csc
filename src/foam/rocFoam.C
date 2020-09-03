#include "rocFoam.H"
#include <vector>

rocFoam::rocFoam()
{}

rocFoam::~rocFoam()
{
    finalizeFoam();
}

int rocFoam::finalizeFoam()
{
    // Delete thing that are allocated here
    if (runTimePtr != nullptr)
    {
      delete runTimePtr;
      runTimePtr = nullptr;
    }

    if (argsPtr != nullptr)
    {
         // Once this is deleted, all information
         // about the openfoam-related stuff is gone.
         // Currently have no better place to delete
         // this. One option is not to delete this.

         delete argsPtr;
         argsPtr = nullptr;
    }

    return 0;
}

int rocFoam::PostProcess(int argc, char *argv[])
#ifdef HAVE_OFE20
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

        if (args.found("list"))
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

        // Externally stored dictionary for functionObjectList
        // if not constructed from runTime
        dictionary functionsDict;

        HashSet<wordRe> selectedFields;

        // Construct functionObjectList
        autoPtr<functionObjectList> functionsPtr
        (
            functionObjectList::New(args, runTime, functionsDict, selectedFields)
        );

        forAll(timeDirs, timei)
        {
            runTime.setTime(timeDirs[timei], timei);

            Info << "Time = " << runTime.timeName() << endl;

            switch (mesh.readUpdate())
            {
                case polyMesh::POINTS_MOVED:
                {
                    functionsPtr->movePoints(mesh);
                    break;
                }
                case polyMesh::TOPO_CHANGE:
                case polyMesh::TOPO_PATCH_CHANGE:
                {
                    mapPolyMesh mpm(mesh);
                    functionsPtr->updateMesh(mpm);
                    break;
                }
                case polyMesh::UNCHANGED:
                {
                    // No additional work
                    break;
                }
                default:
                {
                    FatalErrorIn(args.executable())
                        << "Unhandled enumeration"
                        << abort(FatalError);
                }
            }

            FatalIOError.throwExceptions();

            try
            {
                INCLUDE_FILE(CREATE_FIELDS)

                #ifdef CREATE_FIELDS_2
                INCLUDE_FILE(CREATE_FIELDS_2)
                #endif

                #ifdef CREATE_FIELDS_3
                INCLUDE_FILE(CREATE_FIELDS_3)
                #endif

                functionsPtr->execute();

                // Execute the functionObject 'end()' function for the last time
                if (timei == timeDirs.size()-1)
                {
                    functionsPtr->end();
                }

                // Report to output (avoid overwriting values from simulation)
                profiling::print(Info);
            }
            catch (const IOerror& err)
            {
                Warning<< err << endl;
            }

            Info<< endl;
        }

        Info << "End\n" << endl;

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
#elif defined(HAVE_OF7) || defined(HAVE_OF8)
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
                INCLUDE_FILE(CREATE_FIELDS_2)
                #endif

                #ifdef CREATE_FIELDS_3
                INCLUDE_FILE(CREATE_FIELDS_3)
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

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #undef INCLUDE_FILE
    #undef INCLUDE_FILE2

    #undef CREATE_MESH
    #undef CREATE_FIELDS
    #undef CREATE_CONTROL

    return 0;
}
#endif

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
        "Read control dictionary from specified location"
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
    
#ifdef HAVE_OFE20
    Foam::argList::addOption
    (
        "fields",
        "list",
        "Specify a list of fields to be processed, e.g. '(U T p)'"
    );
#elif defined(HAVE_OF7) || defined(HAVE_OF8)
    Foam::argList::addOption
    (
        "fields",
        "list",
        "Specify a list of fields to be processed, e.g. '(U T p)' - "
        "regular expressions not currently supported"
    );
#endif
    
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

#ifdef HAVE_OFE20
    #include "foamDlOpenLibs.H"
#endif

    return 0;
}

int rocFoam::setRootCaseLists()
{
#ifdef HAVE_OFE20
    //  setRootCaseListOptions.H
    setRootCaseListOptions();
    //--------------------------
    
    //  setRootCase.H
    setRootCase();
    //---------------

    //  setRootCaseListOutput.H
    setRootCaseListOutput();
    //-------------------------
#elif defined(HAVE_OF7) || defined(HAVE_OF8)
    //  listOptions.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^
    listOptions_();
    // ------------------------------------------

    //  setRootCase.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^
    setRootCase();
    // ------------------------------------------

    //  listOutput.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^
    listOutput();
    // ------------------------------------------
#endif
    return 0;
}

#ifdef HAVE_OFE20
int rocFoam::setRootCaseListOptions()
{
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // Declare some "standard" list options
    #ifndef setRootCaseListOptions_H
    #define setRootCaseListOptions_H

    argList::addBoolOption
    (
        "listSwitches",
        "List switches declared in libraries"
        " (see -listUnsetSwitches option)",
        true // advanced
    );
    argList::addBoolOption
    (
        "listRegisteredSwitches",
        "List switches registered for run-time modification"
        " (see -listUnsetSwitches option)",
        true // advanced
    );
    argList::addBoolOption
    (
        "listUnsetSwitches",
        "Modifies switch listing to display values not set in etc/controlDict",
        true // advanced
    );

    #ifdef fvPatchField_H
    argList::addBoolOption
    (
        "listScalarBCs",
        "List scalar field boundary conditions (fvPatchField<scalar>)",
        true // advanced
    );
    argList::addBoolOption
    (
        "listVectorBCs",
        "List vector field boundary conditions (fvPatchField<vector>)",
        true // advanced
    );
    #endif

    #ifdef functionObject_H
    argList::addBoolOption
    (
        "listFunctionObjects",
        "List functionObjects",
        true // advanced
    );
    #endif

    #ifdef fvOption_H
    argList::addBoolOption
    (
        "listFvOptions",
        "List fvOptions",
        true // advanced
    );
    #endif

    #if defined(turbulentTransportModel_H) || defined(turbulentFluidThermoModel_H)
    argList::addBoolOption
    (
        "listTurbulenceModels",
        "List turbulenceModels",
        true // advanced
    );
    #endif


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    #endif
    // ************************************************************************* //
    return 0;
}

int rocFoam::setRootCaseListOutput()
{
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // Process some "standard" list options
    #ifndef setRootCaseListOutput_H
    #define setRootCaseListOutput_H

    {
        Foam::argList &args(*argsPtr);

        //bool listOptions = false;

        if (args.found("listSwitches"))
        {
            debug::listSwitches(args.found("listUnsetSwitches"));
            listOptions = true;
        }

        if (args.found("listRegisteredSwitches"))
        {
            debug::listRegisteredSwitches(args.found("listUnsetSwitches"));
            listOptions = true;
        }

        #ifdef fvPatchField_H
        if (args.found("listScalarBCs"))
        {
            Info<< "scalarBCs"
                << fvPatchField<Foam::scalar>::
                dictionaryConstructorTablePtr_->sortedToc()
                << endl;
            listOptions = true;
        }

        if (args.found("listVectorBCs"))
        {
            Info<< "vectorBCs"
                <<  fvPatchField<Foam::vector>::
                    dictionaryConstructorTablePtr_->sortedToc()
                << endl;
            listOptions = true;
        }
        #endif

        #ifdef functionObject_H
        if (args.found("listFunctionObjects"))
        {
            Info<< "functionObjects"
                << functionObject::dictionaryConstructorTablePtr_->sortedToc()
                << endl;
            listOptions = true;
        }
        #endif

        #ifdef fvOption_H
        if (args.found("listFvOptions"))
        {
            Info<< "fvOptions"
                << fv::option::dictionaryConstructorTablePtr_->sortedToc()
                << endl;
            listOptions = true;
        }
        #endif

        #ifdef turbulentTransportModel_H
        if (args.found("listTurbulenceModels"))
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
        if (args.found("listTurbulenceModels"))
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
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    #endif
    // ************************************************************************* //
   
    return 0;
}

#elif defined(HAVE_OF7) || defined(HAVE_OF8)
int rocFoam::listOptions_()
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

#ifdef HAVE_OF7
    #ifdef turbulentTransportModel_H
    if (args.optionFound("listTurbulenceModels"))
    {
        Info << "Turbulence models"
             << incompressible::turbulenceModel::
                    dictionaryConstructorTablePtr_->sortedToc()
             << endl;
#elif defined(HAVE_OF8)
    #ifdef kinematicMomentumTransportModel_H
    if (args.optionFound("listMomentumTransportModels"))
    {
        Info<< "Turbulence models"
            << incompressible::momentumTransportModel::
                    dictionaryConstructorTablePtr_->sortedToc()
            << endl;
#endif

        Info << "RAS models"
             << incompressible::RASModel::
                    dictionaryConstructorTablePtr_->sortedToc()
             << endl;

        Info << "LES models"
             << incompressible::LESModel::
                    dictionaryConstructorTablePtr_->sortedToc()
             << endl;
        listOptions = true;
    }

#ifdef HAVE_OF7
    #elif defined(turbulentFluidThermoModel_H)
    if (args.optionFound("listTurbulenceModels"))
    {
        Info << "Turbulence models"
             << compressible::turbulenceModel::
                    dictionaryConstructorTablePtr_->sortedToc()
             << endl;
#elif defined(HAVE_OF8)
    #elif defined(fluidThermoMomentumTransportModel_H)
    if (args.optionFound("listMomentumTransportModels"))
    {
        Info<< "Turbulence models"
            << compressible::momentumTransportModel::
                    dictionaryConstructorTablePtr_->sortedToc()
            << endl;
#endif

        Info << "RAS models"
             << compressible::RASModel::
                    dictionaryConstructorTablePtr_->sortedToc()
             << endl;

        Info << "LES models"
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
#endif


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

#ifdef HAVE_OFE20
    adjustTimeStep =
        runTime.controlDict().getOrDefault("adjustTimeStep", false);

    maxCo = runTime.controlDict().getOrDefault<scalar>("maxCo", 1);

    maxDeltaT =
        runTime.controlDict().getOrDefault<scalar>("maxDeltaT", GREAT);
#elif defined(HAVE_OF7) || defined(HAVE_OF8)

    adjustTimeStep =
        runTime.controlDict().lookupOrDefault("adjustTimeStep", false);

    maxCo = runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);

    maxDeltaT =
        runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", great);
#endif

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

    if (adjustTimeStep)
    {
#ifdef HAVE_OFE20
        scalar maxDeltaTFact = maxCo/(CoNum + SMALL);
#elif defined(HAVE_OF7) defined(HAVE_OF8)
        scalar maxDeltaTFact = maxCo/(CoNum + small);
#endif

        scalar deltaTFact = min
        (
            min(maxDeltaTFact, 1.0 + 0.1 * maxDeltaTFact),
            1.2
        );

        double mainDeltaT = min(deltaTFact * runTime.deltaTValue(), maxDeltaT);
        
        runTime.setDeltaT(mainDeltaT);

        Info << "deltaT = " << runTime.deltaTValue() << endl;
    }
    return 0;
}

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
#ifdef HAVE_OFE20
                dimensionedScalar("one", dimless/dimTime, 1),
#elif defined(HAVE_OF7) || defined(HAVE_OF8)
                dimensionedScalar(dimless / dimTime, 1),
#endif
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

#ifdef HAVE_OFE20
    meshPtr = dynamicFvMesh::New((args, runTime));
#elif defined(HAVE_OF7) || defined(HAVE_OF8)
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
#endif

    return 0;
}


