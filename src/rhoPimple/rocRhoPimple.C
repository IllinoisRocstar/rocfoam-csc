#include "rocRhoPimple.H"

//^^^ DEFINITION OF CONSTRUCTORS ^^^^^^^^^^^^^^^^^^^^^^^^^^
rhoPimple::rhoPimple()
{
    solverType = "rocRhoPimple";
}

rhoPimple::rhoPimple(int argc, char *argv[])
{
    solverType = "rocRhoPimple";
    initFOAM(argc, argv);
}
//=========================================================


//^^^ DEFINITION OF COM-RELATED MTHODS ^^^^^^^^^^^^^^^^^^^^
//^^^^^ LOAD MODULES ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
int rhoPimple::loadInternal(const char* name)
{
    load(name);
    return 0;
}

void rhoPimple::load(const char *name)
{
    //  Anouncing default communicator  ^^^^^^^^^
    MPI_Comm tmpComm;
    tmpComm = COM_get_default_communicator();  

    int tmpRank, tmpNProc;
    MPI_Comm_rank(tmpComm, &tmpRank);
    MPI_Comm_size(tmpComm, &tmpNProc);
    
    if (tmpRank == 0)
    {
        std::cout << "rocRhoPimple: Loading rocRhoPimple with name "
                   << name << "." << std::endl;

        std::cout << "rocFoam.load: Rank = " << tmpRank
                  << ", NProc = " << tmpNProc
                  << ", COMM = " << tmpComm << std::endl;

        std::cout << std::endl;
    }

    std::string winName = name;
    int winExist = COM_get_window_handle(winName.c_str());
    if (winExist>0)
    {
        std::cout << "WARNING: Window " << winName << " already exists."
                  << " CSC must create this window name."
                  << std::endl;
        exit(-1);
    }
    else
    {
        COM_new_window(winName, tmpComm);

        Info << "rocFoam.load: Window " << winName
             << " created." << endl;
    }

    // Register object ^^^^^^^^^^^^^^^^^^^^^^^^^^
    rhoPimple *comFoamPtr = new rhoPimple();
    std::string objectName = winName + string(".object");
    COM_new_dataitem(objectName.c_str(), 'w', COM_VOID, 1, "");
    COM_set_object(objectName.c_str(), 0, comFoamPtr);
    COM_window_init_done(winName);
    
    //MPI_Comm_dup(tmpComm, &(comFoamPtr->winComm));
    comFoamPtr->winComm = tmpComm;
    //Foam::PstreamGlobals::MPI_COMM_FOAM = comFoamPtr->winComm;
    MPI_Comm_rank(comFoamPtr->winComm, &(comFoamPtr->ca_myRank));
    MPI_Comm_size(comFoamPtr->winComm, &(comFoamPtr->ca_nProc));
    comFoamPtr->winName = winName;
    comFoamPtr->registerFunctions(winName.c_str());
    //-------------------------------------------

    // Vol window ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    winName = name+string("VOL");
    winExist = COM_get_window_handle(winName.c_str());
    if (winExist>0)
    {
        std::cout << "Window " << winName << " already exists."
                  << " Assure that there is nothing wrong with it."
                  << std::endl;
    }
    else
    {
        COM_new_window(winName, tmpComm);
        COM_window_init_done(winName);

        Info << "rocFoam.load: Window " << winName
             << " created." << endl;
    }
    //-------------------------------------------

    // Surf window ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    winName = name+string("SURF");
    winExist = COM_get_window_handle(winName.c_str());
    if (winExist>0)
    {
        std::cout << "Window " << winName << " already exists."
                  << " Assure that there is nothing wrong with it."
                  << std::endl;
    }
    else
    {
        COM_new_window(winName, tmpComm);
        COM_window_init_done(winName);

        Info << "rocFoam.load: Window " << winName
             << " created." << endl;
    }
    //-------------------------------------------

    return;
}
//-----------------------------------------------

//^^^^^ UNLOAD MODULES ^^^^^^^^^^^^^^^^^^^^^^^^^^
void rhoPimple::unload(const char *name)
{
    Foam::Info << "rocFoam.unload: Unloading rocRhoPimple with name "
               << name << "." << Foam::endl;

    std::string winName = name+std::string("VOL");
    int winExist = COM_get_window_handle(winName.c_str());
    if (winExist>0)
        COM_delete_window(winName);

    winName = name+std::string("SURF");
    winExist = COM_get_window_handle(winName.c_str());
    if (winExist>0)
        COM_delete_window(winName);

    winName = std::string(name);
    winExist = COM_get_window_handle(winName.c_str());
    if (winExist>0)
    {
        comFoam *comFoamPtr = nullptr;
        std::string objectName(winName+".object");
        COM_get_object(objectName.c_str(), 0, &comFoamPtr);

        //comFoamPtr->finalize();
        if (comFoamPtr != nullptr)
            delete comFoamPtr;
        
        COM_delete_window(winName);
    }
}
//-----------------------------------------------
//=========================================================

//^^^ Solver-specific methods ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
int rhoPimple::initFOAM(int argc, char *argv[])
{
    createArgs(argc, argv);

#ifdef HAVE_OFE20
    argList::addNote
    (
        "Transient solver for compressible turbulent flow.\n"
        "With optional mesh motion and mesh topology changes."
    );
#endif

    //  postProcess.H
    PostProcess(argc, argv);
    // --------------

#ifdef HAVE_OFE20
    //  addCheckCaseOptions.H
    addCheckCaseOptions();
    // ----------------------
#endif

    //  setRootCaseLists.H
    setRootCaseLists();
    // -------------------

    //  createTime.H
    createTime();
    // -------------

    //  createDynamicFvMesh.H
    createDynamicFvMesh();
    // ----------------------

    //  createDyMControls.H
    createDyMControls();
    // --------------------

    //  initContinuityErrs.H
    initContinuityErrs();
    // ---------------------

    //  createFields.H
    createFields();
    // ---------------

    //  createFieldRefs.H
    createFieldRefs();
    // ------------------

    //  createRhoUfIfPresent.H
    createRhoUfIfPresent();
    // -----------------------

#if defined(HAVE_OFE20)
    compressible::turbulenceModel& turbulence(*turbulencePtr);
#elif defined(HAVE_OF7)
    compressible::turbulenceModel& turbulence(*turbulencePtr);
#elif defined(HAVE_OF8)
    compressible::momentumTransportModel& turbulence(*turbulencePtr);
#endif

    turbulence.validate();

    if (!LTS)
    {
        //  compressibleCourantNo.H
        compressibleCourantNo();
        // ------------------------

        //  setInitialDeltaT.H
        setInitialDeltaT();
        // -------------------
    }

    Foam::Info << "End of initialization of rhoPimple." << Foam::endl;

    initializeStat = 0;
    return initializeStat;
}

int rhoPimple::createControl()
{
    dynamicFvMesh &mesh(*meshPtr);

    //  createPimpleControl.H
    // pimpleControl pimple(mesh);
    pimplePtr = new pimpleControl(mesh);
    // ----------------------

    return 0;
}

int rhoPimple::createDyMControls()
{
    dynamicFvMesh &mesh(*meshPtr);

    //  createControl.H
    createControl();
    // ----------------
    pimpleControl &pimple(*pimplePtr);

    //  createTimeControls.H
    createTimeControls();
    // ---------------------

#ifdef HAVE_OFE20
    correctPhi = pimple.dict().getOrDefault
    (
        "correctPhi",
        mesh.dynamic()
    );

    checkMeshCourantNo = pimple.dict().getOrDefault
    (
        "checkMeshCourantNo",
        false
    );

    moveMeshOuterCorrectors = pimple.dict().getOrDefault
    (
        "moveMeshOuterCorrectors",
        false
    );
#elif defined(HAVE_OF7) || defined(HAVE_OF8)
    correctPhi = pimple.dict().lookupOrDefault
    (
        "correctPhi",
        mesh.dynamic()
    );

    checkMeshCourantNo = pimple.dict().lookupOrDefault
    (
        "checkMeshCourantNo",
        false
    );

    moveMeshOuterCorrectors = pimple.dict().lookupOrDefault
    (
        "moveMeshOuterCorrectors",
        false
    );
#endif

    return 0;
}

int rhoPimple::initContinuityErrs()
{
    #ifndef initContinuityErrs_H
    #define initContinuityErrs_H

#ifdef HAVE_OFE20
    Foam::Time &runTime(*runTimePtr);
    dynamicFvMesh &mesh(*meshPtr);

    uniformDimensionedScalarField cumulativeContErrIO
    (
        IOobject
        (
            "cumulativeContErr",
            runTime.timeName(),
            "uniform",
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        dimensionedScalar(dimless, Zero)
    );
    cumulativeContErr = cumulativeContErrIO.value();

#elif defined(HAVE_OF7) || defined(HAVE_OF8)

    cumulativeContErr = 0;

#endif

    #endif

    return 0;
}

int rhoPimple::createFields()
{
    Foam::Time &runTime(*runTimePtr);
    Foam::argList &args(*argsPtr);
    dynamicFvMesh &mesh(*meshPtr);
    pimpleControl &pimple(*pimplePtr);

    //  createRDeltaT.H  ^^^^^^^^^^^^^
    createRDeltaT();
    // -------------------------------

    Info << "Reading thermophysical properties\n" << endl;

    pThermoPtr = autoPtr<fluidThermo>(fluidThermo::New(mesh));
    fluidThermo &thermo(*pThermoPtr);

    thermo.validate(args.executable(), "h", "e");
    pPtr = &thermo.p();

    volScalarField &p(*pPtr);
    TPtr = &thermo.T();

    rhoPtr = new volScalarField
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        thermo.rho()
    );
    volScalarField &rho(*rhoPtr);

    Info << "Reading field U\n" << endl;

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

    //  compressibleCreatePhi.H  ^^^^^
    compressibleCreatePhi();
    // -------------------------------
    surfaceScalarField &phi(*phiPtr);

#ifdef HAVE_OFE20
    pressureControlPtr = new pressureControl
    (
        p,
        rho,
        pimple.dict(),
        false
    );
#elif defined(HAVE_OF7)
    pressureControlPtr = new pressureControl
    (
        p,
        rho,
        pimple.dict(),
        false
    );
#elif defined(HAVE_OF8)
    pressureControlPtr = new pressureControl
    (
        p,
        rho,
        pimple.dict(),
        thermo.incompressible()
    );
#endif

    mesh.setFluxRequired(p.name());

    Info << "Creating turbulence model\n" << endl;

#ifdef HAVE_OFE20
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
#elif defined(HAVE_OF7)
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
#elif defined(HAVE_OF8)
    turbulencePtr = autoPtr<compressible::momentumTransportModel>
    (
        compressible::momentumTransportModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );

    Info << "Creating thermophysical transport model\n" << endl;
    thermophysicalTransportPtr = autoPtr<fluidThermophysicalTransportModel> 
    (
        fluidThermophysicalTransportModel::New(turbulencePtr(), thermo)
    );
#endif

#ifdef HAVE_OFE20
    //  createDpdt.H
    createDpdt();
    // -------------

    //  createK.H
    createK();
    // ----------
#elif defined(HAVE_OF7) || defined(HAVE_OF8)
    Info << "Creating field dpdt\n" << endl;
    dpdtPtr = new volScalarField
    (
        IOobject
        (
            "dpdt",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(p.dimensions() / dimTime, 0)
    );

    Info << "Creating field kinetic energy K\n" << endl;
    KPtr = new volScalarField("K", 0.5 * magSqr(U));
#endif


#ifdef HAVE_OF8
    if (initialMassPtr == nullptr)
    {
        initialMassPtr = new dimensionedScalar(fvc::domainIntegrate(rho));
    }
    else
    {
        *initialMassPtr = fvc::domainIntegrate(rho);
    }
#endif

    //  createMRF.H
    createMRF();
    // ------------

#ifdef HAVE_OFE20
    rhoMaxPtr = new dimensionedScalar("rhoMax", dimDensity, GREAT, pimple.dict());
    rhoMinPtr = new dimensionedScalar("rhoMin", dimDensity, Zero, pimple.dict());
#endif

    //  createFvOptions.H
    createFvOptions();
    // ------------------

    return 0;
}

#ifdef HAVE_OFE20
void rhoPimple::createDpdt()
{
    Foam::Time &runTime(*runTimePtr);
    dynamicFvMesh &mesh(*meshPtr);
    volScalarField &p(*pPtr);
    fluidThermo &thermo(*pThermoPtr);
    
    IOobject dpdtHeader
    (
        "dpdt",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    if (mesh.dynamic())
    {
        Info<< "Creating field dpdt for moving meshes\n" << endl;

        // Note
        // - set to READ_IF_PRESENT and AUTO_WRITE to simplify dpdt correction
        //   by meshPhi

        dpdtHeader.readOpt() = IOobject::READ_IF_PRESENT;
        dpdtHeader.writeOpt() = IOobject::AUTO_WRITE;
    }
    else
    {
        Info<< "Creating field dpdt\n" << endl;
    }

    dpdtPtr = new volScalarField(dpdtHeader, fvc::ddt(p));
    volScalarField& dpdt(*dpdtPtr);

    if (!thermo.dpdt())
    {
        dpdt == dimensionedScalar(dpdt.dimensions(), Zero);
        dpdt.writeOpt() = IOobject::NO_WRITE;
    }
}

void rhoPimple::createK()
{
    volVectorField &U(*UPtr);

    Info<< "Creating field kinetic energy K\n" << endl;
    KPtr = new volScalarField("K", 0.5 * magSqr(U));
    volScalarField &K(*KPtr);

    if (U.nOldTimes())
    {
        volVectorField* Uold = &U.oldTime();
        volScalarField* Kold = &K.oldTime();
        *Kold == 0.5*magSqr(*Uold);

        while (Uold->nOldTimes())
        {
            Uold = &Uold->oldTime();
            Kold = &Kold->oldTime();
            *Kold == 0.5*magSqr(*Uold);
        }
    }
}
#endif

int rhoPimple::createMRF()
{
    dynamicFvMesh &mesh(*meshPtr);

    MRFPtr = new IOMRFZoneList(mesh);

    return 0;
}

int rhoPimple::compressibleCreatePhi()
{
    Foam::Time &runTime(*runTimePtr);
    dynamicFvMesh &mesh(*meshPtr);
    volScalarField &rho(*rhoPtr);
    volVectorField &U(*UPtr);

    Info << "Reading/calculating face flux field phi\n" << endl;

    phiPtr = new surfaceScalarField
    (
        IOobject("phi", runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
        linearInterpolate(rho * U) & mesh.Sf()
    );

    return 0;
}

int rhoPimple::createFvOptions()
{
    dynamicFvMesh &mesh(*meshPtr);

    fvOptionsPtr = new Foam::fv::options(mesh);
    Foam::fv::options &fvOptions(*fvOptionsPtr);

    if (!fvOptions.optionList::size())
    {
        Info << "No finite volume options present" << endl;
    }

    return 0;
}

int rhoPimple::createFieldRefs()
{
    fluidThermo &thermo(*pThermoPtr);

    psiPtr = &thermo.psi();

    return 0;
}

int rhoPimple::createRhoUfIfPresent()
{
    Foam::Time &runTime(*runTimePtr);
    dynamicFvMesh &mesh(*meshPtr);
    volScalarField &rho(*rhoPtr);
    volVectorField &U(*UPtr);

    if (mesh.dynamic())
    {
        Info << "Constructing face momentum rhoUf" << endl;

#ifdef HAVE_OFE20
        rhoUfPtr.reset
        (
            new surfaceVectorField
            (
                IOobject
                (
                    "rhoUf",
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                fvc::interpolate(rho*U)
            )
        );
#elif defined(HAVE_OF7) || defined(HAVE_OF8)
        rhoUfPtr = new surfaceVectorField
        (
            IOobject
            (
                "rhoUf",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::interpolate(rho * U)
        );
#endif
    }

    return 0;
}

int rhoPimple::compressibleCourantNo()
{
    Foam::Time &runTime(*runTimePtr);
    dynamicFvMesh &mesh(*meshPtr);
    volScalarField &rho(*rhoPtr);
    surfaceScalarField &phi(*phiPtr);

    CoNum = 0.0;
    meanCoNum = 0.0;

    {
        scalarField sumPhi
        (
            fvc::surfaceSum
                (mag(phi))().primitiveField()
              / rho.primitiveField()
        );

        CoNum = 0.5 * gMax(sumPhi / mesh.V().field()) * runTime.deltaTValue();

        meanCoNum = 0.5 * (gSum(sumPhi) / gSum(mesh.V().field()))
                  * runTime.deltaTValue();
    }

    Info << "Courant Number mean: " << meanCoNum << " max: " << CoNum << endl;

    return 0;
}

int rhoPimple::setInitialDeltaT()
{
    Foam::Time &runTime(*runTimePtr);

    if (adjustTimeStep)
    {
#ifdef HAVE_OFE20
        if ((runTime.timeIndex() == 0) && (CoNum > SMALL))
#elif defined(HAVE_OF7) || defined(HAVE_OF8)
        if ((runTime.timeIndex() == 0) && (CoNum > small))
#endif
        {
            runTime.setDeltaT
            (
                min
                (
                    maxCo*runTime.deltaTValue()/CoNum,
                    min(runTime.deltaTValue(), maxDeltaT)
                )
            );
        }
    }

    return 0;
}

int rhoPimple::loop()
{
    Foam::Time &runTime(*runTimePtr);
    dynamicFvMesh &mesh(*meshPtr);
    pimpleControl &pimple(*pimplePtr);
    volScalarField &rho(*rhoPtr);
    volVectorField &U(*UPtr);
    surfaceScalarField &phi(*phiPtr);
    fluidThermo &thermo(*pThermoPtr);
    IOMRFZoneList &MRF(*MRFPtr);

#ifdef HAVE_OFE20
    compressible::turbulenceModel& turbulence(*turbulencePtr);
#elif defined(HAVE_OF7)
    compressible::turbulenceModel& turbulence(*turbulencePtr);
#elif defined(HAVE_OF8)
    compressible::momentumTransportModel& turbulence(*turbulencePtr);
    fluidThermophysicalTransportModel& thermoTransModel(*thermophysicalTransportPtr);
#endif

    autoPtr<surfaceVectorField> &rhoUf(rhoUfPtr);

    Info << "\nStarting time loop\n" << endl;

#ifdef HAVE_OFE20
    while (runTime.run())
#elif defined(HAVE_OF7)
    while (runTime.run())
#elif defined(HAVE_OF8)
    while (pimple.run(runTime))
#endif
    {
        //  readDyMControls.H
        readDyMControls();
        // ------------------

        // Store divrhoU from the previous mesh so that it can be mapped
        // and used in correctPhi to ensure the corrected phi has the
        // same divergence
        if (correctPhi)
        {
            divrhoUPtr.reset
            (
                new volScalarField
                (
                    "divrhoU",
                    fvc::div(fvc::absolute(phi, rho, U))
                )
            );
        }

        if (LTS)
        {
            //  setRDeltaT.H
            setRDeltaT();
            // -------------
        }
        else
        {
            //  compressibleCourantNo.H
            compressibleCourantNo();
            // ------------------------

            //  setDeltaT.H
            setDeltaT();
            // ------------
        }

        runTime++;

        Info << "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
#ifdef HAVE_OFE20
            if (pimple.firstIter() || moveMeshOuterCorrectors)
#elif defined(HAVE_OF7) || defined(HAVE_OF8)
            if (pimple.firstPimpleIter() || moveMeshOuterCorrectors)
#endif
            {
                // Store momentum to set rhoUf for introduced faces.
                autoPtr<volVectorField> rhoU;
                if (rhoUf.valid())
                {
                    rhoU.reset(new volVectorField("rhoU", rho*U));
                }
                // Do any mesh changes
#ifdef HAVE_OFE20
                mesh.controlledUpdate();
#elif defined(HAVE_OF7) || defined(HAVE_OF8)
                mesh.update();
#endif

                if (mesh.changing())
                {
                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & rhoUf();

                        //  correctPhi.H
                        correctPhi_();
                        // -------------

                        // Make the fluxes relative to the mesh-motion
                        fvc::makeRelative(phi, rho, U);
                    }

                    if (checkMeshCourantNo)
                    {
                        //  meshCourantNo.H
                        meshCourantNo();
                        // ----------------
                    }
                }
            }

#ifdef HAVE_OFE20
            if (pimple.firstIter() && !pimple.SIMPLErho())
#elif defined(HAVE_OF7)
            if (pimple.firstPimpleIter() && !pimple.simpleRho())
#elif defined(HAVE_OF8)
            if
            (
                !mesh.steady()
             && !pimple.simpleRho()
             && pimple.firstPimpleIter()
            )
#endif
            {
                //  rhoEqn.H
                rhoEqn_();
                // ---------
            }

            //  UEqn.H
            UEqn_();
            // -------

            //  EEqn.H
            EEqn_();
            // -------

            // --- Pressure corrector loop
            while (pimple.correct())
#ifdef HAVE_OFE20
            //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            {
                if (pimple.consistent())
                {
                    //  pcEqn.H
                    pcEqn();
                    // --------
                }
                else
                {
                    //  pEqn.H
                    pEqn_();
                    // -------
                }
            }

            if (pimple.turbCorr())
            {
                turbulence.correct();
            }
        }

        rho = thermo.rho();
        //-------------------------------------------
#elif defined(HAVE_OF7)
            //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            {
                if (pimple.consistent())
                {
                    //  pcEqn.H
                    pcEqn();
                    // --------
                }
                else
                {
                    //  pEqn.H
                    pEqn_();
                    // -------
                }
            }

            if (pimple.turbCorr())
            {
                turbulence.correct();
            }
        }

        rho = thermo.rho();
        //-------------------------------------------
#elif defined(HAVE_OF8)
            //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            {
                //  pEqn.H
                pEqn_();
                // --------
            }

            if (pimple.turbCorr())
            {
                turbulence.correct();
                thermoTransModel.correct();
            }
        }

        if (!mesh.steady())
        {
            rho = thermo.rho();
        }
        //-------------------------------------------
#endif

        runTime.write();

#ifdef HAVE_OFE20
        runTime.printExecutionTime(Info);
#elif defined(HAVE_OF7) || defined(HAVE_OF8)
        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl
             << endl;
#endif
    }

    Info << "End\n" << endl;

    loopStat = 0;
    return loopStat;
}


int rhoPimple::step(double* incomingDeltaT, int* gmHandle)
{
    Foam::Time &runTime(*runTimePtr);
    dynamicFvMesh &mesh(*meshPtr);
    pimpleControl &pimple(*pimplePtr);
    volScalarField &rho(*rhoPtr);
    volVectorField &U(*UPtr);
    surfaceScalarField &phi(*phiPtr);
    fluidThermo &thermo(*pThermoPtr);
    IOMRFZoneList &MRF(*MRFPtr);

#ifdef HAVE_OFE20
    compressible::turbulenceModel& turbulence(*turbulencePtr);
#elif defined(HAVE_OF7)
    compressible::turbulenceModel& turbulence(*turbulencePtr);
#elif defined(HAVE_OF8)
    compressible::momentumTransportModel& turbulence(*turbulencePtr);
    fluidThermophysicalTransportModel& thermoTransModel(*thermophysicalTransportPtr);
#endif

    autoPtr<surfaceVectorField> &rhoUf(rhoUfPtr);

    double mandatedTime{0};
    if (incomingDeltaT != nullptr)
    {
        mandatedTime = runTime.value() + *incomingDeltaT;
    }

    int count{0};
    double alpha{0};
    bool continueIter{true};

#ifdef HAVE_OFE20
    while
    (
        runTime.run() &&
        continueIter
    )
#elif defined(HAVE_OF7)
    while
    (
        runTime.run() &&
        continueIter
    )
#elif defined(HAVE_OF8)
    while
    (
        pimple.run(runTime) &&
        continueIter
    )    
#endif
    {
        count++;
        Info << ">>MultiPhysics outer iteration "
             << count << "<<" << endl;

        if (modifiedDeltaT)
        {
            runTime.setDeltaT(unmodifiedDeltaTvalue);

            runTime.controlDict().lookup("adjustTimeStep");
        }

        //  readDyMControls.H
        readDyMControls();
        // ------------------

        // Store divrhoU from the previous mesh so that it can be mapped
        // and used in correctPhi to ensure the corrected phi has the
        // same divergence
        if (correctPhi)
        {
            divrhoUPtr.reset
            (
                new volScalarField
                (
                    "divrhoU",
                    fvc::div(fvc::absolute(phi, rho, U))
                )
            );
        }

        if (LTS)
        {
            //  setRDeltaT.H
            setRDeltaT();
            // -------------
        }
        else //if (count==1)
        {
            //  compressibleCourantNo.H
            compressibleCourantNo();
            // ------------------------

            //  setDeltaT.H
            setDeltaT();
            // ------------

            double flowDeltaT = runTime.deltaTValue();
            double flowCurTime = runTime.value();
            double expectedTime = flowCurTime + flowDeltaT;
            
            if (incomingDeltaT != nullptr)
            {
                modifiedDeltaT = false;

                if (std::abs( expectedTime - mandatedTime ) < 0.0001*flowDeltaT)
                {
                    continueIter = false;
                }
                else if (expectedTime > mandatedTime)
                {
                    double newDeltaT = mandatedTime - flowCurTime;
                    
#ifdef HAVE_OFE20
                    bool adjust{false};
                    runTime.setDeltaT(newDeltaT, adjust);
#elif defined(HAVE_OF7) || defined(HAVE_OF8)
                    runTime.setDeltaTNoAdjust(newDeltaT);
#endif
                    Info << "NewdeltaT according to the Rocstar deltaT = " << newDeltaT
                         << endl;

                    modifiedDeltaT = true;
                    unmodifiedDeltaTvalue = flowDeltaT;

                    continueIter = false;
                }
                else if (expectedTime == mandatedTime)
                {
                    continueIter = false;
                }

                /*
                if (expectedTime >= mandatedTime)
                {
                    continueIter = false;
                }
                */

                const double& incomingDeltaT_{*incomingDeltaT};
                flowDeltaT = runTime.deltaTValue();
                
                alpha +=  flowDeltaT / incomingDeltaT_;

                if (gmHandle != nullptr)
                {
                    if (*gmHandle >= 0)
                    {
                        COM_call_function(*gmHandle, &alpha);

                        updateSurfaceData_incoming(count);

                        Info << " alpha = " << alpha << endl;
                    }
                }
            }
            else
            {
                updateSurfaceData_incoming();
                continueIter = false;
            }

            if (runTime.deltaTValue() < 0)
            {
                Info << "Unphysical deltaT. Exiting the simulation" << endl;
                exit(-1);
            }
        }

        runTime++;

        Info << "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
#ifdef HAVE_OFE20
            if (pimple.firstIter() || moveMeshOuterCorrectors)
#elif defined(HAVE_OF7) || defined(HAVE_OF8)
            if (pimple.firstPimpleIter() || moveMeshOuterCorrectors)
#endif
            {
                // Store momentum to set rhoUf for introduced faces.
                autoPtr<volVectorField> rhoU;
                if (rhoUf.valid())
                {
                    rhoU.reset(new volVectorField("rhoU", rho*U));
                }
                // Do any mesh changes
#ifdef HAVE_OFE20
                mesh.controlledUpdate();
#elif defined(HAVE_OF7) || defined(HAVE_OF8)
                mesh.update();
#endif

                if (mesh.changing())
                {
                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & rhoUf();

                        //  correctPhi.H
                        correctPhi_();
                        // -------------

                        // Make the fluxes relative to the mesh-motion
                        fvc::makeRelative(phi, rho, U);
                    }

                    if (checkMeshCourantNo)
                    {
                        //  meshCourantNo.H
                        meshCourantNo();
                        // ----------------
                    }
                }
            }

#ifdef HAVE_OFE20
            if (pimple.firstIter() && !pimple.SIMPLErho())
#elif defined(HAVE_OF7)
            if (pimple.firstPimpleIter() && !pimple.simpleRho())
#elif defined(HAVE_OF8)
            if
            (
                !mesh.steady()
             && !pimple.simpleRho()
             && pimple.firstPimpleIter()
            )
#endif
            {
                //  rhoEqn.H
                rhoEqn_();
                // ---------
            }

            //  UEqn.H
            UEqn_();
            // -------

            //  EEqn.H
            EEqn_();
            // -------

            // --- Pressure corrector loop
            while (pimple.correct())
#ifdef HAVE_OFE20
            //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            {
                if (pimple.consistent())
                {
                    //  pcEqn.H
                    pcEqn();
                    // --------
                }
                else
                {
                    //  pEqn.H
                    pEqn_();
                    // -------
                }
            }

            if (pimple.turbCorr())
            {
                turbulence.correct();
            }
        }

        rho = thermo.rho();
        //-------------------------------------------
#elif defined(HAVE_OF7)
            //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            {
                if (pimple.consistent())
                {
                    //  pcEqn.H
                    pcEqn();
                    // --------
                }
                else
                {
                    //  pEqn.H
                    pEqn_();
                    // -------
                }
            }

            if (pimple.turbCorr())
            {
                turbulence.correct();
            }
        }

        rho = thermo.rho();
        //-------------------------------------------
#elif defined(HAVE_OF8)
        //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            {
                //  pcEqn.H
                pEqn_();
                // --------
            }

            if (pimple.turbCorr())
            {
                turbulence.correct();
                thermoTransModel.correct();
            }
        }

        if (!mesh.steady())
        {
            rho = thermo.rho();
        }
        //-------------------------------------------
#endif

        runTime.write();

#ifdef HAVE_OFE20
        runTime.printExecutionTime(Info);
#elif defined(HAVE_OF7) || defined(HAVE_OF8)
        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl
             << endl;
#endif
    }

    stepStat = 0;
    return stepStat;
}


int rhoPimple::readDyMControls()
{
    pimpleControl &pimple(*pimplePtr);

    //  readTimeControls.H  ^^^^^^^^^^^^^^^^^^^^^^^
    readTimeControls();
    // --------------------------------------------

    correctPhi = pimple.dict().lookupOrDefault("correctPhi", correctPhi);

    checkMeshCourantNo = pimple.dict().lookupOrDefault
    (
        "checkMeshCourantNo",
        checkMeshCourantNo
    );

    moveMeshOuterCorrectors = pimple.dict().lookupOrDefault
    (
        "moveMeshOuterCorrectors",
        moveMeshOuterCorrectors
    );

    return 0;
}

int rhoPimple::setRDeltaT()
{
    Foam::Time &runTime(*runTimePtr);
    dynamicFvMesh &mesh(*meshPtr);
    pimpleControl &pimple(*pimplePtr);
    volScalarField &rho(*rhoPtr);
    volVectorField &U(*UPtr);
    surfaceScalarField &phi(*phiPtr);
    const volScalarField &psi(*psiPtr);

    volScalarField &rDeltaT = trDeltaT.ref();

    const dictionary &pimpleDict = pimple.dict();

#ifdef HAVE_OFE20
    maxCo = pimpleDict.getOrDefault<scalar>("maxCo", 0.8);

    scalar rDeltaTSmoothingCoeff
    (
        pimpleDict.getOrDefault<scalar>("rDeltaTSmoothingCoeff", 0.02)
    );

    scalar rDeltaTDampingCoeff
    (
        pimpleDict.getOrDefault<scalar>("rDeltaTDampingCoeff", 1.0)
    );

    maxDeltaT = pimpleDict.getOrDefault<scalar>("maxDeltaT", GREAT);

    volScalarField rDeltaT0("rDeltaT0", rDeltaT);

    // Set the reciprocal time-step from the local Courant number
    rDeltaT.ref() = max
    (
        1/dimensionedScalar("maxDeltaT", dimTime, maxDeltaT),
        fvc::surfaceSum(mag(phi))()()
       /((2*maxCo)*mesh.V()*rho())
    );
#elif defined(HAVE_OF7) || defined(HAVE_OF8)
    maxCo = pimpleDict.lookupOrDefault<scalar>("maxCo", 0.8);

    scalar rDeltaTSmoothingCoeff
    (
        pimpleDict.lookupOrDefault<scalar>("rDeltaTSmoothingCoeff", 0.02)
    );

    scalar rDeltaTDampingCoeff
    (
        pimpleDict.lookupOrDefault<scalar>("rDeltaTDampingCoeff", 1.0)
    );

    maxDeltaT = pimpleDict.lookupOrDefault<scalar>("maxDeltaT", great);

    volScalarField rDeltaT0("rDeltaT0", rDeltaT);

    // Set the reciprocal time-step from the local Courant number
    rDeltaT.ref() =
        max(1 / dimensionedScalar(dimTime, maxDeltaT),
            fvc::surfaceSum(mag(phi))()() / ((2 * maxCo) * mesh.V() * rho()));
#endif

    if (pimple.transonic())
    {
        surfaceScalarField phid
        (
            "phid",
            fvc::interpolate(psi) * fvc::flux(U)
        );

        rDeltaT.ref() = max
        (
            rDeltaT(),
            fvc::surfaceSum(mag(phid))()() / ((2 * maxCo) * mesh.V() * psi())
        );
    }

    // Update the boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    Info << "Flow time scale min/max = " << gMin(1 / rDeltaT.primitiveField())
         << ", " << gMax(1 / rDeltaT.primitiveField()) << endl;

    if (rDeltaTSmoothingCoeff < 1.0) {
        fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);
    }

    Info << "Smoothed flow time scale min/max = "
         << gMin(1 / rDeltaT.primitiveField()) << ", "
         << gMax(1 / rDeltaT.primitiveField()) << endl;

    // Limit rate of change of time scale
    // - reduce as much as required
    // - only increase at a fraction of old time scale
    if (rDeltaTDampingCoeff < 1.0 &&
        runTime.timeIndex() > runTime.startTimeIndex() + 1) {
        rDeltaT =
            rDeltaT0 * max(rDeltaT / rDeltaT0, scalar(1) - rDeltaTDampingCoeff);

        Info << "Damped flow time scale min/max = "
             << gMin(1 / rDeltaT.primitiveField()) << ", "
             << gMax(1 / rDeltaT.primitiveField()) << endl;
    }

    return 0;
}

int rhoPimple::correctPhi_()
{
    pimpleControl &pimple(*pimplePtr);
    volScalarField &rho(*rhoPtr);
    volVectorField &U(*UPtr);
    surfaceScalarField &phi(*phiPtr);
    const volScalarField &psi(*psiPtr);
    volScalarField &divrhoU(*divrhoUPtr);
    volScalarField &p(*pPtr);

#ifdef HAVE_OFE20
    CorrectPhi
    (
        U, phi, p, rho, psi,
        dimensionedScalar("rAUf", dimTime, 1),
        divrhoU(), pimple
    );
#elif defined(HAVE_OF7) || defined(HAVE_OF8)
    CorrectPhi
    (
        U, phi, p, rho, psi,
        dimensionedScalar("rAUf", dimTime, 1),
        divrhoU(), pimple, true
    );
#endif

    return 0;
}

int rhoPimple::meshCourantNo()
{
    dynamicFvMesh &mesh(*meshPtr);
    Foam::Time &runTime(*runTimePtr);

    scalar meshCoNum = 0.0;
    scalar meanMeshCoNum = 0.0;

    {
        scalarField sumPhi(fvc::surfaceSum(mag(mesh.phi()))().primitiveField());

        meshCoNum =
            0.5 * gMax(sumPhi / mesh.V().field()) * runTime.deltaTValue();

        meanMeshCoNum = 0.5 * (gSum(sumPhi) / gSum(mesh.V().field())) *
                        runTime.deltaTValue();
    }

    Info << "Mesh Courant Number mean: " << meanMeshCoNum
         << " max: " << meshCoNum << endl;

    return 0;
}

int rhoPimple::rhoEqn_()
{
    volScalarField &rho(*rhoPtr);
    surfaceScalarField &phi(*phiPtr);
    Foam::fv::options &fvOptions(*fvOptionsPtr);

    fvScalarMatrix rhoEqn(fvm::ddt(rho) + fvc::div(phi) == fvOptions(rho));

    fvOptions.constrain(rhoEqn);

    rhoEqn.solve();

    fvOptions.correct(rho);

    return 0;
}

int rhoPimple::UEqn_()
{
    pimpleControl &pimple(*pimplePtr);
    volScalarField &rho(*rhoPtr);
    volVectorField &U(*UPtr);
    surfaceScalarField &phi(*phiPtr);
    Foam::fv::options &fvOptions(*fvOptionsPtr);
    IOMRFZoneList &MRF(*MRFPtr);
#ifdef HAVE_OFE20
    compressible::turbulenceModel& turbulence(*turbulencePtr);
#elif defined(HAVE_OF7)
    compressible::turbulenceModel& turbulence(*turbulencePtr);
#elif defined(HAVE_OF8)
    compressible::momentumTransportModel& turbulence(*turbulencePtr);
#endif
    volScalarField &K(*KPtr);
    volScalarField &p(*pPtr);

    MRF.correctBoundaryVelocity(U);

#ifdef HAVE_OFE20
    tUEqnPtr = tmp<fvVectorMatrix>
    (
        fvm::ddt(rho, U) + fvm::div(phi, U)
      + MRF.DDt(rho, U)
      + turbulence.divDevRhoReff(U)
     == fvOptions(rho, U)
    );
#elif defined(HAVE_OF7)
    tUEqnPtr = tmp<fvVectorMatrix>
    (
        fvm::ddt(rho, U) + fvm::div(phi, U)
      + MRF.DDt(rho, U)
      + turbulence.divDevRhoReff(U)
     == fvOptions(rho, U)
    );
#elif defined(HAVE_OF8)
    tUEqnPtr = tmp<fvVectorMatrix>
    (
        fvm::ddt(rho, U) + fvm::div(phi, U)
      + MRF.DDt(rho, U)
      + turbulence.divDevTau(U)
     == fvOptions(rho, U)
    );
#endif

    tmp<fvVectorMatrix> &tUEqn(tUEqnPtr);

    UEqnPtr = &tUEqn.ref();
    fvVectorMatrix &UEqn(*UEqnPtr);

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));

        fvOptions.correct(U);
        K = 0.5 * magSqr(U);
    }

    return 0;
}

int rhoPimple::EEqn_()
{
    volScalarField &rho(*rhoPtr);
    volVectorField &U(*UPtr);
    surfaceScalarField &phi(*phiPtr);
    Foam::fv::options &fvOptions(*fvOptionsPtr);
    fluidThermo &thermo(*pThermoPtr);
#ifdef HAVE_OFE20
    compressible::turbulenceModel& turbulence(*turbulencePtr);
#elif defined(HAVE_OF7)
    compressible::turbulenceModel& turbulence(*turbulencePtr);
#elif defined(HAVE_OF8)
    //compressible::momentumTransportModel& turbulence(*turbulencePtr);
    const fluidThermophysicalTransportModel& thermoTransModel(*thermophysicalTransportPtr);
#endif
    volScalarField &K(*KPtr);
    volScalarField &dpdt(*dpdtPtr);
    volScalarField &p(*pPtr);

    volScalarField &he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + fvm::div(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      +
        (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi / fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            ) : -dpdt
        )
#ifdef HAVE_OFE20
      - fvm::laplacian(turbulence.alphaEff(), he)
#elif  defined(HAVE_OF7)
      - fvm::laplacian(turbulence.alphaEff(), he)
#elif defined(HAVE_OF8)
      + thermoTransModel.divq(he)
#endif
     == fvOptions(rho, he));

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();

    return 0;
}

#ifdef HAVE_OF7
int rhoPimple::pcEqn()
{
    dynamicFvMesh &mesh(*meshPtr);
    pimpleControl &pimple(*pimplePtr);
    volScalarField &rho(*rhoPtr);
    volVectorField &U(*UPtr);
    surfaceScalarField &phi(*phiPtr);
    pressureControl &pressureControl(*pressureControlPtr);
    Foam::fv::options &fvOptions(*fvOptionsPtr);
    fluidThermo &thermo(*pThermoPtr);
    const volScalarField &psi(*psiPtr);
    IOMRFZoneList &MRF(*MRFPtr);
    autoPtr<surfaceVectorField> &rhoUf(rhoUfPtr);
    volScalarField &K(*KPtr);
    volScalarField &dpdt(*dpdtPtr);
    volScalarField &p(*pPtr);
    fvVectorMatrix &UEqn(*UEqnPtr);
    tmp<fvVectorMatrix> &tUEqn(tUEqnPtr);

    if (!pimple.simpleRho())
    {
        rho = thermo.rho();
    }

    // Thermodynamic density needs to be updated by psi*d(p) after the
    // pressure solution
    const volScalarField psip0(psi * p);

    volScalarField rAU(1.0 / UEqn.A());
    volScalarField rAtU(1.0 / (1.0 / rAU - UEqn.H1()));
    volVectorField HbyA(constrainHbyA(rAU * UEqn.H(), U, p));

    if (pimple.nCorrPiso() <= 1)
    {
        tUEqn.clear();
    }

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
            fvc::interpolate(rho) * fvc::flux(HbyA) 
          + MRF.zeroFilter
            (
                fvc::interpolate(rho * rAU) * fvc::ddtCorr(rho, U, phi, rhoUf))
        )
    );

    fvc::makeRelative(phiHbyA, rho, U);
    MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

    volScalarField rhorAtU("rhorAtU", rho * rAtU);

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, rho, U, phiHbyA, rhorAtU, MRF);

    if (pimple.transonic())
    {
        surfaceScalarField phid
        (
            "phid",
            (fvc::interpolate(psi) / fvc::interpolate(rho)) * phiHbyA
        );

        phiHbyA +=
            fvc::interpolate(rho * (rAtU - rAU))
          * fvc::snGrad(p) * mesh.magSf()
          - fvc::interpolate(psi * p) * phiHbyA / fvc::interpolate(rho);

        HbyA -= (rAU - rAtU) * fvc::grad(p);

        fvScalarMatrix pDDtEqn
        (
            fvc::ddt(rho) + psi * correction(fvm::ddt(p))
          + fvc::div(phiHbyA) + fvm::div(phid, p)
         == fvOptions(psi, p, rho.name())
        );

        while (pimple.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn(pDDtEqn - fvm::laplacian(rhorAtU, p));

            // Relax the pressure equation to ensure diagonal-dominance
            pEqn.relax();

            pEqn.solve();

            if (pimple.finalNonOrthogonalIter())
            {
                phi = phiHbyA + pEqn.flux();
            }
        }
    }
    else
    {
        phiHbyA += fvc::interpolate(rho * (rAtU - rAU))
                 * fvc::snGrad(p) * mesh.magSf();
        HbyA -= (rAU - rAtU) * fvc::grad(p);

        fvScalarMatrix pDDtEqn
        (
            fvc::ddt(rho) + psi * correction(fvm::ddt(p))
          + fvc::div(phiHbyA)
         == fvOptions(psi, p, rho.name())
        );

        while (pimple.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn(pDDtEqn - fvm::laplacian(rhorAtU, p));

            pEqn.solve();

            if (pimple.finalNonOrthogonalIter())
            {
                phi = phiHbyA + pEqn.flux();
            }
        }
    }

    bool limitedp = pressureControl.limit(p);

    // Thermodynamic density update
    thermo.correctRho(psi * p - psip0);

    if (limitedp)
    {
        rho = thermo.rho();
    }

    //  rhoEqn.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    rhoEqn_();
    // ---------------------------------------------

    //  compressibleContinuityErrs.H  ^^^^^^^^^^^^^^
    compressibleContinuityErrs();
    // ---------------------------------------------

    // Explicitly relax pressure for momentum corrector
    p.relax();

    U = HbyA - rAtU * fvc::grad(p);
    U.correctBoundaryConditions();
    fvOptions.correct(U);
    K = 0.5 * magSqr(U);

    if (pimple.simpleRho())
    {
        rho = thermo.rho();
    }

    // Correct rhoUf if the mesh is moving
    fvc::correctRhoUf(rhoUf, rho, U, phi);

    if (thermo.dpdt())
    {
        dpdt = fvc::ddt(p);

        if (mesh.moving())
        {
            dpdt -= fvc::div(fvc::meshPhi(rho, U), p);
        }
    }

    return 0;
}
#elif defined(HAVE_OFE20)
int rhoPimple::pcEqn()
{
    dynamicFvMesh &mesh(*meshPtr);
    pimpleControl &pimple(*pimplePtr);
    volScalarField &rho(*rhoPtr);
    volVectorField &U(*UPtr);
    surfaceScalarField &phi(*phiPtr);
    pressureControl &pressureControl(*pressureControlPtr);
    Foam::fv::options &fvOptions(*fvOptionsPtr);
    fluidThermo &thermo(*pThermoPtr);
    const volScalarField &psi(*psiPtr);
    IOMRFZoneList &MRF(*MRFPtr);
    autoPtr<surfaceVectorField> &rhoUf(rhoUfPtr);
    volScalarField &K(*KPtr);
    volScalarField &dpdt(*dpdtPtr);
    volScalarField &p(*pPtr);
    fvVectorMatrix &UEqn(*UEqnPtr);
    tmp<fvVectorMatrix> &tUEqn(tUEqnPtr);

    if (!pimple.SIMPLErho())
    {
        rho = thermo.rho();
    }

    // Thermodynamic density needs to be updated by psi*d(p) after the
    // pressure solution
    const volScalarField psip0(psi * p);

    volScalarField rAU(1.0 / UEqn.A());
    volScalarField rAtU(1.0 / (1.0 / rAU - UEqn.H1()));
    volVectorField HbyA(constrainHbyA(rAU * UEqn.H(), U, p));

    if (pimple.nCorrPISO() <= 1)
    {
        tUEqn.clear();
    }

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
            fvc::interpolate(rho) * fvc::flux(HbyA) 
          + MRF.zeroFilter
            (
                fvc::interpolate(rho * rAU) * fvc::ddtCorr(rho, U, phi, rhoUf))
        )
    );

    fvc::makeRelative(phiHbyA, rho, U);
    MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

    volScalarField rhorAtU("rhorAtU", rho * rAtU);

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, rho, U, phiHbyA, rhorAtU, MRF);

    if (pimple.transonic())
    {
        surfaceScalarField phid
        (
            "phid",
            (fvc::interpolate(psi) / fvc::interpolate(rho)) * phiHbyA
        );

        phiHbyA +=
            fvc::interpolate(rho * (rAtU - rAU))
          * fvc::snGrad(p) * mesh.magSf()
          - fvc::interpolate(psi * p) * phiHbyA / fvc::interpolate(rho);

        HbyA -= (rAU - rAtU) * fvc::grad(p);

        fvScalarMatrix pDDtEqn
        (
            fvc::ddt(rho) + psi * correction(fvm::ddt(p))
          + fvc::div(phiHbyA) + fvm::div(phid, p)
         == fvOptions(psi, p, rho.name())
        );

        while (pimple.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn(pDDtEqn - fvm::laplacian(rhorAtU, p));

            // Relax the pressure equation to ensure diagonal-dominance
            pEqn.relax();

            pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

            if (pimple.finalNonOrthogonalIter())
            {
                phi = phiHbyA + pEqn.flux();
            }
        }
    }
    else
    {
        phiHbyA += fvc::interpolate(rho * (rAtU - rAU))
                 * fvc::snGrad(p) * mesh.magSf();
        HbyA -= (rAU - rAtU) * fvc::grad(p);

        fvScalarMatrix pDDtEqn
        (
            fvc::ddt(rho) + psi * correction(fvm::ddt(p))
          + fvc::div(phiHbyA)
         == fvOptions(psi, p, rho.name())
        );

        while (pimple.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn(pDDtEqn - fvm::laplacian(rhorAtU, p));

            pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

            if (pimple.finalNonOrthogonalIter())
            {
                phi = phiHbyA + pEqn.flux();
            }
        }
    }

    //  rhoEqn.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    rhoEqn_();
    // ---------------------------------------------

    //  compressibleContinuityErrs.H  ^^^^^^^^^^^^^^
    compressibleContinuityErrs();
    // ---------------------------------------------

    // Explicitly relax pressure for momentum corrector
    p.relax();

    U = HbyA - rAtU * fvc::grad(p);
    U.correctBoundaryConditions();
    fvOptions.correct(U);
    K = 0.5 * magSqr(U);

    if (pressureControl.limit(p))
    {
        p.correctBoundaryConditions();
    }

    thermo.correctRho(psi*p - psip0, *rhoMinPtr, *rhoMaxPtr);
    rho = thermo.rho();

    // Correct rhoUf if the mesh is moving
    fvc::correctRhoUf(rhoUf, rho, U, phi);

    if (thermo.dpdt())
    {
        dpdt = fvc::ddt(p);

        if (mesh.moving())
        {
            dpdt -= fvc::div(fvc::meshPhi(rho, U), p);
        }
    }

    return 0;
}
#endif

int rhoPimple::compressibleContinuityErrs()
{
    volScalarField &rho(*rhoPtr);
    fluidThermo &thermo(*pThermoPtr);

    dimensionedScalar totalMass = fvc::domainIntegrate(rho);

    scalar sumLocalContErr =
        (fvc::domainIntegrate(mag(rho - thermo.rho())) / totalMass).value();

    scalar globalContErr =
        (fvc::domainIntegrate(rho - thermo.rho()) / totalMass).value();

    cumulativeContErr += globalContErr;

    Info << "time step continuity errors : sum local = " << sumLocalContErr
         << ", global = " << globalContErr
         << ", cumulative = " << cumulativeContErr << endl;

    return 0;
}

int rhoPimple::incompressibleContinuityErrs()
{
    dynamicFvMesh &mesh(*meshPtr);
    Foam::Time &runTime(*runTimePtr);
    surfaceScalarField &phi(*phiPtr);

    volScalarField contErr(fvc::div(phi));

    scalar sumLocalContErr = runTime.deltaTValue()*
        mag(contErr)().weightedAverage(mesh.V()).value();

    scalar globalContErr = runTime.deltaTValue()*
        contErr.weightedAverage(mesh.V()).value();
    cumulativeContErr += globalContErr;

    Info<< "time step continuity errors : sum local = " << sumLocalContErr
        << ", global = " << globalContErr
        << ", cumulative = " << cumulativeContErr
        << endl;

    return 0;
}


int rhoPimple::pEqn_()
{
    dynamicFvMesh &mesh(*meshPtr);
    pimpleControl &pimple(*pimplePtr);
    volScalarField &rho(*rhoPtr);
    volVectorField &U(*UPtr);
    surfaceScalarField &phi(*phiPtr);
    pressureControl &pressureControl(*pressureControlPtr);
    Foam::fv::options &fvOptions(*fvOptionsPtr);
    fluidThermo &thermo(*pThermoPtr);
    const volScalarField &psi(*psiPtr);
    IOMRFZoneList &MRF(*MRFPtr);
    volScalarField &K(*KPtr);
    autoPtr<surfaceVectorField> &rhoUf(rhoUfPtr);
    volScalarField &dpdt(*dpdtPtr);
    volScalarField &p(*pPtr);
    fvVectorMatrix &UEqn(*UEqnPtr);
    tmp<fvVectorMatrix> &tUEqn(tUEqnPtr);

#ifdef HAVE_OFE20
    //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    if (!pimple.SIMPLErho())
    {
        rho = thermo.rho();
    }

    // Thermodynamic density needs to be updated by psi*d(p) after the
    // pressure solution
    const volScalarField psip0(psi * p);

    volScalarField rAU(1.0 / UEqn.A());
    surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho * rAU));

    volVectorField HbyA(constrainHbyA(rAU * UEqn.H(), U, p));

    if (pimple.nCorrPISO() <= 1)
    {
        tUEqn.clear();
    }

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::interpolate(rho) * fvc::flux(HbyA)
      + MRF.zeroFilter(rhorAUf * fvc::ddtCorr(rho, U, phi, rhoUf))
    );

    fvc::makeRelative(phiHbyA, rho, U);
    MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, rho, U, phiHbyA, rhorAUf, MRF);

    if (pimple.transonic())
    {
        surfaceScalarField phid
        (
            "phid",
            (fvc::interpolate(psi) / fvc::interpolate(rho)) * phiHbyA
        );

        phiHbyA -= fvc::interpolate(psi * p) * phiHbyA / fvc::interpolate(rho);

        fvScalarMatrix pDDtEqn
        (
            fvc::ddt(rho) + psi * correction(fvm::ddt(p))
          + fvc::div(phiHbyA) + fvm::div(phid, p)
         == fvOptions(psi, p, rho.name())
        );

        while (pimple.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn(pDDtEqn - fvm::laplacian(rhorAUf, p));

            // Relax the pressure equation to ensure diagonal-dominance
            pEqn.relax();

            pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

            if (pimple.finalNonOrthogonalIter())
            {
                phi = phiHbyA + pEqn.flux();
            }
        }
    }
    else
    {
        fvScalarMatrix pDDtEqn
        (
            fvc::ddt(rho) + psi * correction(fvm::ddt(p))
          + fvc::div(phiHbyA)
         == fvOptions(psi, p, rho.name())
        );

        while (pimple.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn(pDDtEqn - fvm::laplacian(rhorAUf, p));
            
            pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

            if (pimple.finalNonOrthogonalIter())
            {
                phi = phiHbyA + pEqn.flux();
            }
        }
    }

    // Explicitly relax pressure for momentum corrector
    p.relax();

    U = HbyA - rAU * fvc::grad(p);
    U.correctBoundaryConditions();
    fvOptions.correct(U);
    K = 0.5 * magSqr(U);

    if (pressureControl.limit(p))
    {
        p.correctBoundaryConditions();
    }

    thermo.correctRho(psi*p - psip0, *rhoMinPtr, *rhoMaxPtr);

    //  rhoEqn.H
    rhoEqn_();
    // ---------

    //  compressibleContinuityErrs.H
    compressibleContinuityErrs();
    // -----------------------------

    rho = thermo.rho();

    // Correct rhoUf if the mesh is moving
    fvc::correctRhoUf(rhoUf, rho, U, phi);

    if (thermo.dpdt())
    {
        dpdt = fvc::ddt(p);

        if (mesh.moving())
        {
            dpdt -= fvc::div(fvc::meshPhi(rho, U), p);
        }
    }
    //-------------------------------------------
#elif defined(HAVE_OF7)
    //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    if (!pimple.simpleRho())
    {
        rho = thermo.rho();
    }

    // Thermodynamic density needs to be updated by psi*d(p) after the
    // pressure solution
    const volScalarField psip0(psi * p);

    volScalarField rAU(1.0 / UEqn.A());
    surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho * rAU));

    volVectorField HbyA(constrainHbyA(rAU * UEqn.H(), U, p));

    if (pimple.nCorrPiso() <= 1)
    {
        tUEqn.clear();
    }

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::interpolate(rho) * fvc::flux(HbyA)
      + MRF.zeroFilter(rhorAUf * fvc::ddtCorr(rho, U, phi, rhoUf))
    );

    fvc::makeRelative(phiHbyA, rho, U);
    MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, rho, U, phiHbyA, rhorAUf, MRF);

    if (pimple.transonic())
    {
        surfaceScalarField phid
        (
            "phid",
            (fvc::interpolate(psi) / fvc::interpolate(rho)) * phiHbyA
        );

        phiHbyA -= fvc::interpolate(psi * p) * phiHbyA / fvc::interpolate(rho);

        fvScalarMatrix pDDtEqn
        (
            fvc::ddt(rho) + psi * correction(fvm::ddt(p))
          + fvc::div(phiHbyA) + fvm::div(phid, p)
         == fvOptions(psi, p, rho.name())
        );

        while (pimple.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn(pDDtEqn - fvm::laplacian(rhorAUf, p));

            // Relax the pressure equation to ensure diagonal-dominance
            pEqn.relax();
            pEqn.solve();

            if (pimple.finalNonOrthogonalIter())
            {
                phi = phiHbyA + pEqn.flux();
            }
        }
    }
    else
    {
        fvScalarMatrix pDDtEqn
        (
            fvc::ddt(rho) + psi * correction(fvm::ddt(p))
          + fvc::div(phiHbyA)
         == fvOptions(psi, p, rho.name())
        );

        while (pimple.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn(pDDtEqn - fvm::laplacian(rhorAUf, p));
            pEqn.solve();

            if (pimple.finalNonOrthogonalIter())
            {
                phi = phiHbyA + pEqn.flux();
            }
        }
    }

    bool limitedp = pressureControl.limit(p);

    // Thermodynamic density update
    thermo.correctRho(psi * p - psip0);

    if (limitedp)
    {
        rho = thermo.rho();
    }

    //  rhoEqn.H
    rhoEqn_();
    // ---------

    //  compressibleContinuityErrs.H
    compressibleContinuityErrs();
    // -----------------------------

    // Explicitly relax pressure for momentum corrector
    p.relax();

    U = HbyA - rAU * fvc::grad(p);
    U.correctBoundaryConditions();
    fvOptions.correct(U);
    K = 0.5 * magSqr(U);

    if (pimple.simpleRho())
    {
        rho = thermo.rho();
    }

    // Correct rhoUf if the mesh is moving
    fvc::correctRhoUf(rhoUf, rho, U, phi);

    if (thermo.dpdt())
    {
        dpdt = fvc::ddt(p);

        if (mesh.moving())
        {
            dpdt -= fvc::div(fvc::meshPhi(rho, U), p);
        }
    }
    //-------------------------------------------
#elif defined(HAVE_OF8)
    //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    const dimensionedScalar& initialMass(*initialMassPtr);

    if ((!mesh.steady() && !pimple.simpleRho()) || pimple.consistent())
    {
        rho = thermo.rho();
    }

    // Thermodynamic density needs to be updated by psi*d(p) after the
    // pressure solution
    const volScalarField psip0(psi*p);

    const volScalarField rAU("rAU", 1.0/UEqn.A());
    const surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));

    tmp<volScalarField> rAtU
    (
        pimple.consistent()
    ? volScalarField::New("rAtU", 1.0/(1.0/rAU - UEqn.H1()))
    : tmp<volScalarField>(nullptr)
    );
    tmp<surfaceScalarField> rhorAtUf
    (
        pimple.consistent()
    ? surfaceScalarField::New("rhoRAtUf", fvc::interpolate(rho*rAtU()))
    : tmp<surfaceScalarField>(nullptr)
    );

    const volScalarField& rAAtU = pimple.consistent() ? rAtU() : rAU;
    const surfaceScalarField& rhorAAtUf =
        pimple.consistent() ? rhorAtUf() : rhorAUf;

    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

    if (pimple.nCorrPiso() <= 1)
    {
        tUEqn.clear();
    }

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::interpolate(rho)*fvc::flux(HbyA)
    + MRF.zeroFilter(rhorAUf*fvc::ddtCorr(rho, U, phi, rhoUf))
    );

    fvc::makeRelative(phiHbyA, rho, U);
    MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

    bool adjustMass = false;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, rho, U, phiHbyA, rhorAAtUf, MRF);

    if (pimple.transonic())
    {
        surfaceScalarField phid
        (
            "phid",
            (fvc::interpolate(psi)/fvc::interpolate(rho))*phiHbyA
        );

        phiHbyA -= fvc::interpolate(psi*p)*phiHbyA/fvc::interpolate(rho);

        if (pimple.consistent())
        {
            phiHbyA += (rhorAAtUf - rhorAUf)*fvc::snGrad(p)*mesh.magSf();
            HbyA += (rAAtU - rAU)*fvc::grad(p);
        }

        fvScalarMatrix pDDtEqn
        (
            fvc::ddt(rho) + psi*correction(fvm::ddt(p))
        + fvc::div(phiHbyA) + fvm::div(phid, p)
        ==
            fvOptions(psi, p, rho.name())
        );

        while (pimple.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn(pDDtEqn - fvm::laplacian(rhorAAtUf, p));

            // Relax the pressure equation to ensure diagonal-dominance
            pEqn.relax();

            pEqn.setReference
            (
                pressureControl.refCell(),
                pressureControl.refValue()
            );

            pEqn.solve();

            if (pimple.finalNonOrthogonalIter())
            {
                phi = phiHbyA + pEqn.flux();
            }
        }
    }
    else
    {
        if (mesh.steady())
        {
            adjustMass = adjustPhi(phiHbyA, U, p);
        }

        if (pimple.consistent())
        {
            phiHbyA += (rhorAAtUf - rhorAUf)*fvc::snGrad(p)*mesh.magSf();
            HbyA += (rAAtU - rAU)*fvc::grad(p);
        }

        fvScalarMatrix pDDtEqn
        (
            fvc::ddt(rho) + psi*correction(fvm::ddt(p))
        + fvc::div(phiHbyA)
        ==
            fvOptions(psi, p, rho.name())
        );

        while (pimple.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn(pDDtEqn - fvm::laplacian(rhorAAtUf, p));

            pEqn.setReference
            (
                pressureControl.refCell(),
                pressureControl.refValue()
            );

            pEqn.solve();

            if (pimple.finalNonOrthogonalIter())
            {
                phi = phiHbyA + pEqn.flux();
            }
        }
    }

    if (mesh.steady())
    {
        incompressibleContinuityErrs();
    }
    else
    {
        const bool limitedp = pressureControl.limit(p);

        // Thermodynamic density update
        thermo.correctRho(psi*p - psip0);

        if (limitedp)
        {
            rho = thermo.rho();
        }

        rhoEqn_();
        compressibleContinuityErrs();
    }

    // Explicitly relax pressure for momentum corrector
    p.relax();

    U = HbyA - rAAtU*fvc::grad(p);
    U.correctBoundaryConditions();
    fvOptions.correct(U);
    K = 0.5*magSqr(U);

    if (mesh.steady())
    {
        pressureControl.limit(p);
    }

    // For steady compressible closed-volume cases adjust the pressure level
    // to obey overall mass continuity
    if (adjustMass && !thermo.incompressible())
    {
        p += (initialMass - fvc::domainIntegrate(thermo.rho()))
            /fvc::domainIntegrate(psi);
        p.correctBoundaryConditions();
    }

    if (mesh.steady() || pimple.simpleRho() || adjustMass)
    {
        rho = thermo.rho();
    }

    // Correct rhoUf if the mesh is moving
    fvc::correctRhoUf(rhoUf, rho, U, phi);

    if ((mesh.steady() || pimple.simpleRho()) && !pimple.transonic())
    {
        rho.relax();
    }

    if (thermo.dpdt())
    {
        dpdt = fvc::ddt(p);

        if (mesh.moving())
        {
            dpdt -= fvc::div(fvc::meshPhi(rho, U), p);
        }
    }
    //-------------------------------------------
#endif

    return 0;
}

rhoPimple::~rhoPimple()
{
   finalizeFoam();

}

int rhoPimple::finalizeFoam()
{
    //delete argsPtr;
    //delete runTimePtr;
    //delete pPtr;
    //delete TPtr;
    //delete psiPtr;
    //delete ePtr;

    if (rhoPtr != nullptr)
    {
        delete rhoPtr;
        rhoPtr = nullptr;
    }
    
    if (UPtr != nullptr)
    {
        delete UPtr;
        UPtr = nullptr;
    }

    if (rhoUPtr != nullptr)
    {
        delete rhoUPtr;
        rhoUPtr = nullptr;
    }

    if (rhoEPtr != nullptr)
    {
        delete rhoEPtr;
        rhoEPtr = nullptr;
    }

    if (phiPtr != nullptr)
    {
        delete phiPtr;
        phiPtr = nullptr;
    }

    if (pimplePtr != nullptr)
    {
        delete pimplePtr;
        pimplePtr = nullptr;
    }

    if (pressureControlPtr != nullptr)
    {
        delete pressureControlPtr;
        pressureControlPtr = nullptr;
    }

    if (dpdtPtr != nullptr)
    {
        delete dpdtPtr;
        dpdtPtr = nullptr;
    }

    if (KPtr != nullptr) 
    {
        delete KPtr;
        KPtr = nullptr;
    }

    if (fvOptionsPtr != nullptr) 
    {
        delete fvOptionsPtr;
        fvOptionsPtr = nullptr;
    }

    if (MRFPtr != nullptr)
    {
        delete MRFPtr;
        MRFPtr = nullptr;
    }

#ifdef HAVE_OF8
    if (initialMassPtr != nullptr)
    {
        delete initialMassPtr;
        initialMassPtr = nullptr;
    }
#endif

#ifdef HAVE_OFE20
    if (rhoMaxPtr != nullptr)
    {
        delete rhoMaxPtr;
        rhoMaxPtr = nullptr;
    }
    if (rhoMinPtr != nullptr)
    {
        delete rhoMinPtr;
        rhoMinPtr = nullptr;
    }
#endif

    //meshPtr.clear();
    turbulencePtr.clear();
#ifdef HAVE_OF8
    thermophysicalTransportPtr.clear();
#endif
    trDeltaT.clear();
    rhoUfPtr.clear();
    //UEqnPtr.clear();
    pThermoPtr.clear();
    divrhoUPtr.clear();
    tUEqnPtr.clear();

    finalizeStat = 0;
    return finalizeStat;
}

//^^^^^ (UN)LOAD METHOD ^^^^^^^^^^^^^^^^^^^^^^^^^
// C/C++ bindings to load rocFoam
extern "C" void rocrhopimple_load_module(const char *name)
{
    rhoPimple::load(name);
}

extern "C" void rocrhopimple_unload_module(const char *name)
{
    rhoPimple::unload(name);
}
//===============================================

double rhoPimple::errorEvaluate(int argc, char *argv[])
{
    createArgs(argc, argv);
    setRootCase();
    
    Foam::argList &args(*argsPtr);
    
    //Commenting the code bellow
    //Not sure what happens
    /*
    if (args.optionFound("list"))
    {
        functionObjectList::list();
        return 0;
    }
    */
    
    createTime();
    createDynamicFvMesh();

    Foam::Time &runTime(*runTimePtr);
    
    Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);
    
    std::vector<volScalarField> rhoVec;
    std::vector<volVectorField> rhoUVec;
    std::vector<volScalarField> UMagVec;
    std::vector<volScalarField> rhoEVec;

    forAll(timeDirs, timei)
    {
        runTime.setTime(timeDirs[timei], timei);

        Info << "iTime = " << timei <<", " << "Time = " << runTime.timeName() << endl;

        FatalIOError.throwExceptions();

        createDyMControls();
        initContinuityErrs();
        createFields();

        volScalarField &rho(*rhoPtr);
        volVectorField &rhoU(*UPtr);
        volScalarField &rhoE(*KPtr);
        
        rhoVec.push_back(rho);
        rhoUVec.push_back(rhoU);
        UMagVec.emplace_back(magSqr(rhoU));
        rhoEVec.push_back(rhoE);

        finalizeFoam();

    }

    Info << "Field vectors created. " << endl;

    double rhoMax{-1.0};
    forAll(rhoVec[0], i)
    {
        rhoVec[0][i]  = std::abs(rhoVec[2][i] - rhoVec[1][i]);

        rhoMax = std::max(rhoVec[0][i], rhoMax);
    }

    double UMagVecMax{-1.0};
    forAll(UMagVec[0], i)
    {
        UMagVec[0][i] = std::abs(UMagVec[2][i] - UMagVec[1][i]);
        
        UMagVecMax = std::max(UMagVec[0][i], UMagVecMax);
    }

    double rhoEVecMax{-1.0};
    forAll(rhoEVec[0], i)
    {
        rhoEVec[0][i] = std::abs(rhoEVec[2][i] - rhoEVec[1][i]);

        rhoEVecMax = std::max(rhoEVec[0][i], rhoEVecMax);
    }

    Info << "Infinity norm(rho)  = " << rhoMax << endl;
    Info << "Infinity norm(UMag) = " << UMagVecMax << endl;
    Info << "Infinity norm(rhoE) = " << rhoEVecMax << endl;

    double maxError = std::max( rhoMax, UMagVecMax );
    testStat = std::max( maxError , rhoEVecMax );
    Info << "Max Infinity norm   = " << testStat << endl;
    
    return testStat;
}
//=========================================================

