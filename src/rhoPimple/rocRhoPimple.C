#include "rocRhoPimple.H"

rocRhoPimple::rocRhoPimple()
    : rocFoam(),
      pimplePtr(NULL),
      pressureControlPtr(NULL),
      dpdtPtr(NULL),
      KPtr(NULL),
      fvOptionsPtr(NULL),
      MRFPtr(NULL),
      UEqnPtr(NULL),
      pThermoPtr(NULL),
      rhoUfPtr(NULL),
      divrhoUPtr(NULL),
      tUEqnPtr(NULL),
      correctPhi(false),
      checkMeshCourantNo(false),
      moveMeshOuterCorrectors(false),
      cumulativeContErr(0.0)
{}

rocRhoPimple::rocRhoPimple(int argc, char *argv[])
    : rocFoam(),
      pimplePtr(NULL),
      pressureControlPtr(NULL),
      dpdtPtr(NULL),
      KPtr(NULL),
      fvOptionsPtr(NULL),
      MRFPtr(NULL),
      UEqnPtr(NULL),
      pThermoPtr(NULL),
      rhoUfPtr(NULL),
      divrhoUPtr(NULL),
      tUEqnPtr(NULL),
      correctPhi(false),
      checkMeshCourantNo(false),
      moveMeshOuterCorrectors(false),
      cumulativeContErr(0.0)
{
    initialize(argc, argv);
}

int rocRhoPimple::initialize(int argc, char *argv[])
{
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

    //  createDyMControls.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    createDyMControls();
    // ---------------------------------------------------

    //  initContinuityErrs.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    initContinuityErrs();
    // ---------------------------------------------------

    //  createFields.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    createFields();
    // ---------------------------------------------------

    //  createFieldRefs.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    createFieldRefs();
    // ---------------------------------------------------

    //  createRhoUfIfPresent.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^
    createRhoUfIfPresent();
    // ---------------------------------------------------

    compressible::turbulenceModel &turbulence(*turbulencePtr);

    turbulence.validate();

    if (!LTS)
    {
        //  compressibleCourantNo.H  ^^^^^^^^^^^^^^^^^^^^^^^^^
        compressibleCourantNo();
        // ---------------------------------------------------

        //  setInitialDeltaT.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        setInitialDeltaT();
        // ---------------------------------------------------
    }

    Foam::Info << "End of initialization." << endl;

    return 0;
}

int rocRhoPimple::createControl()
{
    dynamicFvMesh &mesh(*meshPtr);

    //  createPimpleControl.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^
    // pimpleControl pimple(mesh);
    pimplePtr = new pimpleControl(mesh);
    // ---------------------------------------------------

    return 0;
}

int rocRhoPimple::createDyMControls()
{
    dynamicFvMesh &mesh(*meshPtr);

    //  createControl.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    createControl();
    // ---------------------------------------------------
    pimpleControl &pimple(*pimplePtr);

    //  createTimeControls.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    createTimeControls();
    // ---------------------------------------------------

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

    return 0;
}

int rocRhoPimple::initContinuityErrs()
{

#ifndef initContinuityErrs_H
#define initContinuityErrs_H

    cumulativeContErr = 0;

#endif

    return 0;
}

int rocRhoPimple::createFields()
{
    Foam::Time &runTime(*runTimePtr);
    Foam::argList &args(*argsPtr);
    dynamicFvMesh &mesh(*meshPtr);
    pimpleControl &pimple(*pimplePtr);

    //  createRDeltaT.H  ^^^^^^^^^^^^^
    createRDeltaT();
    // -------------------------------

    Info << "Reading thermophysical properties\n" << endl;

    pThermoPtr = autoPtr<fluidThermo>
    (
        fluidThermo::New(mesh)
    );
    fluidThermo &thermo(*pThermoPtr);

    thermo.validate(args.executable(), "h", "e");

    pPtr = &thermo.p();
    volScalarField &p(*pPtr);

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

    pressureControlPtr = new pressureControl
    (
        p,
        rho,
        pimple.dict(),
        false
    );

    mesh.setFluxRequired(p.name());

    Info << "Creating turbulence model\n" << endl;
    turbulencePtr = autoPtr<compressible::turbulenceModel>
    (
        compressible::turbulenceModel::New
        (
            rho,
            U,
            phi,thermo
        )
    );

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

    //  createMRF.H  ^^^^^^^^^^^^^^^^
    createMRF();
    // -------------------------------

    //  createFvOptions.H  ^^^^^^^^^^^
    createFvOptions();
    // -------------------------------

    return 0;
}

int rocRhoPimple::createMRF()
{
    dynamicFvMesh &mesh(*meshPtr);

    MRFPtr = new IOMRFZoneList(mesh);

    return 0;
}

int rocRhoPimple::compressibleCreatePhi()
{
    Foam::Time &runTime(*runTimePtr);
    dynamicFvMesh &mesh(*meshPtr);
    volScalarField &rho(*rhoPtr);
    volVectorField &U(*UPtr);

    Info << "Reading/calculating face flux field phi\n" << endl;

    phiPtr = new surfaceScalarField
    (
        IOobject("phi", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
        linearInterpolate(rho * U) & mesh.Sf()
    );

    return 0;
}

int rocRhoPimple::createFvOptions()
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

int rocRhoPimple::createFieldRefs()
{
    fluidThermo &thermo(*pThermoPtr);

    psiPtr = &thermo.psi();

    return 0;
}

int rocRhoPimple::createRhoUfIfPresent()
{
    Foam::Time &runTime(*runTimePtr);
    dynamicFvMesh &mesh(*meshPtr);
    volScalarField &rho(*rhoPtr);
    volVectorField &U(*UPtr);

    if (mesh.dynamic())
    {
        Info << "Constructing face momentum rhoUf" << endl;

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
    }

    return 0;
}

int rocRhoPimple::compressibleCourantNo()
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

int rocRhoPimple::setInitialDeltaT()
{
    Foam::Time &runTime(*runTimePtr);

    if (adjustTimeStep)
    {
        if ((runTime.timeIndex() == 0) && (CoNum > small))
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

int rocRhoPimple::loop()
{
    Foam::Time &runTime(*runTimePtr);
    dynamicFvMesh &mesh(*meshPtr);
    pimpleControl &pimple(*pimplePtr);
    volScalarField &rho(*rhoPtr);
    volVectorField &U(*UPtr);
    surfaceScalarField &phi(*phiPtr);
    fluidThermo &thermo(*pThermoPtr);
    IOMRFZoneList &MRF(*MRFPtr);
    compressible::turbulenceModel &turbulence(*turbulencePtr);
    autoPtr<surfaceVectorField> &rhoUf(rhoUfPtr);

    Info << "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        //  readDyMControls.H  ^^^^^^^^^^^
        readDyMControls();
        // -------------------------------

        // Store divrhoU from the previous mesh so that it can be mapped
        // and used in correctPhi to ensure the corrected phi has the
        // same divergence
        if (correctPhi)
        {
            divrhoUPtr = new volScalarField
            (
                "divrhoU",
                fvc::div(fvc::absolute(phi, rho, U))
            );

        }

        if (LTS)
        {
            //  setRDeltaT.H  ^^^^^^^^^^^^^^^
            setRDeltaT();
            // -------------------------------
        }
        else
        {
            //  compressibleCourantNo.H  ^^^^^^^^^^^^^^^^^^^
            compressibleCourantNo();
            // ---------------------------------------------

            //  setDeltaT.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            setDeltaT();
            // ---------------------------------------------
        }

        runTime++;

        Info << "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstPimpleIter() || moveMeshOuterCorrectors)
            {
                // Store momentum to set rhoUf for introduced faces.
                autoPtr<volVectorField> rhoU;
                if (rhoUf.valid())
                {
                    rhoU = new volVectorField("rhoU", rho * U);
                }

                // Do any mesh changes
                mesh.update();

                if (mesh.changing())
                {
                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & rhoUf();

                        //  correctPhi.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                        correctPhi_();
                        // ---------------------------------------------

                        // Make the fluxes relative to the mesh-motion
                        fvc::makeRelative(phi, rho, U);
                    }

                    if (checkMeshCourantNo)
                    {
                        //  meshCourantNo.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^
                        meshCourantNo();
                        // ---------------------------------------------
                    }
                }
            }

            if (pimple.firstPimpleIter() && !pimple.simpleRho())
            {
                //  rhoEqn.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                rhoEqn_();
                // ---------------------------------------------
            }

            //  UEqn.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            UEqn_();
            // -------------------------------------------

            //  EEqn.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            EEqn_();
            // --------------------------------------------

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                if (pimple.consistent())
                {
                    //  pcEqn.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                    pcEqn();
                    // --------------------------------------------
                }
                else
                {
                    //  pEqn.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                    pEqn_();
                    // --------------------------------------------
                }
            }

            if (pimple.turbCorr())
            {
                turbulence.correct();
            }
        }

        rho = thermo.rho();

        runTime.write();

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl
             << endl;
    }

    Info << "End\n" << endl;

    return 0;
}

int rocRhoPimple::readDyMControls()
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

int rocRhoPimple::setRDeltaT()
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

    scalar maxCo(pimpleDict.lookupOrDefault<scalar>("maxCo", 0.8));

    scalar rDeltaTSmoothingCoeff
    (
        pimpleDict.lookupOrDefault<scalar>("rDeltaTSmoothingCoeff", 0.02)
    );

    scalar rDeltaTDampingCoeff
    (
        pimpleDict.lookupOrDefault<scalar>("rDeltaTDampingCoeff", 1.0)
    );

    scalar maxDeltaT
    (
        pimpleDict.lookupOrDefault<scalar>("maxDeltaT", great)
    );

    volScalarField rDeltaT0("rDeltaT0", rDeltaT);

    // Set the reciprocal time-step from the local Courant number
    rDeltaT.ref() =
        max(1 / dimensionedScalar(dimTime, maxDeltaT),
            fvc::surfaceSum(mag(phi))()() / ((2 * maxCo) * mesh.V() * rho()));

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

    // Update tho boundary values of the reciprocal time-step
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

int rocRhoPimple::correctPhi_()
{
    pimpleControl &pimple(*pimplePtr);
    volScalarField &rho(*rhoPtr);
    volVectorField &U(*UPtr);
    surfaceScalarField &phi(*phiPtr);
    const volScalarField &psi(*psiPtr);
    volScalarField &divrhoU(*divrhoUPtr);
    volScalarField &p(*pPtr);

    CorrectPhi
    (
        U, phi, p, rho, psi,
        dimensionedScalar("rAUf", dimTime, 1),
        divrhoU(), pimple, true
    );

    return 0;
}

int rocRhoPimple::meshCourantNo()
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

int rocRhoPimple::rhoEqn_()
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

int rocRhoPimple::UEqn_()
{
    pimpleControl &pimple(*pimplePtr);
    volScalarField &rho(*rhoPtr);
    volVectorField &U(*UPtr);
    surfaceScalarField &phi(*phiPtr);
    Foam::fv::options &fvOptions(*fvOptionsPtr);
    IOMRFZoneList &MRF(*MRFPtr);
    compressible::turbulenceModel &turbulence(*turbulencePtr);
    volScalarField &K(*KPtr);
    volScalarField &p(*pPtr);

    MRF.correctBoundaryVelocity(U);

    tUEqnPtr = tmp<fvVectorMatrix>
    (
        fvm::ddt(rho, U) + fvm::div(phi, U)
      + MRF.DDt(rho, U)
      + turbulence.divDevRhoReff(U)
     == fvOptions(rho, U)
    );


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

int rocRhoPimple::EEqn_()
{
    volScalarField &rho(*rhoPtr);
    volVectorField &U(*UPtr);
    surfaceScalarField &phi(*phiPtr);
    Foam::fv::options &fvOptions(*fvOptionsPtr);
    fluidThermo &thermo(*pThermoPtr);
    compressible::turbulenceModel &turbulence(*turbulencePtr);
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
      - fvm::laplacian(turbulence.alphaEff(), he)
     == fvOptions(rho, he));

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();

    return 0;
}

int rocRhoPimple::pcEqn()
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

int rocRhoPimple::compressibleContinuityErrs()
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

int rocRhoPimple::pEqn_()
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

    //  rhoEqn.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    rhoEqn_();
    // ---------------------------------------------

    //  compressibleContinuityErrs.H  ^^^^^^^^^^^^^^
    compressibleContinuityErrs();
    // ---------------------------------------------

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

    return 0;
}

int rocRhoPimple::readTimeControls()
{
    Foam::Time &runTime(*runTimePtr);

    adjustTimeStep =
        runTime.controlDict().lookupOrDefault("adjustTimeStep", false);

    maxCo = runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);

    maxDeltaT =
        runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", great);

    return 0;
}

int rocRhoPimple::finalize()
{
    delete argsPtr;
    delete runTimePtr;
    delete pPtr;
    delete TPtr;
    delete psiPtr;
    delete ePtr;
    delete rhoPtr;
    delete UPtr;
    delete rhoUPtr;
    delete rhoEPtr;
    delete phiPtr;

    // delete meshPtr;
    // delete turbulencePtr;
    // delete trDeltaT;

    delete pimplePtr;
    delete pressureControlPtr;
    delete dpdtPtr;
    delete KPtr;
    delete fvOptionsPtr;
    delete MRFPtr;
    delete UEqnPtr;

    // delete pThermoPtr;
    // delete rhoUfPtr;
    // delete divrhoUPtr;
    // delete tUEqnPtr;

    return 0;
}
