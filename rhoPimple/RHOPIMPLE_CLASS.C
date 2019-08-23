#include "RHOPIMPLE_CLASS.H"

RHOPIMPLE_CLASS::RHOPIMPLE_CLASS() : FOAM_CLASS()
{}

RHOPIMPLE_CLASS::RHOPIMPLE_CLASS(int argc,char *argv[]) : FOAM_CLASS()
{
   Initialize(argc, argv);
}


int RHOPIMPLE_CLASS::Initialize(int argc,char *argv[]){

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

    Foam::Info << "End of initialization of openFoamPar module." << endl;

    return 0;
}


int RHOPIMPLE_CLASS::createControl()
{

    dynamicFvMesh &mesh(*meshPtr);


    //  createPimpleControl.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^
    pimpleControl pimple(mesh);
    // ---------------------------------------------------

    return 0;
}



int RHOPIMPLE_CLASS::createDyMControls()
{

    dynamicFvMesh &mesh(*meshPtr);


    //  createControl.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    createControl();
    // ---------------------------------------------------

    //  createTimeControls.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    createTimeControls();
    // ---------------------------------------------------

    correctPhi = 
        pimple.dict().lookupOrDefault
        (
            "correctPhi", mesh.dynamic()
        );

    checkMeshCourantNo =
        pimple.dict().lookupOrDefault
        (
            "checkMeshCourantNo", false
        );

    moveMeshOuterCorrectors =
        pimple.dict().lookupOrDefault
        (
            "moveMeshOuterCorrectors", false
        );

    return 0;
}







int RHOPIMPLE_CLASS::initContinuityErrs()
{
    #ifndef initContinuityErrs_H
    #define initContinuityErrs_H

    cumulativeContErr = 0;

    #endif

    return 0;
}

int RHOPIMPLE_CLASS::createFields()
{

    dynamicFvMesh &mesh(*meshPtr);


   //  createRDeltaT.H  ^^^^^^^^^^^^^
   createRDeltaT();
   // -------------------------------
   
    Info<< "Reading thermophysical properties\n" << endl;

    autoPtr<fluidThermo> pThermo
    (
        fluidThermo::New(mesh)
    );
    fluidThermo& thermo = pThermo();
    thermo.validate(args.executable(), "h", "e");

    volScalarField& p = thermo.p();

    volScalarField rho
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

    Info<< "Reading field U\n" << endl;
    volVectorField U
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

    //  compressibleCreatePhi.H  ^^^^^
    compressibleCreatePhi();
    // -------------------------------

    pressureControl pressureControl(p, rho, pimple.dict(), false);

    mesh.setFluxRequired(p.name());

    Info<< "Creating turbulence model\n" << endl;
    autoPtr<compressible::turbulenceModel> turbulence
    (
        compressible::turbulenceModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );

    Info<< "Creating field dpdt\n" << endl;
    volScalarField dpdt
    (
        IOobject
        (
            "dpdt",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(p.dimensions()/dimTime, 0)
    );

    Info<< "Creating field kinetic energy K\n" << endl;
    volScalarField K("K", 0.5*magSqr(U));

   //  createMRF.H  ^^^^^^^^^^^^^^^^
   createMRF();
   // -------------------------------

   //  createFvOptions.H  ^^^^^^^^^^^
   createFvOptions();
   // -------------------------------

    return 0;
}

int RHOPIMPLE_CLASS::createMRF()
{
    IOMRFZoneList MRF(mesh);

    return 0;
}


int RHOPIMPLE_CLASS:compressibleCreatePhi()
{

    dynamicFvMesh &mesh(*meshPtr);



    Info<< "Reading/calculating face flux field phi\n" << endl;

    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(rho*U) & mesh.Sf()
    );

    return 0;
}


int RHOPIMPLE_CLASS::createFvOptions()
{
    fv::options& fvOptions(fv::options::New(mesh));

    if (!fvOptions.optionList::size())
    {
        Info << "No finite volume options present" << endl;
    }

    return 0;
}



int RHOPIMPLE_CLASS::createFieldRefs()
{

    const volScalarField& psi = thermo.psi();

    return 0;
}


int RHOPIMPLE_CLASS::createRhoUfIfPresent()
{

    dynamicFvMesh &mesh(*meshPtr);



    autoPtr<surfaceVectorField> rhoUf;

    if (mesh.dynamic())
    {
        Info<< "Constructing face momentum rhoUf" << endl;

        rhoUf = new surfaceVectorField
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
        );
    }

    return 0;
}


int RHOPIMPLE_CLASS::compressibleCourantNo()
{

    dynamicFvMesh &mesh(*meshPtr);


    CoNum = 0.0;
    meanCoNum = 0.0;

    {
        scalarField sumPhi
        (
            fvc::surfaceSum(mag(phi))().primitiveField()/rho.primitiveField()
        );

        CoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

        meanCoNum =
            0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();
    }

    Info<< "Courant Number mean: " << meanCoNum
        << " max: " << CoNum << endl;

    return 0;
}

int RHOPIMPLE_CLASS::setInitialDeltaT()
{

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


int RHOPIMPLE_CLASS::loop()
{

    dynamicFvMesh &mesh(*meshPtr);


    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {

        //  readDyMControls.H  ^^^^^^^^^^^
        readDyMControls();
        // -------------------------------

        // Store divrhoU from the previous mesh so that it can be mapped
        // and used in correctPhi to ensure the corrected phi has the
        // same divergence
        autoPtr<volScalarField> divrhoU;
        if (correctPhi)
        {
            divrhoU = new volScalarField
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

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstPimpleIter() || moveMeshOuterCorrectors)
            {
                // Store momentum to set rhoUf for introduced faces.
                autoPtr<volVectorField> rhoU;
                if (rhoUf.valid())
                {
                    rhoU = new volVectorField("rhoU", rho*U);
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
                    pEqn();
                    // --------------------------------------------
                }
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        rho = thermo.rho();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;



    return 0;
}

int RHOPIMPLE_CLASS::readDyMControls()
{

    //  readTimeControls.H  ^^^^^^^^^^^^^^^^^^^^^^^
    readTimeControls();
    // --------------------------------------------

    correctPhi = pimple.dict().lookupOrDefault
    (
        "correctPhi",
        correctPhi
    );

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


int RHOPIMPLE_CLASS::setRDeltaT()
{

    dynamicFvMesh &mesh(*meshPtr);


    volScalarField& rDeltaT = trDeltaT.ref();

    const dictionary& pimpleDict = pimple.dict();

    scalar maxCo
    (
        pimpleDict.lookupOrDefault<scalar>("maxCo", 0.8)
    );

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
    rDeltaT.ref() = max
    (
        1/dimensionedScalar(dimTime, maxDeltaT),
        fvc::surfaceSum(mag(phi))()()
       /((2*maxCo)*mesh.V()*rho())
    );

    if (pimple.transonic())
    {
        surfaceScalarField phid
        (
            "phid",
            fvc::interpolate(psi)*fvc::flux(U)
        );

        rDeltaT.ref() = max
        (
            rDeltaT(),
            fvc::surfaceSum(mag(phid))()()
            /((2*maxCo)*mesh.V()*psi())
        );
    }

    // Update tho boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    Info<< "Flow time scale min/max = "
        << gMin(1/rDeltaT.primitiveField())
        << ", " << gMax(1/rDeltaT.primitiveField()) << endl;

    if (rDeltaTSmoothingCoeff < 1.0)
    {
        fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);
    }

    Info<< "Smoothed flow time scale min/max = "
        << gMin(1/rDeltaT.primitiveField())
        << ", " << gMax(1/rDeltaT.primitiveField()) << endl;

    // Limit rate of change of time scale
    // - reduce as much as required
    // - only increase at a fraction of old time scale
    if
    (
        rDeltaTDampingCoeff < 1.0
     && runTime.timeIndex() > runTime.startTimeIndex() + 1
    )
    {
        rDeltaT =
            rDeltaT0
           *max(rDeltaT/rDeltaT0, scalar(1) - rDeltaTDampingCoeff);

        Info<< "Damped flow time scale min/max = "
            << gMin(1/rDeltaT.primitiveField())
            << ", " << gMax(1/rDeltaT.primitiveField()) << endl;
    }


    return 0;
}



int RHOPIMPLE_CLASS::correctPhi_()
{

    CorrectPhi
    (
        U,
        phi,
        p,
        rho,
        psi,
        dimensionedScalar("rAUf", dimTime, 1),
        divrhoU(),
        pimple,
        true
    );

    return 0;
}


int RHOPIMPLE_CLASS::meshCourantNo()
{

    dynamicFvMesh &mesh(*meshPtr);


    scalar meshCoNum = 0.0;
    scalar meanMeshCoNum = 0.0;

    {
        scalarField sumPhi
        (
            fvc::surfaceSum(mag(mesh.phi()))().primitiveField()
        );

        meshCoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

        meanMeshCoNum =
            0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();
    }

    Info<< "Mesh Courant Number mean: " << meanMeshCoNum
        << " max: " << meshCoNum << endl;

    return 0;
}

int RHOPIMPLE_CLASS::rhoEqn_()
{
    fvScalarMatrix rhoEqn
    (
        fvm::ddt(rho)
      + fvc::div(phi)
      ==
        fvOptions(rho)
    );

    fvOptions.constrain(rhoEqn);

    rhoEqn.solve();

    fvOptions.correct(rho);

    return 0;
}

int RHOPIMPLE_CLASS::UEqn_()
{

    MRF.correctBoundaryVelocity(U);

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::ddt(rho, U) + fvm::div(phi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(U)
     ==
        fvOptions(rho, U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));

        fvOptions.correct(U);
        K = 0.5*magSqr(U);
    }

    return 0;
}

int RHOPIMPLE_CLASS::EEqn_()
{

    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + fvm::div(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
          : -dpdt
        )
      - fvm::laplacian(turbulence->alphaEff(), he)
     ==
        fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();

    return 0;
}

int RHOPIMPLE_CLASS::pcEqn()
{

    dynamicFvMesh &mesh(*meshPtr);


    if (!pimple.simpleRho())
    {
        rho = thermo.rho();
    }

    // Thermodynamic density needs to be updated by psi*d(p) after the
    // pressure solution
    const volScalarField psip0(psi*p);

    volScalarField rAU(1.0/UEqn.A());
    volScalarField rAtU(1.0/(1.0/rAU - UEqn.H1()));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

    if (pimple.nCorrPiso() <= 1)
    {
        tUEqn.clear();
    }

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
            fvc::interpolate(rho)*fvc::flux(HbyA)
          + MRF.zeroFilter
            (
                fvc::interpolate(rho*rAU)*fvc::ddtCorr(rho, U, phi, rhoUf)
            )
        )
    );

    fvc::makeRelative(phiHbyA, rho, U);
    MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

    volScalarField rhorAtU("rhorAtU", rho*rAtU);

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, rho, U, phiHbyA, rhorAtU, MRF);

    if (pimple.transonic())
    {
        surfaceScalarField phid
        (
            "phid",
            (fvc::interpolate(psi)/fvc::interpolate(rho))*phiHbyA
        );

        phiHbyA +=
            fvc::interpolate(rho*(rAtU - rAU))*fvc::snGrad(p)*mesh.magSf()
          - fvc::interpolate(psi*p)*phiHbyA/fvc::interpolate(rho);

        HbyA -= (rAU - rAtU)*fvc::grad(p);

        fvScalarMatrix pDDtEqn
        (
            fvc::ddt(rho) + psi*correction(fvm::ddt(p))
          + fvc::div(phiHbyA) + fvm::div(phid, p)
         ==
            fvOptions(psi, p, rho.name())
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
        phiHbyA += fvc::interpolate(rho*(rAtU - rAU))*fvc::snGrad(p)*mesh.magSf();
        HbyA -= (rAU - rAtU)*fvc::grad(p);

        fvScalarMatrix pDDtEqn
        (
            fvc::ddt(rho) + psi*correction(fvm::ddt(p))
          + fvc::div(phiHbyA)
         ==
            fvOptions(psi, p, rho.name())
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
    thermo.correctRho(psi*p - psip0);

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

    U = HbyA - rAtU*fvc::grad(p);
    U.correctBoundaryConditions();
    fvOptions.correct(U);
    K = 0.5*magSqr(U);

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


int RHOPIMPLE_CLASS::compressibleContinuityErrs()
{

    dimensionedScalar totalMass = fvc::domainIntegrate(rho);

    scalar sumLocalContErr =
        (fvc::domainIntegrate(mag(rho - thermo.rho()))/totalMass).value();

    scalar globalContErr =
        (fvc::domainIntegrate(rho - thermo.rho())/totalMass).value();

    cumulativeContErr += globalContErr;

    Info<< "time step continuity errors : sum local = " << sumLocalContErr
        << ", global = " << globalContErr
        << ", cumulative = " << cumulativeContErr
        << endl;


    return 0;
}

int pEqn_()
{

    dynamicFvMesh &mesh(*meshPtr);


    if (!pimple.simpleRho())
    {
        rho = thermo.rho();
    }

    // Thermodynamic density needs to be updated by psi*d(p) after the
    // pressure solution
    const volScalarField psip0(psi*p);

    volScalarField rAU(1.0/UEqn.A());
    surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
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

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, rho, U, phiHbyA, rhorAUf, MRF);

    if (pimple.transonic())
    {
        surfaceScalarField phid
        (
            "phid",
            (fvc::interpolate(psi)/fvc::interpolate(rho))*phiHbyA
        );

        phiHbyA -= fvc::interpolate(psi*p)*phiHbyA/fvc::interpolate(rho);

        fvScalarMatrix pDDtEqn
        (
            fvc::ddt(rho) + psi*correction(fvm::ddt(p))
          + fvc::div(phiHbyA) + fvm::div(phid, p)
         ==
            fvOptions(psi, p, rho.name())
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
            fvc::ddt(rho) + psi*correction(fvm::ddt(p))
          + fvc::div(phiHbyA)
         ==
            fvOptions(psi, p, rho.name())
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
    thermo.correctRho(psi*p - psip0);

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

    U = HbyA - rAU*fvc::grad(p);
    U.correctBoundaryConditions();
    fvOptions.correct(U);
    K = 0.5*magSqr(U);

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

int RHOPIMPLE_CLASS::readTimeControls()
{

    adjustTimeStep =
        runTime.controlDict().lookupOrDefault("adjustTimeStep", false);

    maxCo =
        runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);

    maxDeltaT =
        runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", FOAM_CLASS::great);

    return 0;
}

