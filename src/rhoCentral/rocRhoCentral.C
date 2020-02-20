#include "rocRhoCentral.H"

//^^^ DEFINITION OF CONSTRUCTORS ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
rhoCentral::rhoCentral()
    : posPtr(NULL),
      negPtr(NULL),
      amaxSfPtr(NULL),
      pThermoPtr(NULL),
      fluxScheme(""),
      inviscid(false)
{
    solverType = const_cast<char *>("rocRhoCentral");
};

rhoCentral::rhoCentral(int argc, char *argv[])
    : posPtr(NULL),
      negPtr(NULL),
      amaxSfPtr(NULL),
      pThermoPtr(NULL),
      fluxScheme(""),
      inviscid(false)
{
    solverType = const_cast<char *>("rocRhoCentral");
    initialize(argc, argv);
}
//===================================================================


//^^^ DEFINITION OF COM-RELATED MTHODS ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//^^^^^ LOAD MODULES ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
void rhoCentral::load(const char *name)
{
    //  Anouncing default communicator  ^^^^^^^^^^^^^^^^^^^
    MPI_Comm tmpComm;
    tmpComm = COM_get_default_communicator();  

    int tmpRank, tmpNProc;
    MPI_Comm_rank(tmpComm, &tmpRank);
    MPI_Comm_size(tmpComm, &tmpNProc);
    
    if (tmpRank == 0)
    {
        std::cout << "rocFoam.load: Loading rhoCentral with name "
                   << name << "." << std::endl;

        std::cout << "rocFoam.load: Rank = " << tmpRank
                  << ", NProc = " << tmpNProc
                  << ", COMM = " << tmpComm << std::endl;

        std::cout << std::endl;
    }


    std::string volName = name+string("VOL");
    std::string srfName = name+string("SRF");


    //  Register module with COM ^^^^^^^^^^^^^^^^^^^^^^^^^^
    rhoCentral *comFoamPtr = new rhoCentral();

    //COM_new_window(name, MPI_COMM_NULL);
    COM_new_window(volName, tmpComm);
    //COM_new_window(srfName, tmpComm);

    comFoamPtr->winVolName = volName;
    comFoamPtr->winSrfName = srfName;

    //MPI_Comm_dup(tmpComm, &(comFoamPtr->winComm));
    comFoamPtr->winComm = tmpComm;
    
    //Foam::PstreamGlobals::MPI_comFoam_to_openFoam = comFoamPtr->winComm;
    Foam::PstreamGlobals::MPI_COMM_FOAM = comFoamPtr->winComm;
    
    MPI_Comm_rank(comFoamPtr->winComm, &(comFoamPtr->winRank));
    MPI_Comm_size(comFoamPtr->winComm, &(comFoamPtr->winNProc));

    std::string objectName = volName + string(".object");

    COM_new_dataitem(objectName.c_str(), 'w', COM_VOID, 1, "");

    COM_set_object(objectName.c_str(), 0, comFoamPtr);

    COM_window_init_done(volName);
    //COM_window_init_done(srfName);  

    comFoamPtr->registerFunctions(volName.c_str());

    return;
}
//---------------------------------------------------------


//^^^^^ UNLOAD MODULES ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
void rhoCentral::unload(const char *name)
{
    Foam::Info << "rocFoam.unload: Unloading rocRhoCentral with name "
               << name << "." << Foam::endl;

    comFoam *comFoamPtr = NULL;

    std::string volName = name+string("VOL");
    std::string srfName = name+string("SRF");

    std::string objectName(volName+".object");

    COM_get_object(objectName.c_str(), 0, &comFoamPtr);

    //comFoamPtr->finalize();
    delete comFoamPtr;

    COM_delete_window(std::string(volName));
    //COM_delete_window(std::string(srfName));
}
//---------------------------------------------------------

//===================================================================


//^^^ DEFINITION OF OPENFOAM-RELATED MTHODS ^^^^^^^^^^^^^^^^^^^^^^^^^
int rhoCentral::initialize(int argc, char *argv[])
{
#define NO_CONTROL

    // Not quite sure where this line should be
    createArgs(argc, argv);

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

    Foam::Info << "End of initialization of rhoCentral." << Foam::endl;

    // Setting comFoam variables that will be registered ^^
    //   with COM. This are needed for flow control when
    //   solver runs step-by-step.
    Foam::Time &runTime(*runTimePtr);
    winTime = runTime.value();
    winDeltaT = runTime.deltaTValue();
    winRun = static_cast<int>(runTime.run());
    //-----------------------------------------------------

    initializeStat = 0;
    return initializeStat;
}


int rhoCentral::loop()
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
    Foam::psiThermo &thermo(*pThermoPtr);
    compressible::turbulenceModel &turbulence(*turbulencePtr);

    dimensionedScalar v_zero("v_zero", dimVolume / dimTime, 0.0);

    // Courant numbers used to adjust the time-step
    // scalar CoNum = 0.0;
    // scalar meanCoNum = 0.0;

    Info << "\nStarting time loop\n" << endl;

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

        volScalarField rPsi("rPsi", 1.0 / psi);
        surfaceScalarField rPsi_pos(interpolate(rPsi, pos, T.name()));
        surfaceScalarField rPsi_neg(interpolate(rPsi, neg, T.name()));

        surfaceScalarField e_pos(interpolate(e, pos, T.name()));
        surfaceScalarField e_neg(interpolate(e, neg, T.name()));

        surfaceVectorField U_pos("U_pos", rhoU_pos / rho_pos);
        surfaceVectorField U_neg("U_neg", rhoU_neg / rho_neg);

        surfaceScalarField p_pos("p_pos", rho_pos * rPsi_pos);
        surfaceScalarField p_neg("p_neg", rho_neg * rPsi_neg);

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

        surfaceScalarField a_pos
        (
            "a_pos",
            ap / (ap - am)
        );

        // surfaceScalarField amaxSf("amaxSf", max(mag(am), mag(ap)));
        if (amaxSfPtr == NULL)
        {
            amaxSfPtr = new surfaceScalarField("amaxSf", max(mag(am), mag(ap)));
        }
        
        surfaceScalarField &amaxSf(*amaxSfPtr);

        amaxSf = max(mag(am), mag(ap));

        surfaceScalarField aSf("aSf", am * a_pos);

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

        Info << "Time = " << runTime.timeName() << nl << endl;

        phi = aphiv_pos * rho_pos + aphiv_neg * rho_neg;

        surfaceVectorField phiUp
        (
            (aphiv_pos * rhoU_pos + aphiv_neg * rhoU_neg) +
            (a_pos * p_pos + a_neg * p_neg) * mesh.Sf()
        );

        surfaceScalarField phiEp
        (
            "phiEp",
            aphiv_pos * (rho_pos * (e_pos + 0.5 * magSqr(U_pos)) + p_pos) +
            aphiv_neg * (rho_neg * (e_neg + 0.5 * magSqr(U_neg)) + p_neg) +
            aSf * p_pos - aSf * p_neg
        );

        // Make flux for pressure-work absolute
        if (mesh.moving())
        {
            phiEp += mesh.phi() * (a_pos * p_pos + a_neg * p_neg);
        }

        volScalarField muEff("muEff", turbulence.muEff());
        volTensorField tauMC("tauMC", muEff * dev2(Foam::T(fvc::grad(U))));

        // --- Solve density
        solve(fvm::ddt(rho) + fvc::div(phi));

        // --- Solve momentum
        solve(fvm::ddt(rhoU) + fvc::div(phiUp));

        U.ref() = rhoU() / rho();

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
                fvc::interpolate(muEff) * mesh.magSf() * fvc::snGrad(U) +
                fvc::dotInterpolate(mesh.Sf(), tauMC)
            ) & (a_pos * U_pos + a_neg * U_neg)
        );

        solve
        (
            fvm::ddt(rhoE)
          + fvc::div(phiEp)
          - fvc::div(sigmaDotU)
        );


        e = rhoE / rho - 0.5 * magSqr(U);
        e.correctBoundaryConditions();
        thermo.correct();
        rhoE.boundaryFieldRef() ==
            rho.boundaryField() *
                (e.boundaryField() + 0.5 * magSqr(U.boundaryField()));

        if (!inviscid)
        {
            solve
            (
                fvm::ddt(rho, e) - fvc::ddt(rho, e) -
                fvm::laplacian(turbulence.alphaEff(), e)
            );
            thermo.correct();
            rhoE = rho * (e + 0.5 * magSqr(U));
        }

        p.ref() = rho() / psi();

        p.correctBoundaryConditions();
        rho.boundaryFieldRef() == psi.boundaryField() * p.boundaryField();

        turbulence.correct();

        runTime.write();

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl
             << endl;
    }

    Info << "End\n" << endl;

    loopStat = 0;
    return loopStat;
}



int rhoCentral::step()
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
    Foam::psiThermo &thermo(*pThermoPtr);
    compressible::turbulenceModel &turbulence(*turbulencePtr);

    dimensionedScalar v_zero("v_zero", dimVolume / dimTime, 0.0);

    // Courant numbers used to adjust the time-step
    // scalar CoNum = 0.0;
    // scalar meanCoNum = 0.0;

    //Info << "\nStarting time loop\n" << endl;

    //while (runTime.run())
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

        volScalarField rPsi("rPsi", 1.0 / psi);
        surfaceScalarField rPsi_pos(interpolate(rPsi, pos, T.name()));
        surfaceScalarField rPsi_neg(interpolate(rPsi, neg, T.name()));

        surfaceScalarField e_pos(interpolate(e, pos, T.name()));
        surfaceScalarField e_neg(interpolate(e, neg, T.name()));

        surfaceVectorField U_pos("U_pos", rhoU_pos / rho_pos);
        surfaceVectorField U_neg("U_neg", rhoU_neg / rho_neg);

        surfaceScalarField p_pos("p_pos", rho_pos * rPsi_pos);
        surfaceScalarField p_neg("p_neg", rho_neg * rPsi_neg);

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

        surfaceScalarField a_pos
        (
            "a_pos",
            ap / (ap - am)
        );

        // surfaceScalarField amaxSf("amaxSf", max(mag(am), mag(ap)));
        if (amaxSfPtr == NULL)
        {
            amaxSfPtr = new surfaceScalarField("amaxSf", max(mag(am), mag(ap)));
        }
        
        surfaceScalarField &amaxSf(*amaxSfPtr);

        amaxSf = max(mag(am), mag(ap));

        surfaceScalarField aSf("aSf", am * a_pos);

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

        Info << "Time = " << runTime.timeName() << nl << endl;

        phi = aphiv_pos * rho_pos + aphiv_neg * rho_neg;

        surfaceVectorField phiUp
        (
            (aphiv_pos * rhoU_pos + aphiv_neg * rhoU_neg) +
            (a_pos * p_pos + a_neg * p_neg) * mesh.Sf()
        );

        surfaceScalarField phiEp
        (
            "phiEp",
            aphiv_pos * (rho_pos * (e_pos + 0.5 * magSqr(U_pos)) + p_pos) +
            aphiv_neg * (rho_neg * (e_neg + 0.5 * magSqr(U_neg)) + p_neg) +
            aSf * p_pos - aSf * p_neg
        );

        // Make flux for pressure-work absolute
        if (mesh.moving())
        {
            phiEp += mesh.phi() * (a_pos * p_pos + a_neg * p_neg);
        }

        volScalarField muEff("muEff", turbulence.muEff());
        volTensorField tauMC("tauMC", muEff * dev2(Foam::T(fvc::grad(U))));

        // --- Solve density
        solve(fvm::ddt(rho) + fvc::div(phi));

        // --- Solve momentum
        solve(fvm::ddt(rhoU) + fvc::div(phiUp));

        U.ref() = rhoU() / rho();

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
                fvc::interpolate(muEff) * mesh.magSf() * fvc::snGrad(U) +
                fvc::dotInterpolate(mesh.Sf(), tauMC)
            ) & (a_pos * U_pos + a_neg * U_neg)
        );

        solve
        (
            fvm::ddt(rhoE)
          + fvc::div(phiEp)
          - fvc::div(sigmaDotU)
        );


        e = rhoE / rho - 0.5 * magSqr(U);
        e.correctBoundaryConditions();
        thermo.correct();
        rhoE.boundaryFieldRef() ==
            rho.boundaryField() *
                (e.boundaryField() + 0.5 * magSqr(U.boundaryField()));

        if (!inviscid)
        {
            solve
            (
                fvm::ddt(rho, e) - fvc::ddt(rho, e) -
                fvm::laplacian(turbulence.alphaEff(), e)
            );
            thermo.correct();
            rhoE = rho * (e + 0.5 * magSqr(U));
        }

        p.ref() = rho() / psi();

        p.correctBoundaryConditions();
        rho.boundaryFieldRef() == psi.boundaryField() * p.boundaryField();

        turbulence.correct();

        runTime.write();

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl
             << endl;
    }

    //Info << "End\n" << endl;

    // Setting comFoam variables that will be registered ^^
    //   with COM. This are needed for flow control when
    //   solver runs step-by-step.
    winTime = runTime.value();
    winDeltaT = runTime.deltaTValue();
    winRun = static_cast<int>(runTime.run());
    //-----------------------------------------------------

    stepStat = 0;
    return stepStat;
}

int rhoCentral::createFields()
{
    Foam::Time &runTime(*runTimePtr);
    dynamicFvMesh &mesh(*meshPtr);

    //  createRDeltaT.H  ^^^^^^^^^^^^^
    createRDeltaT();
    // -------------------------------

    Info << "Reading thermophysical properties\n" << endl;

    pThermoPtr = autoPtr<psiThermo>(psiThermo::New(mesh));
    Foam::psiThermo &thermo(*pThermoPtr);

    ePtr = &thermo.he();
    volScalarField &e(*ePtr);

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
        rho * (e + 0.5 * magSqr(U))
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

    phiPtr = new surfaceScalarField
    (
        "phi",
        fvc::flux(rhoU)
    );
    surfaceScalarField &phi(*phiPtr);

    Info << "Creating turbulence model\n" << endl;

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

int rhoCentral::createFieldRefs()
{
    Foam::psiThermo &thermo(*pThermoPtr);

    /*
    volScalarField& p = thermo.p();
    const volScalarField& T = thermo.T();
    const volScalarField& psi = thermo.psi(); */
    const volScalarField &mu = thermo.mu();

    pPtr = &thermo.p();
    TPtr = &thermo.T();
    psiPtr = &thermo.psi();

    // bool inviscid(true);
    inviscid = true;
    if (max(mu.primitiveField()) > 0.0)
    {
        inviscid = false;
    }

    return 0;
}

int rhoCentral::readFluxScheme()
{
    dynamicFvMesh &mesh(*meshPtr);

    // word fluxScheme("Kurganov");
    word fluxScheme("Kurganov");
    if (mesh.schemesDict().readIfPresent("fluxScheme", fluxScheme))
    {
        if ((fluxScheme == "Tadmor") || (fluxScheme == "Kurganov"))
        {
            Info << "fluxScheme: " << fluxScheme << endl;
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


int rhoCentral::centralCourantNo()
{
    Foam::Time &runTime(*runTimePtr);
    dynamicFvMesh &mesh(*meshPtr);
    surfaceScalarField &amaxSf(*amaxSfPtr);

    if (mesh.nInternalFaces())
    {
        scalarField sumAmaxSf(fvc::surfaceSum(amaxSf)().primitiveField());

        CoNum =
            0.5 * gMax(sumAmaxSf / mesh.V().field()) * runTime.deltaTValue();

        meanCoNum = 0.5 * (gSum(sumAmaxSf) / gSum(mesh.V().field())) *
                    runTime.deltaTValue();
    }

    Info << "Mean and max Courant Numbers = " << meanCoNum << " " << CoNum
         << endl;

    return 0;
}

int rhoCentral::setRDeltaT()
{
    Foam::Time &runTime(*runTimePtr);
    dynamicFvMesh &mesh(*meshPtr);
    surfaceScalarField &amaxSf(*amaxSfPtr);

    volScalarField &rDeltaT = trDeltaT.ref();

    scalar rDeltaTSmoothingCoeff
    (
        runTime.controlDict().lookupOrDefault<scalar>
        (
            "rDeltaTSmoothingCoeff",
            0.02
        )
    );

    // Set the reciprocal time-step from the local Courant number
    rDeltaT.ref() = max(1 / dimensionedScalar(dimTime, maxDeltaT),
                        fvc::surfaceSum(amaxSf)()() / ((2 * maxCo) * mesh.V()));

    // Update tho boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);

    Info << "Flow time scale min/max = "
         << gMin(1 / rDeltaT.primitiveField())
         << ", " << gMax(1 / rDeltaT.primitiveField()) << endl;

    return 0;
}

rhoCentral::~rhoCentral()
{
    finalizeStat = finalize();
}

int rhoCentral::finalize()
{
    // Delete thing that are allocated here
    if (posPtr != NULL) {delete posPtr; posPtr = NULL;}
    if (negPtr != NULL) {delete negPtr; negPtr = NULL;}
    if (amaxSfPtr != NULL) {delete amaxSfPtr; amaxSfPtr = NULL;}
    if (UPtr != NULL) {delete UPtr; UPtr=NULL;}
    if (rhoPtr != NULL) {delete rhoPtr; rhoPtr=NULL;}
    if (rhoUPtr != NULL) {delete rhoUPtr; rhoUPtr=NULL;}
    if (rhoEPtr != NULL) {delete rhoEPtr; rhoEPtr=NULL;}
    if (phiPtr != NULL) {delete phiPtr; phiPtr=NULL;}

    //delete argsPtr; Let it be the last thing to delete in the
    //                parrent class:rocFoam

    //delete pThermoPtr;
    //delete runTimePtr;
    //delete ePtr; a pointer
    //delete pPtr; a pointer
    //delete TPtr; a pointer
    //delete psiPtr; a pointer

    //delete meshPtr;
    //delete turbulencePtr;
    //delete trDeltaT;

    finalizeStat = 0;
    
    return finalizeStat;
}


double rhoCentral::errorEvaluate(int argc, char *argv[])
{
    createArgs(argc, argv);
    setRootCase();
    
    Foam::argList &args(*argsPtr);
    
    if (args.optionFound("list"))
    {
        functionObjectList::list();
        return 0;
    }
    
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

        createFields();

        volScalarField &rho(*rhoPtr);
        volVectorField &rhoU(*rhoUPtr);
        volScalarField &rhoE(*rhoEPtr);
        
        rhoVec.push_back(rho);
        rhoUVec.push_back(rhoU);
        UMagVec.push_back(magSqr(rhoU));
        rhoEVec.push_back(rhoE);
    }

    Info << "Field vectors created. " << endl;

    rhoVec[0]  = mag(rhoVec[2]  - rhoVec[1]);
    UMagVec[0] = mag(UMagVec[2]-UMagVec[1]);
    rhoEVec[0] = mag(rhoEVec[2] - rhoEVec[1]);
    
    Info << "Infinity norm = " << max(rhoVec[0]) << endl;
    Info << "Infinity norm = " << max(UMagVec[0]) << endl;
    Info << "Infinity norm = " << max(rhoEVec[0]) << endl;

    double maxError = std::max( max(rhoVec[0]).value(), max(UMagVec[0]).value() );
    testStat = std::max( maxError , max(rhoEVec[0]).value() );
    
    return testStat;
}

//===================================================================
