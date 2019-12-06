#include "rocRhoCentral.H"

int main(int argc, char *argv[])
{
    rhoCentral rocFoamCentral(argc, argv);

    //rocFoamCentral.overrideTimeStep = true;
    //rocFoamCentral.overrideDeltaT = 0.0007;

    //rocFoamCentral.loop();

    Info << "\nStarting time loop\n" << endl;
    while (rocFoamCentral.runTimePtr->run())
    {
        rocFoamCentral.step();
    }
    Info << "End\n" << endl;
    rocFoamCentral.stepStat = 0;

    //rocFoamCentral.finalize();

    return 0;
}
