#include "rocRhoPimple.H"

int main(int argc, char *argv[])
{
    rhoPimple rocFoamPimple(argc, argv);


    Info << "\nStarting time loop\n" << endl;
    while (rocFoamPimple.runTimePtr->run())
    {
        rocFoamPimple.step();
    }
    Info << "End\n" << endl;
    rocFoamPimple.stepStat = 0;

    //rocFoamPimple.loop();
    //rocFoamPimple.finalize();

    return 0;
}
