#include "rocRhoPimple.H"

int main(int argc, char *argv[])
{
    rhoPimple rocFoamPimple(argc, argv);


    Info << "\nStarting time loop\n" << endl;

#ifdef HAVE_OFE20
    while (rocFoamPimple.runTimePtr->run())
#elif defined(HAVE_OF7)
    while (rocFoamPimple.runTimePtr->run())
#elif defined(HAVE_OF8)
    Foam::Time &runTime(*(rocFoamPimple.runTimePtr));
    pimpleControl& pimCtrl(rocFoamPimple.getPimpleControl());

    while ( pimCtrl.run(runTime) )
#endif
    {
        rocFoamPimple.step();
    }
    Info << "End\n" << endl;
    rocFoamPimple.stepStat = 0;

    //rocFoamPimple.loop();
    //rocFoamPimple.finalize();

    return 0;
}
