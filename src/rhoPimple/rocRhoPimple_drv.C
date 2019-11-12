#include "rocRhoPimple.H"

int main(int argc, char *argv[])
{
    rhoPimple rocFoamPimple(argc, argv);

    rocFoamPimple.loop();
    //rocFoamPimple.finalize();

    return 0;
}
