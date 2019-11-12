#include "comRhoPimple.H"

int main(int argc, char *argv[])
{
    comRhoPimpleModule rocFoamPimple(argc, argv);

    rocFoamPimple.loop();
    //rocFoamPimple.finalize();

    return 0;
}
