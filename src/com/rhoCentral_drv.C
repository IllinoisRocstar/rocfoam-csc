#include "comRhoCentral.H"

int main(int argc, char *argv[])
{
    comRhoCentralModule rocFoamCentral(argc, argv);

    rocFoamCentral.loop();
    //rocFoamCentral.finalize();

    return 0;
}
