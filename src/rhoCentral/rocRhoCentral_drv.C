#include "rocRhoCentral.H"

int main(int argc, char *argv[])
{
    rhoCentral rocFoamCentral(argc, argv);

    rocFoamCentral.loop();
    //rocFoamCentral.finalize();

    return 0;
}
