#include "rocRhoCentral.H"

int main(int argc, char *argv[])
{
    rocRhoCentral rocFoamCentral(argc, argv);

    rocFoamCentral.loop();
    rocFoamCentral.finalize();

    return 0;
}
