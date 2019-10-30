#include "comFoam.H"

int main(int argc, char *argv[])
{
 
    comFoamModule  comFoam(argc, argv)
    
    //comFoam.comStatus();
    comFoam.comLoop();
    comFoam.comFinalize();

    return 0;
}
