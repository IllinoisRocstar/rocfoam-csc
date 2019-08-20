#include "RHOCENTRAL_CLASS.H"

int main(int argc, char *argv[])
{
   RHOCENTRAL_CLASS rocFoamCentral(argc, argv);

   rocFoamCentral.loop();

   return 0;
}
