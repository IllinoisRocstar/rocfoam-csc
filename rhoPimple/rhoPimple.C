#include "RHOPIMPLE_CLASS.H"

int main(int argc, char *argv[])
{
   RHOPIMPLE_CLASS rocRhoPimple(argc, argv);

   rocRhoPimple.loop();
   rocRhoPimple.finalize();

   return 0;
}
