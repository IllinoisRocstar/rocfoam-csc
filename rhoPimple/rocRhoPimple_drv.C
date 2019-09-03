#include "rocRhoPimple.H"

int main(int argc, char *argv[]) {
    rocRhoPimple rocRhoPimple(argc, argv);

    rocRhoPimple.loop();
    rocRhoPimple.finalize();

    return 0;
}

