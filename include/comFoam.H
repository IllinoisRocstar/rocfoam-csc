#ifndef COMFOAM_
#define COMFOAM_

#include "com.h"
#include "com_devel.hpp"
#include "rocRhoCentral.H"

class rocFoamModule : public rocRhoCentral, public COM_Object
{
public:
    int initializeAll(int *pargc, void **pargv, int *verbIn);


private:
    std::string my_window_name; /// Tracks *this* window name.

protected:
    // ====== Data to register with IMPACT =======
    //  Surface Mesh
    /*
    int numPointsSurface;
    int numElementsSurface;
    std::vector<double> surfaceCoordinates;
    std::vector<unsigned int> surfaceConnectivity;
    std::map<int, int> surfaceNodeMap;
    std::map<int, int> interfaceToFoamNodeMap;
    std::map<int, int> foamGlobalToPatchNodeMap;
    std::vector<double> surfacePressure;
    std::vector<double> surfaceTraction;
    std::vector<double> surfaceDisplacement;
    std::vector<double> solidDisplacement;
    std::vector<double> time;
    std::vector<double> deltaTime;
    int verbosity;
    */

    // function status variables
    /*
    std::vector<int> initStatus;
    std::vector<int> runStatus;
    */

    // timing variables
    // std::vector<double> endTime;

    // ====== OpenFOAM stuff ========
    //
    //    General, init stuff
    // Foam::argList *argsPtr;
    // Foam::Time *runTimePtr;
    //
    //    Fluid data
    /*
    autoPtr<dynamicFvMesh> meshPtr;
    IOdictionary *transportPropertiesPtr;
    dimensionedScalar *nuPtr;
    dimensionedScalar *rhoFluidPtr;
    volScalarField *pPtr;
    volVectorField *UPtr;
    surfaceScalarField *phiPtr;
    label pRefCell;
    scalar pRefValue;
    */   
}
