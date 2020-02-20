#include "comFoam.H"

comFoam::comFoam()
    : nPoints(0),
      nFaces(0),
      nCells(0),
      nPatches(0),
      winVolName(""),
      winSrfName(""),
      solverType(""),
      winComm(NULL),
      winNProc(0),
      winRank(0),
      winTime(0.0),
      winDeltaT(0.0),
      winRun(1)
{};

comFoam::comFoam(int *pargc, void **pargv, const char *name)
    : nPoints(0),
      nFaces(0),
      nCells(0),
      nPatches(0),
      winVolName(""),
      winSrfName(""),
      solverType(""),
      winComm(NULL),
      winNProc(0),
      winRank(0),
      winTime(0.0),
      winDeltaT(0.0),
      winRun(1)
{
    flowInit(pargc, pargv, name);
}

//^^^ DEFINITION OF COM-RELATED MTHODS ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
int comFoam::flowInit(int *pargc, void **pargv, const char *name)
{
    MPI_Comm tmpComm = COM_get_default_communicator();  

    int tmpRank;
    MPI_Comm_rank(tmpComm, &tmpRank);
    
    if (tmpRank == 0)
    {
        std::cout << "rocFoam.flowInit: Initializing flow solver."
                  << std::endl;
    }

std::cout<<__LINE__<< std::endl;

    //  OpenFOAM initializer ^^^^^^^^^^^^^^^^^^^^
    comFoam *comFoamPtr = NULL;

std::cout<<__LINE__<<std::endl;

    std::string volName = name+string("VOL");

std::cout<<__LINE__<<std::endl;

    std::string objectName = volName+string(".object");

std::cout<<__LINE__<<std::endl;

    COM_get_object(objectName.c_str(), 0, &comFoamPtr);


std::cout<<__LINE__<<std::endl;

    int argc = *pargc;
    char** argv = reinterpret_cast<char**>(pargv);

    comFoamPtr->initialize(argc, argv);

Foam::Info<<__LINE__<<endl;

    //  Other initializations ^^^^^^^^^^^^^^^^^^^
    // extractData can be called here, or in
    // rocFoam driver with comExtractData

    comFoamPtr->extractData();

Foam::Info<<__LINE__<<endl;

    return 0;
}


int comFoam::flowLoop(const char *name)
{
    Foam::Info << "rocFoam.flowLoop: Iterating flow solver." << Foam::endl;

    //  Call the flow iterator ^^^^^^^^^^^^^^^^^^
    comFoam *comFoamPtr = NULL;
    std::string volName = name+string("VOL");
    std::string objectName = volName+string(".object");
    COM_get_object(objectName.c_str(), 0, &comFoamPtr);

    comFoamPtr->loop();
    
    return 0;
}

int comFoam::flowStep(const char *name)
{

    Foam::Info << "rocFoam.flowStep: Stepping flow solver." << Foam::endl;

    //  Call the flow iterator ^^^^^^^^^^^^^^^^^^
    comFoam *comFoamPtr = NULL;
    std::string volName = name+string("VOL");
    std::string objectName = volName+string(".object");
    COM_get_object(objectName.c_str(), 0, &comFoamPtr);

    comFoamPtr->step();
    
    return 0;
}

int comFoam::flowExtractData(const char *name)
{

    Foam::Info << "rocFoam.extractData: Extracting flow data." << Foam::endl;

    //  Call the flow iterator ^^^^^^^^^^^^^^^^^^
    comFoam *comFoamPtr = NULL;
    std::string volName = name+string("VOL");
    std::string objectName = volName+string(".object");
    COM_get_object(objectName.c_str(), 0, &comFoamPtr);

    comFoamPtr->extractData();
    
    return 0;
}

//^^^^^ REGISTER FUNCTIONS ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
int comFoam::registerFunctions(const char *name)
{
    //  Anouncing default communicator  ^^^^^^^^^^^^^^^^^^^
    MPI_Comm tmpComm = COM_get_default_communicator();  

    int tmpRank;
    MPI_Comm_rank(tmpComm, &tmpRank);
    
    if (tmpRank == 0)
    {
        std::cout << "rocFoam.flowRegister: Registering flow functions with name "
                  << name << std::endl;
    }
    
    //  Register module with COM ^^^^^^^^^^^^^^^^^^^^^^^^^^
    comFoam *comFoamPtr = NULL;

    //std::string name="ROCFOAM";
    std::string objectName = name+string(".object");
    COM_get_object(objectName.c_str(), 0, &comFoamPtr);

    /// Register functions ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    std::vector<COM_Type> types(13,COM_INT);

    types[0] = COM_RAWDATA;
    types[2] = COM_VOID;

    std::string functionName = name+string(".flowInit");
    COM_set_member_function
    (
        functionName.c_str(),
        reinterpret_cast<Member_func_ptr>(&comFoam::flowInit),
        objectName.c_str(),
        "biii",
        &types[0]
    );

    functionName = name+string(".flowLoop");
    COM_set_member_function
    (
        functionName.c_str(),
        reinterpret_cast<Member_func_ptr>(&comFoam::flowLoop),
        objectName.c_str(),
        "bi",
        &types[0]
    );

    functionName = name+string(".flowStep");
    COM_set_member_function
    (
        functionName.c_str(),
        reinterpret_cast<Member_func_ptr>(&comFoam::flowStep),
        objectName.c_str(),
        "bi",
        &types[0]
    );

    functionName = name+string(".flowExtractData");
    COM_set_member_function
    (
        functionName.c_str(),
        reinterpret_cast<Member_func_ptr>(&comFoam::flowExtractData),
        objectName.c_str(),
        "bi",
        &types[0]
    );

    //COM_set_member_function
    //(
    //    (name + string(".flowFin")).c_str(),
    //    reinterpret_cast<Member_func_ptr>(&rhoCentral::flowFin),
    //    objectName.c_str(), "b", &types[0]
    //);

    //  Registering data of this module to COM ^^^^^^^^^^^^
    std::string dataName = name+string(".winNProc");
    COM_new_dataitem( dataName.c_str(), 'w', COM_INT, 1, "");
    COM_set_size(     dataName.c_str(), 0, 1);
    COM_set_array(    dataName.c_str(), 0, &(comFoamPtr->winNProc));

    dataName = name+string(".winTime");
    COM_new_dataitem( dataName.c_str(), 'w', COM_DOUBLE, 1, "");
    COM_set_size(     dataName.c_str(), 0, 1);
    COM_set_array(    dataName.c_str(), 0, &(comFoamPtr->winTime));

    dataName = name+string(".winDeltaT");
    COM_new_dataitem( dataName.c_str(), 'w', COM_DOUBLE, 1, "");
    COM_set_size(     dataName.c_str(), 0, 1);
    COM_set_array(    dataName.c_str(), 0, &(comFoamPtr->winDeltaT) );

    dataName = name+string(".winRun");
    COM_new_dataitem( dataName.c_str(), 'w', COM_INT, 1, "");
    COM_set_size(     dataName.c_str(), 0, 1);
    COM_set_array(    dataName.c_str(), 0, &(comFoamPtr->winRun));

    COM_window_init_done(name); 

    return 0;
}
//---------------------------------------------------------

//===================================================================


int comFoam::extractData()
{

    const dynamicFvMesh  &mesh(*meshPtr);
    const Foam::Time     &runTime(*runTimePtr);
    const volScalarField &p(*pPtr);
    const volVectorField &U(*UPtr);
    const volScalarField &T(*TPtr);
    const volScalarField &rho(*rhoPtr);

    // Mesh and conmnectivities ^^^^^^^^^^^^^^^^^
    const pointField &points     = mesh.points();
    const cellList   &cells      = mesh.cells();
    const faceList   &faces      = mesh.faces();
    const labelList  &faceOwner  = mesh.faceOwner();
    const labelList  &faceNeighb = mesh.faceNeighbour();    
    const polyBoundaryMesh &patches = mesh.boundaryMesh();
    //-------------------------------------------

    // Temporary Vectors ^^^^^^^^^^^^^^^^^^^^^^^^
    std::vector< std::vector<double> > vecVecTmpDouble;
    std::vector<double> vecTmpDouble;
    
    std::vector<int> vecTmpInt;
    //-------------------------------------------

    nPoints = mesh.nPoints();
    nFaces  = mesh.nFaces();
    nCells  = mesh.nCells();
    nPatches = patches.size();

    //  Registering node data of this module to COM ^^^^^^^


    //dynamicFvMesh &mesh(*meshPtr);

    /*
    {
        const pointField &points = mesh.points();

        nComponents = 3;
        nNodes = mesh.nPoints();
        
        int nTotal = nComponents * nNodes;
        comXYZ = new double[nTotal];

        forAll(points, i)
        {
            for(int j=0; j< nComponents; j++)
            *(comXYZ + i * nComponents) = points[i][j];
        }
    }
    */
    

    // Field Vectors ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    vecPoints.clear();
    vecFaceToPointConn.clear();
    vecOwners.clear();
    vecNeighbs.clear();
    
    vecFieldVel.clear();
    vecFieldRho.clear();
    vecFieldP.clear();
    vecFieldT.clear();
    //-------------------------------------------

    // Patch Vectors ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    patchGlobalPointIndex.clear();
    patchLocalConnectivity.clear();

    vecPatchVel.clear();
    vecPatchRho.clear();
    vecPatchP.clear();
    vecPatchT.clear();

    vecPatchName.clear();
    vecPatchType.clear();
    vecPatchInGroup.clear();

    vecPatchStart.clear();
    vecPatchNFaces.clear();
    vecPatchNPoints.clear();
    //-------------------------------------------

    // Save point positions ^^^^^^^^^^^^^^^^^^^^^
    forAll(points, i)
    {
        vecTmpDouble.clear();
        for(int j=0; j<nComponents; j++)
        {
            vecTmpDouble.push_back(points[i][j]);
        }
        vecPoints.push_back(vecTmpDouble);
    }

    // Save faceToCell and  faceToPoint conn ^^^^
    forAll(faces, i)
    {
        const labelList &pointsList = faces[i];
        
        vecOwners.push_back(faceOwner[i]);
        vecNeighbs.push_back(faceNeighb[i]);

        vecTmpInt.clear();
        forAll(pointsList, j)
        {
            vecTmpInt.push_back(pointsList[j]);
        }
        
        vecFaceToPointConn.push_back(vecTmpInt);
    }
    //-------------------------------------------

    // Save cell-centered field data ^^^^^^^^^^^^
    forAll(cells, i)
    {
        vecTmpDouble.clear();
        for(int j=0; j<nComponents; j++)
        {
            vecTmpDouble.push_back(U[i].component(j));
        }
        vecFieldVel.push_back(vecTmpDouble);

        vecFieldRho.push_back(rho[i]);
        vecFieldP.push_back(p[i]);
        vecFieldT.push_back(T[i]);
    }
    //-------------------------------------------

    // Save Patch faceToNode connectivity  ^^^^^^
    // Local and Global
    forAll(patches, ipatch)
    {
        const polyPatch &patch = patches[ipatch];

        const word &patchName = patch.name();
        const word &patchType = patch.type();
        const wordList &patchInGroup = patch.inGroups();

        const label &patchStart = patch.start();
        const int &patchSize = patch.size();

        vecPatchName   .push_back(patchName);
        vecPatchType   .push_back(patchType);
        vecPatchInGroup.push_back(patchInGroup);
        vecPatchStart  .push_back(patchStart);
        vecPatchNFaces .push_back(patchSize);

        std::vector<int> vecTmpInt_patchGlobalPointIndex;
        std::vector<std::vector<int>> vecTmpInt_patchLocalConnectivity;

        for(int iface=0; iface<patchSize; iface++)
        {
            const label &face = patchStart + iface;
            const labelList &pointsList = faces[face];

            vecTmpInt.clear();
            forAll(pointsList, ipoint)
            {
                const int &point = pointsList[ipoint];
                std::vector<int>::iterator index = std::find
                                  (
                                    vecTmpInt_patchGlobalPointIndex.begin(),
                                    vecTmpInt_patchGlobalPointIndex.end(),
                                    point
                                  );

                if  ( index == vecTmpInt_patchGlobalPointIndex.end() )
                {
                    vecTmpInt_patchGlobalPointIndex.push_back(point);

                    index = vecTmpInt_patchGlobalPointIndex.end() - 1;
                }

                int indexVal = std::distance(vecTmpInt_patchGlobalPointIndex.begin(), index);
                vecTmpInt.push_back( indexVal );
            }
            vecTmpInt_patchLocalConnectivity.push_back(vecTmpInt);
        }

        patchGlobalPointIndex.push_back(vecTmpInt_patchGlobalPointIndex);
        patchLocalConnectivity.push_back(vecTmpInt_patchLocalConnectivity);

        vecPatchNPoints.push_back(vecTmpInt_patchGlobalPointIndex.size()) ;
    }
    //-------------------------------------------

    //  Save Flow Quantities on Patches ^^^^^^^^^
    vecPatchVel.clear();
    vecPatchRho.clear();
    vecPatchP.clear();
    vecPatchT.clear();
    for(int ipatch=0; ipatch<nPatches; ipatch++)
    {
        const int &nFaces = vecPatchNFaces[ipatch];

        //  Patch Face Velocity
        vecVecTmpDouble.clear();
        if (vecPatchType[ipatch] != string("empty"))
        {
            for(int iface=0; iface<nFaces; iface++)
            {
                vecTmpDouble.clear();
                for(int k=0; k<nComponents; k++)
                {
                    vecTmpDouble.push_back( U.boundaryField()[ipatch][iface].component(k) );
                }
                vecVecTmpDouble.push_back(vecTmpDouble);
            }
        }
        vecPatchVel.push_back(vecVecTmpDouble);

        //  Patch Face Density
        vecTmpDouble.clear();
        if (vecPatchType[ipatch] != string("empty"))
        {
            for(int iface=0; iface<nFaces; iface++)
            {
                vecTmpDouble.push_back( rho.boundaryField()[ipatch][iface] );
            }
        }
        vecPatchRho.push_back(vecTmpDouble);

        //  Patch Face Pressure
        vecTmpDouble.clear();
        if (vecPatchType[ipatch] != string("empty"))
        {
            for(int iface=0; iface<nFaces; iface++)
            {
                vecTmpDouble.push_back( p.boundaryField()[ipatch][iface] );
            }
        }
        vecPatchP.push_back(vecTmpDouble);

        //  Patch Face Temperature
        vecTmpDouble.clear();
        if (vecPatchType[ipatch] != string("empty"))
        {
            for(int iface=0; iface<nFaces; iface++)
            {
                vecTmpDouble.push_back( T.boundaryField()[ipatch][iface] );
            }
        }
        vecPatchT.push_back(vecTmpDouble);
    }
    //-------------------------------------------

for(int ipatch=0; ipatch<nPatches; ipatch++)
{
    Foam::Info << " Patch = " << vecPatchName[ipatch] << " patchID = " << ipatch << endl;
    int faceStart = vecPatchStart[ipatch];
    int nFaces =  vecPatchVel[ipatch].size(); // vecPatchNFaces[ipatch];
    
    Foam::Info << " Global Face Indices = ";
    for(int iface=0; iface<nFaces; iface++)
    {
        Foam::Info << faceStart+iface << " ";
    }
    Foam::Info << endl;
    
    for(int iface=0; iface<nFaces; iface++)
    {
        Foam::Info  << "     "
                    << "faceID = " << iface
                    << " V = " << vecPatchVel[ipatch][iface][0]
                    << ", " << vecPatchVel[ipatch][iface][1] 
                    << ", " << vecPatchVel[ipatch][iface][2]
                    << " Rho = " << vecPatchRho[ipatch][iface]
                    << " P = " << vecPatchP[ipatch][iface] 
                    << " T = " << vecPatchT[ipatch][iface]
                    << endl;
    }
    Foam::Info << endl;
}
        
Foam::Info << endl;



    return 0;
}



