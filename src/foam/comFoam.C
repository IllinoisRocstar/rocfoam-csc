#include "comFoam.H"

comFoam::comFoam()
    : winNameVol(""),
      solverType(""),
      winComm(NULL),
      winNProc(0),
      winRank(0),
      winTime(0.0),
      winDeltaT(0.0),
      winRun(1)
{};

comFoam::comFoam(int *pargc, void **pargv, int *verbIn)
    : winNameVol(""),
      solverType(""),
      winComm(NULL),
      winNProc(0),
      winRank(0),
      winTime(0.0),
      winDeltaT(0.0),
      winRun(1)
{
    flowInit(pargc, pargv, verbIn);
}

//^^^ DEFINITION OF COM-RELATED MTHODS ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
int comFoam::flowInit(int *pargc, void **pargv, int *verbIn)
{
    MPI_Comm tmpComm = COM_get_default_communicator();  

    int tmpRank;
    MPI_Comm_rank(tmpComm, &tmpRank);
    
    if (tmpRank == 0)
    {
        std::cout << "rocFoam.flowInit: Initializing flow solver."
                  << std::endl;
    }

    //  OpenFOAM initializer ^^^^^^^^^^^^^^^^^^^^
    comFoam *comFoamPtr = NULL;
    std::string name="ROCFOAM";
    std::string globalName(name+".global");
    COM_get_object(globalName.c_str(), 0, &comFoamPtr);

    int argc = *pargc;
    char** argv = reinterpret_cast<char**>(pargv);

    comFoamPtr->initialize(argc, argv);
    
    //  Other initializations ^^^^^^^^^^^^^^^^^^^

    return 0;
}


int comFoam::flowLoop()
{

    Foam::Info << "rocFoam.flowLoop: Iterating flow solver." << Foam::endl;

    //  Call the flow iterator ^^^^^^^^^^^^^^^^^^
    comFoam *comFoamPtr = NULL;

    std::string name="ROCFOAM";
    std::string globalName(name+".global");
    COM_get_object(globalName.c_str(), 0, &comFoamPtr);

    comFoamPtr->loop();
    
    return 0;
}

int comFoam::flowStep()
{

    Foam::Info << "rocFoam.flowStep: Stepping flow solver." << Foam::endl;

    //  Call the flow iterator ^^^^^^^^^^^^^^^^^^
    comFoam *comFoamPtr = NULL;

    std::string name="ROCFOAM";
    std::string globalName(name+".global");
    COM_get_object(globalName.c_str(), 0, &comFoamPtr);

    comFoamPtr->step();
    
    return 0;
}

//^^^^^ REGISTER FUNCTIONS ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
int comFoam::flowRegister()
{
    //  Anouncing default communicator  ^^^^^^^^^^^^^^^^^^^
    MPI_Comm tmpComm = COM_get_default_communicator();  

    int tmpRank;
    MPI_Comm_rank(tmpComm, &tmpRank);
    
    if (tmpRank == 0)
    {
        std::cout << "rocFoam.flowRegister: Registering flow functions."
                  << std::endl << std::endl;
    }
    
    //  Register module with COM ^^^^^^^^^^^^^^^^^^^^^^^^^^
    comFoam *comFoamPtr = NULL;

    std::string name="ROCFOAM";
    std::string globalName(name+".global");
    COM_get_object(globalName.c_str(), 0, &comFoamPtr);

    /// Register functions ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    std::vector<COM_Type> types(13,COM_INT);

    types[0] = COM_RAWDATA;
    types[2] = COM_VOID;

    COM_set_member_function
    (
        (name + string(".flowInit")).c_str(),
        reinterpret_cast<Member_func_ptr>(&comFoam::flowInit),
        globalName.c_str(), "biii", &types[0]
    );


    COM_set_member_function
    (
        (name + string(".flowLoop")).c_str(),
        reinterpret_cast<Member_func_ptr>(&comFoam::flowLoop),
        globalName.c_str(), "b", &types[0]
    );

    COM_set_member_function
    (
        (name + string(".flowStep")).c_str(),
        reinterpret_cast<Member_func_ptr>(&comFoam::flowStep),
        globalName.c_str(), "b", &types[0]
    );

    //COM_set_member_function
    //(
    //    (name + string(".flowFin")).c_str(),
    //    reinterpret_cast<Member_func_ptr>(&rhoCentral::flowFin),
    //    globalName.c_str(), "b", &types[0]
    //);

    //  Registering data of this module to COM ^^^^^^^^^^^^
    std::string dataName="";

    dataName = name+string(".winNProc");
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
    const pointField &points    = mesh.points();
    const faceList  &faces      = mesh.faces();
    const labelList &faceOwner  = mesh.faceOwner();
    const labelList &faceNeighb = mesh.faceNeighbour();    
    const polyBoundaryMesh &patches = mesh.boundaryMesh();

    int nPoints  = mesh.nPoints();
    int nFaces   = faces.size();
    int nPatches = patches.size();
    //-------------------------------------------



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
    std::vector< std::vector<double> > vecPoints;
    std::vector< std::vector<int> > vecFaceToPointConn;
    std::vector<int> vecOwners;
    std::vector<int> vecNeighbs;
    


    std::vector< std::vector<double> > vecFieldVel;
    std::vector<double> vecFieldRho;
    std::vector<double> vecFieldP;
    std::vector<double> vecFieldT;
    //-------------------------------------------
    
    // Temporary Vectors ^^^^^^^^^^^^^^^^^^^^^^^^
    std::vector<double> vecTmpDouble;
    
    std::vector<int> vecTmpInt;
    //-------------------------------------------


    // Save point positions and filed data ^^^^^^
    forAll(points, i)
    {
        vecTmpDouble.clear();
        for(int j=0; j<nComponents; j++)
        {
            vecTmpDouble.push_back(points[i][j]);
        }
        vecPoints.push_back(vecTmpDouble);
        
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

    // Save faceToCell and  faceToPoint conn ^^^^
    forAll(faces, i)
    {
        const labelList &points = faces[i];
        
        vecOwners.push_back(faceOwner[i]);
        vecNeighbs.push_back(faceNeighb[i]);

        vecTmpInt.clear();
        forAll(points, j)
        {
            vecTmpInt.push_back(points[j]);
        }
        
        vecFaceToPointConn.push_back(vecTmpInt);
    }
    //-------------------------------------------


        
        
        

        std::vector< std::vector< std::vector<double> >> patchVel;
        std::vector< std::vector<double> > patchRho;
        std::vector< std::vector<double> > patchP;
        std::vector< std::vector<double> > patchT;

        std::vector<std::string> vecPatchName;
        std::vector<std::string> vecPatchType;
        // Note: I need to change this
        std::vector< wordList > vecPatchInGroup;

        std::vector<int> vecPatchStart;
        std::vector<int> vecPatchNFaces;
        std::vector<int> vecPatchNPoints;

        


        std::vector< std::vector<int> > patchGlobalPointIndex;
        std::vector< std::vector< std::vector<int> >> patchLocalConnectivity;


        forAll(patches, ipatch)
        {

//if (ipatch >= 1) break;

            const polyPatch &patch = patches[ipatch];

            const word &patchName = patch.name();
            const word &patchType = patch.type();
            const wordList &patchInGroup = patch.inGroups();

            const label& patchStart = patch.start();
            const int& patchSize = patch.size();

            vecPatchName   .push_back(patchName);
            vecPatchType   .push_back(patchType);
            vecPatchInGroup.push_back(patchInGroup);
            vecPatchStart  .push_back(patchStart);
            vecPatchNFaces .push_back(patchSize);

            //std::vector<int> vecTmpInt_globalFaceIndex;
            //std::vector<int> vecTmpInt_localFaceIndex;
            
            std::vector<int> vecTmpInt_patchGlobalPointIndex;

            //std::vector<int> VecTmpInt_globalConnectivity;
            std::vector<std::vector<int>> vecTmpInt_patchLocalConnectivity;

            Foam::Info << " Patch " << patchName << " patchStart = " << patchStart << " patchSize = " << patchSize << endl;

            //VecTmpInt_globalFaceIndex.clear();

            vecTmpInt_patchGlobalPointIndex.clear();
            vecTmpInt_patchLocalConnectivity.clear();
            for(int iface=0; iface<patchSize; iface++) //iface<5; iface++) 
            {
                const label& face = patchStart + iface; //patch.whichFace(patchStart + iFace);
                
                //vecTmpInt_localFaceIndex.push_back(iface);
                //vecTmpInt_globalFaceIndex.push_back(face);

                const labelList &points = faces[face];
                int nPoints = points.size();

//Foam::Info << "      face " << " local-to-global " << iface << endl << "       ";

                
                std::vector<int> vecTmpInt;
                //vecTmpInt.clear();
                forAll(points, ipoint)
                {
                    const int& point = points[ipoint];

//Foam::Info << point << " ";


                    //VecTmpInt_globalConnectivity.push_back( point );

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


//Foam::Info << endl;

            }

            patchGlobalPointIndex.push_back(vecTmpInt_patchGlobalPointIndex);
            patchLocalConnectivity.push_back(vecTmpInt_patchLocalConnectivity);

            vecPatchNPoints.push_back(vecTmpInt_patchGlobalPointIndex.size()) ;
        }


Foam::Info << endl;

        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {

if (vecPatchName[ipatch] != string("outlet")) continue;

            Foam::Info << " Patch " << vecPatchName[ipatch]
                       << " patchStart = " << vecPatchStart[ipatch]
                       << " nFaces = " << vecPatchNFaces[ipatch]
                       << " nPoints = " << vecPatchNPoints[ipatch]
                       << endl;
        }



        //  Save Flow Quantities
        std::vector< std::vector<double> > vecVecTmpDouble;
        std::vector<double> vecTmpDouble;
        

        patchVel.clear();
        patchRho.clear();
        patchP.clear();
        patchT.clear();
        for(int ipatch=0; ipatch<nPatches ; ipatch++)
        {

            int faceStart = vecPatchStart[ipatch];
            int nFaces = vecPatchNFaces[ipatch];

            vecVecTmpDouble.clear();
            if (vecPatchType[ipatch] != string("empty"))
            {
                for(int iface=0; iface<nFaces; iface++)
                {
Foam::Info << ipatch << iface;

                    vecTmpDouble.clear();
                    for(int k=0; k<nComponents; k++)
                    {
                        vecTmpDouble.push_back( U.boundaryField()[ipatch][iface].component(k) );
                    }
                    vecVecTmpDouble.push_back(vecTmpDouble);

Foam::Info << " " << ipatch << iface << endl;

                }
            }
            patchVel.push_back(vecVecTmpDouble);



            //  Patch Face Density
            vecTmpDouble.clear();
            if (vecPatchType[ipatch] != string("empty"))
            {
                for(int iface=0; iface<nFaces; iface++)
                {
                    vecTmpDouble.push_back( rho.boundaryField()[ipatch][iface] );
                }
            }
            patchRho.push_back(vecTmpDouble);

            //  Patch Face Pressure
            vecTmpDouble.clear();
            if (vecPatchType[ipatch] != string("empty"))
            {
                for(int iface=0; iface<nFaces; iface++)
                {
                    vecTmpDouble.push_back( p.boundaryField()[ipatch][iface] );
                }
            }
            patchP.push_back(vecTmpDouble);


            //  Patch Face Temperature
            vecTmpDouble.clear();
            if (vecPatchType[ipatch] != string("empty"))
            {
                for(int iface=0; iface<nFaces; iface++)
                {
                    vecTmpDouble.push_back( T.boundaryField()[ipatch][iface] );
                }
            }
            patchT.push_back(vecTmpDouble);



        }


        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {
            Foam::Info << " Patch = " << vecPatchName[ipatch] << " patchID = " << ipatch << endl;

            int faceStart = vecPatchStart[ipatch];
            int nFaces =  patchVel[ipatch].size(); // vecPatchNFaces[ipatch];
            
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
                            << " V = " << patchVel[ipatch][iface][0]
                            << ", " << patchVel[ipatch][iface][1] 
                            << ", " << patchVel[ipatch][iface][2]
                            << " Rho = " << patchRho[ipatch][iface]
                            << " P = " << patchP[ipatch][iface] 
                            << " T = " << patchT[ipatch][iface]
                            << endl;
            }
            Foam::Info << endl;
        }
        
Foam::Info << endl;






/*
        for(int ipatch=0; ipatch<nPatches; ipatch++)
        {
            int npoints = vecPatchNPoints[ipatch];

            //  Patch Nodeal Velocity
            vecVecTmpDouble.clear();
            for(int ipoint=0; ipoint<npoints; ipoint++)
            {
                int index = patchGlobalPointIndex[ipatch][ipoint];
            
                vecTmpDouble.clear();
                for(int k=0; k<nComponents; k++)
                {
                    vecTmpDouble.push_back(U[index].component(k));
                }
                vecVecTmpDouble.push_back(vecTmpDouble);
            }
            patchVel.push_back(vecVecTmpDouble);

            //  Patch Nodeal Density
            vecTmpDouble.clear();
            for(int ipoint=0; ipoint<npoints; ipoint++)
            {
                int index = patchGlobalPointIndex[ipatch][ipoint];
                vecTmpDouble.push_back(rho[index]);
            }
            patchRho.push_back(vecTmpDouble);

            //  Patch Nodeal Pressure
            vecTmpDouble.clear();
            for(int ipoint=0; ipoint<npoints; ipoint++)
            {
                int index = patchGlobalPointIndex[ipatch][ipoint];
                vecTmpDouble.push_back(p[index]);
            }
            patchP.push_back(vecTmpDouble);

            //  Patch Nodeal Temperature
            vecTmpDouble.clear();
            for(int ipoint=0; ipoint<npoints; ipoint++)
            {
                int index = patchGlobalPointIndex[ipatch][ipoint];
                vecTmpDouble.push_back(T[index]);
            }
            patchT.push_back(vecTmpDouble);
        }

        for( std::vector<std::vector<int>>::iterator itVecVec=patchGlobalPointIndex.begin(); itVecVec != patchGlobalPointIndex.end(); itVecVec++)
        {
            int ipatch = std::distance( patchGlobalPointIndex.begin(), itVecVec);

if (vecPatchName[ipatch] != string("outlet")) continue;

            Foam::Info << " Patch " << vecPatchName[ipatch] << " local-to-global " << ipatch << endl;

            int npoints = vecPatchNPoints[ipatch];
            Foam::Info << " Global Node Indices = ";
            for(int ipoint=0; ipoint<npoints; ipoint++)
            {
                Foam::Info << patchGlobalPointIndex[ipatch][ipoint] << " ";
            }
            Foam::Info << endl;

            
            for( std::vector<int>::iterator itVec=itVecVec->begin(); itVec != itVecVec->end(); itVec++)
            {
                int ipoint = std::distance( itVecVec->begin(), itVec);
            
            
                Foam::Info  << "     "
                           << "PointID = " << *itVec
                           << " V = " << patchVel[ipatch][ipoint][1]
                              << ", " << patchVel[ipatch][ipoint][2] 
                              << ", " << patchVel[ipatch][ipoint][3]
                           << " Rho = " << patchRho[ipatch][ipoint]
                           << " P = " << patchP[ipatch][ipoint] 
                           << " T = " << patchT[ipatch][ipoint]  << endl;
            }
            Foam::Info << endl;
        }
        
        
*/
Foam::Info << endl;

/*
        for (
                std::vector<std::vector<std::vector<int>>>::iterator
                    itVecVecVec=patchLocalConnectivity.begin();
                itVecVecVec != patchLocalConnectivity.end();
                itVecVecVec++
            )
        {
            int index = std::distance( 
                                        patchLocalConnectivity.begin(),
                                        itVecVecVec
                                     );

            Foam::Info << " Patch " << vecPatchName[index] << " connectivity " << endl;

            for( std::vector<std::vector<int>>::iterator itVecVec=itVecVecVec->begin(); itVecVec != itVecVecVec->end(); itVecVec++)
            {
                int faceindex = std::distance( itVecVecVec->begin(), itVecVec);
                Foam::Info << "      face " << " local-to-global " << faceindex << endl << "       ";
                
                for( std::vector<int>::iterator itVec=itVecVec->begin(); itVec != itVecVec->end(); itVec++)
                {
                    Foam::Info << *itVec << " ";
                }
                Foam::Info << endl;
            }
        }
        
        Foam::Info << endl;
*/


/*
            
            //  Patch Nodeal Velocity
            vecVecTmp.clear();
            for(int j=0; j<patchSize; j++)
            {
                int index = patchStart+j;

                vecTmp.clear();
                for(int k=0; k<nComponents; k++)
                {
                    vecTmp.push_back(U[index].component(k));
                }

                vecVecTmp.push_back(vecTmp);
            }
            vecPatchVel.push_back(vecVecTmp);



            //  Patch Nodeal Density
            vecTmp.clear();
            for(int j=0; j<patchSize; j++)
            {
                int index=patchStart+j;
                vecTmp.push_back(rho[index]);
            }
            vecPatchRho.push_back(vecTmp);

            //  Patch Nodeal Pressure
            vecTmp.clear();
            for(int j=0; j<patchSize; j++)
            {
                int index=patchStart+j;
                vecTmp.push_back(p[index]);
            }
            vecPatchP.push_back(vecTmp);

            //  Patch Nodeal Temperature
            vecTmp.clear();
            for(int j=0; j<patchSize; j++)
            {
                int index=patchStart+j;
                vecTmp.push_back(T[index]);
            }
            vecPatchT.push_back(vecTmp);


            //Foam::Info << i << " " << patchName << " " 
            //           << patchType << " " << patchInGroup 
            //           << " " << patchStart << " " << patchSize << endl;
        }
        
        Foam::Info <<  " SIZE OF T ==== " << vecPatchT.size() << endl;


        for(std::vector< std::vector<double> >::iterator vvI=vecPatchRho.begin(); vvI != vecPatchRho.end(); vvI++)
        {
            Foam::Info << " Rho of partch " << std::distance(vecPatchRho.begin(), vvI) << " = ";
            for(std::vector<double>::iterator vI=vvI->begin(); vI != vvI->begin()+5; vI++)
            {
                    Foam::Info << *vI << " ";
            }
            Foam::Info << endl;;
        }
*/



    //dataName = name+string(".nc");
    //COM_set_size( dataName.c_str(), 0, nNodes);
    //COM_set_array(dataName.c_str(), 0, comXYZ, 3);
        
    

    /*
    const faceList &faces = mesh.faces();
    int nFaces = faces.size();
    comFaces  = new int[nFaces];
    comOwner  = new int[nFaces];
    comNeighb = new int[nFaces];
    */




/*

    const cellList &cells = mesh.cells();
    const faceList &faces = mesh.faces();
    int nCells = cells.size();

    std::vector<int> vecPoints;
    std::vector< std::vector<int> > vecFaces;
    std::vector< std::vector< std::vector<int> > > vecCells;


Foam::Info << __LINE__ << endl;


    forAll(cells, i)
    {

        vecFaces.clear();

        const faceList &faces = mesh.faces();
        forAll(faces, j)
        {



            const labelList &points = faces[j];
            
            vecPoints.clear();
            forAll(points, k)
            {            
                vecPoints.push_back( points[k] );
            }
            vecFaces.push_back( vecPoints );

        }

//        vecCells.push_back( vecFaces );


    }

Foam::Info << __LINE__ << endl;

*/

//    std::vector< std::vector< std::vector<int> > >::iterator vvvIt = vecCells.begin();

/*    for (std::vector< std::vector<int> >::iterator vvIt=vvvIt->begin(); vvIt != vvvIt->end(); vvIt++)
    {
        for (std::vector<int>::iterator vIt=vvIt->begin(); vIt != vvIt->end(); vIt++)
        {
            Foam::Info << "Ponits Vector = " << *vIt;
        }
        Foam::Info << Foam::endl;
    }
*/


    //std::cout << "Cell Vector   = " << *vvvIt << endl;
    //Foam::Info << "Face Vector   = " << *vvIt << Foam::endl;
    //Foam::Info << "Ponits Vector = " << *vIt << Foam::endl;

    return 0;
}



