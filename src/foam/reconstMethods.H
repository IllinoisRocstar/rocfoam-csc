void comFoam::copyWindow(const char *name1, const char *name2)
{
    std::string winName1 = name1; //+string("VOL");
    std::string winName2 = name2; //+string("VOL");
    Info << endl
         << "rocFoam.copyWindow: window names "
         << winName1.c_str() << " & "
         << winName2.c_str() << endl;


    // Get window information ^^^^^^^^^^^^^^^^^^^
    int nItemps;
    std::string datNames;
    COM_get_dataitems(winName1.c_str(), &nItemps, datNames);

    Info << "rocFoam.copyWindow: nItemps in "
         << "winName1 = "
         << nItemps
         << endl;
    // ------------------------------------------

    //  Copy pane nodes & connectivities ^^^^^^^^
    int nPanes;
    int* paneList;
    COM_get_panes(winName1.c_str(), &nPanes, &paneList);
    Info << "rocFoam.copyWindow: Number of Panes in "
         << "winName1 = "
         << nPanes << endl;

    for(int ipane=0; ipane<nPanes; ipane++)
    {
        COM_clone_dataitem((winName2+".mesh").c_str(),
                           (winName1+".mesh").c_str(),
                           ipane);
    }
    // ------------------------------------------

    //  Window data ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    std::istringstream issDataNames(datNames);
    std::vector<std::string> vecDataNames;
    for (int item=0; item<nItemps; ++item)
    {
        std::string dataName="";

        issDataNames >> dataName;

        if (dataName != string("object"))
        {
            vecDataNames.push_back(dataName);

            std::string dataName1 = winName1+string(".")+dataName;
            std::string dataName2 = winName2+string(".")+dataName;

            COM_clone_dataitem(dataName2.c_str(),
                               dataName1.c_str());

            //COM_clone_dataitem(winName2.c_str(), winName1.c_str(), 0, dataName.c_str(), 0);
            //COM_clone_dataitem(winName2, winName1, 0, dataName);

            Info << "rocFoam.copyWindow: dataName["
                 << item << "] = " << dataName
                 << ", copied from " << winName1.c_str()
                 << " to " << winName2.c_str() << endl;
         }
         else
         {
                Info << "rocFoam.copyWindow: dataName["
                 << item << "] = " << dataName
                 << ", skipped " << endl;
         }
         
    }
    Info << endl;
    // ------------------------------------------

    COM_window_init_done(winName2); 
}


bool comFoam::nameExists(const std::vector<std::string>& dataItemNames,
                          const std::string& dataName)
{
    bool status = (std::find(dataItemNames.begin(),
                   dataItemNames.end(), dataName)
                   != dataItemNames.end());

    if (! status)
    {
        std::cout << std::endl << " Warning: dataItemName " << dataName
                  << " does not exist." << std::endl << std::endl;
    }

    return status;
}


int comFoam::reconstMesh(int argc, char **argv, const char* name)
{
    // Form point data container ^^^^^^^^^^^^^^^^
    Foam::pointField of_Points(*ca_nPoints);

    for(int ipoint=0; ipoint<*ca_nPoints; ipoint++)
    {
        
        Foam::point pointTmp;
        for(int jcomp=0; jcomp<nComponents; jcomp++)
        {
            pointTmp[jcomp] = *(ca_Points+ipoint*nComponents+jcomp);

        }
        of_Points[ipoint] = pointTmp;
    }
    //-------------------------------------------


    // Form face data container ^^^^^^^^^^^^^^^^^
    Foam::faceList of_Faces(*ca_nFaces);
    
    Foam::labelList of_Owner(*ca_nFaces);
    Foam::labelList of_Neighb(*ca_nFaces);
    
    int faceIndex = 0;
    for(int itype=0; itype<*ca_faceToPointConn_types; itype++)
    {
        int nPoints = ca_faceToPointConn_map[itype];
        int nFaces = ca_faceToPointConn_size[itype];
        Foam::labelList pointsList(nPoints);
        for(int iface=0; iface<nFaces; iface++)
        {

            int faceID = ca_faceToFaceMap[faceIndex];
            for(int ipoint=0; ipoint<nPoints; ipoint++)
            {
                int index = ipoint+iface*nPoints;
                pointsList[ipoint] = ca_faceToPointConn[itype][index];
            }            
            
            of_Faces[faceID] = face(pointsList);

            of_Owner[faceID] = ca_faceOwner[faceIndex];
            of_Neighb[faceID] = ca_faceNeighb[faceIndex];
            
            faceIndex++;
        }
    }
    //-------------------------------------------
/*for(int iface=0; iface<*ca_nFaces; iface++)*/
/*{*/
/*    Info << of_Owner[iface] << " " << of_Neighb[iface] << endl;*/
/*}*/

Info << " HELLO 0" << endl;


    createArgs(argc, argv);

Info << " HELLO 1" << endl;

    PostProcess(argc, argv);

Info << " HELLO 10" << endl;

    //  setRootCaseLists.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    setRootCaseLists();
    // ---------------------------------------------------

Info << " HELLO 2" << endl;

    //  createTime.H  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    createTime();
    // ---------------------------------------------------

Info << " HELLO 3" << endl;

    setRootCase();

    Foam::Time &runTime(*runTimePtr);
    Info << "Create mesh for time = " << runTime.timeName() << nl << endl;

    const bool syncPar = true;

Info << " HELLO 4" << endl;
    
    Foam::polyMesh meshTmp = Foam::polyMesh
    (
        IOobject
        (
            dynamicFvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            IOobject::NO_READ
        ),
        Foam::move(of_Points),
        Foam::move(of_Faces),
        Foam::move(of_Owner),
        Foam::move(of_Neighb),
        syncPar
    );


    meshPtr = dynamicFvMesh::New(meshTmp);
    dynamicFvMesh &mesh(*meshPtr);

Info << " HELLO 5" << endl;

    int nPatches = *ca_nPatches;


    Foam::polyBoundaryMesh of_polyBoundaryMesh
    (
        IOobject
        (
            "polyBoundaryMesh",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        nPatches
    );


Info << " HELLO 6" << endl;

    std::vector<polyPatch> of_vecPatches;
    
    List<polyPatch*> of_Patches(nPatches);
    List<polyPatch*> of_Patches_(nPatches);

    for(int ipatch=0; ipatch<nPatches; ipatch++)
    {
        word patchName = patchNameStr[ipatch].c_str();
        word patchType = patchTypeStr[ipatch].c_str();
                
        label patchSize = *ca_patchSize[ipatch];
        label patchStart = *ca_patchStart[ipatch];
        label patchIndex = ipatch;
        
        
        Foam::polyPatch of_Patch
        (
            patchName,
            patchSize,
            patchStart,
            patchIndex,
            of_polyBoundaryMesh,
            patchType
        );
        

        //of_vecPatches.push_back(of_Patch);


        //const polyPatch& tmpPatch=of_vecPatches[ipatch];

        //of_Patches[ipatch] = tmpPatch.ptr();

        
        /*of_Patches[ipatch] = new polyPatch
        (
            patchName,
            patchSize,
            patchStart,
            patchIndex,
            of_polyBoundaryMesh,
            patchType
        );
        
        //of_Patches_[ipatch] = of_Patches[ipatch];
        */
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    for(int ipatch=0; ipatch<nPatches; ipatch++)
    {
        Info << " patchName =  " << patches[ipatch].name() << endl;
        Info << " patchType =  " << patches[ipatch].type() << endl;
        Info << " patchStrat =  " << patches[ipatch].start() << endl;
        Info << " patchSize =  " << patches[ipatch].size() << endl;
    }


    return 0;
}


