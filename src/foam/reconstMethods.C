void comFoam::copyWindow(const string& winName1, const string& winName2)
{
    std::cout << "rocFoam.copyWindow: window names "
              << winName1.c_str() << " & "
              << winName2.c_str() << std::endl;

    // Get window information ^^^^^^^^^^^^^^^^^^^
    int nItemps;
    std::string datNames;
    COM_get_dataitems(winName1.c_str(), &nItemps, datNames);

    std::cout << "rocFoam.copyWindow: nItemps in "
              << "winName1 = "
              << nItemps
              << std::endl;
    // ------------------------------------------

    //  Copy pane nodes & connectivities ^^^^^^^^
    int nPanes;
    int* paneList;
    COM_get_panes(winName1.c_str(), &nPanes, &paneList);
    std::cout << "rocFoam.copyWindow: Number of Panes in window "
              << winName1 << " = "
              << nPanes << std::endl;

    //for(int ipane=0; ipane<nPanes; ipane++)
    //{
        COM_clone_dataitem((winName2+".mesh").c_str(),
                           (winName1+".mesh").c_str(),
                           0);
    //}
    // ------------------------------------------

    //  Window data ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    std::istringstream issDataNames(datNames);
    std::vector<std::string> vecDataNames;
    int fileCount = 0;
    for (int item=0; item<nItemps; ++item)
    {
        std::string dataName="";

        issDataNames >> dataName;

        if (dataName == std::string("object"))
        {
            std::cout << "rocFoam.copyWindow: Cannot copy dataName" 
                      << " \"object\" in window " << winName1
                      << std::endl;
            exit(-1);
        }
        else
        {
            vecDataNames.push_back(dataName);

            std::string dataName1 = winName1+std::string(".")+dataName;
            std::string dataName2 = winName2+std::string(".")+dataName;

            COM_clone_dataitem(dataName2.c_str(),
                               dataName1.c_str(),
                               0);

            std::string subName = dataName.substr(0,4);
            if (subName != "file")
            {
                std::cout << "rocFoam.copyWindow: dataName["
                          << item << "] = " << dataName
                          << " copied from " << winName1.c_str()
                          << " to " << winName2.c_str()
                          << std::endl;
            }
            else
            {
                fileCount++;
            }
         }
    }

    if (fileCount>0)
    {
        std::cout << "rocFoam.copyWindow: Total number of "
        << fileCount << " file-related items copied from "
        << winName1.c_str()
        << " to " << winName2.c_str() << std::endl;
    }
    
    std::cout << std::endl;
    // ------------------------------------------

    COM_window_init_done(winName2);
}


bool comFoam::nameExists(const std::vector<std::string>& dataItemNames,
                          const std::string& dataName)
{
    bool status = (std::find(dataItemNames.begin(),
                   dataItemNames.end(), dataName)
                   != dataItemNames.end());

    if (!status)
    {
        std::cout << ">>DataItemName " << dataName
                  << " does not exist.<<" << std::endl;
    }

    return status;
}

/*
int comFoam::reconstDynamicFvMesh()
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

    // Creating mesh objectt ^^^^^^^^^^^^^^^^^^^^
    Foam::Time &runTime(*runTimePtr);
    Info << "Create mesh for time = " << runTime.timeName() << nl << endl;

    meshPtr = dynamicFvMesh::New
    (
        IOobject
        (
            fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            IOobject::READ_IF_PRESENT
        ),
        Foam::move(of_Points),
        Foam::move(of_Faces),
        Foam::move(of_Owner),
        Foam::move(of_Neighb)
    );
    dynamicFvMesh& mesh(*meshPtr);
    //-------------------------------------------

    // Creating patches ^^^^^^^^^^^^^^^^^^^^^^^^^
    int nPatches = *ca_nPatches;
    Foam::polyBoundaryMesh of_polyBoundaryMesh
    (
        IOobject
        (
            "polyBoundaryMesh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        nPatches
    );
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
    }

    List<polyPatch*> of_Patches(*ca_nPatches);
    std::vector<polyPatch> of_vecPatches;
    for(int ipatch=0; ipatch<nPatches; ipatch++)
    {
        word patchName = patchNameStr[ipatch].c_str();
        word patchType = patchTypeStr[ipatch].c_str();
                
        label patchSize = *ca_patchSize[ipatch];
        label patchStart = *ca_patchStart[ipatch];
        label patchIndex = ipatch;

        of_Patches[ipatch] = new polyPatch
        (
            patchName,
            patchSize,
            patchStart,
            patchIndex,
            of_polyBoundaryMesh,
            patchType
        );
    }
    mesh.addPatches(of_Patches);
    //-------------------------------------------

    return 0;
}
*/
