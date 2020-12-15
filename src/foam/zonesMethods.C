#include "comFoam.H"

int comFoam::createZonesData()
{
    const dynamicFvMesh& mesh(*meshPtr);

    // cellZoness ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    const cellZoneMesh& cellZones = mesh.cellZones();
    const int& nCellZones = cellZones.size();
    
    ca_nCellZones = new int(nCellZones);

    if (nCellZones>0)
    {
        if (cellZonesTypeStr == nullptr)
            cellZonesTypeStr = new std::string[nCellZones]{};
        
        if (cellZonesNameStr == nullptr)
            cellZonesNameStr = new std::string[nCellZones]{};

        int maxTypeLength{0};
        int maxNameLength{0};
        forAll(cellZones, izone)
        {
            cellZonesTypeStr[izone] = cellZones[izone].type();
            std::string strTmp(cellZonesTypeStr[izone]);
            maxTypeLength = std::max( maxTypeLength, static_cast<int>(strTmp.length())+1 );

            cellZonesNameStr[izone] = cellZones[izone].name();
            strTmp = cellZonesNameStr[izone];
            maxNameLength = std::max( maxNameLength, static_cast<int>(strTmp.length())+1 );
        }

        ca_cellZonesTypeMaxLength = new int(maxTypeLength);
        ca_cellZonesNameMaxLength = new int(maxNameLength);

        // zone Types; CHAR type
        ca_cellZonesType = new char[nCellZones * maxTypeLength]{};
        forAll(cellZones, izone)
        {
            int startIndex = izone * maxTypeLength;
            int endIndex = startIndex + cellZonesTypeStr[izone].length();
            int count{0};
            for (int i=startIndex; i<endIndex; i++)
            {
                ca_cellZonesType[i] = cellZonesTypeStr[izone][count];
                count++;
            }
            ca_cellZonesType[endIndex] = '\0';

            for (int i=endIndex+1; i<startIndex+maxTypeLength; i++)
            {
                ca_cellZonesType[i] = ' ';
            }
        }

        // name Names; CHAR type
        ca_cellZonesName = new char[nCellZones * maxNameLength]{};
        forAll(cellZones, izone)
        {
            int startIndex = izone * maxNameLength;
            int endIndex = startIndex + cellZonesNameStr[izone].length();
            int count{0};
            for (int i=startIndex; i<endIndex; i++)
            {
                ca_cellZonesName[i] = cellZonesNameStr[izone][count];
                count++;
            }
            ca_cellZonesName[endIndex] = '\0';

            for (int i=endIndex+1; i<startIndex+maxNameLength; i++)
            {
                ca_cellZonesName[i] = ' ';
            }
        }

        ca_cellZonesCount = new int[nCellZones]{};
        ca_cellZonesList = new int*[nCellZones]{};

        forAll(cellZones, izone)
        {
            const labelList& cellList = cellZones[izone];
            int nCells = cellList.size();
            ca_cellZonesCount[izone] = nCells;

            if (nCells>0)
            {
                ca_cellZonesList[izone] = new int[nCells]{};
                
                forAll(cellList, icell)
                {
                    ca_cellZonesList[izone][icell] = cellList[icell];
                }
            }
        }
    }
    //-------------------------------------------

    // faceZoness ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    const faceZoneMesh& faceZones = mesh.faceZones();
    const int& nFaceZones = faceZones.size();
    
    ca_nFaceZones = new int(nFaceZones);

    if (nFaceZones>0)
    {
        if (faceZonesTypeStr == nullptr)
            faceZonesTypeStr = new std::string[nFaceZones]{};
        
        if (faceZonesNameStr == nullptr)
            faceZonesNameStr = new std::string[nFaceZones]{};

        int maxTypeLength{0};
        int maxNameLength{0};
        forAll(faceZones, izone)
        {
            faceZonesTypeStr[izone] = faceZones[izone].type();
            std::string strTmp(faceZonesTypeStr[izone]);
            maxTypeLength = std::max( maxTypeLength, static_cast<int>(strTmp.length())+1 );

            faceZonesNameStr[izone] = faceZones[izone].name();
            strTmp = faceZonesNameStr[izone];
            maxNameLength = std::max( maxNameLength, static_cast<int>(strTmp.length())+1 );
        }

        ca_faceZonesTypeMaxLength = new int(maxTypeLength);
        ca_faceZonesNameMaxLength = new int(maxNameLength);

        // faceZones Type; CHAR type
        ca_faceZonesType = new char[nFaceZones * maxTypeLength]{};
        forAll(faceZones, izone)
        {
            int startIndex = izone * maxTypeLength;
            int endIndex = startIndex + faceZonesTypeStr[izone].length();
            int count{0};
            for (int i=startIndex; i<endIndex; i++)
            {
                ca_faceZonesType[i] = faceZonesTypeStr[izone][count];
                count++;
            }
            ca_faceZonesType[endIndex] = '\0';

            for (int i=endIndex+1; i<startIndex+maxTypeLength; i++)
            {
                ca_faceZonesType[i] = ' ';
            }
        }

        // faceZones Name; CHAR type
        ca_faceZonesName = new char[nFaceZones * maxNameLength]{};
        forAll(faceZones, izone)
        {
            int startIndex = izone * maxNameLength;
            int endIndex = startIndex + faceZonesNameStr[izone].length();
            int count{0};
            for (int i=startIndex; i<endIndex; i++)
            {
                ca_faceZonesName[i] = faceZonesNameStr[izone][count];
                count++;
            }
            ca_faceZonesName[endIndex] = '\0';

            for (int i=endIndex+1; i<startIndex+maxNameLength; i++)
            {
                ca_faceZonesName[i] = ' ';
            }
        }

        ca_faceZonesCount = new int[nFaceZones]{};
        ca_faceZonesList = new int*[nFaceZones]{};
        ca_faceZonesFlipMap = new int*[nFaceZones]{};

        forAll(faceZones, izone)
        {
            const labelList& faceList = faceZones[izone];
            const boolList&  flipMap = faceZones[izone].flipMap();

            int nFaces = faceList.size();
            ca_faceZonesCount[izone] = nFaces;

            if (nFaces>0)
            {
                ca_faceZonesList[izone] = new int[nFaces]{};
                ca_faceZonesFlipMap[izone] = new int[nFaces]{};
                
                forAll(faceList, iface)
                {
                    ca_faceZonesList[izone][iface] = faceList[iface];
                    ca_faceZonesFlipMap[izone][iface] = flipMap[iface];
                }
            }
        }
    }
    //-------------------------------------------

    // pointZoness ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    const pointZoneMesh& pointZones = mesh.pointZones();
    const int& nPointZones = pointZones.size();
    
    ca_nPointZones = new int(nPointZones);

    if (nPointZones>0)
    {
        if (pointZonesTypeStr == nullptr)
            pointZonesTypeStr = new std::string[nPointZones]{};
        
        if (pointZonesNameStr == nullptr)
            pointZonesNameStr = new std::string[nPointZones]{};

        int maxTypeLength{0};
        int maxNameLength{0};
        forAll(pointZones, izone)
        {
            pointZonesTypeStr[izone] = pointZones[izone].type();
            std::string strTmp(pointZonesTypeStr[izone]);
            maxTypeLength = std::max( maxTypeLength, static_cast<int>(strTmp.length())+1 );

            pointZonesNameStr[izone] = pointZones[izone].name();
            strTmp = pointZonesNameStr[izone];
            maxNameLength = std::max( maxNameLength, static_cast<int>(strTmp.length())+1 );
        }

        ca_pointZonesTypeMaxLength = new int(maxTypeLength);
        ca_pointZonesNameMaxLength = new int(maxNameLength);

        // pointZones Type; CHAR type
        ca_pointZonesType = new char[nPointZones * maxTypeLength]{};
        forAll(pointZones, izone)
        {
            int startIndex = izone * maxTypeLength;
            int endIndex = startIndex + pointZonesTypeStr[izone].length();
            int count{0};
            for (int i=startIndex; i<endIndex; i++)
            {
                ca_pointZonesType[i] = pointZonesTypeStr[izone][count];
                count++;
            }
            ca_pointZonesType[endIndex] = '\0';

            for (int i=endIndex+1; i<startIndex+maxTypeLength; i++)
            {
                ca_pointZonesType[i] = ' ';
            }
        }

        // pointZones Name; CHAR type
        ca_pointZonesName = new char[nPointZones * maxNameLength]{};
        forAll(pointZones, izone)
        {
            int startIndex = izone * maxNameLength;
            int endIndex = startIndex + pointZonesNameStr[izone].length();
            int count{0};
            for (int i=startIndex; i<endIndex; i++)
            {
                ca_pointZonesName[i] = pointZonesNameStr[izone][count];
                count++;
            }
            ca_pointZonesName[endIndex] = '\0';

            for (int i=endIndex+1; i<startIndex+maxNameLength; i++)
            {
                ca_pointZonesName[i] = ' ';
            }
        }

        ca_pointZonesCount = new int[nPointZones]{};
        ca_pointZonesList = new int*[nPointZones]{};

        forAll(pointZones, izone)
        {
            const labelList& pointList = pointZones[izone];
            int nPoints = pointList.size();
            ca_pointZonesCount[izone] = nPoints;

            if (nPoints>0)
            {
                ca_pointZonesList[izone] = new int[nPoints]{};
                
                forAll(pointList, ipoint)
                {
                    ca_pointZonesList[izone][ipoint] = pointList[ipoint];
                }
            }
        }
    }
    //-------------------------------------------


    return 0;
}

int comFoam::registerZonesData(const char *name)
{
    std::string volName = name+std::string("VOL");
    
    std::stringstream output{};
    output << "rocFoam.registerZonesData: "
         << "Registering zones data with window name "
         << volName;
    verbose_message(output.str(), true);

    // grid and field data
    int paneID = Pstream::myProcNo()+1;// Use this paneID for volume connectivity

    output = std::stringstream{};
    output << "procID = " << ca_myRank
         << ", paneID = " << paneID
         << " ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^";
    verbose_message(output.str(), true);

    // cellZones ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    std::string dataName = volName+std::string(".nCellZones");
    COM_new_dataitem( dataName, 'p', COM_INT, 1, "");
    COM_set_size(     dataName, paneID, 1);
    COM_set_array(    dataName, paneID, ca_nCellZones);
    output = std::stringstream{};
    output << "  " << dataName.c_str() << " registered.";
    verbose_message(output.str(), true);

    int nCellZones = *ca_nCellZones;
    if (nCellZones>0)
    {
        int cellZonesTypeTotalSize{ *ca_cellZonesTypeMaxLength * nCellZones};
        dataName = volName+std::string(".cellZonesType");
        COM_new_dataitem( dataName, 'p', COM_CHAR, 1, "");
        COM_set_size( dataName, paneID, cellZonesTypeTotalSize);
        COM_set_array(dataName, paneID, ca_cellZonesType);
        output = std::stringstream{};
        output << "  " << dataName.c_str() << " registered.";
        verbose_message(output.str(), true);

        int cellZonesNameTotalSize{ *ca_cellZonesNameMaxLength * nCellZones};
        dataName = volName+std::string(".cellZonesName");
        COM_new_dataitem( dataName, 'p', COM_CHAR, 1, "");
        COM_set_size( dataName, paneID, cellZonesNameTotalSize);
        COM_set_array(dataName, paneID, ca_cellZonesName);
        output = std::stringstream{};
        output << "  " << dataName.c_str() << " registered.";
        verbose_message(output.str(), true);

        dataName = volName+std::string(".cellZonesTypeMaxLength");
        COM_new_dataitem( dataName, 'p', COM_INT, 1, "");
        COM_set_size(     dataName, paneID, 1);
        COM_set_array(    dataName, paneID, ca_cellZonesTypeMaxLength);
        output = std::stringstream{};
        output << "  " << dataName.c_str() << " registered.";
        verbose_message(output.str(), true);

        dataName = volName+std::string(".cellZonesNameMaxLength");
        COM_new_dataitem( dataName, 'p', COM_INT, 1, "");
        COM_set_size(     dataName, paneID, 1);
        COM_set_array(    dataName, paneID, ca_cellZonesNameMaxLength);
        output = std::stringstream{};
        output << "  " << dataName.c_str() << " registered.";
        verbose_message(output.str(), true);

        dataName = volName+std::string(".cellZonesCount");
        COM_new_dataitem( dataName, 'p', COM_INT, 1, "");
        COM_set_size(     dataName, paneID, nCellZones);
        COM_set_array(    dataName, paneID, ca_cellZonesCount);
        output = std::stringstream{};
        output << "  " << dataName.c_str() << " registered.";
        verbose_message(output.str(), true);

        for (int izone=0; izone<nCellZones; izone++)
        {
            if (ca_cellZonesCount[izone]>0)
            {
                dataName = volName+".cellZonesList"+std::to_string(izone);
                COM_new_dataitem( dataName, 'p', COM_INT, 1, "");
                COM_set_size(     dataName, paneID, ca_cellZonesCount[izone]);
                COM_set_array(    dataName, paneID, ca_cellZonesList[izone]);
                output = std::stringstream{};
                output << "  " << dataName.c_str() << " registered.";
                verbose_message(output.str(), true);
            }
        }
    }
    //-------------------------------------------

    // faceZones ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    dataName = volName+std::string(".nFaceZones");
    COM_new_dataitem( dataName, 'p', COM_INT, 1, "");
    COM_set_size(     dataName, paneID, 1);
    COM_set_array(    dataName, paneID, ca_nFaceZones);
    output = std::stringstream{};
    output << "  " << dataName.c_str() << " registered.";
    verbose_message(output.str(), true);

    int nFaceZones = *ca_nFaceZones;
    if (nFaceZones>0)
    {
        int faceZonesTypeTotalSize{ *ca_faceZonesTypeMaxLength * nFaceZones};
        dataName = volName+std::string(".faceZonesType");
        COM_new_dataitem( dataName, 'p', COM_CHAR, 1, "");
        COM_set_size( dataName, paneID,faceZonesTypeTotalSize);
        COM_set_array(dataName, paneID,ca_faceZonesType);
        output = std::stringstream{};
        output << "  " << dataName.c_str() << " registered.";
        verbose_message(output.str(), true);

        int faceZonesNameTotalSize{ *ca_faceZonesNameMaxLength * nFaceZones};
        dataName = volName+std::string(".faceZonesName");
        COM_new_dataitem( dataName, 'p', COM_CHAR, 1, "");
        COM_set_size( dataName, paneID,faceZonesNameTotalSize);
        COM_set_array(dataName, paneID,ca_faceZonesName);
        output = std::stringstream{};
        output << "  " << dataName.c_str() << " registered.";
        verbose_message(output.str(), true);

        dataName = volName+std::string(".faceZonesTypeMaxLength");
        COM_new_dataitem( dataName, 'p', COM_INT, 1, "");
        COM_set_size(     dataName, paneID, 1);
        COM_set_array(    dataName, paneID, ca_faceZonesTypeMaxLength);
        output = std::stringstream{};
        output << "  " << dataName.c_str() << " registered.";
        verbose_message(output.str(), true);

        dataName = volName+std::string(".faceZonesNameMaxLength");
        COM_new_dataitem( dataName, 'p', COM_INT, 1, "");
        COM_set_size(     dataName, paneID, 1);
        COM_set_array(    dataName, paneID, ca_faceZonesNameMaxLength);
        output = std::stringstream{};
        output << "  " << dataName.c_str() << " registered.";
        verbose_message(output.str(), true);

        dataName = volName+std::string(".faceZonesCount");
        COM_new_dataitem( dataName, 'p', COM_INT, 1, "");
        COM_set_size(     dataName, paneID, nFaceZones);
        COM_set_array(    dataName, paneID, ca_faceZonesCount);
        output = std::stringstream{};
        output << "  " << dataName.c_str() << " registered.";
        verbose_message(output.str(), true);

        for (int izone=0; izone<nFaceZones; izone++)
        {
            if (ca_faceZonesCount[izone]>0)
            {
                dataName = volName+".faceZonesList"+std::to_string(izone);
                COM_new_dataitem( dataName, 'p', COM_INT, 1, "");
                COM_set_size(     dataName, paneID, ca_faceZonesCount[izone]);
                COM_set_array(    dataName, paneID, ca_faceZonesList[izone]);
                output = std::stringstream{};
                output << "  " << dataName.c_str() << " registered.";
                verbose_message(output.str(), true);

                dataName = volName+".faceZonesFlipMap"+std::to_string(izone);
                COM_new_dataitem( dataName, 'p', COM_INT, 1, "");
                COM_set_size(     dataName, paneID, ca_faceZonesCount[izone]);
                COM_set_array(    dataName, paneID, ca_faceZonesFlipMap[izone]);
                output = std::stringstream{};
                output << "  " << dataName.c_str() << " registered.";
                verbose_message(output.str(), true);
            }
        }
    }
    //-------------------------------------------

    // pointZones ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    dataName = volName+std::string(".nPointZones");
    COM_new_dataitem( dataName, 'p', COM_INT, 1, "");
    COM_set_size(     dataName, paneID, 1);
    COM_set_array(    dataName, paneID, ca_nPointZones);
    output = std::stringstream{};
    output << "  " << dataName.c_str() << " registered.";
    verbose_message(output.str(), true);

    int nPointZones = *ca_nPointZones;
    if (nPointZones>0)
    {
        int pointZonesTypeTotalSize{ *ca_pointZonesTypeMaxLength * nPointZones};
        dataName = volName+std::string(".pointZonesType");
        COM_new_dataitem( dataName, 'p', COM_CHAR, 1, "");
        COM_set_size( dataName, paneID,pointZonesTypeTotalSize);
        COM_set_array(dataName, paneID,ca_pointZonesType);
        output = std::stringstream{};
        output << "  " << dataName.c_str() << " registered.";
        verbose_message(output.str(), true);

        int pointZonesNameTotalSize{ *ca_pointZonesNameMaxLength * nPointZones};
        dataName = volName+std::string(".pointZonesName");
        COM_new_dataitem( dataName, 'p', COM_CHAR, 1, "");
        COM_set_size( dataName, paneID,pointZonesNameTotalSize);
        COM_set_array(dataName, paneID,ca_pointZonesName);
        output = std::stringstream{};
        output << "  " << dataName.c_str() << " registered.";
        verbose_message(output.str(), true);

        dataName = volName+std::string(".pointZonesTypeMaxLength");
        COM_new_dataitem( dataName, 'p', COM_INT, 1, "");
        COM_set_size(     dataName, paneID, 1);
        COM_set_array(    dataName, paneID, ca_pointZonesTypeMaxLength);
        output = std::stringstream{};
        output << "  " << dataName.c_str() << " registered.";
        verbose_message(output.str(), true);

        dataName = volName+std::string(".pointZonesNameMaxLength");
        COM_new_dataitem( dataName, 'p', COM_INT, 1, "");
        COM_set_size(     dataName, paneID, 1);
        COM_set_array(    dataName, paneID, ca_pointZonesNameMaxLength);
        output = std::stringstream{};
        output << "  " << dataName.c_str() << " registered.";
        verbose_message(output.str(), true);

        dataName = volName+std::string(".pointZonesCount");
        COM_new_dataitem( dataName, 'p', COM_INT, 1, "");
        COM_set_size(     dataName, paneID, nPointZones);
        COM_set_array(    dataName, paneID, ca_pointZonesCount);
        output = std::stringstream{};
        output << "  " << dataName.c_str() << " registered.";
        verbose_message(output.str(), true);

        for (int izone=0; izone<nPointZones; izone++)
        {
            if (ca_cellZonesCount[izone]>0)
            {
                dataName = volName+".pointZonesList"+std::to_string(izone);
                COM_new_dataitem( dataName, 'p', COM_INT, 1, "");
                COM_set_size(     dataName, paneID, ca_pointZonesCount[izone]);
                COM_set_array(    dataName, paneID, ca_pointZonesList[izone]);
                output = std::stringstream{};
                output << "  " << dataName.c_str() << " registered.";
                verbose_message(output.str(), true);
            }
        }
    }
    //-------------------------------------------

    return 0;
}

int comFoam::reconstZonesData(const char *name)
{
    std::string volName = name+std::string("VOL");

    std::stringstream output{};
    output << "rocFoam.reconstZonesData, procID = "
           << ca_myRank
           << ", Retreiving zones data form window "
           << volName << ".";
    verbose_message(output.str(), true);

    std::string regNames;
    int numDataItems=0;
    COM_get_dataitems(volName.c_str(), &numDataItems, regNames);

    std::vector<std::string> dataItemNames;
    dataItemNames.clear();
    std::istringstream Istr(regNames);
    for (int i=0; i<numDataItems; ++i)
    {
        std::string nameTmp;
        Istr >> nameTmp;

        std::string subName1 = nameTmp.substr(0,13);
        std::string subName2 = nameTmp.substr(0,14);
        std::string subName3 = nameTmp.substr(0,16);

        if
        (
            nameTmp  == "nCellZones" ||
            nameTmp  == "cellZonesTypeMaxLength" ||
            nameTmp  == "cellZonesNameMaxLength" ||
            nameTmp  == "cellZonesType" ||
            nameTmp  == "cellZonesName" ||
            nameTmp  == "cellZonesCount" ||
            subName1 == "cellZonesList" ||

            nameTmp  == "nFaceZones" ||
            nameTmp  == "faceZonesTypeMaxLength" ||
            nameTmp  == "faceZonesNameMaxLength" ||
            nameTmp  == "faceZonesType" ||
            nameTmp  == "faceZonesName" ||
            nameTmp  == "faceZonesCount" ||
            subName1 == "faceZonesList" ||
            subName3 == "faceZonesFlipMap" ||

            nameTmp  == "nPointZones" ||
            nameTmp  == "pointZonesTypeMaxLength" ||
            nameTmp  == "pointZonesNameMaxLength" ||
            nameTmp  == "pointZonesType" ||
            nameTmp  == "pointZonesName" ||
            nameTmp  == "pointZonesCount" ||
            subName2 == "pointZonesList"
        )
        {
            dataItemNames.push_back(nameTmp);

            output = std::stringstream{};
            output << "  DataItem[" << i << "] = " << nameTmp;
            verbose_message(output.str(), true);
        }
    }
    output = std::stringstream{};
    output << "  Number of items = " << dataItemNames.size()
           << std::endl;
    verbose_message(output.str(), true);

    int paneID = ca_myRank+1;// Use this paneID for volume connectivity
    output = std::stringstream{};
    output << "  procID = " << ca_myRank
         << ", paneID = " << paneID
         << " ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^";
    verbose_message(output.str(), true);

    std::string dataName{};
    std::string regName{};

    // cellZone data ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    dataName = std::string("nCellZones");
    if (nameExists(dataItemNames, dataName))
    {
        regName = volName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_nCellZones);
        output = std::stringstream{};
        output << "  " << dataName.c_str() << " = " << *ca_nCellZones;
        verbose_message(output.str(), true);
    }
    else
    {
        ca_nCellZones = new int(0);
    }
    int nCellZones{*ca_nCellZones};

    if (nCellZones>0)
    {
        // cellZonesType
        dataName = std::string("cellZonesTypeMaxLength");
        if (nameExists(dataItemNames, dataName))
        {
            regName = volName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_cellZonesTypeMaxLength);
            output = std::stringstream{};
            output << "  " << dataName.c_str() << " = " << *ca_cellZonesTypeMaxLength;
            verbose_message(output.str(), true);
        }

        int nComp{};
        dataName = std::string("cellZonesType");
        if(nameExists(dataItemNames, dataName))
        {
            regName = volName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_cellZonesType);
            COM_get_size(regName.c_str(), paneID, &nComp);

            int cellZonesTypeTotalSize{nCellZones * *ca_cellZonesTypeMaxLength};
            compareWarningExit(cellZonesTypeTotalSize, nComp,
                                "cellZonesTypeTotalSize", "nComp");

            cellZonesTypeStr = new std::string[nCellZones]{};
            for (int i=0; i<nCellZones; i++)
            {
                int startIndex = i * *ca_cellZonesTypeMaxLength;
                char* charTmp = &ca_cellZonesType[startIndex];
                cellZonesTypeStr[i] = charTmp;
                
                output = std::stringstream{};
                output << "  " << dataName.c_str() << "[" << i << "] = "
                    << charTmp << ", " << cellZonesTypeStr[i];
                verbose_message(output.str(), true);
            }
        }

        // cellZonesName
        dataName = std::string("cellZonesNameMaxLength");
        if (nameExists(dataItemNames, dataName))
        {
            regName = volName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_cellZonesNameMaxLength);
            output = std::stringstream{};
            output << "  " << dataName.c_str() << " = " << *ca_cellZonesNameMaxLength;
            verbose_message(output.str(), true);
        }

        dataName = std::string("cellZonesName");
        if(nameExists(dataItemNames, dataName))
        {
            regName = volName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_cellZonesName);
            COM_get_size(regName.c_str(), paneID, &nComp);

            int cellZonesNameTotalSize{nCellZones * *ca_cellZonesNameMaxLength};
            compareWarningExit(cellZonesNameTotalSize, nComp,
                                "cellZonesNameTotalSize", "nComp");

            cellZonesNameStr = new std::string[nCellZones]{};
            for (int i=0; i<nCellZones; i++)
            {
                int startIndex = i * *ca_cellZonesNameMaxLength;
                char* charTmp = &ca_cellZonesName[startIndex];
                cellZonesNameStr[i] = charTmp;
                
                output = std::stringstream{};
                output << "  " << dataName.c_str() << "[" << i << "] = "
                    << charTmp << ", " << cellZonesNameStr[i];
                verbose_message(output.str(), true);
            }
        }

        dataName = std::string("cellZonesCount");
        if(nameExists(dataItemNames, dataName))
        {
            regName = volName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_cellZonesCount);
            COM_get_size(regName.c_str(), paneID, &nComp);
            for(int icomp=0; icomp<nComp; icomp++)
            {
                output = std::stringstream{};
                output << "  " << dataName.c_str() << "[" << icomp << "] = "
                    << ca_cellZonesCount[icomp];
                verbose_message(output.str(), true);
            }
            compareWarningExit(nCellZones, nComp,
                                "nCellZones", "cellZonesCount");
        }

        ca_cellZonesList = new int*[nCellZones]{};
        for (int izone=0; izone<nCellZones; izone++)
        {
            dataName = "cellZonesList"+std::to_string(izone);
            if(nameExists(dataItemNames, dataName))
            {
                regName = volName+".cellZonesList"+std::to_string(izone);
                COM_get_array(regName.c_str(), paneID, &ca_cellZonesList[izone]);
                COM_get_size(regName.c_str(), paneID, &nComp);
                output = std::stringstream{};
                output << "  " << dataName.c_str() << " has " << nComp
                    << " elements";
                verbose_message(output.str(), true);
            }
            else
            {
                output = std::stringstream{};
                output << "  " << dataName.c_str() << " is a nullptr.";
                verbose_message(output.str(), true);
            }
        }
    }
    //-------------------------------------------

    // faceZone data ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    dataName = std::string("nFaceZones");
    if (nameExists(dataItemNames, dataName))
    {
        regName = volName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_nFaceZones);
        output = std::stringstream{};
        output << "  " << dataName.c_str() << " = " << *ca_nFaceZones;
        verbose_message(output.str(), true);
    }
    else
    {
        ca_nFaceZones = new int(0);
    }
    int nFaceZones{*ca_nFaceZones};

    if (nFaceZones>0)
    {
        // faceZonesType
        dataName = std::string("faceZonesTypeMaxLength");
        if (nameExists(dataItemNames, dataName))
        {
            regName = volName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_faceZonesTypeMaxLength);
            output = std::stringstream{};
            output << "  " << dataName.c_str() << " = " << *ca_faceZonesTypeMaxLength;
            verbose_message(output.str(), true);
        }

        int nComp{};
        dataName = std::string("faceZonesType");
        if(nameExists(dataItemNames, dataName))
        {
            regName = volName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_faceZonesType);
            COM_get_size(regName.c_str(), paneID, &nComp);

            int faceZonesTypeTotalSize{nFaceZones * *ca_faceZonesTypeMaxLength};
            compareWarningExit(faceZonesTypeTotalSize, nComp,
                                "faceZonesTypeTotalSize", "nComp");

            faceZonesTypeStr = new std::string[nFaceZones]{};
            for (int i=0; i<nFaceZones; i++)
            {
                int startIndex = i * *ca_faceZonesTypeMaxLength;
                char* charTmp = &ca_faceZonesType[startIndex];
                faceZonesTypeStr[i] = charTmp;
                
                output = std::stringstream{};
                output << "  " << dataName.c_str() << "[" << i << "] = "
                    << charTmp << ", " << faceZonesTypeStr[i];
                verbose_message(output.str(), true);
            }
        }

        // faceZonesName
        dataName = std::string("faceZonesNameMaxLength");
        if (nameExists(dataItemNames, dataName))
        {
            regName = volName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_faceZonesNameMaxLength);
            output = std::stringstream{};
            output << "  " << dataName.c_str() << " = " << *ca_faceZonesNameMaxLength;
            verbose_message(output.str(), true);
        }

        dataName = std::string("faceZonesName");
        if(nameExists(dataItemNames, dataName))
        {
            regName = volName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_faceZonesName);
            COM_get_size(regName.c_str(), paneID, &nComp);

            int faceZonesNameTotalSize{nFaceZones * *ca_faceZonesNameMaxLength};
            compareWarningExit(faceZonesNameTotalSize, nComp,
                                "faceZonesNameTotalSize", "nComp");

            faceZonesNameStr = new std::string[nFaceZones]{};
            for (int i=0; i<nFaceZones; i++)
            {
                int startIndex = i * *ca_faceZonesNameMaxLength;
                char* charTmp = &ca_faceZonesName[startIndex];
                faceZonesNameStr[i] = charTmp;
                
                output = std::stringstream{};
                output << "  " << dataName.c_str() << "[" << i << "] = "
                    << charTmp << ", " << faceZonesNameStr[i];
                verbose_message(output.str(), true);
            }
        }

        dataName = std::string("faceZonesCount");
        if(nameExists(dataItemNames, dataName))
        {
            regName = volName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_faceZonesCount);
            COM_get_size(regName.c_str(), paneID, &nComp);
            for(int icomp=0; icomp<nComp; icomp++)
            {
                output = std::stringstream{};
                output << "  " << dataName.c_str() << "[" << icomp << "] = "
                    << ca_faceZonesCount[icomp];
                verbose_message(output.str(), true);
            }
            compareWarningExit(nFaceZones, nComp,
                                "nFaceZones", "faceZonesCount");
        }

        ca_faceZonesList = new int*[nFaceZones]{};
        for (int izone=0; izone<nFaceZones; izone++)
        {
            dataName = "faceZonesList"+std::to_string(izone);
            if(nameExists(dataItemNames, dataName))
            {
                regName = volName+".faceZonesList"+std::to_string(izone);
                COM_get_array(regName.c_str(), paneID, &ca_faceZonesList[izone]);
                COM_get_size(regName.c_str(), paneID, &nComp);
                output = std::stringstream{};
                output << "  " << dataName.c_str() << " has " << nComp
                    << " elements";
                verbose_message(output.str(), true);
            }
            else
            {
                output = std::stringstream{};
                output << "  " << dataName.c_str() << " is a nullptr.";
                verbose_message(output.str(), true);
            }
        }

        ca_faceZonesFlipMap = new int*[nFaceZones]{};
        for (int izone=0; izone<nFaceZones; izone++)
        {
            dataName = "faceZonesFlipMap"+std::to_string(izone);
            if(nameExists(dataItemNames, dataName))
            {
                regName = volName+".faceZonesFlipMap"+std::to_string(izone);
                COM_get_array(regName.c_str(), paneID, &ca_faceZonesFlipMap[izone]);
                COM_get_size(regName.c_str(), paneID, &nComp);
                output = std::stringstream{};
                output << "  " << dataName.c_str() << " has " << nComp
                    << " elements";
                verbose_message(output.str(), true);
            }
            else
            {
                output = std::stringstream{};
                output << "  " << dataName.c_str() << " is a nullptr.";
                verbose_message(output.str(), true);
            }
        }
    }
    //-------------------------------------------

    // pointZone data ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    dataName = std::string("nPointZones");
    if (nameExists(dataItemNames, dataName))
    {
        regName = volName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_nPointZones);
        output = std::stringstream{};
        output << "  " << dataName.c_str() << " = " << *ca_nPointZones;
        verbose_message(output.str(), true);
    }
    else
    {
        ca_nPointZones = new int(0);
    }
    int nPointZones{*ca_nPointZones};

    if (nPointZones>0)
    {
        // pointZonesType
        dataName = std::string("pointZonesTypeMaxLength");
        if (nameExists(dataItemNames, dataName))
        {
            regName = volName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_pointZonesTypeMaxLength);
            output = std::stringstream{};
            output << "  " << dataName.c_str() << " = " << *ca_pointZonesTypeMaxLength;
            verbose_message(output.str(), true);
        }

        int nComp{};
        dataName = std::string("pointZonesType");
        if(nameExists(dataItemNames, dataName))
        {
            regName = volName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_pointZonesType);
            COM_get_size(regName.c_str(), paneID, &nComp);

            int pointZonesTypeTotalSize{nPointZones * *ca_pointZonesTypeMaxLength};
            compareWarningExit(pointZonesTypeTotalSize, nComp,
                                "pointZonesTypeTotalSize", "nComp");

            pointZonesTypeStr = new std::string[nPointZones]{};
            for (int i=0; i<nPointZones; i++)
            {
                int startIndex = i * *ca_pointZonesTypeMaxLength;
                char* charTmp = &ca_pointZonesType[startIndex];
                pointZonesTypeStr[i] = charTmp;
                
                output = std::stringstream{};
                output << "  " << dataName.c_str() << "[" << i << "] = "
                    << charTmp << ", " << pointZonesTypeStr[i];
                verbose_message(output.str(), true);
            }
        }

        // pointZonesName
        dataName = std::string("pointZonesNameMaxLength");
        if (nameExists(dataItemNames, dataName))
        {
            regName = volName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_pointZonesNameMaxLength);
            output = std::stringstream{};
            output << "  " << dataName.c_str() << " = " << *ca_pointZonesNameMaxLength;
            verbose_message(output.str(), true);
        }

        dataName = std::string("pointZonesName");
        if(nameExists(dataItemNames, dataName))
        {
            regName = volName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_pointZonesName);
            COM_get_size(regName.c_str(), paneID, &nComp);

            int pointZonesNameTotalSize{nPointZones * *ca_pointZonesNameMaxLength};
            compareWarningExit(pointZonesNameTotalSize, nComp,
                                "pointZonesNameTotalSize", "nComp");

            pointZonesNameStr = new std::string[nPointZones]{};
            for (int i=0; i<nPointZones; i++)
            {
                int startIndex = i * *ca_pointZonesNameMaxLength;
                char* charTmp = &ca_pointZonesName[startIndex];
                pointZonesNameStr[i] = charTmp;
                
                output = std::stringstream{};
                output << "  " << dataName.c_str() << "[" << i << "] = "
                    << charTmp << ", " << pointZonesNameStr[i];
                verbose_message(output.str(), true);
            }
        }

        dataName = std::string("pointZonesCount");
        if(nameExists(dataItemNames, dataName))
        {
            regName = volName+std::string(".")+dataName;
            COM_get_array(regName.c_str(), paneID, &ca_pointZonesCount);
            COM_get_size(regName.c_str(), paneID, &nComp);
            for(int icomp=0; icomp<nComp; icomp++)
            {
                output = std::stringstream{};
                output << "  " << dataName.c_str() << "[" << icomp << "] = "
                    << ca_pointZonesCount[icomp];
                verbose_message(output.str(), true);
            }
            compareWarningExit(nPointZones, nComp,
                                "nPointZones", "pointZonesCount");
        }

        ca_pointZonesList = new int*[nPointZones]{};
        for (int izone=0; izone<nPointZones; izone++)
        {
            dataName = "pointZonesList"+std::to_string(izone);
            if(nameExists(dataItemNames, dataName))
            {
                regName = volName+".pointZonesList"+std::to_string(izone);
                COM_get_array(regName.c_str(), paneID, &ca_pointZonesList[izone]);
                COM_get_size(regName.c_str(), paneID, &nComp);
                output = std::stringstream{};
                output << "  " << dataName.c_str() << " has " << nComp
                    << " elements";
                verbose_message(output.str(), true);
            }
            else
            {
                output = std::stringstream{};
                output << "  " << dataName.c_str() << " is a nullptr.";
                verbose_message(output.str(), true);
            }
        }
    }
    //-------------------------------------------

    return 0;
}

int comFoam::deleteZonesData()
{
    // cellZones ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    if (cellZonesTypeStr != nullptr)
    {
        delete [] cellZonesTypeStr;
        cellZonesTypeStr = nullptr;
    }

    if (cellZonesNameStr != nullptr)
    {
        delete [] cellZonesNameStr;
        cellZonesNameStr = nullptr;
    }

    if (ca_cellZonesCount != nullptr)
    {
        delete [] ca_cellZonesCount;
        ca_cellZonesCount = nullptr;
    }

    if (ca_cellZonesList != nullptr)
    {
        for(int i=0; i<*ca_nCellZones; i++)
        {
            if (ca_cellZonesList[i] != nullptr)
            {
                delete [] ca_cellZonesList[i];
                ca_cellZonesList[i] = nullptr;
            }
        }

        delete [] ca_cellZonesList;
        ca_cellZonesList = nullptr;
    }

    if (ca_cellZonesType != nullptr)
    {
        delete [] ca_cellZonesType;
        ca_cellZonesType = nullptr;
    }

    if (ca_cellZonesName != nullptr)
    {
        delete [] ca_cellZonesName;
        ca_cellZonesName = nullptr;
    }

    if (ca_cellZonesTypeMaxLength != nullptr)
    {
        delete ca_cellZonesTypeMaxLength;
        ca_cellZonesTypeMaxLength = nullptr;
    }

    if (ca_cellZonesNameMaxLength != nullptr)
    {
        delete ca_cellZonesNameMaxLength;
        ca_cellZonesNameMaxLength = nullptr;
    }

    if (ca_nCellZones != nullptr)
    {
        delete ca_nCellZones;
        ca_nCellZones = nullptr;
    }
    //-------------------------------------------

    // faceZones ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    if (faceZonesTypeStr != nullptr)
    {
        delete [] faceZonesTypeStr;
        faceZonesTypeStr = nullptr;
    }

    if (faceZonesNameStr != nullptr)
    {
        delete [] faceZonesNameStr;
        faceZonesNameStr = nullptr;
    }

    if (ca_faceZonesCount != nullptr)
    {
        delete [] ca_faceZonesCount;
        ca_faceZonesCount = nullptr;
    }

    if (ca_faceZonesList != nullptr)
    {
        for(int i=0; i<*ca_nFaceZones; i++)
        {
            if (ca_faceZonesList[i] != nullptr)
            {
                delete [] ca_faceZonesList[i];
                ca_faceZonesList[i] = nullptr;
            }
        }

        delete [] ca_faceZonesList;
        ca_faceZonesList = nullptr;
    }

    if (ca_faceZonesFlipMap != nullptr)
    {
        for(int i=0; i<*ca_nFaceZones; i++)
        {
            if (ca_faceZonesFlipMap[i] != nullptr)
            {
                delete [] ca_faceZonesFlipMap[i];
                ca_faceZonesFlipMap[i] = nullptr;
            }
        }

        delete [] ca_faceZonesFlipMap;
        ca_faceZonesFlipMap = nullptr;
    }

    if (ca_faceZonesType != nullptr)
    {
        delete [] ca_faceZonesType;
        ca_faceZonesType = nullptr;
    }

    if (ca_faceZonesName != nullptr)
    {
        delete [] ca_faceZonesName;
        ca_faceZonesName = nullptr;
    }

    if (ca_faceZonesTypeMaxLength != nullptr)
    {
        delete ca_faceZonesTypeMaxLength;
        ca_faceZonesTypeMaxLength = nullptr;
    }

    if (ca_faceZonesNameMaxLength != nullptr)
    {
        delete ca_faceZonesNameMaxLength;
        ca_faceZonesNameMaxLength = nullptr;
    }

    if (ca_nFaceZones != nullptr)
    {
        delete ca_nFaceZones;
        ca_nFaceZones = nullptr;
    }
    //-------------------------------------------

    // pointZones ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    if (pointZonesTypeStr != nullptr)
    {
        delete [] pointZonesTypeStr;
        pointZonesTypeStr = nullptr;
    }

    if (pointZonesNameStr != nullptr)
    {
        delete [] pointZonesNameStr;
        pointZonesNameStr = nullptr;
    }

    if (ca_pointZonesCount != nullptr)
    {
        delete [] ca_pointZonesCount;
        ca_pointZonesCount = nullptr;
    }

    if (ca_pointZonesList != nullptr)
    {
        for(int i=0; i<*ca_nPointZones; i++)
        {
            if (ca_pointZonesList[i] != nullptr)
            {
                delete [] ca_pointZonesList[i];
                ca_pointZonesList[i] = nullptr;
            }
        }

        delete [] ca_pointZonesList;
        ca_pointZonesList = nullptr;
    }

    if (ca_pointZonesType != nullptr)
    {
        delete [] ca_pointZonesType;
        ca_pointZonesType = nullptr;
    }

    if (ca_pointZonesName != nullptr)
    {
        delete [] ca_pointZonesName;
        ca_pointZonesName = nullptr;
    }

    if (ca_pointZonesTypeMaxLength != nullptr)
    {
        delete ca_pointZonesTypeMaxLength;
        ca_pointZonesTypeMaxLength = nullptr;
    }

    if (ca_pointZonesNameMaxLength != nullptr)
    {
        delete ca_pointZonesNameMaxLength;
        ca_pointZonesNameMaxLength = nullptr;
    }

    if (ca_nPointZones != nullptr)
    {
        delete ca_nPointZones;
        ca_nPointZones = nullptr;
    }
    //-------------------------------------------

    return 0;
}