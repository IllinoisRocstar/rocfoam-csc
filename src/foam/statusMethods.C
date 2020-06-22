// Setting comFoam variables that will be registered
//   with COM. These are needed for flow control when
//   solver runs step-by-step.

int comFoam::createStatusData()
{
    if (ca_runStat == nullptr)
        ca_runStat = new int{0};

    if (ca_time == nullptr)
        ca_time = new double{0};

    if (ca_deltaT == nullptr)
        ca_deltaT = new double{0};

    if (ca_deltaT0 == nullptr)
        ca_deltaT0 = new double{0};

    if (ca_timeIndex == nullptr)
        ca_timeIndex = new int{0};

    if (ca_movingMesh == nullptr)
        ca_movingMesh = new int{0};

    if (ca_dynamicFvMesh == nullptr)
        ca_dynamicFvMesh = new char[genCharSize]{'0'};

    if (ca_timeName == nullptr)
        ca_timeName = new char[genCharSize]{'0'};

    return 0;
}

int comFoam::updateStatusData()
{
    const Foam::Time &runTime(*runTimePtr);
    const dynamicFvMesh& mesh(*meshPtr);

    if (ca_runStat != nullptr)
        *ca_runStat = static_cast<int>(runTime.run());

    if (ca_time != nullptr)
        *ca_time = runTime.value();

    if (ca_deltaT != nullptr)
        *ca_deltaT = runTime.deltaTValue();

    if (ca_deltaT0 != nullptr)
        *ca_deltaT0 = runTime.deltaT0Value();

    if (ca_timeIndex != nullptr)
        *ca_timeIndex = runTime.timeIndex();

    if (ca_movingMesh != nullptr)
        *ca_movingMesh = mesh.moving();

    if (ca_dynamicFvMesh != nullptr)
    {
        for (int i=0; i<genCharSize; i++)
            ca_dynamicFvMesh[i] = ' ';
            

        IOdictionary dynamicMeshDict
        (
            IOobject
            (
                "dynamicMeshDict",
                runTime.constant(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            )
        );
        //word dynamicFvMesh1_ = dynamicMeshDict.lookup("dynamicFvMesh");
        
        std::string dynamicFvMesh2_; // = dynamicFvMesh1_;
        int length = dynamicFvMesh2_.length()+1;
        if (length > genCharSize)
        {
            std::cout << "Warning:: genCharSize is not big enough,"
                      << " genCharSize = " << genCharSize
                      << " dynamicFvMesh2_.size = "
                      << length
                      << std::endl;
        }
        std::strcpy(ca_dynamicFvMesh, dynamicFvMesh2_.c_str());
        //( std::string(ca_dynamicFvMesh) == "dynamicMotionSolverFvMesh")
    }
    
    if (ca_timeName != nullptr)
    {
        for (int i=0; i<genCharSize; i++)
            ca_timeName[i] = ' ';

        std::string timeNameStr = runTime.timeName();
        int length = timeNameStr.length()+1;
        if (length > genCharSize)
        {
            std::cout << "Warning:: genCharSize is not big enough,"
                      << " genCharSize = " << genCharSize
                      << " timeName.size = "
                      << length
                      << std::endl;
        }
        std::strcpy(ca_timeName, timeNameStr.c_str());
    }

    return 0;
}


int comFoam::registerStatusData(const char *name)
{
    std::string volName = name+std::string("VOL");
    
    Info << "rocFoam.registerStatusData: "
         << "Registering status data with name "
         << volName
         << endl;

    // Genral data registered with window
    // Solver data    
    std::string dataName = volName+std::string(".time");
    COM_new_dataitem( dataName, 'w', COM_DOUBLE, 1, "");
    COM_set_size(     dataName, 0, 1);
    COM_set_array(    dataName, 0, ca_time);
    Info << dataName << " registered." << endl;

    dataName = volName+std::string(".timeName");
    COM_new_dataitem( dataName, 'w', COM_CHAR, 1, "");
    COM_set_size(     dataName, 0, genCharSize);
    COM_set_array(    dataName, 0, ca_timeName);
    Info << dataName << " registered." << endl;

    dataName = volName+std::string(".deltaT");
    COM_new_dataitem( dataName, 'w', COM_DOUBLE, 1, "");
    COM_set_size(     dataName, 0, 1);
    COM_set_array(    dataName, 0, ca_deltaT);
    Info << dataName << " registered." << endl;

    dataName = volName+std::string(".deltaT0");
    COM_new_dataitem( dataName, 'w', COM_DOUBLE, 1, "");
    COM_set_size(     dataName, 0, 1);
    COM_set_array(    dataName, 0, ca_deltaT0);
    Info << dataName << " registered." << endl;

    dataName = volName+std::string(".timeIndex");
    COM_new_dataitem( dataName, 'w', COM_INT, 1, "");
    COM_set_size(     dataName, 0, 1);
    COM_set_array(    dataName, 0, ca_timeIndex);
    Info << dataName << " registered." << endl;

    dataName = volName+std::string(".runStat");
    COM_new_dataitem( dataName, 'w', COM_INT, 1, "");
    COM_set_size(     dataName, 0, 1);
    COM_set_array(    dataName, 0, ca_runStat);
    Info << dataName << " registered." << endl;

    dataName = volName+std::string(".movingMesh");
    COM_new_dataitem( dataName, 'w', COM_INT, 1, "");
    COM_set_size(     dataName, 0, 1);
    COM_set_array(    dataName, 0, ca_movingMesh);
    Info << dataName << " registered." << endl;

    dataName = volName+std::string(".dynamicFvMesh");
    COM_new_dataitem( dataName, 'w', COM_CHAR, 1, "");
    COM_set_size(     dataName, 0, genCharSize);
    COM_set_array(    dataName, 0, ca_dynamicFvMesh);
    Info << dataName << " registered." << endl;

    return 0;
}

int comFoam::reconstStatusData(const char *name)
{
    std::string volName = name+std::string("VOL");

    std::cout << "rocFoam.reconstructStatusData, procID = "
              << Pstream::myProcNo()
              << ", Retreiving surface data form window "
              << volName << "."
              << std::endl;

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

        std::string subName = nameTmp.substr(0,4);
        if (nameTmp == "time" ||
            nameTmp == "timeIndex" ||
            nameTmp == "deltaT" ||
            nameTmp == "deltaT0" ||
            nameTmp == "timeName" ||
            nameTmp == "runStat" ||
            nameTmp == "movingMesh" ||
            nameTmp == "dynamicFvMesh")
        {
            dataItemNames.push_back(nameTmp);
            std::cout << "  DataItem[" << i << "] = " << nameTmp << std::endl;
        }
    }
    std::cout << "  Bumber of items = " << dataItemNames.size()
              << std::endl << std::endl;

    // Flow stat data ^^^^^^^^^^^^^^^^^^^^^^^^^^^
    std::string dataName = std::string("time");
    nameExists(dataItemNames, dataName);
    std::string regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), 0, &ca_time);
    std::cout << "    " << dataName.c_str() << " = " << *ca_time << std::endl;

    dataName = std::string("timeIndex");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), 0, &ca_timeIndex);
    std::cout << "    " << dataName.c_str() << " = " << *ca_timeIndex << std::endl;

    dataName = std::string("deltaT");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), 0, &ca_deltaT);
    std::cout << "    " << dataName.c_str() << " = " << *ca_deltaT << std::endl;

    dataName = std::string("deltaT0");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), 0, &ca_deltaT0);
    std::cout << "    " << dataName.c_str() << " = " << *ca_deltaT0 << std::endl;

    int nComp;
    dataName = std::string("timeName");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), 0, &ca_timeName);
    COM_get_size(regName.c_str(), 0, &nComp);
    std::cout << "    " << dataName.c_str() << " = " << ca_timeName << std::endl;

    dataName = std::string("runStat");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), 0, &ca_runStat);
    std::cout << "    " << dataName.c_str() << " = " << *ca_runStat << std::endl;

    dataName = std::string("movingMesh");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), 0, &ca_movingMesh);
    std::cout << "    " << dataName.c_str() << " = " << *ca_movingMesh << std::endl;

    dataName = std::string("dynamicFvMesh");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), 0, &ca_dynamicFvMesh);
    COM_get_size(regName.c_str(), 0, &nComp);
    std::cout << "    " << dataName.c_str() << " = " << ca_dynamicFvMesh << std::endl;
    //-------------------------------------------

    return 0;
}


int comFoam::deleteStatusData()
{
    if (ca_runStat != nullptr)
    {
        delete ca_runStat;
        ca_runStat = nullptr;
    }
    if (ca_time != nullptr)
    {
        delete ca_time;
        ca_time = nullptr;
    }
    if (ca_deltaT != nullptr)
    {
        delete ca_deltaT;
        ca_deltaT = nullptr;
    }
    if (ca_deltaT0 != nullptr)
    {
        delete ca_deltaT0;
        ca_deltaT0 = nullptr;
    }
    if (ca_timeIndex != nullptr)
    {
        delete ca_timeIndex;
        ca_timeIndex = nullptr;
    }
    if (ca_movingMesh != nullptr)
    {
        delete ca_movingMesh;
        ca_movingMesh = nullptr;
    }

    if (ca_dynamicFvMesh != nullptr)
    {
        delete ca_dynamicFvMesh;
        ca_dynamicFvMesh = nullptr;
    }

    if (ca_timeName != nullptr)
    {
        delete ca_timeName;
        ca_timeName = nullptr;
    }

    return 0;
}
