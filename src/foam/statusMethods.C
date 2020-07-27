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

    if (ca_isDynamicFvMesh == nullptr)
        ca_isDynamicFvMesh = new int{0};

    if (ca_dynamicSolverType == nullptr)
        ca_dynamicSolverType = new char[genCharSize]{' '};

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

    if (ca_isDynamicFvMesh != nullptr)
        ca_isDynamicFvMesh = new int(0);

    if (ca_dynamicSolverType != nullptr)
        ca_dynamicSolverType = new char[genCharSize]{' '};
    
    dictionary meshDict(mesh.dynamicMeshDict());
    if (meshDict.found("dynamicFvMesh"))
    {
       *ca_isDynamicFvMesh = 1;
       
        if (meshDict.found("solver"))
        {
            word solver = meshDict.lookup("solver");

            //solver = "displacementLaplacian"
            std::string solver_ = solver;


            if (solver_.length() >= genCharSize)
            {
                std::cout << "Warning:: genCharSize is not big enough,"
                          << " genCharSize = " << genCharSize
                          << " solver_.length() = "
                          << solver_.length()
                          << std::endl;
                exit(-1);
            }

            for (size_t i=0; i<solver_.length(); i++)
            {
                ca_dynamicSolverType[i] = solver_[i];
            }
            ca_dynamicSolverType[ solver_.length() ] = '\0';

            for (size_t i=solver_.length()+1; i<genCharSize; i++)
            {
                ca_dynamicSolverType[i] = ' ';
            }
            //std::strcpy(ca_dynamicSolverType, solver_.c_str());
        }
    }

   
    if (ca_timeName != nullptr)
    {
        for (size_t i=0; i<genCharSize; i++)
            ca_timeName[i] = ' ';

        std::string timeNameStr = runTime.timeName();
        size_t length = timeNameStr.length()+1;
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

    dataName = volName+std::string(".isDynamicFvMesh");
    COM_new_dataitem( dataName, 'w', COM_INT, 1, "");
    COM_set_size(     dataName, 0, 1);
    COM_set_array(    dataName, 0, ca_isDynamicFvMesh);
    Info << dataName << " registered." << endl;

    dataName = volName+std::string(".dynamicSolverType");
    COM_new_dataitem( dataName, 'w', COM_CHAR, 1, "");
    COM_set_size(     dataName, 0, genCharSize);
    COM_set_array(    dataName, 0, ca_dynamicSolverType);
    Info << dataName << " registered." << endl;

    COM_window_init_done(volName);

    return 0;
}

int comFoam::reconstStatusData(const char *name)
{
    std::string volName = name+std::string("VOL");

    std::cout << "rocFoam.reconstructStatusData, procID = "
              << ca_myRank
              << ", Retreiving status data form window "
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

    dataName = std::string("isDynamicFvMesh");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), 0, &ca_isDynamicFvMesh);
    std::cout << "    " << dataName.c_str() << " = " << *ca_isDynamicFvMesh << std::endl;

    dataName = std::string("dynamicSolverType");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), 0, &ca_dynamicSolverType);
    COM_get_size(regName.c_str(), 0, &nComp);
    std::cout << "    " << dataName.c_str() << " = " << ca_dynamicSolverType << std::endl;
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
    if (ca_isDynamicFvMesh != nullptr)
    {
        delete ca_isDynamicFvMesh;
        ca_isDynamicFvMesh = nullptr;
    }

    if (ca_dynamicSolverType != nullptr)
    {
        delete [] ca_dynamicSolverType;
        ca_dynamicSolverType = nullptr;
    }

    if (ca_timeName != nullptr)
    {
        delete [] ca_timeName;
        ca_timeName = nullptr;
    }

    return 0;
}
