#include "comFoam.H"

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

    if (ca_dynamicFvMeshType == nullptr)
        ca_dynamicFvMeshType = new char[genCharSize]{' '};

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

    if (ca_dynamicFvMeshType != nullptr)
    {
        for (size_t i=0; i<genCharSize; i++)
        {
            ca_dynamicFvMeshType[i] = ' ';
        }
    }

    if (ca_dynamicSolverType != nullptr)
    {
        for (size_t i=0; i<genCharSize; i++)
        {
            ca_dynamicSolverType[i] = ' ';
        }
    }
    
#ifdef HAVE_OFE20
    dictionary meshDict
               (
                    IOdictionary
                    (
                        IOobject
                        (
                            "dynamicMeshDict",
                            runTime.constant(),
                            mesh,
                            IOobject::READ_IF_PRESENT,
                            IOobject::NO_WRITE,
                            false
                        )
                    )
               );
#elif defined(HAVE_OF7) || defined(HAVE_OF8)
    dictionary meshDict(mesh.dynamicMeshDict());
#endif

    bool foundSolver{false};
    readMeshDict(meshDict, foundSolver);

    if (ca_timeName != nullptr)
    {
        for (size_t i=0; i<genCharSize; i++)
            ca_timeName[i] = ' ';

        std::string timeNameStr = runTime.timeName();
        size_t length = timeNameStr.length()+1;
        if (length > genCharSize)
        {
            WarningInFunction
                << "Warning:: genCharSize is not big enough,"
                << " genCharSize = " << genCharSize
                << " timeName.size = "
                << length
                << endl;
        }
        std::strcpy(ca_timeName, timeNameStr.c_str());
    }

    return 0;
}


int comFoam::registerStatusData(const char *name)
{
    std::string volName = name+std::string("VOL");
    
    std::stringstream output{};
    output << "rocFoam.registerStatusData: "
         << "Registering status data with name "
         << volName;
    verbose_message(output.str(), true);


    // Genral data registered with window
    // Solver data    
    std::string dataName = volName+std::string(".time");
    COM_new_dataitem( dataName, 'w', COM_DOUBLE, 1, "");
    COM_set_size(     dataName, 0, 1);
    COM_set_array(    dataName, 0, ca_time);
    output.str("");  //output = std::stringstream{};
    output << dataName << " registered.";
    verbose_message(output.str(), true);

    dataName = volName+std::string(".timeName");
    COM_new_dataitem( dataName, 'w', COM_CHAR, 1, "");
    COM_set_size(     dataName, 0, genCharSize);
    COM_set_array(    dataName, 0, ca_timeName);
    output.str("");  //output = std::stringstream{};
    output << dataName << " registered.";
    verbose_message(output.str(), true);

    dataName = volName+std::string(".deltaT");
    COM_new_dataitem( dataName, 'w', COM_DOUBLE, 1, "");
    COM_set_size(     dataName, 0, 1);
    COM_set_array(    dataName, 0, ca_deltaT);
    output.str("");  //output = std::stringstream{};
    output << dataName << " registered.";
    verbose_message(output.str(), true);

    dataName = volName+std::string(".deltaT0");
    COM_new_dataitem( dataName, 'w', COM_DOUBLE, 1, "");
    COM_set_size(     dataName, 0, 1);
    COM_set_array(    dataName, 0, ca_deltaT0);
    output.str("");  //output = std::stringstream{};
    output << dataName << " registered.";
    verbose_message(output.str(), true);

    dataName = volName+std::string(".timeIndex");
    COM_new_dataitem( dataName, 'w', COM_INT, 1, "");
    COM_set_size(     dataName, 0, 1);
    COM_set_array(    dataName, 0, ca_timeIndex);
    output.str("");  //output = std::stringstream{};
    output << dataName << " registered.";
    verbose_message(output.str(), true);

    dataName = volName+std::string(".runStat");
    COM_new_dataitem( dataName, 'w', COM_INT, 1, "");
    COM_set_size(     dataName, 0, 1);
    COM_set_array(    dataName, 0, ca_runStat);
    output.str("");  //output = std::stringstream{};
    output << dataName << " registered.";
    verbose_message(output.str(), true);

    dataName = volName+std::string(".dynamicFvMeshType");
    COM_new_dataitem( dataName, 'w', COM_CHAR, 1, "");
    COM_set_size(     dataName, 0, genCharSize);
    COM_set_array(    dataName, 0, ca_dynamicFvMeshType);
    output.str("");  //output = std::stringstream{};
    output << dataName << " registered.";
    verbose_message(output.str(), true);

    dataName = volName+std::string(".isDynamicFvMesh");
    COM_new_dataitem( dataName, 'w', COM_INT, 1, "");
    COM_set_size(     dataName, 0, 1);
    COM_set_array(    dataName, 0, ca_isDynamicFvMesh);
    output.str("");  //output = std::stringstream{};
    output << dataName << " registered.";
    verbose_message(output.str(), true);

    dataName = volName+std::string(".dynamicSolverType");
    COM_new_dataitem( dataName, 'w', COM_CHAR, 1, "");
    COM_set_size(     dataName, 0, genCharSize);
    COM_set_array(    dataName, 0, ca_dynamicSolverType);
    output.str("");  //output = std::stringstream{};
    output << dataName << " registered.";
    verbose_message(output.str(), true);

    COM_window_init_done(volName);

    return 0;
}

int comFoam::reconstStatusData(const char *name)
{
    std::string volName = name+std::string("VOL");

    std::stringstream output{};
    output << "rocFoam.reconstructStatusData, procID = "
           << ca_myRank
           << ", Retreiving status data form window "
           << volName << ".";
    verbose_message(output.str(),true);

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
            nameTmp == "dynamicFvMesh" ||
            nameTmp == "dynamicSolverType"
            )
        {
            dataItemNames.push_back(nameTmp);

            output.str("");  //output = std::stringstream{};
            output << "  DataItem[" << i << "] = " << nameTmp;
            verbose_message(output.str(), true);
        }
    }
    output.str("");  //output = std::stringstream{};
    output << "  Bumber of items = " << dataItemNames.size()
           << std::endl;
    verbose_message(output.str(), true);

    // Flow stat data ^^^^^^^^^^^^^^^^^^^^^^^^^^^
    std::string dataName = std::string("time");
    nameExists(dataItemNames, dataName);
    std::string regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), 0, &ca_time);
    output.str("");  //output = std::stringstream{};
    output << "    " << dataName.c_str() << " = " << *ca_time;
    verbose_message(output.str(), true);

    dataName = std::string("timeIndex");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), 0, &ca_timeIndex);
    output.str("");  //output = std::stringstream{};
    output << "    " << dataName.c_str() << " = " << *ca_timeIndex;
    verbose_message(output.str(), true);

    dataName = std::string("deltaT");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), 0, &ca_deltaT);
    output.str("");  //output = std::stringstream{};
    output << "    " << dataName.c_str() << " = " << *ca_deltaT;
    verbose_message(output.str(), true);

    dataName = std::string("deltaT0");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), 0, &ca_deltaT0);
    output.str("");  //output = std::stringstream{};
    output << "    " << dataName.c_str() << " = " << *ca_deltaT0;
    verbose_message(output.str(), true);

    int nComp;
    dataName = std::string("timeName");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), 0, &ca_timeName);
    COM_get_size(regName.c_str(), 0, &nComp);
    output.str("");  //output = std::stringstream{};
    output << "    " << dataName.c_str() << " = " << ca_timeName;
    verbose_message(output.str(), true);

    dataName = std::string("runStat");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), 0, &ca_runStat);
    output.str("");  //output = std::stringstream{};
    output << "    " << dataName.c_str() << " = " << *ca_runStat;
    verbose_message(output.str(), true);

    dataName = std::string("isDynamicFvMesh");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), 0, &ca_isDynamicFvMesh);
    output.str("");  //output = std::stringstream{};
    output << "    " << dataName.c_str() << " = " << *ca_isDynamicFvMesh;
    verbose_message(output.str(), true);

    dataName = std::string("dynamicFvMeshType");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), 0, &ca_dynamicFvMeshType);
    COM_get_size(regName.c_str(), 0, &nComp);
    output.str("");  //output = std::stringstream{};
    output << "    " << dataName.c_str() << " = " << ca_dynamicFvMeshType;
    verbose_message(output.str(), true);

    dataName = std::string("dynamicSolverType");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), 0, &ca_dynamicSolverType);
    COM_get_size(regName.c_str(), 0, &nComp);
    output.str("");  //output = std::stringstream{};
    output << "    " << dataName.c_str() << " = " << ca_dynamicSolverType;
    verbose_message(output.str(), true);
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

    if (ca_dynamicFvMeshType != nullptr)
    {
        delete [] ca_dynamicFvMeshType;
        ca_dynamicFvMeshType = nullptr;
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
