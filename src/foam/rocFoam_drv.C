//#include "comLoadUnload.H"
#include "rocRhoPimple.H"
#include "rocRhoCentral.H"

COM_EXTERN_MODULE(rocfoam);
COM_EXTERN_MODULE(Rocout);
COM_EXTERN_MODULE(Rocin);
using namespace COM;

MPI_Comm masterComm;
MPI_Comm newComm;

int myRank{0};
int nProcs{1};
bool runParallel{false};
bool rocstarStyle{false};

char *solverType;
std::string status;

//  Status Variables ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
// volume window

char getDataItemLoc;
COM_Type getDataItemType;
int numElementNodes;
std::string getDataItemUnits;

int timeArrayLength;

//  Function Handlers ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
std::vector<std::string> winNames;
std::vector<int> flowInitHandle;
std::vector<int> flowStatHandle;
std::vector<int> flowLoopHandle;
std::vector<int> flowStepHandle;
std::vector<int> flowFinHandle;
std::vector<int> flowRestartInitHandle;

// Rocstar functions handls
std::vector<int> initializeHandle;
std::vector<int> update_solutionHandle;
std::vector<int> finalizeHandle;

// Registered veriables with COM ^^^^^^^^^^^^^^^^^^^^^^^^^^
int *fluidRun;
double *fluidTime;
double *fluidDeltaT;

// Declaration of functions ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
int comDrvInit(int argc, char *argv[]);
int comGetFunctionHandles(const char *name);
//int flowReconst(const char *name);
int comGetVolDataItems(const char *name);
int comGetRunStatItems(const char *name);
int comDrvStat(const char *name);
int comDrvLoop(const char *name);
int comDrvStep(const char *name);
int comDrvFin(const char *name);

int comDrvStart(int argc, char *argv[]);
int comDrvRestart(int argc, char *argv[]);
int comDrvOutput(std::string winName, std::string suffix = "");

int comDrvRestart_Rocstar();
int comDrvStep_Rocstar(const char *name);
int comDrvFin_Rocstar(const char *name);

int main(int argc, char *argv[])
{

    comDrvInit(argc, argv);
    
    if (status == "preprocess")
    {
        comDrvStart(argc, argv);
        comDrvFin(winNames[0].c_str());
    }
    else if (status == "production")
    {
        if (!rocstarStyle)
        {
            comDrvRestart(argc, argv);
            std::string lookUpWindow = winNames[0]+string("VOL");
            comGetRunStatItems(lookUpWindow.c_str());
            comDrvStep(winNames[0].c_str());
            comDrvFin(winNames[0].c_str());
        }
        else
        {
            comDrvRestart_Rocstar();

            std::string lookUpWindow = winNames[0]+string("VOL");
            comGetRunStatItems(lookUpWindow.c_str());

            comDrvStep_Rocstar(winNames[0].c_str());
            comDrvFin_Rocstar(winNames[0].c_str());
        }
    }
    else
    {
        std::cout << " Warning: status = " << status
                  << " is unknown." << std::endl;
    }

    return 0;
}

int comDrvInit(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    masterComm = MPI_COMM_WORLD;

    MPI_Comm_rank(masterComm, &myRank);
    MPI_Comm_size(masterComm, &nProcs);

    if (myRank==0)
    {
        std::cout << "rocFoam.main: Setting up communicator..."
                  << std::endl;

        std::cout << "rocFoam.main:Rank " << myRank
                  << ", NProc = " << nProcs
                  << ", COMM = " << masterComm
                  << std::endl;

        std::cout << std::endl;
    }

    COM_init(&argc, &argv);

    // A new communicator can be generated for
    //   openfoam solver
    newComm = masterComm;
    COM_set_default_communicator(newComm);

    // Run in parallel mode?
    runParallel = false;
    solverType = const_cast<char *>("rocRhoPimple");
    
    std::string arg;
    std::stringstream ss;
    if (argc > 1)
    {
        for (int i=1; i<argc; ++i)
        {
            ss.clear();
            ss.str("");
            ss << argv[i];

            if (ss.str() == "-parallel")
            {
                runParallel = true;
            }
            /*
            else if (ss.str() == "-rocRhoCentral")
            {
                solverType = const_cast<char *>("rocRhoCentral");
            }
            else if (ss.str() == "-rocRhoPimple")
            {
                solverType = const_cast<char *>("rocRhoPimple");
            }
            */
            else if (ss.str() == "-preprocess")
            {
                status = "preprocess";
            }
            else if (ss.str() == "-production")
            {
                status = "production";
            }
            else if (ss.str() == "-rocstar")
            {
                rocstarStyle = true;
            }
            /* else
            {
                if (myRank==0)
                {
                    std::cout << "rocFoam.main: Unknown argumnet"
                              << ss.str() << std::endl;
                }
                return -1;
            } */
        }
    }

    if (status != "preprocess" &&
        status != "production")
    {
        std::cout << "WARNING: simulation status unknown." << std::endl;
        exit(-1);
    }

    if (solverType != const_cast<char *>("rocRhoCentral") &&
        solverType != const_cast<char *>("rocRhoPimple"))
    {
        std::cout << "WARNING: solverType unknown." << std::endl;
        exit(-1);
    }


    if (runParallel && nProcs > 1)
    {
        if (myRank==0)
        {
            std::cout << "rocFoam.main: Running in PRALLEL with solver " 
                      << solverType << "." << std::endl;
        }
    }
    else if (!runParallel && nProcs > 1)
    {
        if (myRank==0)
        {
            std::cout << "rocFoam.main: NProc>1 detected for a serial job."
                      << std::endl;
            return -1;
        }    
    
    }
    else
    {
        runParallel = false;

        if (myRank==0)
        {
            std::cout << "rocFoam.main: Running in SERIAL with solver " 
                      << solverType << "." << std::endl;
        }
    }
    if (myRank==0) std::cout << std::endl;

    return 0;
}

int comDrvStart(int argc, char *argv[])
{
    std::string winName{"ROCFOAM"};
    winNames.push_back(winName);
    
    //rocfoam_load_module(winName.c_str(), solverType);
    COM_LOAD_MODULE_STATIC_DYNAMIC(rocfoam, winName.c_str());
    comGetFunctionHandles(winName.c_str());

    int myArgc = 1;
    char *myArgv[2];
    myArgv[0] = argv[0];
    myArgv[1] = nullptr;

    if (runParallel)
    {
        myArgc = 2;
        myArgv[1] = const_cast<char *>("-parallel");
    }

    //  Fluid initializer ^^^^^^^^^^^^^^^^^^^^^^^
    COM_call_function(flowInitHandle[0],
                      &myArgc, &myArgv,
                      winName.c_str());

    comDrvOutput(winName);

    return 0;
}

int comDrvOutput(std::string winName_, std::string suffix_)
{
    std::string winName{winName_};
    std::cout << "Writing data windows with the base name "
              << winName << std::endl;
    COM_LOAD_MODULE_STATIC_DYNAMIC(SimOUT, "OUT");
    int OUT_write = COM_get_function_handle("OUT.write_dataitem");
    for (int count=0; count<2; count++)
    {
        std::string strTmp;
        if (count == 0)
        {
            strTmp = "VOL";
        }
        else if (count == 1)
        {
            strTmp = "SURF";
        }

        std::string lookUpWindow = winName+strTmp;
        std::string targetName{lookUpWindow+suffix_};
        
        std::string whatToWrite = lookUpWindow+std::string(".all");
        int whatToWriteHandle = COM_get_dataitem_handle(whatToWrite.c_str());

        std::string path = std::string("./")
                           + winName+"/";
        if (count==0)
        {
            if (myRank == 0)
            {
                std::string a = "rmdir /s /q " + path;
                system(a.c_str());
            }
            MPI_Barrier(newComm);
        }

        std::string fullPath = path + targetName+std::string("_");
        std::cout << "  Target window = " << lookUpWindow << std::endl;
        std::cout << "  What to write = " << whatToWrite << std::endl;
        std::cout << "  Path = " << fullPath << std::endl;

        char* outputPath = new char[40]{};
        std::strcpy(outputPath, fullPath.c_str());

        char* material = new char[40]{};
        std::strcpy(material, lookUpWindow.c_str());

        char* timeName = new char[40]{'0'};

        char* mfile_pre{nullptr};

        COM_call_function
        (
            OUT_write,
            outputPath,
            &whatToWriteHandle,
            material,
            timeName,
            mfile_pre,
            &newComm
        );

        delete [] outputPath;
        outputPath = nullptr;

        delete [] material;
        material = nullptr;

        delete [] timeName;
        timeName = nullptr;

        std::cout << "File "
                  << fullPath << " created."
                  << std::endl;

        // Creating text files for SimIO using input files
        std::string attrFile;
        if (count == 0)
        {
            attrFile = "fluid_in_00.000000.txt";
        }
        else if (count == 1)
        {
            attrFile = "ifluid_in_00.000000.txt";
        }
        fullPath = path + attrFile;

        for (int iproc=0; iproc<nProcs; iproc++)
        {
            MPI_Barrier(newComm);
            if(iproc != myRank)
                continue;

            std::stringstream intToOs;
            intToOs << std::setw(4) << std::setfill('0');
            intToOs << myRank;
            std::string outputFile = targetName +
                                     std::string("_") +
                                     intToOs.str() + "*";

            std::string content = "@Proc: "+std::to_string(myRank)+"\n";
            content +="@Files: " + outputFile + "\n";

            int nPanes;
            int* paneList;
            COM_get_panes(lookUpWindow.c_str(), &nPanes, &paneList);

            if (nPanes>0)
            {
                content +="@Panes: ";
                for (int ipane=0; ipane<nPanes; ipane++)
                {
                    content += std::to_string(paneList[ipane]) + " ";
                }
                content += "\n\n";
            }

            std::ofstream outpuFile;
            if (iproc==0)
            {
                outpuFile.open(fullPath, std::ofstream::trunc);

                std::cout << "File "
                          << fullPath << " created."
                          << std::endl;
            }
            else
            {
                outpuFile.open(fullPath, std::ofstream::app);

                std::cout << "File "
                          << fullPath << " appended."
                          << std::endl;
            }
            
            outpuFile << content;
            outpuFile.close();
        }
    }

    COM_UNLOAD_MODULE_STATIC_DYNAMIC(SimOUT, "OUT");
    std::cout << "Unloaded SIMOUT" << std::endl;

    return 0;
}

int comDrvRestart(int argc, char *argv[])
{
    std::string winNameOld = "ROCFOAM0";
    std::string winName = "ROCFOAM";

    //  Fluid initializer ^^^^^^^^^^^^^^^^^^^^^^^
    std::cout << "Reading data windows" << std::endl;
    COM_LOAD_MODULE_STATIC_DYNAMIC(SimIN, "IN");
    //int IN_read = COM_get_function_handle("IN.read_window");
    int IN_read = COM_get_function_handle("IN.read_by_control_file");
    for (int count=0; count<2; count++)
    {
        MPI_Barrier(newComm);

        std::string strTmp;
        if (count == 0 )
        {
            strTmp = "VOL";
        }
        else if (count == 1 )
        {
            strTmp = "SURF";
        }
        
        std::string targetName; // = winName+strTmp+std::string("_");
        if (count == 0 )
        {
            targetName = "fluid_in_00.000000.txt";
        }
        else if (count == 1 )
        {
            targetName = "ifluid_in_00.000000.txt";
        }
        
        std::string lookUpWindow = winNameOld+strTmp;
        std::string path = std::string("./")
                          + winName+"/";

        std::string fullPath = path+targetName;
                  
//        std::ostringstream intToOs;
//        intToOs << std::setw(4) << std::setfill('0');
//        intToOs << myRank;
//        std::string whatToRead = fullPath + intToOs.str() + "*";

        std::string whatToRead = fullPath; //+ "*";

        std::cout << "Proc " << myRank << ": "
                  << "Begin reading "
                  << whatToRead << "."
                  << std::endl;

        COM_call_function
        (
            IN_read,
            whatToRead.c_str(),
            lookUpWindow.c_str(),
            &newComm
        );

        std::cout << "Proc " << myRank << ": "
                  << "Finished reading "
                  << whatToRead << "."
                  << std::endl;
    }
    COM_UNLOAD_MODULE_STATIC_DYNAMIC(SimIN, "IN");
    std::cout << "Unloaded SimIN" << std::endl;

    winNames.push_back(winName);

    COM_LOAD_MODULE_STATIC_DYNAMIC(rocfoam, winName.c_str());
    comGetFunctionHandles(winName.c_str());

    comFoam::copyWindow((winNameOld+"VOL").c_str(),
                        (winName+"VOL").c_str());
    comFoam::copyWindow((winNameOld+"SURF").c_str(),
                        (winName+"SURF").c_str());
    
    int myArgc = 1;
    char *myArgv[2];
    myArgv[0] = argv[0];
    myArgv[1] = nullptr;

    if (runParallel)
    {
        myArgc = 2;
        myArgv[1] = const_cast<char *>("-parallel");
    }

    COM_call_function(flowRestartInitHandle[0],
                      &myArgc, &myArgv,
                      winName.c_str());

    return 0;
}

int comDrvRestart_Rocstar()
{
    std::string winNameOld = "ROCFOAM";
    std::for_each(winNameOld.begin(), winNameOld.end(), [](char & c) {
        c = ::tolower(c);
    });

    std::string winName = "ROCFOAM";
    std::for_each(winName.begin(), winName.end(), [](char & c) {
        c = ::toupper(c);
    });
    
    //  Fluid initializer ^^^^^^^^^^^^^^^^^^^^^^^
    std::cout << "Reading data window with the base name "
              << winNameOld << std::endl;
    COM_LOAD_MODULE_STATIC_DYNAMIC(SimIN, "IN");
    //int IN_read = COM_get_function_handle("IN.read_window");
    int IN_read = COM_get_function_handle("IN.read_by_control_file");
    for (int count=0; count<2; count++)
    {
        std::string strTmp;
        if (count == 0 )
        {
            strTmp = "VOL";
        }
        else if (count == 1 )
        {
            strTmp = "SURF";
        }
        
        std::string targetName;
        if (count == 0 )
        {
            targetName = "fluid_in_00.000000.txt";
        }
        else if (count == 1 )
        {
            targetName = "ifluid_in_00.000000.txt";
        }
        
        std::string lookUpWindow = winNameOld+strTmp;
        std::string path = std::string("./")
                          + winName+"/";

        std::string fullPath = path+targetName;
        std::cout << "Reading file " << fullPath << std::endl;

        std::string whatToRead = fullPath;

        COM_call_function
        (
            IN_read,
            whatToRead.c_str(),
            lookUpWindow.c_str(),
            &newComm
        );

        std::cout << "Finished reading "
                  << whatToRead << " window."
                  << std::endl;
    }
    COM_UNLOAD_MODULE_STATIC_DYNAMIC(SimIN, "IN");
    std::cout << "Unloaded SimIN" << std::endl;

    winNames.push_back(winName);
    
    COM_LOAD_MODULE_STATIC_DYNAMIC(rocfoam, winName.c_str());
    comGetFunctionHandles(winName.c_str());

    double initTime = 0;
    int initHndl = -1;
    int obtHndl = -1;
    std::string volName = winNameOld+"VOL";
    std::string surfName = winNameOld+"SURF";

    COM_call_function
    (
        initializeHandle[0],
        &initTime,
        &newComm,
        &initHndl,
        surfName.c_str(),
        volName.c_str(),
        &obtHndl
    );

    COM_delete_window(volName);
    COM_delete_window(surfName);
    
    return 0;
}


int comGetFunctionHandles(const char *name)
{
    //  Get the handle for the initialize function ^^^^^^^^
    std::vector<std::string>::iterator location = std::find(
            winNames.begin(), winNames.end(), string(name));
    int index = std::distance(winNames.begin(), location);
    std::string winName = winNames[index];

    std::string functionName = winName+string(".flowInit");
    int intTmp = COM_get_function_handle(functionName.c_str());
    flowInitHandle.push_back(intTmp);
    if (intTmp <= 0)
    { // fail
        std::cout << "Could not get handle for "
                  << functionName.c_str()
                  << std::endl;
        exit(-1);
    }
    else
    {
        if (myRank==0)
        {
            std::cout << "Acquired a handle for "
                  << functionName.c_str()
                  << ", handle = " << intTmp
                  << std::endl;
        }
    }

    //  Get the handle for the loop function ^^^^^^^^^^^^^^
    functionName = winName+string(".flowLoop");
    intTmp = COM_get_function_handle(functionName.c_str());
    flowLoopHandle.push_back(intTmp);
    if (intTmp <= 0)
    { // fail
        std::cout << "Could not get handle for "
                  << functionName.c_str()
                  << std::endl;
        exit(-1);
    }
    else
    {
        if (myRank==0)
        {
            std::cout << "Acquired a handle for "
                  << functionName.c_str()
                  << ", handle = " << intTmp
                  << std::endl;
        }
    }

    //  Get the handle for the step function ^^^^^^^^^^^^^^
    functionName = winName+string(".flowStep");
    intTmp = COM_get_function_handle(functionName.c_str());
    flowStepHandle.push_back(intTmp);
    if (intTmp <= 0)
    { // fail
        std::cout << "Could not get handle for "
                  << functionName.c_str()
                  << std::endl;
        exit(-1);
    }
    else
    {
        if (myRank==0)
        {
            std::cout << "Acquired a handle for "
                  << functionName.c_str()
                  << ", handle = " << intTmp
                  << std::endl;
        }
    }

    //  Get the handle for the step function ^^^^^^^^^^^^^^
    functionName = winName+string(".flowRestartInit");
    intTmp = COM_get_function_handle(functionName.c_str());
    flowRestartInitHandle.push_back(intTmp);

    if (intTmp <= 0)
    { // fail
        std::cout << "Could not get handle for "
                  << functionName.c_str()
                  << std::endl;
        exit(-1);
    }
    else
    {
        if (myRank==0)
        {
            std::cout << "Acquired a handle for "
                  << functionName.c_str()
                  << ", handle = " << intTmp
                  << std::endl;
        }
    }

    functionName = winName+string(".initialize");
    intTmp = COM_get_function_handle(functionName.c_str());
    initializeHandle.push_back(intTmp);
    if (intTmp <= 0)
    { // fail
        std::cout << "Could not get handle for "
                  << functionName.c_str()
                  << std::endl;
        exit(-1);
    }
    else
    {
        if (myRank==0)
        {
            std::cout << "Acquired a handle for "
                  << functionName.c_str()
                  << ", handle = " << intTmp
                  << std::endl;
        }
    }

    functionName = winName+string(".update_solution");
    intTmp = COM_get_function_handle(functionName.c_str());
    update_solutionHandle.push_back(intTmp);
    if (intTmp <= 0)
    { // fail
        std::cout << "Could not get handle for "
                  << functionName.c_str()
                  << std::endl;
        exit(-1);
    }
    else
    {
        if (myRank==0)
        {
            std::cout << "Acquired a handle for "
                  << functionName.c_str()
                  << ", handle = " << intTmp
                  << std::endl;
        }
    }

    functionName = winName+string(".finalize");
    intTmp = COM_get_function_handle(functionName.c_str());
    finalizeHandle.push_back(intTmp);
    if (intTmp <= 0)
    { // fail
        std::cout << "Could not get handle for "
                  << functionName.c_str()
                  << std::endl;
        exit(-1);
    }
    else
    {
        if (myRank==0)
        {
            std::cout << "Acquired a handle for "
                  << functionName.c_str()
                  << ", handle = " << intTmp
                  << std::endl;
        }
    }

    return 0;
}

int comGetRunStatItems(const char *name)
{
    std::string volName = name;

    int numDataItems=0;
    std::vector<std::string> dataItemNames;

         
    std::string output;
    COM_get_dataitems(volName.c_str(), &numDataItems, output);
    Info << "  numDataItems = " << numDataItems << endl;

    std::istringstream Istr(output);
    dataItemNames.clear();
    for (int i=0; i<numDataItems; ++i)
    {
        std::string nameTmp;
        Istr >> nameTmp;
        dataItemNames.push_back(nameTmp);
        //Info << "  DataItem[" << i << "] = " << nameTmp << endl;
    }
    //Info << endl;


    std::string dataName = string("time");
    std::string regName = volName+string(".")+dataName;
    bool ifCorrect = (std::find(dataItemNames.begin(),
                                dataItemNames.end(), dataName)
                                != dataItemNames.end());
    if (ifCorrect)
    {
        COM_get_array(regName.c_str(), 0, &fluidTime);
        Info << "  " << dataName.c_str() << " = " << *fluidTime << endl;
    }

    dataName = string("deltaT");
    regName = volName+string(".")+dataName;
    ifCorrect = (std::find(dataItemNames.begin(),
                                dataItemNames.end(), dataName)
                                != dataItemNames.end());
    if (ifCorrect)
    {
        COM_get_array(regName.c_str(), 0, &fluidDeltaT);
        Info << "  " << dataName.c_str() << " = " << *fluidDeltaT << endl;
    }

    dataName = string("runStat");
    regName = volName+string(".")+dataName;
    ifCorrect = (std::find(dataItemNames.begin(),
                                dataItemNames.end(), dataName)
                                != dataItemNames.end());
    if (ifCorrect)
    {
        COM_get_array(regName.c_str(), 0, &fluidRun);
        Info << "  " << dataName.c_str() << " = " << *fluidRun << endl;
    }
    
    return 0;
}

int comDrvLoop(char* const name)
{
    //  Call the flow iterator ^^^^^^^^^^^^^^^^^^
    std::vector<std::string>::iterator location = std::find(
            winNames.begin(), winNames.end(), string(name));
    int index = std::distance(winNames.begin(), location);
    //std::string volName = winNames[index]+string("VOL");


    COM_call_function(flowLoopHandle[index], name);
    
    return 0;
}

int comDrvStep(const char* name)
{
    std::vector<std::string>::iterator location = std::find(
            winNames.begin(), winNames.end(), string(name));
    int index = std::distance(winNames.begin(), location);
    //std::string volName = winNames[index]+string("VOL");

    //  Call the flow stepper ^^^^^^^^^^^^^^^^^^
    Info << "\nStarting time loop\n" << endl;

    while (*fluidRun)
    {
        COM_call_function(flowStepHandle[index]);
    }

    Info << "End\n" << endl;

    Info << "rocFoam.main: Stepping status = " 
         << *fluidRun
         << endl;

    Info << "rocFoam.main: Simulation time = "
         << *fluidTime
         << endl;

    Info << "rocFoam.main: Simulation Time step = "
         << *fluidDeltaT
         << endl;

    return 0;
}

int comDrvStep_Rocstar(const char* name)
{
    std::vector<std::string>::iterator location = std::find(
            winNames.begin(), winNames.end(), string(name));
    int index = std::distance(winNames.begin(), location);
    //std::string volName = winNames[index]+string("VOL");

    //  Call the flow stepper ^^^^^^^^^^^^^^^^^^
    Info << "\nStarting time loop\n" << endl;

    while (*fluidRun)
    {
        double currentTime = *fluidTime;
        double timeStep{1}; //*fluidDeltaT;
        int bcHandle{-1};
        int gmHandle{-1};

        COM_call_function
        (
            update_solutionHandle[index],
            &currentTime,
            &timeStep,
            &bcHandle,
            &gmHandle
        );
    }

    Info << "End\n" << endl;

    Info << "rocFoam.main: Stepping status = " 
         << *fluidRun
         << endl;

    Info << "rocFoam.main: Simulation time = "
         << *fluidTime
         << endl;

    Info << "rocFoam.main: Simulation Time step = "
         << *fluidDeltaT
         << endl;

    return 0;
}

int comDrvFin(const char* name)
{
    //  Call the flow unloader ^^^^^^^^^^^^^^^^^^
    //COM_UNLOAD_MODULE_STATIC_DYNAMIC(rocfoam, "ROCFOAM");
    //rocfoam_unload_module(name, solverType);
    COM_UNLOAD_MODULE_STATIC_DYNAMIC(rocfoam, name);

    COM_set_default_communicator(masterComm);
    
    COM_finalize();

    MPI_Barrier(masterComm);
    MPI_Finalize();
    
    return 0;
}

int comDrvFin_Rocstar(const char* name)
{
    //  Call the flow unloader ^^^^^^^^^^^^^^^^^^
    std::vector<std::string>::iterator location = std::find(
            winNames.begin(), winNames.end(), string(name));
    int index = std::distance(winNames.begin(), location);

    COM_call_function( finalizeHandle[index] );

    COM_UNLOAD_MODULE_STATIC_DYNAMIC(rocfoam, name);
    
    COM_finalize();

    MPI_Barrier(masterComm);
    MPI_Finalize();
    
    return 0;
}

