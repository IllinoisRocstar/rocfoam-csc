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

int comDrvRestart_Rocstar();
int comDrvStep_Rocstar(const char *name);

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
    //solverType = const_cast<char *>("rocRhoCentral");
    
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
            else if (ss.str() == "-rocRhoCentral")
            {
                solverType = const_cast<char *>("rocRhoCentral");
            }
            else if (ss.str() == "-rocRhoPimple")
            {
                solverType = const_cast<char *>("rocRhoPimple");
            }
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

    std::cout << "Writing data windows" << std::endl;
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
        std::string whatToWrite = lookUpWindow+std::string(".all");
        int whatToWriteHandle = COM_get_dataitem_handle(whatToWrite.c_str());

        std::string pathTmp = std::string("./")
                            + winName+"/"
                            + lookUpWindow+std::string("_");

        char* outputPath = new char[40]{};
        std::strcpy(outputPath, pathTmp.c_str());

        char* material = new char[40]{};
        std::strcpy(material, lookUpWindow.c_str());

        char* timeName = new char[40]{};

        COM_call_function
        (
            OUT_write,
            outputPath,
            &whatToWriteHandle,
            material,
            timeName
        );

        delete [] outputPath;
        outputPath = nullptr;

        delete [] material;
        material = nullptr;

        delete [] timeName;
        timeName = nullptr;

        std::cout << "Finished writting "
                  << strTmp << " window."
                  << std::endl;
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
    int IN_read = COM_get_function_handle("IN.read_window");
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
        
        
        std::string lookUpWindow = winNameOld+strTmp;

        std::string pathTmp = std::string("./")
                            + winName+"/"
                            + winName+strTmp+std::string("_");

        std::cout << "Reading file " << pathTmp << std::endl;

        std::string whatToRead = pathTmp+"*";

        COM_call_function
        (
            IN_read,
            whatToRead.c_str(),
            lookUpWindow.c_str()
        );

        std::cout << "Finished reading "
                  << pathTmp << " window."
                  << std::endl;
    }
    COM_UNLOAD_MODULE_STATIC_DYNAMIC(SimIN, "IN");
    std::cout << "Unloaded SimIN" << std::endl;

    winNames.push_back(winName);
    //rocfoam_load_module(winName.c_str(), solverType);
    COM_LOAD_MODULE_STATIC_DYNAMIC(rocfoam, winName.c_str());
    comGetFunctionHandles(winName.c_str());

    comFoam::copyWindow((winNameOld+"VOL").c_str(),
                        (winName+"VOL").c_str());
    comFoam::copyWindow((winNameOld+"SURF").c_str(),
                        (winName+"SURF").c_str());
    
    //rocfoam_unload_module(winNameOld.c_str(), solverType);
    //COM_UNLOAD_MODULE_STATIC_DYNAMIC(rocfoam, winNameOld.c_str());

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
    std::string winNameOld = "ROCFOAM0";
    std::string winName = "ROCFOAM";
    //  Fluid initializer ^^^^^^^^^^^^^^^^^^^^^^^
    std::cout << "Reading data windows" << std::endl;
    COM_LOAD_MODULE_STATIC_DYNAMIC(SimIN, "IN");
    int IN_read = COM_get_function_handle("IN.read_window");
    for (int count=0; count<2; count++)
    {
        std::string strTmp;
        if (count == 0 )
        {
            strTmp = "_vol";
        }
        else if (count == 1 )
        {
            strTmp = "_srf";
        }
        
        
        std::string lookUpWindow = winNameOld+strTmp;

        std::string pathTmp = std::string("./")
                            + winName+"/"
                            + winName+strTmp+std::string("_");

        std::cout << "Reading file " << pathTmp << std::endl;

        std::string whatToRead = pathTmp+"*";

        COM_call_function
        (
            IN_read,
            whatToRead.c_str(),
            lookUpWindow.c_str()
        );

        std::cout << "Finished reading "
                  << pathTmp << " window."
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
    std::string volName = winName+"_vol";
    std::string surfName = winName+"_srf";


    COM_call_function
    (
        initializeHandle[0],
        initTime,
        newComm,
        initHndl,
        surfName.c_str(),
        volName.c_str(),
        obtHndl
    );

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

int comGetVolDataItems(const char *name)
{
    std::string volName = name;

    int numDataItems=0;
    std::vector<std::string> dataItemNames;
    int numPanes;
    int* paneList;

    double *volCoord;
    int volNumNodes = 0;
    int nComp = 0;

    int numCells = 0;
    double* cellVel;
    double* cellPres;
    double* cellTemp;
    double* cellRho;

    int* Conn;
    int numConn;
    int numElem;

    int* Owner;
    int* Neighb;
    int numFaces;

    std::vector<std::string> connNames;
    Info << endl
         << "rocFoam.main: Retreiving data form window "
         << volName << "."
         << endl;
         
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
        Info << "  DataItem[" << i << "] = " << nameTmp << endl;
    }
    Info << endl;

    // Volume data ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    std::string dataName = string("nProc");
    std::string regName = volName+string(".")+dataName;
    bool ifCorrect = (std::find(dataItemNames.begin(),
                                dataItemNames.end(), dataName)
                                != dataItemNames.end());
    if (ifCorrect)
    {
        int *nProcReg;
        COM_get_array(regName.c_str(), 0, &nProcReg);
        Info << "  " << dataName.c_str() << " = " << *nProcReg << endl;
    }

    dataName = string("time");
    regName = volName+string(".")+dataName;
    ifCorrect = (std::find(dataItemNames.begin(),
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

    dataName = string("nPoints");
    regName = volName+string(".")+dataName;
    ifCorrect = (std::find(dataItemNames.begin(),
                                dataItemNames.end(), dataName)
                                != dataItemNames.end());
    if (ifCorrect)
    {
        int *nPoints;
        COM_get_array(regName.c_str(), 0, &nPoints);
        Info << "  " << dataName.c_str() << " = " << *nPoints << endl;
    }

    dataName = string("nCells");
    regName = volName+string(".")+dataName;
    ifCorrect = (std::find(dataItemNames.begin(),
                                dataItemNames.end(), dataName)
                                != dataItemNames.end());
    if (ifCorrect)
    {
        int *nCells;
        COM_get_array(regName.c_str(), 0, &nCells);
        Info << "  " << dataName.c_str() << " = " << *nCells << endl;
    }

    dataName = string("cellToPointConn_types");
    regName = volName+string(".")+dataName;
    ifCorrect = (std::find(dataItemNames.begin(),
                                dataItemNames.end(), dataName)
                                != dataItemNames.end());
    if (ifCorrect)
    {
        int *cellToPointConn_types;
        COM_get_array(regName.c_str(), 0, &cellToPointConn_types);
        Info << "  " << dataName.c_str() << " = " << *cellToPointConn_types << endl;
    }

    dataName = string("cellToPointConn_map");
    regName = volName+string(".")+dataName;
    ifCorrect = (std::find(dataItemNames.begin(),
                                dataItemNames.end(), dataName)
                                != dataItemNames.end());
    if (ifCorrect)
    {
        int *cellToPointConn_map;
        COM_get_array(regName.c_str(), 0, &cellToPointConn_map);
        COM_get_size(regName.c_str(), 0, &nComp);
        //Info << "  " << dataName.c_str() << " size = " << nComp << endl;
        for(int icomp=0; icomp<nComp; icomp++)
        {
            Info << "  " << dataName.c_str() << "[" << icomp << "] = "
                 << cellToPointConn_map[icomp] << endl;
        }
    }

    dataName = string("cellToPointConn_size");
    regName = volName+string(".")+dataName;
    ifCorrect = (std::find(dataItemNames.begin(),
                                dataItemNames.end(), dataName)
                                != dataItemNames.end());
    if (ifCorrect)
    {
        int *cellToPointConn_size;
        COM_get_array(regName.c_str(), 0, &cellToPointConn_size);
        COM_get_size(regName.c_str(), 0, &nComp);
        //Info << "  " << dataName.c_str() << " size = " << nComp << endl;
        for(int icomp=0; icomp<nComp; icomp++)
        {
            Info << "  " << dataName.c_str() << "[" << icomp << "] = "
                 << cellToPointConn_size[icomp] << endl;
        }
    }
    //-------------------------------------------

    // Face data ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    int* nFaces;
    dataName = string("nFaces");
    regName = volName+string(".")+dataName;
    ifCorrect = (std::find(dataItemNames.begin(),
                                dataItemNames.end(), dataName)
                                != dataItemNames.end());
    if (ifCorrect)
    {
        COM_get_array(regName.c_str(), 0, &nFaces);
        Info << "  " << dataName.c_str() << " = " << *nFaces << endl;
    }

    dataName = string("faceToPointConn_types");
    regName = volName+string(".")+dataName;
    ifCorrect = (std::find(dataItemNames.begin(),
                                dataItemNames.end(), dataName)
                                != dataItemNames.end());
    if (ifCorrect)
    {
        int *faceToPointConn_types;
        COM_get_array(regName.c_str(), 0, &faceToPointConn_types);
        Info << "  " << dataName.c_str() << " = " << *faceToPointConn_types << endl;
    }


    dataName = string("faceToPointConn_map");
    regName = volName+string(".")+dataName;
    ifCorrect = (std::find(dataItemNames.begin(),
                                dataItemNames.end(), dataName)
                                != dataItemNames.end());
    if (ifCorrect)
    {
        int *faceToPointConn_map;
        COM_get_array(regName.c_str(), 0, &faceToPointConn_map);
        COM_get_size(regName.c_str(), 0, &nComp);
        //Info << "  " << dataName.c_str() << " size = " << nComp << endl;
        for(int icomp=0; icomp<nComp; icomp++)
        {
            Info << "  " << dataName.c_str() << "[" << icomp << "] = "
                 << faceToPointConn_map[icomp] << endl;
        }
    }

    dataName = string("faceToPointConn_size");
    regName = volName+string(".")+dataName;
    ifCorrect = (std::find(dataItemNames.begin(),
                                dataItemNames.end(), dataName)
                                != dataItemNames.end());
    if (ifCorrect)
    {
        int *faceToPointConn_size;
        COM_get_array(regName.c_str(), 0, &faceToPointConn_size);
        COM_get_size(regName.c_str(), 0, &nComp);
        //Info << "  " << dataName.c_str() << " size = " << nComp << endl;
        for(int icomp=0; icomp<nComp; icomp++)
        {
            Info << "  " << dataName.c_str() << "[" << icomp << "] = "
                 << faceToPointConn_size[icomp] << endl;
        }
    }
    //-------------------------------------------
    
    // Surface data ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    int* nPatches;
    dataName = string("nPatches");
    regName = volName+string(".")+dataName;
    ifCorrect = (std::find(dataItemNames.begin(),
                                dataItemNames.end(), dataName)
                                != dataItemNames.end());
    if (ifCorrect)
    {
        COM_get_array(regName.c_str(), 0, &nPatches);
        Info << "  " << dataName.c_str() << " = " << *nPatches << endl;
    }

    //  List of panes in this window ^^^^^^^^^^^^
    COM_get_panes(volName.c_str(), &numPanes, &paneList);
    Info << "  Number of Panes = " << numPanes << endl;

    for (int i=0; i<numPanes; ++i) 
    {
        Info << "  Pane[" << i
             <<"]^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;

        dataName = string("patchName");
        regName = volName+string(".")+dataName;
        ifCorrect = (std::find(dataItemNames.begin(),
                                    dataItemNames.end(), dataName)
                                    != dataItemNames.end());
        if (ifCorrect)
        {
            char* patchName;
        
            COM_get_array(regName.c_str(), paneList[i], &patchName);
            COM_get_size(regName.c_str(), paneList[i], &nComp);

            std::istringstream iss(patchName);
            
            std::string strPatchName;
            iss >> strPatchName;
            strPatchName.resize(nComp);
            Info << "    " << dataName.c_str() << " = " << strPatchName.c_str() << endl;
        }

        dataName = string("patchType");
        regName = volName+string(".")+dataName;
        ifCorrect = (std::find(dataItemNames.begin(),
                                    dataItemNames.end(), dataName)
                                    != dataItemNames.end());
        if (ifCorrect)
        {
            char* patchType;
        
            COM_get_array(regName.c_str(), paneList[i], &patchType);
            COM_get_size(regName.c_str(), paneList[i], &nComp);

            std::istringstream iss(patchType);
            
            std::string strPatchType;
            iss >> strPatchType;
            strPatchType.resize(nComp);
            Info << "    " << dataName.c_str() << " = " << strPatchType.c_str() << endl;
        }

        dataName = string("patchStart");
        regName = volName+string(".")+dataName;
        ifCorrect = (std::find(dataItemNames.begin(),
                                    dataItemNames.end(), dataName)
                                    != dataItemNames.end());
        if (ifCorrect)
        {
            int *patchStart;
            COM_get_array(regName.c_str(), paneList[i], &patchStart);
            Info << "    " << dataName.c_str() << " = " << *patchStart << endl;
        }

        dataName = string("patchSize");
        regName = volName+string(".")+dataName;
        ifCorrect = (std::find(dataItemNames.begin(),
                                    dataItemNames.end(), dataName)
                                    != dataItemNames.end());
        if (ifCorrect)
        {
            int *patchSize;
            COM_get_array(regName.c_str(), paneList[i], &patchSize);
            Info << "    " << dataName.c_str() << " = " << *patchSize << endl;
        }

        dataName = string("patchPointToPointMap_size");
        regName = volName+string(".")+dataName;
        ifCorrect = (std::find(dataItemNames.begin(),
                                    dataItemNames.end(), dataName)
                                    != dataItemNames.end());
        if (ifCorrect)
        {
            int *patchPointToPointMap_size;
            COM_get_array(regName.c_str(), paneList[i], &patchPointToPointMap_size);
            Info << "    " << dataName.c_str() << " = " << *patchPointToPointMap_size << endl;
        }

        dataName = string("patchFaceToPointConn_types");
        regName = volName+string(".")+dataName;
        ifCorrect = (std::find(dataItemNames.begin(),
                                    dataItemNames.end(), dataName)
                                    != dataItemNames.end());
        if (ifCorrect)
        {
            int *patchFaceToPointConn_types;
            COM_get_array(regName.c_str(), paneList[i], &patchFaceToPointConn_types);
            Info << "    " << dataName.c_str() << " = " << *patchFaceToPointConn_types << endl;
        }


        dataName = string("patchFaceToPointConn_map");
        regName = volName+string(".")+dataName;
        ifCorrect = (std::find(dataItemNames.begin(),
                                    dataItemNames.end(), dataName)
                                    != dataItemNames.end());
        if (ifCorrect)
        {
            int *patchFaceToPointConn_map;
            COM_get_array(regName.c_str(), paneList[i], &patchFaceToPointConn_map);
            COM_get_size(regName.c_str(), paneList[i], &nComp);
            //Info << "    " << dataName.c_str() << " size = " << nComp << endl;
            for(int icomp=0; icomp<nComp; icomp++)
            {
                Info << "    " << dataName.c_str() << "[" << icomp << "] = "
                     << patchFaceToPointConn_map[icomp] << endl;
            }
        }

        dataName = string("patchFaceToPointConn_size");
        regName = volName+string(".")+dataName;
        ifCorrect = (std::find(dataItemNames.begin(),
                                    dataItemNames.end(), dataName)
                                    != dataItemNames.end());
        if (ifCorrect)
        {
            int *patchFaceToPointConn_size;
            COM_get_array(regName.c_str(), paneList[i], &patchFaceToPointConn_size);
            COM_get_size(regName.c_str(), paneList[i], &nComp);
            //Info << "    " << dataName.c_str() << " size = " << nComp << endl;
            for(int icomp=0; icomp<nComp; icomp++)
            {
                Info << "    " << dataName.c_str() << "[" << icomp << "] = "
                     << patchFaceToPointConn_size[icomp] << endl;
            }
        }





        // Point and connectivity stuff ^^^^^^^^^
        dataName = string("nc");
        regName = volName+string(".")+dataName;
        COM_get_array(regName.c_str(), paneList[i], &volCoord, &nComp);
        COM_get_size(regName.c_str(), paneList[i], &volNumNodes);
        Info << "    " << dataName.c_str() << " points = " << volNumNodes
             << ", components = " << nComp << endl;
//        for(int ipoint=0; ipoint<volNumNodes; ipoint++)
//        {
//          Info << "Node " << ipoint << " volCoordinates = ";
//          for(int icomp=0; icomp<nComp; icomp++)
//          {
//              Info << *(volCoord+ipoint*nComp+icomp) << " ";
//          }
//          Info << endl;
//        }

        std::string connNames;
        COM_get_connectivities(volName.c_str(), paneList[i], &numConn, connNames);
        std::istringstream connISS(connNames);
        //Info << "    Connectivity types =" << numConn << endl;

        for (int j=0; j<numConn; ++j)
        {
            std::string connName;
            connISS >> connName;
            //connNames.push_back(connName);

            dataName = volName+string(".")+connName;
            COM_get_array(dataName.c_str(), paneList[i], &Conn, &nComp);
            COM_get_size(dataName.c_str(), paneList[i], &numElem);

            Info << "    Connectivity["
                 << j << "] = " << connName
                 << ", elements = " << numElem
                 << ", components =" << nComp << endl;

//            for(int icell=0; icell<numElem; icell++)
//            {
//                Info << "    Element " << icell << " connectivities = ";
//                for(int icomp=0; icomp<nComp; icomp++)
//                {
//                    Info << *(Conn+icell*nComp+icomp) << " ";
//                }
//                Info << endl;
//            }

        }
        //---------------------------------------



        // Mapping data
        dataName = string("cellToCellMap");
        regName = volName+string(".")+dataName;
        ifCorrect = (std::find(dataItemNames.begin(),
                                    dataItemNames.end(), dataName)
                                    != dataItemNames.end());
        if (ifCorrect)
        {
            int* cellToCellMap;
            
            COM_get_array(regName.c_str(), paneList[i], &cellToCellMap, &nComp);
            COM_get_size(regName.c_str(), paneList[i], &numCells);
            Info << "    " << dataName.c_str() << " elements = " << numCells
                 << ", components = " << nComp << endl;
//            for(int icell=0; icell<numCells; icell++)
//            {
//                Info << "Cell " << icell << " velocity = ";
//                for(int icomp=0; icomp<nComp; icomp++)
//                {
//                    Info << *(cellVel+icell*nComp+icomp) << " ";
//                }
//                Info << endl;
//            }
        }

        dataName = string("faceToFaceMap");
        regName = volName+string(".")+dataName;
        ifCorrect = (std::find(dataItemNames.begin(),
                                    dataItemNames.end(), dataName)
                                    != dataItemNames.end());
        if (ifCorrect)
        {
            int* faceToFaceMap;
            
            COM_get_array(regName.c_str(), paneList[i], &faceToFaceMap, &nComp);
            COM_get_size(regName.c_str(), paneList[i], &numCells);
            Info << "    " << dataName.c_str() << " elements = " << numCells
                 << ", components = " << nComp << endl;
//            for(int icell=0; icell<numCells; icell++)
//            {
//                Info << "Cell " << icell << " velocity = ";
//                for(int icomp=0; icomp<nComp; icomp++)
//                {
//                    Info << *(cellVel+icell*nComp+icomp) << " ";
//                }
//                Info << endl;
//            }
        }


        dataName = string("patchPointToPointMap");
        regName = volName+string(".")+dataName;
        ifCorrect = (std::find(dataItemNames.begin(),
                                    dataItemNames.end(), dataName)
                                    != dataItemNames.end());
        if (ifCorrect)
        {
            int* patchPointToPointMap;
            
            COM_get_array(regName.c_str(), paneList[i], &patchPointToPointMap, &nComp);
            COM_get_size(regName.c_str(), paneList[i], &numCells);
            Info << "    " << dataName.c_str() << " elements = " << numCells
                 << ", components = " << nComp << endl;
//            for(int icell=0; icell<numCells; icell++)
//            {
//                Info << "Cell " << icell << " velocity = ";
//                for(int icomp=0; icomp<nComp; icomp++)
//                {
//                    Info << *(cellVel+icell*nComp+icomp) << " ";
//                }
//                Info << endl;
//            }
        }


        dataName = string("patchFaceToFaceMap");
        regName = volName+string(".")+dataName;
        ifCorrect = (std::find(dataItemNames.begin(),
                                    dataItemNames.end(), dataName)
                                    != dataItemNames.end());
        if (ifCorrect)
        {
            int* patchFaceToFaceMap;
            
            COM_get_array(regName.c_str(), paneList[i], &patchFaceToFaceMap, &nComp);
            COM_get_size(regName.c_str(), paneList[i], &numCells);
            Info << "    " << dataName.c_str() << " elements = " << numCells
                 << ", components = " << nComp << endl;
//            for(int icell=0; icell<numCells; icell++)
//            {
//                Info << "Cell " << icell << " velocity = ";
//                for(int icomp=0; icomp<nComp; icomp++)
//                {
//                    Info << *(cellVel+icell*nComp+icomp) << " ";
//                }
//                Info << endl;
//            }
        }
        //---------------------------------------



        // Field data ^^^^^^^^^^^^^^^^^^^^^^^^^^^
        dataName = string("vel");
        regName = volName+string(".")+dataName;
        ifCorrect = (std::find(dataItemNames.begin(),
                                    dataItemNames.end(), dataName)
                                    != dataItemNames.end());
        if (ifCorrect)
        {
            COM_get_array(regName.c_str(), paneList[i], &cellVel, &nComp);
            COM_get_size(regName.c_str(), paneList[i], &numCells);
            Info << "    " << dataName.c_str() << " elements = " << numCells
                 << ", components = " << nComp << endl;
//            for(int icell=0; icell<numCells; icell++)
//            {
//                Info << "Cell " << icell << " velocity = ";
//                for(int icomp=0; icomp<nComp; icomp++)
//                {
//                    Info << *(cellVel+icell*nComp+icomp) << " ";
//                }
//                Info << endl;
//            }
        }
    
        dataName = string("pres");
        regName = volName+string(".")+dataName;
        ifCorrect = (std::find(dataItemNames.begin(),
                                    dataItemNames.end(), dataName)
                                    != dataItemNames.end());
        if (ifCorrect)
        {
            COM_get_array(regName.c_str(), paneList[i], &cellPres, &nComp);
            COM_get_size(regName.c_str(), paneList[i], &numCells);
            Info << "    " << dataName.c_str() << " elements = " << numCells
                 << ", components = " << nComp << endl;
//            for(int icell=0; icell<numCells; icell++)
//            {
//                Info << "Cell " << icell << " pressure = ";
//                for(int icomp=0; icomp<nComp; icomp++)
//                {
//                    Info << *(cellPres+icell*nComp+icomp) << " ";
//                }
//                Info << endl;
//            }
        }

    
        dataName = string("temp");
        regName = volName+string(".")+dataName;
        ifCorrect = (std::find(dataItemNames.begin(),
                                    dataItemNames.end(), dataName)
                                    != dataItemNames.end());
        if (ifCorrect)
        {
            COM_get_array(regName.c_str(), paneList[i], &cellTemp, &nComp);
            COM_get_size(regName.c_str(), paneList[i], &numCells);
            Info << "    " << dataName.c_str() << " elements = " << numCells
                 << ", components = " << nComp << endl;
//            for(int icell=0; icell<numCells; icell++)
//            {
//                Info << "Cell " << icell << " temperature = ";
//                for(int icomp=0; icomp<nComp; icomp++)
//                {
//                    Info << *(cellTemp+icell*nComp+icomp) << " ";
//                }
//                Info << endl;
//            }
        }


        dataName = string("rho");
        regName = volName+string(".")+dataName;
        ifCorrect = (std::find(dataItemNames.begin(),
                                    dataItemNames.end(), dataName)
                                    != dataItemNames.end());
        if (ifCorrect)
        {
            COM_get_array(regName.c_str(), paneList[i], &cellRho, &nComp);
            COM_get_size(regName.c_str(), paneList[i], &numCells);
            Info << "    " << dataName.c_str() << " elements = " << numCells
                 << ", components = " << nComp << endl;
//            for(int icell=0; icell<numCells; icell++)
//            {
//                Info << "Cell " << icell << " density = ";
//                for(int icomp=0; icomp<nComp; icomp++)
//                {
//                    Info << *(cellRho+icell*nComp+icomp) << " ";
//                }
//                Info << endl;
//            }
        }

        dataName = string("owner");
        regName = volName+string(".")+dataName;
        ifCorrect = (std::find(dataItemNames.begin(),
                                    dataItemNames.end(), dataName)
                                    != dataItemNames.end());
        if (ifCorrect)
        {
            COM_get_array(regName.c_str(), paneList[i], &Owner, &nComp);
            COM_get_size(regName.c_str(), paneList[i], &numFaces);
            Info << "    " << dataName.c_str() << " elements = " << numFaces
                 << ", components = " << nComp << endl;
//            for(int iface=0; iface<numFaces; iface++)
//            {
//                Info << "Face " << iface << " owner = ";
//                for(int icomp=0; icomp<nComp; icomp++)
//                {
//                    Info << *(Owner+iface*nComp+icomp) << " ";
//                }
//                Info << endl;
//            }
        }

        dataName = string("neighbor");
        regName = volName+string(".")+dataName;
        ifCorrect = (std::find(dataItemNames.begin(),
                                    dataItemNames.end(), dataName)
                                    != dataItemNames.end());
        if (ifCorrect)
        {
            COM_get_array(regName.c_str(), paneList[i], &Neighb, &nComp);
            COM_get_size(regName.c_str(), paneList[i], &numFaces);
            Info << "    " << dataName.c_str() << " elements = " << numFaces
                 << ", components = " << nComp << endl;
//            for(int iface=0; iface<numFaces; iface++)
//            {
//                Info << "Face " << iface << " neighbor = ";
//                for(int icomp=0; icomp<nComp; icomp++)
//                {
//                    Info << *(Neighb+iface*nComp+icomp) << " ";
//                }
//                Info << endl;
//            }
        }

        Info << "  --------------------------------------------------"
             << endl << endl;
    }

    Info << "----------------------------------------------------"
         << endl;

    COM_free_buffer(&paneList);

    return 0;
}

/*
int flowReconst(const char* name)
{
    //  Call the flow iterator ^^^^^^^^^^^^^^^^^^
    std::vector<std::string>::iterator location = std::find(
            winNames.begin(), winNames.end(), string(name));
    int index = std::distance(winNames.begin(), location);
    //std::string volName = winNames[index]+string("VOL");

    COM_call_function(flowReconstDataHandle[index], name);
    
    return 0;
}
*/

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
        double currentTime = *fluidRun;
        double timeStep{10}; //*fluidDeltaT;
        int handle{-1};

        COM_call_function
        (
            update_solutionHandle[index],
            currentTime,
            timeStep,
            handle
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

