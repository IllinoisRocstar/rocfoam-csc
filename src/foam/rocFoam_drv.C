#include "comLoadUnload.H"

COM_EXTERN_MODULE(Rocout);
COM_EXTERN_MODULE(Rocin);


MPI_Comm masterComm;
MPI_Comm newComm;

int masterRank;
int masterNProc;
bool runParallel;

char *solverType;

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

int main(int argc, char *argv[])
{

/*
for (;;)
{
std::string strTmp;
std::cin >> strTmp;
std::cout << "strTmp1 = " << strTmp << std::endl;
strTmp = comFoam::removeTrailZero(strTmp);
std::cout << "strTmp2 = " << strTmp << std::endl;
std::cin.get();
}
*/



    comDrvInit(argc, argv);
    //comDrvStat(const char *name);
    //comDrvLoop(const char *name);
    
    std::string lookUpWindow1 = winNames[1]+string("VOL");
    comGetRunStatItems(lookUpWindow1.c_str());
    comDrvStep(winNames[1].c_str());

    //comDrvFin(winNames[0].c_str());

    rocfoam_unload_module(winNames[0].c_str(), solverType);
    rocfoam_unload_module(winNames[1].c_str(), solverType);


    return 0;
}

int comDrvInit(int argc, char *argv[])
{

    MPI_Init(&argc, &argv);
    masterComm = MPI_COMM_WORLD;

    MPI_Comm_rank(masterComm, &masterRank);
    MPI_Comm_size(masterComm, &masterNProc);


    if (masterRank==0)
    {
        std::cout << "rocFoam.main: Setting up communicator..."
                  << std::endl;

        std::cout << "rocFoam.main:Rank " << masterRank
                  << ", NProc = " << masterNProc
                  << ", COMM = " << masterComm
                  << std::endl;

        std::cout << std::endl;
    }

    COM_init(&argc, &argv);

    // A new communicator can be generated for
    //   openfoam solver
    newComm = masterComm;

    // Run in parallel mode?
    runParallel = false;
    solverType = const_cast<char *>("rocRhoCentral");
    
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
            /* else
            {
                if (masterRank==0)
                {
                    std::cout << "rocFoam.main: Unknown argumnet"
                              << ss.str() << std::endl;
                }
                return -1;
            } */
        }
    }

    if (runParallel && masterNProc > 1)
    {
        if (masterRank==0)
        {
            std::cout << "rocFoam.main: Running in PRALLEL with solver " 
                      << solverType << "." << std::endl;
        }
    }
    else if (!runParallel && masterNProc > 1)
    {
        if (masterRank==0)
        {
            std::cout << "rocFoam.main: NProc>1 detected for a serial job."
                      << std::endl;
            return -1;
        }    
    
    }
    else
    {
        runParallel = false;

        if (masterRank==0)
        {
            std::cout << "rocFoam.main: Running in SERIAL with solver " 
                      << solverType << "." << std::endl;
        }
    }
    if (masterRank==0) std::cout << std::endl;

    //  Setting the defual communicator. Is it needed?
    COM_set_default_communicator(newComm);

    winNames.push_back("ROCFOAM");
    
    rocfoam_load_module(winNames[0].c_str(), solverType);
    comGetFunctionHandles(winNames[0].c_str());

    //  Make a dummy argc/argv for OpenFOAM. ^^^^
    //  No options passed from the command
    //  line will be used by the driver

    int myArgc = 1;
    char *myArgv[2];
    myArgv[0] = argv[0];
    myArgv[1] = nullptr;

    if (runParallel)
    {
        myArgc = 2;
        myArgv[1] = const_cast<char *>("-parallel");
    }

    //int verb=3;

    //  Fluid initializer ^^^^^^^^^^^^^^^^^^^^^^^
    COM_call_function(flowInitHandle[0],
                      &myArgc, &myArgv,
                      winNames[0].c_str());

    std::string lookUpWindow1 = "";
    std::string lookUpWindow2 = "";


    //lookUpWindow1 = winNames[0]+string("VOL");
    //comGetVolDataItems(lookUpWindow1.c_str());
    //lookUpWindow1 = winNames[0]+string("SURF");
    //comGetVolDataItems(lookUpWindow1.c_str());


std::cout << "Writing the window" << std::endl;

COM_LOAD_MODULE_STATIC_DYNAMIC(SimOUT, "OUT");
int OUT_write = COM_get_function_handle("OUT.write_dataitem");

lookUpWindow1 = winNames[0]+string("VOL");
std::string whatToWrite = lookUpWindow1+std::string(".mesh");
int whatToWriteHandle = COM_get_dataitem_handle(whatToWrite.c_str());

char* outputPath = new char[40]{' '};
std::string strTmp = std::string("./")
                   + winNames[0]+"/"
                    +lookUpWindow1+std::string("_");
std::strcpy(outputPath, strTmp.c_str());

char* material = new char[40]{' '};
std::strcpy(material, lookUpWindow1.c_str());

char* timeName = new char[40]{' '};

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

std::cout << "Finished writting surface window" << std::endl;
COM_UNLOAD_MODULE_STATIC_DYNAMIC(SimOUT, "OUT");
std::cout << "Unloaded SIMOUT" << std::endl;

std::cin.get();


    lookUpWindow2 = "ROCFOAM1";
    winNames.push_back(lookUpWindow2);
    rocfoam_load_module(lookUpWindow2.c_str(), solverType);
    comGetFunctionHandles(lookUpWindow2.c_str());

    lookUpWindow1 = winNames[0];
    lookUpWindow2 = winNames[1];
    comFoam::copyWindow(lookUpWindow1.c_str(), lookUpWindow2.c_str());

    COM_call_function(flowRestartInitHandle[1],
                      &myArgc, &myArgv,
                      winNames[1].c_str());
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
        return -2;
    }
    else
    {
        if (masterRank==0)
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
        return -2;
    }
    else
    {
        if (masterRank==0)
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
        return -2;
    }
    else
    {
        if (masterRank==0)
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
        return -2;
    }
    else
    {
        if (masterRank==0)
        {
            std::cout << "Acquired a handle for "
                  << functionName.c_str()
                  << std::endl;
        }
    }

    //  Get the handle for the finalize function ^^^^^^^^^^
    /*int flowFinHandle = COM_get_function_handle("ROCFOAM.flowFin");
    if (flowFinHandle <= 0)
    { // fail
        std::cout << "rocFoam.main: Could not get handle for finalize."
                  << std::endl;
        return -2;
    }
    else
    {
        std::cout << "rocFoam.main: Acquired a handle for finLauncher."
                  << std::endl;    
    }*/

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

int comDrvFin(const char* name)
{
    //  Call the flow unloader ^^^^^^^^^^^^^^^^^^
    //COM_UNLOAD_MODULE_STATIC_DYNAMIC(rocfoam, "ROCFOAM");
    rocfoam_unload_module(name, solverType);

    COM_set_default_communicator(masterComm);
    
    COM_finalize();

    MPI_Barrier(masterComm);
    MPI_Finalize();
    
    return 0;
}    

