#include "comLoadUnload.H"

MPI_Comm masterComm;
MPI_Comm newComm;

int masterRank;
int masterNProc;
bool runParallel;

char *solverType;

//  Status Variables ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
int numDataItems=0;
std::vector<std::string> dataItemNames;
int numPanes;
int* paneList;

double *Coord;
int numNodes = 0;
int nComp = 0;

int numCells = 0;
double* cellVel;
double* cellPres;
double* cellTemp;
double* cellRho;


int *Conn8;
int numConn;
int numElem;

std::vector<std::string> connNames;

char getDataItemLoc;
COM_Type getDataItemType;
int numElementNodes;
std::string getDataItemUnits;

int timeArrayLength;

//  Function Handlers ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
int flowInitHandle;
int flowStatHandle;
int flowLoopHandle;
int flowStepHandle;
int flowFinHandle;
int flowExtractDataHandle;

// Registered veriables with COM ^^^^^^^^^^^^^^^^^^^^^^^^^^
int *fluidRun;
double *fluidTime;
double *fluidDeltaT;

// Declaration of functions ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
int comDrvInit(int argc, char *argv[]);
int comGetFunctionHandles();
int comExtractData();
int comGetDataItems(const char *name);
int comDrvStat();
int comDrvLoop();
int comDrvStep();
int comDrvFin();

int main(int argc, char *argv[])
{
    comDrvInit(argc, argv);
    //comDrvStat();
    //comDrvLoop();
    comDrvStep();
    comDrvFin();

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
    
    comfoam_load_module("ROCFOAM", solverType);

    comGetFunctionHandles();

    //  Make a dummy argc/argv for OpenFOAM. ^^^^
    //  No options passed from the command
    //  line will be used by the driver

    int myArgc = 1;
    char *myArgv[2];
    myArgv[0] = argv[0];
    myArgv[1] = NULL;
    myArgv[1] = NULL;

    if (runParallel)
    {
        myArgc = 2;
        myArgv[1] = const_cast<char *>("-parallel");
    }

    //int verb=3;

    //  Fluid initializer ^^^^^^^^^^^^^^^^^^^^^^^
    COM_call_function(flowInitHandle, &myArgc, &myArgv, "ROCFOAM");

    comGetDataItems("ROCFOAM");
  
    return 0;
}

int comGetFunctionHandles()
{
    //  Get the handle for the initialize function ^^^^^^^^
    std::string name = "ROCFOAM"+string("VOL");

    std::string functionName = name+string(".flowInit");
    flowInitHandle = COM_get_function_handle(functionName.c_str());
    if (flowInitHandle <= 0)
    { // fail
        std::cout << "rocFoam.main: Could not get handle for initialize."
                  << std::endl;
        return -2;
    }
    else
    {
        if (masterRank==0)
        {
            std::cout << "rocFoam.main: Acquired a handle for initialize."
                      << std::endl;    
        }
    }

    //  Get the handle for the loop function ^^^^^^^^^^^^^^
    functionName = name+string(".flowLoop");
    flowLoopHandle = COM_get_function_handle(functionName.c_str());
    if (flowLoopHandle <= 0)
    { // fail
        std::cout << "rocFoam.main: Could not get handle for loop."
                  << std::endl;
        return -2;
    }
    else
    {
        if (masterRank==0)
        {
            std::cout << "rocFoam.main: Acquired a handle for loop."
                      << std::endl;
        }
    }

    //  Get the handle for the step function ^^^^^^^^^^^^^^
    functionName = name+string(".flowStep");
    flowStepHandle = COM_get_function_handle(functionName.c_str());
    if (flowStepHandle <= 0)
    { // fail
        std::cout << "rocFoam.main: Could not get handle for step."
                  << std::endl;
        return -2;
    }
    else
    {
        if (masterRank==0)
        {
            std::cout << "rocFoam.main: Acquired a handle for step."
                      << std::endl;
        }
    }


    //  Get the handle for the step function ^^^^^^^^^^^^^^
    functionName = name+string(".flowExtractData");
    flowExtractDataHandle =
        COM_get_function_handle(functionName.c_str());

    if (flowExtractDataHandle <= 0)
    { // fail
        std::cout << "rocFoam.main: Could not get handle for flowExtractData."
                  << std::endl;
        return -2;
    }
    else
    {
        if (masterRank==0)
        {
            std::cout << "rocFoam.main: Acquired a handle for flowExtractData."
                      << std::endl;
            std::cout << std::endl;
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

int comGetDataItems(const char *name)
{
    std::string volName = name+string("VOL");

    std::string output;
    COM_get_dataitems(volName.c_str(), &numDataItems, output);
    std::istringstream Istr(output);

    Info << "rocFoam.main: numDataItems "
         << numDataItems
         << endl;

    for (int i=0; i<numDataItems; ++i)
    {
        std::string nameTmp;
        Istr >> nameTmp;
        dataItemNames.push_back(nameTmp);
        std::cout << "rocFoam.main: DataItem # "
                  << i << ": " << nameTmp << std::endl;
    }
    Info << endl;

    int *nProcReg;
    std::string dataName = volName+string(".winNProc");
    COM_get_array(dataName.c_str(), 0, &nProcReg);
    Info << "rocFoam.main: The ommunicator uses "
         << *nProcReg << " prcesses"
         << endl;

    dataName = volName+string(".winRun");
    COM_get_array(dataName.c_str(), 0, &fluidRun);
    Info << "rocFoam.main: Stepping status = " 
         << *fluidRun
         << endl;

    dataName = volName+string(".winTime");
    COM_get_array(dataName.c_str(), 0, &fluidTime);
    Info << "rocFoam.main: Simulation time = "
         << *fluidTime
         << endl;

    dataName = volName+string(".winDeltaT");
    COM_get_array(dataName.c_str(), 0, &fluidDeltaT);
    Info << "rocFoam.main: Simulation Time step = "
         << *fluidDeltaT
         << endl;

    dataName = volName+string(".nc");
    COM_get_array(dataName.c_str(), 1, &Coord, &nComp);
    COM_get_size(dataName.c_str(), 1, &numNodes);
    
    Info << "rocFoam.main: Domain has "
         << numNodes << " points and "
         << nComp << " components."
         << endl;
//    for(int ipoint=0; ipoint<numNodes; ipoint++)
//    {
//        Info << "Node " << ipoint << " coordinates = ";
//        for(int icomp=0; icomp<nComp; icomp++)
//        {
//            Info << *(Coord+ipoint*nComp+icomp) << " ";
//        }
//        Info << endl;
//    }


    dataName = volName+string(".:q8");
    COM_get_array(dataName.c_str(), 1, &Conn8, &nComp);
    COM_get_size(dataName.c_str(), 1, &numElem);
    
    Info << "rocFoam.main: Domain has "
         << numElem << " connectivities and "
         << nComp << " components."
         << endl;

/*
    for(int icell=0; icell<numElem; icell++)
    {
        Info << "Cell " << icell << " connectivities = ";
        for(int icomp=0; icomp<nComp; icomp++)
        {
            Info << *(Conn8+icell*nComp+icomp) << " ";
        }
        Info << endl;
    }
*/


    dataName = volName+string(".vel");
    COM_get_array(dataName.c_str(), 1, &cellVel, &nComp);
    COM_get_size(dataName.c_str(), 1, &numCells);

    Info << "rocFoam.main: Velocity has "
         << numCells << " cells and "
         << nComp << " components."
         << endl;
//    for(int icell=0; icell<numCells; icell++)
//    {
//        Info << "Cell " << icell << " velocity = ";
//        for(int icomp=0; icomp<nComp; icomp++)
//        {
//            Info << *(cellVel+icell*nComp+icomp) << " ";
//        }
//        Info << endl;
//    }
    
    dataName = volName+string(".pres");
    COM_get_array(dataName.c_str(), 1, &cellPres, &nComp);
    COM_get_size(dataName.c_str(), 1, &numCells);

    Info << "rocFoam.main: Pressure has "
         << numCells << " cells and "
         << nComp << " components."
         << endl;
    for(int icell=0; icell<numCells; icell++)
    {
        Info << "Cell " << icell << " pressure = ";
        for(int icomp=0; icomp<nComp; icomp++)
        {
            Info << *(cellPres+icell*nComp+icomp) << " ";
        }
        Info << endl;
    }

    
    dataName = volName+string(".temp");
    COM_get_array(dataName.c_str(), 1, &cellTemp);
    COM_get_size(dataName.c_str(), 1, &numCells);

    Info << "rocFoam.main: Temperature has "
         << numCells << " cells and "
         << nComp << " components."
         << endl;
//    for(int icell=0; icell<numCells; icell++)
//    {
//        Info << "Cell " << icell << " temperature = ";
//        for(int icomp=0; icomp<nComp; icomp++)
//        {
//            Info << *(cellTemp+icell*nComp+icomp) << " ";
//        }
//        Info << endl;
//    }


    dataName = volName+string(".rho");
    COM_get_array(dataName.c_str(), 1, &cellRho);
    COM_get_size(dataName.c_str(), 1, &numCells);

    Info << "rocFoam.main: Density has "
         << numCells << " cells and "
         << nComp << " components."
         << endl;
//    for(int icell=0; icell<numCells; icell++)
//    {
//        Info << "Cell " << icell << " density = ";
//        for(int icomp=0; icomp<nComp; icomp++)
//        {
//            Info << *(cellRho+icell*nComp+icomp) << " ";
//        }
//        Info << endl;
//    }

    return 0;
}

int comExtractData()
{
    //  Call the flow iterator ^^^^^^^^^^^^^^^^^^
    COM_call_function(flowExtractDataHandle, "ROCFOAM");
    
    return 0;
}


int comDrvLoop()
{
    //  Call the flow iterator ^^^^^^^^^^^^^^^^^^
    COM_call_function(flowLoopHandle, "ROCFOAM");
    
    return 0;
}

int comDrvStep()
{
    //  Call the flow stepper ^^^^^^^^^^^^^^^^^^

    Info << "\nStarting time loop\n" << endl;

    //while (*fluidRun)
    {
        COM_call_function(flowStepHandle, "ROCFOAM");
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

int comDrvFin()
{
    //  Call the flow unloader ^^^^^^^^^^^^^^^^^^
    //COM_UNLOAD_MODULE_STATIC_DYNAMIC(comfoam, "ROCFOAM");
    std::string name = "ROCFOAM";
    comfoam_unload_module(name.c_str(), solverType);

    COM_set_default_communicator(masterComm);
    
    COM_finalize();

    MPI_Barrier(masterComm);
    MPI_Finalize();
    
    return 0;
}    

int comDrvStat()
{
    //  Get information about what was ^^^^^^^^^^
    //  registered in this window  
    std::string name = "ROCFOAM"+string("VOL");

    std::string output;
    COM_get_dataitems(name.c_str(), &numDataItems, output);
    std::istringstream Istr(output);

    std::cout << "rocFoam.main: numDataItems "
              << numDataItems << std::endl;

    for (int i=0; i<numDataItems; ++i)
    {
        std::string name;
        Istr >> name;
        dataItemNames.push_back(name);
        std::cout << "rocFoam.main: DataItem # "
                  << i << ": " << name << std::endl;
    }

    //  List of panes in this window ^^^^^^^^^^^^
    COM_get_panes(name.c_str(), &numPanes, &paneList);
    std::cout << "rocFoam.main: Number of Panes "
              << numPanes << std::endl;

    for (int i=0; i<numPanes; ++i) 
    {
        std::cout << "rocFoam.main: Pane ID # "
                  << i+1 << "="<< paneList[i] << std::endl;
    }

    //  Only one pane for serial runs ^^^^^^^^^^^
    int pane = paneList[0];

    //  Get for grid coordinates ^^^^^^^^^^^^^^^^
    std::string dataName = name+string(".nc");
    COM_get_array(dataName.c_str(), pane, &Coord);

    //  Check for expected number of nodes ^^^^^^
    COM_get_size(dataName.c_str(), pane, &numNodes);

    //  Get connectivity tables for panes ^^^^^^^
    std::string stringNames;
    COM_get_connectivities(name.c_str(), pane, &numConn, stringNames);
    std::istringstream ConnISS(stringNames);

    for (int i=0; i<numConn; ++i)
    {
        std::string nameTmp;
        ConnISS >> nameTmp;
        connNames.push_back(nameTmp);
        std::cout << "rocFoam.main: Connectivity Table # "
                  << i+1 << ": " << nameTmp << std::endl;
    }

    //  Number of nodes per element ^^^^^^^^^^^^^
    dataName = name+connNames[0];
    std::string fullConnName = dataName.c_str();
    COM_get_dataitem
    (
        fullConnName,
        &getDataItemLoc,
        &getDataItemType, 
        &numElementNodes,
        &getDataItemUnits
    );

    std::cout << "rocFoam.main: getDataItemLoc "
              << getDataItemLoc << std::endl;

    std::cout << "rocFoam.main: getDataItemType "
              << getDataItemType << std::endl;

    std::cout << "rocFoam.main: numElementNodes "
              << numElementNodes << std::endl;

    std::cout << "rocFoam.main: getDataItemUnits "
              << getDataItemUnits << std::endl;
    
    COM_get_array(fullConnName.c_str(), pane, &Conn8);
    COM_get_size(fullConnName, pane, &numElem);


    std::cout << "rocFoam.main: Conn numElem "
              << numElem << std::endl;

    //  Put elements into a vector so we can ^^^^
    //  build the solver agent 
    std::vector<unsigned int> connVector;
    for (int i=0; i<numElem; ++i)
    {
        for (int j=0; j<numElementNodes; ++j)
        {
            connVector.push_back((Conn8[i*numElementNodes+j]));
        }
    }

    //  Get non-mesh data items ^^^^^^^^^^^^^^^^^
    dataName = name+string(".time");
    
    COM_get_dataitem
    (
        dataName,
        &getDataItemLoc,
        &getDataItemType, 
        &timeArrayLength,
        &getDataItemUnits
    );

    std::cout << "rocFoam.main: timeArrayLength "
              << timeArrayLength << std::endl;

    return 0;
}
