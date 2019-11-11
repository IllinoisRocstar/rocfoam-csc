#include "comFoam.H"

COM_EXTERN_MODULE(comfoam);

MPI_Comm masterComm;
MPI_Comm newComm;

int masterRank;
int masterNProc;
bool runParallel;

std::string solverType;

//  Status Variables ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
int numDataItems=0;
std::vector<std::string> dataItemNames;
int numPanes;
int* paneList;

double* Coord;
int numNodes = 0;

int numConn;
int numElem;
int* Conn;
std::vector<std::string> connNames;

char getDataItemLoc;
COM_Type getDataItemType;
int numElementNodes;
std::string getDataItemUnits;

int timeArrayLength;

//  Function Handlers ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
int flowInitHandle;
int flowStatHandle;
int flowLoopHandle;
int flowFinHandle;

void comDrvInit(int argc, char *argv[]);
void comDrvStat();
void comDrvRun();
void comDrvFin();

int main(int argc, char *argv[])
{
    comDrvInit(argc,argv);
    //comDrvStat();
    comDrvRun();
    comDrvFin();

    return 0;
}


void comDrvInit(int argc, char *argv[])
{

    MPI_Init(&argc, &argv);
    masterComm = MPI_COMM_WORLD;

    MPI_Comm_rank(masterComm, &masterRank);
    MPI_Comm_size(masterComm, &masterNProc);


    if (masterRank==0)
    {
        std::cout << "comFoam.main: Setting up communicator..."
                  << std::endl;

        std::cout << "comFoam.main:Rank " << masterRank
                  << ", NPROC = " << masterNProc << ", "
                  << ", Init COMM = " << masterComm
                  << std::endl;
    }

    COM_init(&argc, &argv);

    // A new communicator can be generated for
    //   openfoam solver
    newComm = masterComm;

    // Run in parallel mode?
    runParallel = false;
    solverType = string("rocRhoCentral");
    

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
                solverType = string("rocRhoCentral");
            }
            else
            {
                if (masterRank==0)
                {
                    std::cout << "comFoam.Main: Unknown argumnet"
                              << ss.str() << std::endl;
                }
                throw -1;
            }

        }
    }

    if (runParallel && masterNProc > 1)
    {
        if (masterRank==0)
        {
            std::cout << "comFoam.Main: Running in PRALLEL with solver " 
                      << solverType << "." << std::endl;
        }
    }
    else
    {
        runParallel = false;

        if (masterRank==0)
        {
            std::cout << "comFoam.Main: Running in SERIAL with solver " 
                      << solverType << "." << std::endl;
        }
    }
    if (masterRank==0) std::cout << std::endl;

    if (!runParallel && masterNProc > 1)
    {
        if (masterRank==0)
        {
            std::cout << "comFoam.Main: NProc>1 detected for a serial job."
                      << std::endl;
            throw -1;
        }    
    
    }


    //  Setting the defual communicator. Is it needed?
    COM_set_default_communicator(newComm);
    if (masterRank==0)
    {
        std::cout << "comFoam.Main: New COMM = "
                  << newComm << std::endl;
        std::cout << std::endl;
    }

    COM_LOAD_MODULE_STATIC_DYNAMIC(comfoam, "CFModule");   

    // getting number of processes  
    if (masterRank==0)
    {
        int *nProcReg;
        COM_get_array("CFModule.winNProc", 0, &nProcReg);
        std::cout << "The communicator registered in OFModule uses "
                  << *nProcReg << " prcesses"
                  << std::endl;

        std::cout << std::endl;
    }

    //  Get the handle for the initialize function ^^^^^^^^
    flowInitHandle = COM_get_function_handle("CFModule.flowInit");
    if (flowInitHandle <= 0)
    { // fail
        std::cout << "comFoam.main: Could not get handle for initialize."
                  << std::endl;
        throw -2;
    }
    else
    {
        std::cout << "comFoam.main: Acquired a handle for initialize."
                  << std::endl;    

        std::cout << std::endl;
    }

    //  Get the handle for the loop function ^^^^^^^^^^^^^^
    flowLoopHandle = COM_get_function_handle("CFModule.flowLoop");
    if (flowLoopHandle <= 0)
    { // fail
        std::cout << "comFoam.main: Could not get handle for loop."
                  << std::endl;
        throw -2;
    }
    else
    {
        std::cout << "comFoam.main: Acquired a handle for loop."
                  << std::endl;    

        std::cout << std::endl;
    }

    //  Get the handle for the finalize function ^^^^^^^^^^
    /*int flowFinHandle = COM_get_function_handle("CFModule.flowFin");
    if (flowFinHandle <= 0)
    { // fail
        std::cout << "comFoam.main: Could not get handle for finalize."
                  << std::endl;
        throw -2;
    }
    else
    {
        std::cout << "comFoam.main: Acquired a handle for finLauncher."
                  << std::endl;    
    }*/

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

    int verb=3;

    //  Fluid initializer ^^^^^^^^^^^^^^^^^^^^^^^
    COM_call_function(flowInitHandle, &myArgc, &myArgv, &verb);
    
    return;
}


void comDrvStat()
{
    //  Get information about what was ^^^^^^^^^^
    //  registered in this window  
    std::string output;
    COM_get_dataitems("CFModule", &numDataItems, output);
    std::istringstream Istr(output);

    std::cout << "comFoam.main: numDataItems "
              << numDataItems << std::endl;

    for (int i=0; i<numDataItems; ++i)
    {
        std::string name;
        Istr >> name;
        dataItemNames.push_back(name);
        std::cout << "comFoam.main: DataItem # "
                  << i << ": " << name << std::endl;
    }

    //  List of panes in this window ^^^^^^^^^^^^
    COM_get_panes("CFModule", &numPanes, &paneList);
    std::cout << "comFoam.main: Number of Panes "
              << numPanes << std::endl;

    for (int i=0; i<numPanes; ++i) 
    {
        std::cout << "comFoam.main: Pane ID # "
                  << i+1 << "="<< paneList[i] << std::endl;
    }

    //  Only one pane for serial runs ^^^^^^^^^^^
    int pane = paneList[0];

    //  Get for grid coordinates ^^^^^^^^^^^^^^^^
    COM_get_array("CFModule.nc", pane, &Coord);

    //  Check for expected number of nodes ^^^^^^
    COM_get_size("CFModule.nc", pane, &numNodes);

    //  Get connectivity tables for panes ^^^^^^^
    std::string stringNames;
    COM_get_connectivities("CFModule", pane, &numConn, stringNames);
    std::istringstream ConnISS(stringNames);

    for (int i=0; i<numConn; ++i)
    {
        std::string name;
        ConnISS >> name;
        connNames.push_back(name);
        std::cout << "comFoam.main: Connectivity Table # "
                  << i+1 << ": " << name << std::endl;
    }

    //  Number of nodes per element ^^^^^^^^^^^^^
    std::string fullConnName("CFModule."+connNames[0]);
    COM_get_dataitem
    (
        fullConnName,
        &getDataItemLoc,
        &getDataItemType, 
        &numElementNodes,
        &getDataItemUnits
    );

    std::cout << "comFoam.main: getDataItemLoc "
              << getDataItemLoc << std::endl;

    std::cout << "comFoam.main: getDataItemType "
              << getDataItemType << std::endl;

    std::cout << "comFoam.main: numElementNodes "
              << numElementNodes << std::endl;

    std::cout << "comFoam.main: getDataItemUnits "
              << getDataItemUnits << std::endl;
    
    COM_get_array(fullConnName.c_str(), pane, &Conn);
    COM_get_size(fullConnName, pane, &numElem);


    std::cout << "comFoam.main: Conn numElem "
              << numElem << std::endl;

    //  Put elements into a vector so we can ^^^^
    //  build the solver agent 
    std::vector<unsigned int> connVector;
    for (int i=0; i<numElem; ++i)
    {
        for (int j=0; j<numElementNodes; ++j)
        {
            connVector.push_back((Conn[i*numElementNodes+j]));
        }
    }

    //  Get non-mesh data items ^^^^^^^^^^^^^^^^^
    std::string name("CFModule.time");
    COM_get_dataitem
    (
        name,
        &getDataItemLoc,
        &getDataItemType, 
        &timeArrayLength,
        &getDataItemUnits
    );

    std::cout << "comFoam.main: timeArrayLength "
              << timeArrayLength << std::endl;

    return;
}



void comDrvRun()
{
    //  Call the flow iterator ^^^^^^^^^^^^^^^^^^
    COM_call_function(flowLoopHandle);
    
    return;
}


void comDrvFin()
{
    //  Call the flow unloader ^^^^^^^^^^^^^^^^^^
    COM_UNLOAD_MODULE_STATIC_DYNAMIC(comfoam, "CFModule");

    COM_set_default_communicator(masterComm);
    
    COM_finalize();

    MPI_Barrier(masterComm);
    MPI_Finalize();
}    


