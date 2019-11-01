#include "comFoam.H"

COM_EXTERN_MODULE(comfoam);

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

std::cout << "HERE 1" << std::endl;

    comDrvInit(argc,argv);
    
//    comDrvStat();

    comDrvRun();

    comDrvFin();

     return 0;
}


void comDrvInit(int argc, char *argv[])
{

    COM_init(&argc, &argv);

    COM_LOAD_MODULE_STATIC_DYNAMIC(comfoam, "CFModule");   

    //  Get the handle for the initialize function ^^^^^^^^
    flowInitHandle = COM_get_function_handle("CFModule.flowInit");
    if (flowInitHandle <= 0)
    { // fail
        std::cout << "CFModuleDriver:main: Could not get handle for initialize."
                  << std::endl;
        throw -1;
    }
    else
    {
        std::cout << "CFModuleDriver:main: Acquired a handle for initialize."
                  << std::endl;    
    }

    //  Get the handle for the loop function ^^^^^^^^^^^^^^
    flowLoopHandle = COM_get_function_handle("CFModule.flowLoop");
    if (flowLoopHandle <= 0)
    { // fail
        std::cout << "CFModuleDriver:main: Could not get handle for loop."
                  << std::endl;
        throw -1;
    }
    else
    {
        std::cout << "CFModuleDriver:main: Acquired a handle for loop."
                  << std::endl;    
    }

    //  Get the handle for the finalize function ^^^^^^^^^^
    /*int flowFinHandle = COM_get_function_handle("CFModule.flowFin");
    if (flowFinHandle <= 0)
    { // fail
        std::cout << "CFModuleDriver:main: Could not get handle for finalize."
                  << std::endl;
        throw -1;
    }
    else
    {
        std::cout << "CFModuleDriver:main: Acquired a handle for finLauncher."
                  << std::endl;    
    }*/

    //  Make a dummy argc/argv for OpenFOAM. ^^^^
    //  No options passed from the command
    //  line will be used by the driver
    int myArgc = 1;
    char *myArgv[2];
    myArgv[0] = argv[0];
    myArgv[1] = NULL;
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

    std::cout << "CFModuleDriver:main: numDataItems "
              << numDataItems << std::endl;

int j=200;
std::cout << "HERE " << j++ << std::endl;


    for (int i=0; i<numDataItems; ++i)
    {
        std::string name;
        Istr >> name;
        dataItemNames.push_back(name);
        std::cout << "CFModuleDriver:main: DataItem # "
                  << i << ": " << name << std::endl;
    }


std::cout << "HERE " << j++ << std::endl;

    //  List of panes in this window ^^^^^^^^^^^^
    COM_get_panes("CFModule", &numPanes, &paneList);
    std::cout << "CFModuleDriver:main: Number of Panes "
              << numPanes << std::endl;

    for (int i=0; i<numPanes; ++i) 
    {
        std::cout << "CFModuleDriver:main: Pane ID # "
                  << i+1 << "="<< paneList[i] << std::endl;
    }


if (paneList==NULL)
{
std::cout << "Not Allocated " << j++ << std::endl;
}
else
{
std::cout << "Not Allocated " << j++ << std::endl;
}



std::cout << "HERE " << j++ << std::endl;

    //  Only one pane for serial runs ^^^^^^^^^^^
    int pane = paneList[0];

std::cout << "HERE " << j++ << std::endl;

    //  Get for grid coordinates ^^^^^^^^^^^^^^^^
    COM_get_array("CFModule.nc", pane, &Coord);

std::cout << "HERE " << j++ << std::endl;

    //  Check for expected number of nodes ^^^^^^
    COM_get_size("CFModule.nc", pane, &numNodes);

std::cout << "HERE " << j++ << std::endl;

    //  Get connectivity tables for panes ^^^^^^^
    std::string stringNames;
    COM_get_connectivities("CFModule", pane, &numConn, stringNames);
    std::istringstream ConnISS(stringNames);

std::cout << "HERE " << j++ << std::endl;

    
    for (int i=0; i<numConn; ++i)
    {
        std::string name;
        ConnISS >> name;
        connNames.push_back(name);
        std::cout << "CFModuleDriver:main: Connectivity Table # "
                  << i+1 << ": " << name << std::endl;
    }

std::cout << "HERE " << j++ << std::endl;

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

std::cout << "HERE " << j++ << std::endl;


    std::cout << "CFModuleDriver:main: getDataItemLoc "
              << getDataItemLoc << std::endl;

    std::cout << "CFModuleDriver:main: getDataItemType "
              << getDataItemType << std::endl;

    std::cout << "CFModuleDriver:main: numElementNodes "
              << numElementNodes << std::endl;

    std::cout << "CFModuleDriver:main: getDataItemUnits "
              << getDataItemUnits << std::endl;

std::cout << "HERE " << j++ << std::endl;

    
    COM_get_array(fullConnName.c_str(), pane, &Conn);
    COM_get_size(fullConnName, pane, &numElem);


    std::cout << "CFModuleDriver:main: Conn numElem "
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

    std::cout << "CFModuleDriver:main: timeArrayLength "
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

    COM_UNLOAD_MODULE_STATIC_DYNAMIC(comfoam, "CFModule");
    
    COM_finalize();

}    


