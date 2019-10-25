#include "com.h"
#inlcude "comFoam.H"

COM_EXTERN_MODULE(RocFoam);


int main(int argc, char *argv[])
{
    COM_init(&argc, &argv);

    COM_LOAD_MODULE_STATIC_DYNAMIC( rocFoamModule, "RFModule");   

    //  Get the handle for the initialize function ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    int rfInitHandle = COM_get_function_handle("RFModule.initializeAll");
    if(rfInitHandle <= 0)
    { // fail
        std::cout << "RFModuleDriver:main: Could not get handle for init function."
                  << std::endl;
        throw -1;
    }
    else
    {
        std::cout << "RFModuleDriver:main: Acquired a handle for init function."
                  << std::endl;    
    }

    int rfLoopHandle = COM_get_function_handle("RFModule.loop");
    { // fail
        std::cout << "RFModuleDriver:main: Could not get handle for loop function."
                  << std::endl;
        throw -1;
    }
    else
    {
        std::cout << "RFModuleDriver:main: Acquired a handle for loop function."
                  << std::endl;    
    }

    int rfFinHandle = COM_get_function_handle("RFModule.finlaize");
    { // fail
        std::cout << "RFModuleDriver:main: Could not get handle for finalize function."
                  << std::endl;
        throw -1;
    }
    else
    {
        std::cout << "RFModuleDriver:main: Acquired a handle for finalize function."
                  << std::endl;    
    }

    //  Make a dummy argc/argv for OpenFOAM. ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    //  No options passed from the command
    //  line will be used by the driver
    int myArgc = 1;
    char *myArgv[2];
    myArgv[0] = argv[0];
    myArgv[1] = NULL;
    int verb=3;

    //  Call the flow initializer  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    COM_call_function(rfInitHandle, &myArgc, &myArgv, &verb);

    //  Get information about what was registered in this window  ^^^^^^^^^^^^^
    int numDataItems=0;
    std::string output;
    COM_get_dataitems("RFModule", &numDataItems, output);
    std::istringstream Istr(output);
    std::vector<std::string> dataItemNames;

    for (int i=0; i<numDataItems; ++i)
    {
        std::string name;
        Istr >> name;
        dataItemNames.push_back(name);
        std::cout << "RFModuleDriver:main: DataItem # "
                  << i << ": " << name << std::endl;
    }

    //  List of panes in this window ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    int numPanes;
    int* paneList;
    COM_get_panes("RFModule", &numPanes, &paneList);
    std::cout << "RFModuleDriver:main: Number of Panes "
              << numPanes << std::endl;

    for (int i=0; i<numPanes; ++i) 
    {
        std::cout << "RFModuleDriver:main: Pane ID # "
                  << i+1 << "="<< paneList[i] << std::endl;
    }

    //  Only one pane for serial runs ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    int pane = paneList[0];

    double* Coord;
    COM_get_array("RFModule.nc", pane, &Coord);

    //  Check for expected number of nodes ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    int numNodes = 0;
    COM_get_size("RFModule.nc", pane, &numNodes);

    //  Get connectivity tables for panes ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    int numConn;
    std::string stringNames;
    COM_get_connectivities("RFModule", pane, &numConn, stringNames);
    std::istringstream ConnISS(stringNames);
    std::vector<std::string> connNames;

    for (int i=0; i<numConn; ++i)
    {
        std::string name;
        ConnISS >> name;
        connNames.push_back(name);
        std::cout << "RFModuleDriver:main: Connectivity Table # "
                  << i+1 << ": " << name << std::endl;
    }

    //  Number of nodes per element ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    char getDataItemLoc;
    COM_Type getDataItemType;
    int numElementNodes;
    std::string getDataItemUnits;
    std::string fullConnName("RFModule."+connNames[0]);
    COM_get_dataitem
    (
        fullConnName,
        &getDataItemLoc,
        &getDataItemType, 
        &numElementNodes,
        &getDataItemUnits
    );

    std::cout << "RFModuleDriver:main: getDataItemLoc "
              << getDataItemLoc << std::endl;

    std::cout << "RFModuleDriver:main: getDataItemType "
              << getDataItemType << std::endl;

    std::cout << "RFModuleDriver:main: numElementNodes "
              << numElementNodes << std::endl;

    std::cout << "RFModuleDriver:main: getDataItemUnits "
              << getDataItemUnits << std::endl;

    int* Conn;
    int numElem;
    COM_get_array(fullConnName.c_str(), pane, &Conn);
    COM_get_size(fullConnName, pane, &numElem);


    std::cout << "RFModuleDriver:main: Conn numElem "
              << numElem << std::endl;

    //  Put elements into a vector so we can build the solver agent ^^^^^^^^^^^
    std::vector<unsigned int> connVector;
    for (int i=0; i<numElem; ++i)
    {
        for (int j=0; j<numElementNodes; ++j)
        {
            connVector.push_back((Conn[i*numElementNodes+j]));
        }
    }

    //  Get non-mesh data items ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    int arrayLength;
    std::string name("RFModule.time");
    COM_get_dataitem
    (
        name,
        &getDataItemLoc,
        &getDataItemType, 
        &arrayLength,
        &getDataItemUnits
    );

    std::cout << "RFModuleDriver:main: timeArrayLength "
              << arrayLength << std::endl;

    double *time;
    COM_get_array(name.c_str(), pane, &time);

    name = "RFModule.endTime";
    COM_get_dataitem
    (
        name,
        &getDataItemLoc,
        &getDataItemType, 
        &arrayLength,
        &getDataItemUnits
    );

    double *endTime;
    COM_get_array
    (
        name.c_str(),
        pane,
        &endTime
    );

    std::cout << "RFModuleDriver:main: endTime "
              << *endTime << std::endl;


    //  Call the flow iterator ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    //COM_call_function(rfLoopHandle);
    
    //  Call the flow finalize function to deallocate pointers ^^^^^^^^^^^^^^^^
    COM_call_function(rfFinHandle);

    //  Unload rocFoam  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    COM_UNLOAD_MODULE_STATIC_DYNAMIC(rocFoamModule, "RFModule");
    COM_finalize();
    
    return 0;
}
