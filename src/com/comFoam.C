#include "comFoam.H"

COM_EXTERN_MODULE(comFoam);

comFoamModule::comFoamModule()
{
    std::cout << "comFoam cunstructor = "
              << test_bool << " " << listOptions_ << std::endl;
}

comFoamModule::comFoamModule(int argc, char *argv[])
{
    std::cout << "comFoam cunstructor = "
              << test_bool << " " << listOptions_ << std::endl;

    comInitialize(argc, argv);
}



//^^^^^ INITIALIZE MODULES ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
int comFoamModule::comInitialize(int argc, char *argv[])
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
    flowFinHandle = COM_get_function_handle("CFModule.flowFin");
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
    }


    /*
    flowStatHandle = COM_get_function_handle("CFModule.statFlow");
    if (flowStatHandle <= 0)
    { // fail
        std::cout << "CFModuleDriver:main: Could not get handle for status."
                  << std::endl;
        throw -1;
    }
    else
    {
        std::cout << "CFModuleDriver:main: Acquired a handle for status."
                  << std::endl;    
    }
    */


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
}


int comFoamModule::flowInit(int *pargc, void **pargv, int *verbIn)
{
    int argc = *pargc;
    char** argv = (char**)(pargv);
    verbosity = *verbIn;

    std::cout << "OF Module verbosity = " << verbosity << std::endl;
    
    //  OpenFOAM initializer ^^^^^^^^^^^^^^^^^^^^
    initialize(argc, argv);
    
    //  Other initializations ^^^^^^^^^^^^^^^^^^^

    return 0;
}
//-----------------------------------------------------------------------------


//^^^^^ LOOP MODULES ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
int comFoamModule::comLoop()
{
    //  Call the flow iterator ^^^^^^^^^^^^^^^^^^
    COM_call_function(flowLoopHandle);
    
    return 0;
}

int comFoamModule::flowLoop()
{
    //  Call the flow iterator ^^^^^^^^^^^^^^^^^^
    loop();
    
    return 0;
}
//-----------------------------------------------------------------------------



//^^^^^ FINALIZE MODULES ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
int comFoamModule::comFinalize()
{
    //  Unload comFoam  ^^^^^^^^^^^^^^^^^^^^^^^^^
    COM_UNLOAD_MODULE_STATIC_DYNAMIC(comfoam, "CFModule");
}

void comFoamModule::Unload(const std::string &name)
{
    std::cout << "RFModule.Unload: Unloading comFoamModule with name "
              << name << "." << std::endl;
              
    //comFoamModule *comFoamPtr = NULL;
    //std::string globalName(name+".global");

    //COM_get_object(globalName.c_str(), 0, &comFoamPtr);

    //COM_assertion_msg(comFoamPtr->validate_object()==0, "Invalid object");
    //delete comFoamPtr;

    //int &flowFinHandle = comFoamModule::flowFinHandle
    int flowFinHandle_;
    flowFinHandle_= comFoamModule::flowFinHandle;

    
    // Deallocate flow poniters ^^^^^^^^^^^^^^^^^
    COM_call_function(*flowFinHandle_);
    
    // Deallocate other poniters ^^^^^^^^^^^^^^^^


    COM_delete_window(std::string(name));
    COM_finalize();
}

int comFoamModule::flowFin()
{
    //  Call the flow iterator ^^^^^^^^^^^^^^^^^^
    finalize();
    
    return 0;
}
//-----------------------------------------------------------------------------


void comFoamModule::Load(const char *name)
{
    std::cout << "RFModule.Load: Loading comFoamModule with name "
              << name << "." << std::endl;

    //  Register module with COM ^^^^^^^^^^^
    //comFoamModule *comFoamPtr = new comFoamModule();

//    std::string &windowName = comFoamModule::windowName;
//    std::string &globalName = comFoamModule::globalName;

    std::string *(comFoamModule::windowName);
    std::string *(comFoamModule::globalName);


    COM_new_window(name, MPI_COMM_NULL);
    windowName = name;

    globalName = name + string(".global");

    COM_new_dataitem(globalName.c_str(), 'w', COM_VOID, 1, "");

    COM_set_object(globalName.c_str(), 0, comFoam);


    /// Register functions
    std::vector<COM_Type> types(13,COM_INT);

    types[0] = COM_RAWDATA;
    types[2] = COM_VOID;

    COM_set_member_function
    (
        (name + string(".flowInit")).c_str(),
        (Member_func_ptr)(&comFoamModule::flowInit),
        globalName.c_str(), "biii", &types[0]
    );


    COM_set_member_function
    (
        (name + string(".flowStatus")).c_str(),
        (Member_func_ptr)(&comFoamModule::flowStatus),
        globalName.c_str(), "b", &types[0]
    );

    COM_set_member_function
    (
        (name + string(".flowLoop")).c_str(),
        (Member_func_ptr)(&comFoamModule::flowLoop),
        globalName.c_str(), "b", &types[0]
    );

    COM_set_member_function
    (
        (name + string(".flowFin")).c_str(),
        (Member_func_ptr)(&comFoamModule::flowFin),
        globalName.c_str(), "b", &types[0]
    );

    COM_window_init_done(name); 

    return;
}



/// @brief C/C++ bindings to load IcoFoamModule
extern "C" void comfoam_load_module(const char *name)
{
  comFoamModule::Load(name);
}

/// @brief C/C++ bindings to unload IcoFoamModule
extern "C" void comfoam_unload_module(const char *name)
{
  comFoamModule::Unload(name);
}









int comFoamModule::comStatus()
{
    COM_call_function(flowStatHandle);
    
    return 0;
}


int comFoamModule::flowStat()
{

    //  Get information about what was ^^^^^^^^^^
    //  registered in this window  
    std::string output;
    COM_get_dataitems("CFModule", &numDataItems, output);
    std::istringstream Istr(output);

    std::cout << "CFModuleDriver:main: numDataItems "
              << numDataItems << std::endl;

    for (int i=0; i<numDataItems; ++i)
    {
        std::string name;
        Istr >> name;
        dataItemNames.push_back(name);
        std::cout << "CFModuleDriver:main: DataItem # "
                  << i << ": " << name << std::endl;
    }

    //  List of panes in this window ^^^^^^^^^^^^
    COM_get_panes("CFModule", &numPanes, &paneList);
    std::cout << "CFModuleDriver:main: Number of Panes "
              << numPanes << std::endl;

    for (int i=0; i<numPanes; ++i) 
    {
        std::cout << "CFModuleDriver:main: Pane ID # "
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
        std::cout << "CFModuleDriver:main: Connectivity Table # "
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

    std::cout << "CFModuleDriver:main: getDataItemLoc "
              << getDataItemLoc << std::endl;

    std::cout << "CFModuleDriver:main: getDataItemType "
              << getDataItemType << std::endl;

    std::cout << "CFModuleDriver:main: numElementNodes "
              << numElementNodes << std::endl;

    std::cout << "CFModuleDriver:main: getDataItemUnits "
              << getDataItemUnits << std::endl;

    
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

    return 0;
}


