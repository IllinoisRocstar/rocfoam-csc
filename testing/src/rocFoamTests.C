#include "gtest.h"

// Global variables used to pass arguments to the tests
char **ARGV;
int ARGC;

// Test fixture which tests COM. Derived from the Google test primer
class rocFoamTest : public ::testing::Test
{
protected:
    void SetUp() {};

    void TearDown() {};
        
    MPI_Comm masterComm;
    MPI_Comm newComm;

    int masterRank;
    int masterNProc;
    bool runParallel;

    char *solverType;

    //  Function Handlers ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    int flowInitHandle;
    int flowStatHandle;
    int flowLoopHandle;
    int flowFinHandle;


    //  Function stats ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    int comDrvInitStat = -1;
    int comDrvRunStat = -1;
    int comDrvFinStat = -1;

    //  Function definitions ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
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
        

        //std::string arg;
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
                    throw -1;
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

        if (!runParallel && masterNProc > 1)
        {
            if (masterRank==0)
            {
                std::cout << "rocFoam.main: NProc>1 detected for a serial job."
                          << std::endl;
                throw -1;
            }    
        
        }


        //  Setting the defual communicator. Is it needed?
        COM_set_default_communicator(newComm);
        
        comfoam_load_module("ROCFOAM", solverType);

        // getting number of processes  
        if (masterRank==0)
        {
            int *nProcReg;
            COM_get_array("ROCFOAM.winNProc", 0, &nProcReg);
            std::cout << "The communicator registered in OFModule uses "
                      << *nProcReg << " prcesses"
                      << std::endl;

            std::cout << std::endl;
        }

        //  Get the handle for the initialize function ^^^^^^^^
        flowInitHandle = COM_get_function_handle("ROCFOAM.flowInit");
        if (flowInitHandle <= 0)
        { // fail
            std::cout << "rocFoam.main: Could not get handle for initialize."
                      << std::endl;
            throw -2;
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
        flowLoopHandle = COM_get_function_handle("ROCFOAM.flowLoop");
        if (flowLoopHandle <= 0)
        { // fail
            std::cout << "rocFoam.main: Could not get handle for loop."
                      << std::endl;
            throw -2;
        }
        else
        {
            if (masterRank==0)
            {

                std::cout << "rocFoam.main: Acquired a handle for loop."
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
            throw -2;
        }
        else
        {
            std::cout << "rocFoam.main: Acquired a handle for finLauncher."
                      << std::endl;    
        }*/

        //  Make a dummy argc/argv for OpenFOAM. ^^^^
        //  No options passed from the command
        //  line will be used by the driver

        int verb=3;
        int myArgc = 0;
        char *myArgv[argc];
        
        for (int i=0; i<argc; i++)
        {
            myArgv[i] = NULL;

            ss.clear();
            ss.str("");
            ss << argv[i];

            
            if (ss.str() != "-"+string(solverType))
            {
                myArgv[myArgc] = argv[i];
                myArgc ++;
            }
        }

        //  Fluid initializer ^^^^^^^^^^^^^^^^^^^^^^^
        COM_call_function(flowInitHandle, &myArgc, &myArgv, &verb);
        
        return 0;
    }


    int comDrvRun()
    {
        //  Call the flow iterator ^^^^^^^^^^^^^^^^^^
        COM_call_function(flowLoopHandle);
        
        return 0;
    }


    int comDrvFin()
    {
        //  Call the flow unloader ^^^^^^^^^^^^^^^^^^
        //COM_UNLOAD_MODULE_STATIC_DYNAMIC(comfoam, "ROCFOAM");
        comfoam_unload_module("ROCFOAM", solverType);

        COM_set_default_communicator(masterComm);
        
        COM_finalize();

        MPI_Barrier(masterComm);
        MPI_Finalize();
        
        return 0;
    }  

};

TEST_F(rocFoamTest, rhoFoam)
{
    comDrvInitStat = comDrvInit(ARGC, ARGV);
    EXPECT_EQ(comDrvInitStat, 0) << "Testing rocFoam: initilize unsuccessful"
                                 << std::endl;

    comDrvRunStat = comDrvRun();
    EXPECT_EQ(comDrvRunStat, 0) << "Testing rocFoam: loop unsuccessful"
                                << std::endl;

    comDrvFinStat = comDrvFin();
    EXPECT_EQ(comDrvFinStat, 0) << "Testing rocFoam: finalize unsuccessful"
                                << std::endl;
}

int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    ARGC = argc;
    ARGV = argv;
    return RUN_ALL_TESTS();
}
