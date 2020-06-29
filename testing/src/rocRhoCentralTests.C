#include "rocRhoCentral.H"
#include "gtest.h"

// Global variables used to pass arguments to the tests
char** ARGV;
int ARGC;

// Test fixture which tests COM. Derived from the Google test primer
class rocRhoCentralTest : public ::testing::Test
{
protected:
    void SetUp() {};

    void TearDown() {};
    
    rhoCentral *rocFoam;
};

TEST_F(rocRhoCentralTest, rhoCentral)
{
    rocFoam = new rhoCentral;
    rocFoam->initFOAM(ARGC, ARGV);
    EXPECT_EQ(rocFoam->initializeStat, 0) << "Testing rocRhoCentral: initilize unsuccessful"
                                          << std::endl;

    rocFoam->loop();
    EXPECT_EQ(rocFoam->loopStat, 0) << "Testing rocRhoCentral: loop unsuccessful"
                                    << std::endl;

    rocFoam->finalize();
    EXPECT_EQ(rocFoam->finalizeStat, 0) << "Testing rocRhoCentral: finalize unsuccessful"
                                        << std::endl;

    delete rocFoam; rocFoam=NULL;
    EXPECT_EQ(rocFoam, nullptr)  << "Testing rocRhoCentral: deleting unsuccessful"
                                 << std::endl;
}

int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    ARGC = argc;
    ARGV = argv;
    return RUN_ALL_TESTS();
}
