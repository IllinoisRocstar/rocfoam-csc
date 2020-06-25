#include "rocRhoPimple.H"
#include "gtest.h"

// Global variables used to pass arguments to the tests
char** ARGV;
int ARGC;

// Test fixture which tests COM. Derived from the Google test primer
class rocRhoPimpleTest : public ::testing::Test
{
protected:
    void SetUp() {};

    void TearDown() {};
    
    rhoPimple *rocFoam;
};

TEST_F(rocRhoPimpleTest, rhoPimple)
{
    rocFoam = new rhoPimple;
    rocFoam->initFOAM(ARGC, ARGV);
    EXPECT_EQ(rocFoam->initializeStat, 0) << "Testing rocRhoPimple: initilize unsuccessful"
                                          << std::endl;

    rocFoam->loop();
    EXPECT_EQ(rocFoam->loopStat, 0) << "Testing rocRhoPimple: loop unsuccessful"
                                    << std::endl;

    rocFoam->finalize();
    EXPECT_EQ(rocFoam->finalizeStat, 0) << "Testing rocRhoPimple: finalize unsuccessful"
                                        << std::endl;

    delete rocFoam; rocFoam=NULL;
    EXPECT_EQ(rocFoam, nullptr)  << "Testing rocRhoPimple: deleting unsuccessful"
                                 << std::endl;
}

int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    ARGC = argc;
    ARGV = argv;
    return RUN_ALL_TESTS();
}
