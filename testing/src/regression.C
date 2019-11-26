#include "rocRhoCentral.H"

int main(int argc, char *argv[])
{

    char *solverType;

    std::stringstream ss;
    std::stringstream caseI;
    std::stringstream caseII;


    std::string word;


std::cout << __LINE__ << std::endl;

    if (argc > 1)
    {
        for (int i=1; i<argc; ++i)
        {
            ss.clear();
            ss.str("");
            ss << argv[i];
            
            //ss >> word;
            //word = ss.str();

//std::cout << __LINE__ << " " << word << " " << word.length()+1 << std::endl;

            
            if (ss.str() == "-rocRhoCentral")
            {
                solverType = const_cast<char *>("rocRhoCentral");
                
//                char *tmpChar = new char[word.length()+1];
//std::cout << __LINE__ << " " << tmpChar << std::endl;
//                strcpy(tmpChar, word.c_str());
//std::cout << __LINE__ << " " << tmpChar << " " << *tmpChar << std::endl;

                //solverType = *tmpChar;
                
//                char *solverType1 = new char(*tmpChar);
//std::cout << __LINE__ << " " << *solverType1 << std::endl;
//                delete [] tmpChar;
            }
            else if (ss.str() == "-rocRhoPimple")
            {
                solverType = const_cast<char *>("rocRhoPimple");
            }
            else if (ss.str() == "-caseI")
            {
                caseI << argv[i+1] ;
            }
            else if (ss.str() == "-caseII")
            {
                caseII << argv[i+1] ;
            }
        }
    }


    int argc1 = 3;
    char *argv1[argc1];
    argv1[0] = solverType; //argv[0];
    argv1[1] = const_cast<char *>("-case");


    word = caseI.str();
    int length = word.length();
    char *tmpCaseI = new char[length+1];

    strcpy(tmpCaseI, word.c_str());
    argv1[2] = tmpCaseI;


    int argc2 = 3;
    char *argv2[argc2];
    argv2[0] = solverType; //argv[0];
    argv2[1] = const_cast<char *>("-case");


    word = caseII.str();
    length = word.length();
    char *tmpCaseII = new char[length+1];

    strcpy(tmpCaseII, word.c_str());
    argv2[2] = tmpCaseII;



    rhoCentral rocFoam1(argc1, argv1);
    
    rhoCentral rocFoam2(argc2, argv2);



    delete [] tmpCaseI;
    delete [] tmpCaseII;




return 0;

    return 0;
}
