#ifndef FILECONTAINER_
#define FILECONTAINER_

#include <iostream>

//using namespace std;

class fileContainer
{
public:
    std::string name;
    std::string path;
    std::string content;
    int size;

public:
    fileContainer();
    ~fileContainer();





private:
    int initSet();
};

#endif
