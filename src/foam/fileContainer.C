#include "fileContainer.H"

fileContainer::fileContainer()
{
    initSet();
}

fileContainer::~fileContainer(){};

int fileContainer::initSet()
{
    name = "";
    path = "";
    content = "";
    size = 0;
    
    return 0;
}
