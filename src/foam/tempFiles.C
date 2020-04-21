int comFoam::readInitFiles(const std::string& rootAddr)
{
    std::string fullAddr=rootAddr;
    std::vector<fileContainer> vecFile;
    int fileCount=0;
    
    readRecursive(rootAddr, fullAddr, vecFile, fileCount);

    ca_nFiles      = new int(fileCount);
    ca_fileSize    = new int[*ca_nFiles];
    ca_fileName    = new char*[*ca_nFiles];
    ca_filePath    = new char*[*ca_nFiles];
    ca_fileContent = new char*[*ca_nFiles];

    for(int ifile=0; ifile<fileCount; ifile++)
    {
        ca_fileSize[ifile] = vecFile[ifile].size; //This does not
                                                  //include the Null char

        std::string tmpStr = vecFile[ifile].name;
        ca_fileName[ifile] = new char [tmpStr.length()+1];
        std::strcpy(ca_fileName[ifile], tmpStr.c_str());

        tmpStr = vecFile[ifile].path;
        ca_filePath[ifile] = new char [tmpStr.length()+1];
        std::strcpy(ca_filePath[ifile], tmpStr.c_str());

        tmpStr = vecFile[ifile].content;
        ca_fileContent[ifile] = new char [tmpStr.length()+1];
        std::strcpy(ca_fileContent[ifile], tmpStr.c_str());
    }

    return 0;
}

int comFoam::createInitFiles(const std::string& rootAddr)
{
    //std::string rootAddr="./fluidTmp";
    namespace BF = boost::filesystem;

    std::vector<fileContainer> vecFile;
   
    for(int ifile=0; ifile<*ca_nFiles; ifile++)
    {
        fileContainer tmpFile;
        tmpFile.size = ca_fileSize[ifile];
        tmpFile.name = std::string(ca_fileName[ifile]);
        tmpFile.path = std::string(ca_filePath[ifile]);
        tmpFile.content = std::string(ca_fileContent[ifile]);

        vecFile.push_back(tmpFile);
    }

    for(int ifile=0; ifile<*ca_nFiles; ifile++)
    {
        std::string localDir = vecFile[ifile].path;
        std::string fullDir = rootAddr+"/"+localDir;
        std::string fullAddr = fullDir+"/"+vecFile[ifile].name;
        
        BF::path fullPath(fullDir);
        BF::create_directories(fullPath);

        std::ofstream outpuFile;
        outpuFile.open(fullAddr);
        
        outpuFile << vecFile[ifile].content;
        outpuFile.close();
    }

    return 0;
}

int comFoam::deleteInitFiles(const std::string& addr)
{
    namespace BF = boost::filesystem;

    BF::path rootPath(addr);
    BF::remove_all(rootPath);

    return 0;
}


int comFoam::registerInitFiles(const char *name)
{
    std::cout << "rocFoam.registerInitFiles: "
               << "Registering flow data with name "
               << name
               << std::endl;

    std::cout << "Files^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
              << std::endl;

    std::string volName = name+std::string("VOL");

    // Genral file data ^^^^^^^^^^^^^^^^^^^^^^^^
    std::string dataName = volName+std::string(".nFiles");
    COM_new_dataitem( dataName, 'w', COM_INT, 1, "");
    COM_set_size(     dataName, 0, 1);
    COM_set_array(    dataName, 0, ca_nFiles);
    std::cout << dataName << " registered." << std::endl;
    
    dataName = volName+std::string(".fileSize");
    COM_new_dataitem( dataName, 'w', COM_INT, 1, "");
    COM_set_size(     dataName, 0, *ca_nFiles);
    COM_set_array(    dataName, 0, ca_fileSize);
    //std::cout << dataName << " registered." << std::endl;

    // Register file paths
    for(int ifile=0; ifile<*ca_nFiles; ifile++)
    {
        std::string charToStr = std::string(ca_filePath[ifile]);
        int charSize = charToStr.length()+1;

        dataName = volName+".filePath"+std::to_string(ifile);
        COM_new_dataitem( dataName, 'w', COM_CHAR, 1, "");
        COM_set_size( dataName, 0, charSize);
        COM_set_array(dataName, 0, ca_filePath[ifile]);
        //std::cout << dataName << " registered." << std::endl;
    }

    // Register file names
    for(int ifile=0; ifile<*ca_nFiles; ifile++)
    {
        std::string charToStr = std::string(ca_fileName[ifile]);
        int charSize = charToStr.length()+1;

        dataName = volName+".fileName"+std::to_string(ifile);
        COM_new_dataitem( dataName, 'w', COM_CHAR, 1, "");
        COM_set_size( dataName, 0, charSize);
        COM_set_array(dataName, 0, ca_fileName[ifile]);
        //std::cout << "  File "
        //          << ca_filePath[ifile]<<"/"<<ca_fileName[ifile]
        //          << " registered."
        //          << std::endl;
    }

    // Register file contents
    for(int ifile=0; ifile<*ca_nFiles; ifile++)
    {
        std::string charToStr = std::string(ca_fileContent[ifile]);
        int charSize = charToStr.length()+1;

        dataName = volName+".fileContent"+std::to_string(ifile);
        COM_new_dataitem( dataName, 'w', COM_CHAR, 1, "");
        COM_set_size( dataName, 0, charSize);
        COM_set_array(dataName, 0, ca_fileContent[ifile]);
        //std::cout << dataName << " registered." << std::endl;
    }

    std::cout << "----------------------------------------------------"
              << std::endl;

    return 0;
}

int comFoam::reconstCaInitFiles(const char *name, const std::string& rootAddr)
{
    /*
    std::cout << "rocFoam.reconstructInitFiles: "
               << "Registering flow data with name "
               << name
               << std::endl;
   */
    std::string volName = name+std::string("VOL");

    std::string regNames;
    int numDataItems=0;
    COM_get_dataitems(volName.c_str(), &numDataItems, regNames);
    //std::cout << "  numDataItems = " << numDataItems << std::endl;

    std::vector<std::string> dataItemNames;
    dataItemNames.clear();
    std::istringstream Istr(regNames);
    for (int i=0; i<numDataItems; ++i)
    {
        std::string nameTmp;
        Istr >> nameTmp;
        dataItemNames.push_back(nameTmp);
        
        std::string subName = nameTmp.substr(0,4);
        if (subName == "file" || subName == "nFil")
        {
            dataItemNames.push_back(nameTmp);
            //Info << "  DataItem[" << i << "] = " << nameTmp << endl;
        }
    }

    // Genral file data ^^^^^^^^^^^^^^^^^^^^^^^^
    std::string dataName = std::string("nFiles");
    nameExists(dataItemNames, dataName);
    std::string regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), 0, &ca_nFiles);
    std::cout << "  " << dataName.c_str() << " = " << *ca_nFiles << std::endl;

    dataName = std::string("fileSize");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    int nComp;
    COM_get_array(regName.c_str(), 0, &ca_fileSize);
    COM_get_size(regName.c_str(), 0, &nComp);

    std::vector<fileContainer> vecFile;

    ca_fileName    = new char*[*ca_nFiles];
    ca_filePath    = new char*[*ca_nFiles];
    ca_fileContent = new char*[*ca_nFiles];
    //-------------------------------------------

    // Read file-specific data ^^^^^^^^^^^^^^^^^^
    for(int ifile=0; ifile<*ca_nFiles; ifile++)
    {
        dataName = std::string("fileName")+std::to_string(ifile);
        nameExists(dataItemNames, dataName);
        regName = volName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), 0, &ca_fileName[ifile]);
        COM_get_size(regName.c_str(), 0, &nComp);

        dataName = std::string("filePath")+std::to_string(ifile);
        nameExists(dataItemNames, dataName);
        regName = volName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), 0, &ca_filePath[ifile]);
        COM_get_size(regName.c_str(), 0, &nComp);

        dataName = std::string("fileContent")+std::to_string(ifile);
        nameExists(dataItemNames, dataName);
        regName = volName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), 0, &ca_fileContent[ifile]);
        COM_get_size(regName.c_str(), 0, &nComp);

        if (nComp != ca_fileSize[ifile]+1)
        {
            std::cout << " Warning: Registered and retrieved sizes of file"
                      << ca_fileName[ifile] << " do not match:"
                      << " fileSize = " << ca_fileSize[ifile]+1
                      << " nComp = " << nComp << std::endl;
        }
    }

    // Register file status ^^^^^^^^^^^^^^^^^^^^^
    for(int ifile=0; ifile<*ca_nFiles; ifile++)
    {
        //std::cout << "    file[" << ifile << "] = "
        //                  << vecFile[ifile].path.c_str() << "/"
        //                  << vecFile[ifile].name.c_str() << std::endl;

        /*
        std::cout << "    fileName[" << ifile << "] = "
                  << vecFile[ifile].name.c_str() << std::endl;

        std::cout << "    filePath[" << ifile << "] = "
                  << vecFile[ifile].path.c_str() << std::endl;

        std::cout << "    fileSize[" << ifile << "] = "
                  << vecFile[ifile].size << std::endl;

        std::cout << "    fileContent[" << ifile << "] = \n"
                  << vecFile[ifile].content << std::endl;
        */
    }
    //-------------------------------------------


    createInitFiles(rootAddr);

    return 0;
}


int comFoam::readRecursive
(
    const std::string& rootAddr,
    std::string fullAddr,
    std::vector<fileContainer>& vecFile,
    int& fileCount
)
{
    namespace BF = boost::filesystem;

    // Separate the local work directory from the global
    BF::path rootPath;
    rootPath = BF::initial_path();

    int firstElemet = 0;
    for (auto it : rootPath)
    {
        //std::cout << firstElemet << " " << it <<std::endl;
        firstElemet++;
    }

    // Walk through the full address
    BF::path fullPath(fullAddr);
    try
    {
        if (BF::exists(fullPath))
        {
            if (BF::is_regular_file(fullPath))
            {
                // If a regular file, separate the local address
                BF::path localPath("");
                int count=0;
                for (auto it : fullPath)
                {
                    if (count >= firstElemet)
                    {
                        localPath = localPath /= it;
                    }
                    count++;
                }

                BF::path tmpPath;
                tmpPath = localPath.parent_path();
                std::string localAddr = tmpPath.string();

                tmpPath = localPath.filename();
                std::string fileName = tmpPath.string();
                
                std::ifstream inputFile(localPath.string());
                if (!inputFile.good())
                {
                    std::cout << "Warnning: possibly " 
                              << localPath.string() << " does not exist."
                              << std::endl;
                }
                else
                {
                    inputFile.seekg (0, inputFile.end);
                    int size = inputFile.tellg();
                    inputFile.seekg (0, inputFile.beg);

                    char* content = new char[size];
                    inputFile.read(content, size);
                    
                    fileContainer tmpFile;
                    tmpFile.name = fileName;
                    tmpFile.path = localAddr;
                    tmpFile.size = size;
                    tmpFile.content = string(content);
                    tmpFile.content.resize(size);
                    
                    //std::cout << fileCount << " Reading from file "
                    //          << localAddr+"/"+fileName
                    //          << std::endl;


                    vecFile.push_back(tmpFile);

                    //std::cout << vecFile[fileCount].name << std::endl;
                    //std::cout << vecFile[fileCount].path << std::endl;
                    //std::cout << vecFile[fileCount].size << std::endl;
                    //std::cout << vecFile[fileCount].content << std::endl;
                }
                fileCount++;
            }
            else if (BF::is_directory(fullPath))
            {
                //std::string fileName = BF::canonical(fullPath).string();
                //std::cout << fileName << " is a directory containing:\n";
                for (auto x : BF::directory_iterator(fullPath))
                {
                    BF::path newPath=x.path();
                    std::string neWaddr = BF::canonical(newPath).string();
                     
                    //std::cout << "    " << fileName << '\n';
                    readRecursive(rootAddr, neWaddr, vecFile, fileCount);
                }
            }
            else
            {
                std::string fileName = BF::canonical(fullPath).string();
                std::cout << fileName
                          << " exists, but is not a regular file or directory"
                          << std::endl;
            }
        }
        else
        {
            std::cout << fullPath << " does not exist" << std::endl;
        }
    }

    catch (const BF::filesystem_error& ex)
    {
        std::cout << ex.what() << std::endl;
    }
    
    return 0;
}


int comFoam::deleteFilesData()
{
    if (ca_fileContent != NULL)
    {
        for (int ifile=0; ifile<*ca_nFiles; ifile++)
        {
            if (ca_fileContent[ifile] != NULL)
            {
                delete [] ca_fileContent[ifile];
            }
        }
        delete [] ca_fileContent;
        ca_fileContent = NULL;
    }

    if (ca_filePath != NULL)
    {
        for (int ifile=0; ifile<*ca_nFiles; ifile++)
        {
            if (ca_filePath[ifile] != NULL)
            {
                delete [] ca_filePath[ifile];
            }
        }
        delete [] ca_filePath;
        ca_filePath = NULL;
    }

    if (ca_fileName != NULL)
    {
        for (int ifile=0; ifile<*ca_nFiles; ifile++)
        {
            if (ca_fileName[ifile] != NULL)
            {
                delete [] ca_fileName[ifile];
            }
        }
        delete [] ca_fileName;
        ca_fileName = NULL;
    }

    if (ca_fileSize != NULL)
    {
        delete [] ca_fileSize;
        ca_fileSize = NULL;
    }

    if (ca_nFiles != NULL)
    {
        delete [] ca_nFiles;
        ca_nFiles = NULL;
    }

    return 0;
}






