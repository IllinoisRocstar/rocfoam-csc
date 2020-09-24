int comFoam::readFilesData(const std::string& rootAddr)
{
    std::string fullAddr=rootAddr;
    std::string locaParAddr = "";
    
    int myRank = 0;
    if (Pstream::parRun())
    {
        myRank = ca_myRank;
        locaParAddr = "processor"+std::to_string(myRank);
        //fullAddr = rootAddr+locaParAddr;
    }

    std::vector<fileContainer> vecFile;
    int fileCount=0;
    
    readRecursive(locaParAddr, fullAddr, vecFile, fileCount);

    ca_nFiles      = new int(fileCount);
    ca_fileSize    = new int[*ca_nFiles];
    ca_fileName    = new char*[*ca_nFiles];
    ca_filePath    = new char*[*ca_nFiles];
    ca_fileContent = new char*[*ca_nFiles];

    for(int ifile=0; ifile<fileCount; ifile++)
    {
        ca_fileSize[ifile] = vecFile[ifile].size; //This does not
                                                  //include the nullptr char

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

int comFoam::createFilesData()
{
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

    createSysConstFiles(tmpFluidDir, vecFile);
    createFieldFiles(tmpFluidDir, vecFile);
    createUniformTimeFile(tmpFluidDir);
    createPointsFile(tmpFluidDir);
    createOwnerFile(tmpFluidDir);
    createNeighborFile(tmpFluidDir);
    createFacesFile(tmpFluidDir);
    createBoundaryFile(tmpFluidDir, vecFile);
    createConnectivityFiles(tmpFluidDir, vecFile);

    return 0;
}


int comFoam::createSysConstFiles
            (
                const std::string& rootAddr,
                const std::vector<fileContainer>& vecFile
            )
{
    namespace BF = boost::filesystem;
    for(int ifile=0; ifile<*ca_nFiles; ifile++)
    {
        std::string fileName = vecFile[ifile].name;
        std::string localDir = vecFile[ifile].path;
        std::string fullDir = rootAddr+"/"+localDir;
        std::string fullAddr = fullDir+"/"+fileName;

        if (localDir=="system" ||
            localDir=="constant" )
        {
            BF::path fullPath(fullDir);
            BF::create_directories(fullPath);

            std::ofstream outpuFile;
            outpuFile.open(fullAddr);
            
            outpuFile << vecFile[ifile].content;
            outpuFile.close();

            std::cout << "    " << fullAddr
                      << " created." << std::endl;
        }
    }
    
    

    return 0;
}

int comFoam::createFieldFiles
            (
                const std::string& rootAddr,
                const std::vector<fileContainer>& vecFile
            )
{
    namespace BF = boost::filesystem;
    
    int nTotal=*ca_nFiles;
    //if (ca_Rho != nullptr)
         nTotal++;
    //if (ca_Phi != nullptr)
        nTotal++;
    //if (ca_RhoUf != nullptr)
        nTotal++;
    //if (ca_Disp != nullptr)
        nTotal++;
    
    for (int ifile=0; ifile<nTotal; ifile++)
    {
        // rho: ifile == *ca_nFiles
        // phi: ifile == *ca_nFiles+1
    
        std::string locaParAddr = "";
        if (ca_nProc>1)
        {
            locaParAddr = "processor"+std::to_string(ca_myRank)+"/";
        }

        std::string fileName = "";
        std::string localDir = "";
        std::string content = "";
        if (ifile<*ca_nFiles)
        {
            fileName = vecFile[ifile].name;
            localDir = vecFile[ifile].path;
            content = vecFile[ifile].content;
        }
        else if (ifile==*ca_nFiles && ca_Rho != nullptr)
        {
            fileName = "rho";
            localDir = locaParAddr+"0";
            content = createBaseFile
                      (
                        "rho",
                        "vol",
                        "Scalar",
                        "[1 -3 0 0 0 0 0]"
                      );
        }
        else if (ifile==*ca_nFiles+1 && ca_Phi != nullptr)
        {
            fileName = "phi";
            localDir = locaParAddr+"0";
            content = createBaseFile
                      (
                        "phi",
                        "surface",
                        "Scalar",
                        "[1 0 -1 0 0 0 0]"
                      );
        }
        else if (ifile==*ca_nFiles+2 && ca_RhoUf != nullptr)
        {
            fileName = "rhoUf";
            localDir = locaParAddr+"0";
            content = createBaseFile
                      (
                        "rhoUf",
                        "surface",
                        "Vector",
                        "[1 -2 -1 0 0 0 0]"
                      );
        }
        else if (ifile==*ca_nFiles+3 && ca_Disp != nullptr)
        {
            fileName = "pointDisplacementNew";
            localDir = locaParAddr+"0";
            content = createBaseFile
                      (
                        "pointDisplacement",
                        "point",
                        "Vector",
                        "[0 1 0 0 0 0 0]"
                      );
        }

        std::string fullDir = rootAddr+"/"+localDir;
        std::string fullAddr = fullDir+"/"+fileName;

        //if (localDir=="0")
        if (localDir==locaParAddr+"0")
        {
            std::string locationStr = "location";
            size_t locationStart = findWord(fullAddr, content, locationStr);
            if (locationStart !=0 && locationStart != std::string::npos)
            {
                std::string scStr = ";";
                size_t locationEnd = findChar(fullAddr, content, scStr, locationStart);
                size_t length = locationEnd-locationStart;
                content.erase(locationStart, length);
                std::string newStr  = "location    \"";
                newStr += std::string(ca_timeName);
                newStr += "\"";
                content.insert(locationStart, newStr);
            }

            if (fileName == "U" ||
                fileName == "p" ||
                fileName == "T" ||
                fileName == "rho" ||
                fileName == "phi" ||
                fileName == "rhoUf" ||
                fileName == "alphat" ||
                fileName == "epsilon" ||
                fileName == "omega" ||
                fileName == "k" ||
                fileName == "nut" ||
                fileName == "pointDisplacement" ||
                fileName == "pointDisplacementNew"
               )
            {   

                std::string intFieldStr = "internalField";
                std::string scStr = ";";
                size_t intFieldStart= findWordOnly(fullAddr, content, intFieldStr);
                size_t intFieldEnd  = findChar(fullAddr, content, scStr, intFieldStart);
                size_t length = intFieldEnd-intFieldStart;
                content.erase(intFieldStart, length);

                std::string newStr = "internalField   nonuniform List";
                if (fileName == "U" && ca_Vel!= nullptr)
                {
                    newStr += "<vector>\n";
                    int ncells = *ca_nCells;
                    newStr += std::to_string(ncells);
                    newStr += "\n(\n";

                    for(int icell=0; icell<ncells; icell++)
                    {
                        int cellIndex = ca_cellToCellMap_inverse[icell];
                        newStr += "(";
                        for(int jcomp=0; jcomp<nComponents; jcomp++)
                        {
                            int localComp = jcomp + cellIndex*nComponents;
                            
                            std::ostringstream doubleToOs;
                            doubleToOs << std::scientific << std::setprecision(IODigits);
                            doubleToOs << ca_Vel[localComp];

                            newStr += removeTrailZero(doubleToOs.str()); //doubleToOs.str();
                            if (jcomp<nComponents-1)
                                newStr += " ";
                        }
                        newStr += ")\n";
                    }
                    newStr += ")\n";
                    
                    content.insert(intFieldStart, newStr);
                }
                else if (fileName == "p" && ca_P!= nullptr)
                {
                    newStr += "<scalar>\n";
                    int ncells = *ca_nCells;
                    newStr += std::to_string(ncells);
                    newStr += "\n(\n";
                    
                    for(int icell=0; icell<ncells; icell++)
                    {
                        int cellIndex = ca_cellToCellMap_inverse[icell];

                        std::ostringstream doubleToOs;                            
                        doubleToOs << std::scientific << std::setprecision(IODigits);
                        doubleToOs << ca_P[cellIndex];

                        newStr += removeTrailZero(doubleToOs.str()); //doubleToOs.str();
                        newStr += "\n";
                    }
                    newStr += ")\n";
                    
                    content.insert(intFieldStart, newStr);
                }
                else if (fileName == "T" && ca_T!= nullptr)
                {
                    newStr += "<scalar>\n";
                    int ncells = *ca_nCells;
                    newStr += std::to_string(ncells);
                    newStr += "\n(\n";
                    
                    for(int icell=0; icell<ncells; icell++)
                    {
                        int cellIndex = ca_cellToCellMap_inverse[icell];

                        std::ostringstream doubleToOs;                            
                        doubleToOs << std::scientific << std::setprecision(IODigits);
                        doubleToOs << ca_T[cellIndex];

                        newStr += removeTrailZero(doubleToOs.str()); //doubleToOs.str();
                        newStr += "\n";
                    }
                    newStr += ")\n";
                    
                    content.insert(intFieldStart, newStr);
                }
                else if (fileName == "rho" && ca_Rho!= nullptr)
                {
                    newStr += "<scalar>\n";
                    int ncells = *ca_nCells;
                    newStr += std::to_string(ncells);
                    newStr += "\n(\n";
                    
                    for(int icell=0; icell<ncells; icell++)
                    {
                        int cellIndex = ca_cellToCellMap_inverse[icell];

                        std::ostringstream doubleToOs;                            
                        doubleToOs << std::scientific << std::setprecision(IODigits);
                        doubleToOs << ca_Rho[cellIndex];

                        newStr += removeTrailZero(doubleToOs.str()); //doubleToOs.str();
                        newStr += "\n";
                    }
                    newStr += ")\n";
                    
                    content.insert(intFieldStart, newStr);
                }
                else if (fileName == "phi" && ca_Phi != nullptr)
                {
                    newStr += "<scalar>\n";

                    int procStartIndex{0};
                    for (int iproc=0; iproc<ca_myRank; iproc++)
                    {
                        procStartIndex += ca_nPatches[iproc];
                    }
                    int index  = procStartIndex + 0; //ipatch;
                    int nfaces = ca_patchStart[index];

                    newStr += std::to_string(nfaces);
                    newStr += "\n(\n";

                    for(int iface=0; iface<nfaces; iface++)
                    {
                        int faceIndex = ca_faceToFaceMap_inverse[iface];

                        std::ostringstream doubleToOs;                            
                        doubleToOs << std::scientific << std::setprecision(IODigits);
                        doubleToOs << ca_Phi[faceIndex];

                        newStr += removeTrailZero(doubleToOs.str()); //doubleToOs.str();
                        newStr += "\n";
                    }
                    newStr += ")\n";
                    
                    content.insert(intFieldStart, newStr);
                }
                else if (fileName == "rhoUf" && ca_RhoUf != nullptr)
                {
                    newStr += "<vector>\n";

                    int procStartIndex{0};
                    for (int iproc=0; iproc<ca_myRank; iproc++)
                    {
                        procStartIndex += ca_nPatches[iproc];
                    }
                    int index  = procStartIndex + 0; //ipatch;
                    int nfaces = ca_patchStart[index];

                    newStr += std::to_string(nfaces);
                    newStr += "\n(\n";

                    for(int iface=0; iface<nfaces; iface++)
                    {
                        int faceIndex = ca_faceToFaceMap_inverse[iface];
                        newStr += "(";

                        for(int jcomp=0; jcomp<nComponents; jcomp++)
                        {
                            int localComp = jcomp + faceIndex*nComponents;
                            
                            std::ostringstream doubleToOs;
                            doubleToOs << std::scientific << std::setprecision(IODigits);
                            doubleToOs << ca_RhoUf[localComp];

                            newStr += removeTrailZero(doubleToOs.str()); //doubleToOs.str();
                            if (jcomp<nComponents-1)
                                newStr += " ";
                        }
                        newStr += ")\n";
                    }
                    newStr += ")\n";

                    content.insert(intFieldStart, newStr);
                }
                // Turbulence data ^^^^^^^^^^^^^^^^^^
                else if (fileName == "alphat" && ca_AlphaT!= nullptr)
                {
                    newStr += "<scalar>\n";
                    int ncells = *ca_nCells;
                    newStr += std::to_string(ncells);
                    newStr += "\n(\n";
                    
                    for(int icell=0; icell<ncells; icell++)
                    {
                        int cellIndex = ca_cellToCellMap_inverse[icell];

                        std::ostringstream doubleToOs;                            
                        doubleToOs << std::scientific << std::setprecision(IODigits);
                        doubleToOs << ca_AlphaT[cellIndex];

                        newStr += removeTrailZero(doubleToOs.str()); //doubleToOs.str();
                        newStr += "\n";
                    }
                    newStr += ")\n";
                    
                    content.insert(intFieldStart, newStr);
                }
                else if (fileName == "epsilon" && ca_Epsilon!= nullptr)
                {
                    newStr += "<scalar>\n";
                    int ncells = *ca_nCells;
                    newStr += std::to_string(ncells);
                    newStr += "\n(\n";
                    
                    for(int icell=0; icell<ncells; icell++)
                    {
                        int cellIndex = ca_cellToCellMap_inverse[icell];

                        std::ostringstream doubleToOs;                            
                        doubleToOs << std::scientific << std::setprecision(IODigits);
                        doubleToOs << ca_Epsilon[cellIndex];

                        newStr += removeTrailZero(doubleToOs.str()); //doubleToOs.str();
                        newStr += "\n";
                    }
                    newStr += ")\n";
                    
                    content.insert(intFieldStart, newStr);
                }
                else if (fileName == "omega" && ca_Omega!= nullptr)
                {
                    newStr += "<scalar>\n";
                    int ncells = *ca_nCells;
                    newStr += std::to_string(ncells);
                    newStr += "\n(\n";
                    
                    for(int icell=0; icell<ncells; icell++)
                    {
                        int cellIndex = ca_cellToCellMap_inverse[icell];

                        std::ostringstream doubleToOs;                            
                        doubleToOs << std::scientific << std::setprecision(IODigits);
                        doubleToOs << ca_Omega[cellIndex];

                        newStr += removeTrailZero(doubleToOs.str()); //doubleToOs.str();
                        newStr += "\n";
                    }
                    newStr += ")\n";
                    
                    content.insert(intFieldStart, newStr);
                }
                else if (fileName == "k" && ca_K!= nullptr)
                {
                    newStr += "<scalar>\n";
                    int ncells = *ca_nCells;
                    newStr += std::to_string(ncells);
                    newStr += "\n(\n";
                    
                    for(int icell=0; icell<ncells; icell++)
                    {
                        int cellIndex = ca_cellToCellMap_inverse[icell];

                        std::ostringstream doubleToOs;                            
                        doubleToOs << std::scientific << std::setprecision(IODigits);
                        doubleToOs << ca_K[cellIndex];

                        newStr += removeTrailZero(doubleToOs.str()); //doubleToOs.str();
                        newStr += "\n";
                    }
                    newStr += ")\n";
                    
                    content.insert(intFieldStart, newStr);
                }
                else if (fileName == "nut" && ca_NuT!= nullptr)
                {
                    newStr += "<scalar>\n";
                    int ncells = *ca_nCells;
                    newStr += std::to_string(ncells);
                    newStr += "\n(\n";
                    
                    for(int icell=0; icell<ncells; icell++)
                    {
                        int cellIndex = ca_cellToCellMap_inverse[icell];

                        std::ostringstream doubleToOs;                            
                        doubleToOs << std::scientific << std::setprecision(IODigits);
                        doubleToOs << ca_NuT[cellIndex];

                        newStr += removeTrailZero(doubleToOs.str()); //doubleToOs.str();
                        newStr += "\n";
                    }
                    newStr += ")\n";
                    
                    content.insert(intFieldStart, newStr);
                }
                else if ((fileName == "pointDisplacement" ||
                         fileName == "pointDisplacementNew")
                         && ca_Disp != nullptr)
                {
                    newStr += "<vector>\n";
                    int npoints = *ca_nPoints;
                    newStr += std::to_string(npoints);
                    newStr += "\n(\n";

                    for(int ipoint=0; ipoint<npoints; ipoint++)
                    {
                        newStr += "(";
                        for(int jcomp=0; jcomp<nComponents; jcomp++)
                        {
                            int localComp = jcomp + ipoint*nComponents;
                            
                            std::ostringstream doubleToOs;
                            doubleToOs << std::scientific << std::setprecision(IODigits);
                            doubleToOs << ca_Disp[localComp];

                            newStr += removeTrailZero(doubleToOs.str()); //doubleToOs.str();
                            if (jcomp<nComponents-1)
                                newStr += " ";
                        }
                        newStr += ")\n";
                    }
                    newStr += ")\n";
                    
                    content.insert(intFieldStart, newStr);
                }
                //-----------------------------------
            }
            // Files that are the in folder "0"" and
            // not being taken care of yet. Need to revisit this later
            else if (fileName == "cellToRegion")
            {}
            //-----------------------------------
            else
            {
                std::cout << "========== WARNING ===============" << std::endl
                          << "Solution field for file " 
                          << fileName << " is not found." << std::endl;
                continue;
            }
            
            for(int ipatch=0; ipatch<ca_nPatches[ca_myRank]; ipatch++)
            {


                if (fileName == "U" ||
                    fileName == "p" ||
                    fileName == "T" ||
                    fileName == "rho" ||
                    fileName == "phi" ||
                    fileName == "rhoUf" ||
                    fileName == "alphat" ||
                    fileName == "epsilon" ||
                    fileName == "omega" ||
                    fileName == "k" ||
                    fileName == "nut"||
                    fileName == "pointDisplacement" ||
                    fileName == "pointDisplacementNew"
                   )
                {
                    int procStartIndex{0};
                    for (int iproc=0; iproc<ca_myRank; iproc++)
                    {
                        procStartIndex += ca_nPatches[iproc];
                    }
                    int index  = procStartIndex + ipatch;
                    int nfaces = ca_patchSize[index];

                    if (nfaces<=0)
                    {
                        continue;
                    }
                    
                    std::string patchName = patchNameStr[ipatch];
                    size_t patchStart = findWordOnly(fullAddr, content, patchName);
                    
                    std::string brckOpen = "{";
                    std::string brckClose = "}";
                    size_t brckStart = findChar(fullAddr, content, brckOpen, patchStart);
                    size_t brckEnd   = findChar(fullAddr, content, brckClose, brckStart);
                    
                    std::string typeStr = "type";
                    std::string scStr = ";";
                    size_t typeStart = findWordOnly(fullAddr, content, typeStr, brckStart, brckEnd);
                    size_t typeEnd   = findChar(fullAddr, content, scStr, typeStart);

                    size_t length = typeEnd-typeStart;
                    std::string subType = content.substr(typeStart, length);
                    std::stringstream iStr(subType);

                    std::string typeName;
                    while(iStr >> typeName)
                    {
                    }
                    if (typeName=="fixedValue"       ||
                        typeName=="symmetryPlane"    ||
                        typeName=="slip"             ||
                        typeName=="noSlip"           ||
                        typeName=="zeroGradient"     ||
                        typeName=="totalTemperature" ||
                        typeName=="totalPressure"    ||
                        typeName=="empty"
                       ) continue;
                    
                    std::string valueStr = "value";
                    scStr = ";";
                    
                    size_t valueStart = 0;
                    if (fileName == "pointDisplacement" ||
                        fileName == "pointDisplacementNew")
                    {
                        valueStart = findWord(fullAddr,
                                             content,
                                             valueStr,
                                             brckStart,
                                             brckEnd);

                        if (valueStart<=brckStart || brckEnd<=valueStart)
                            continue;
                    }
                    else
                    {
                        valueStart = findWordOnly(fullAddr,
                                                 content,
                                                 valueStr,
                                                 brckStart,
                                                 brckEnd);
                    }
                    
                    size_t valueEnd   = findChar(fullAddr, content, scStr, valueStart);

                    length = valueEnd-valueStart;
                    content.erase(valueStart, length);

                    std::string newStr = "value           nonuniform List";

                    if (fileName == "U" && ca_patchVel!= nullptr)
                    {
                        newStr += "<vector>\n";
                        newStr += std::to_string(nfaces);
                        newStr += "\n(\n";

                        for(int iface=0; iface<nfaces; iface++)
                        {
                            int faceIndex = ca_patchFaceToFaceMap_inverse[ipatch][iface];

                            newStr += "(";
                            for(int jcomp=0; jcomp<nComponents; jcomp++)
                            {
                                int localComp = jcomp + faceIndex*nComponents;
                                
                                std::ostringstream doubleToOs;                            
                                doubleToOs << std::scientific << std::setprecision(IODigits);
                                doubleToOs << ca_patchVel[ipatch][localComp];
                                newStr += removeTrailZero(doubleToOs.str());
                                if (jcomp<nComponents-1)
                                    newStr += " ";
                            }
                            newStr += ")\n";
                        }
                        newStr += ")\n";
                    }
                    else if (fileName == "p" && ca_patchP!= nullptr)
                    {
                        newStr += "<scalar>\n";
                        //int nfaces = ca_patchSize[ipatch];
                        newStr += std::to_string(nfaces);
                        newStr += "\n(\n";

                        for(int iface=0; iface<nfaces; iface++)
                        {
                            int faceIndex = ca_patchFaceToFaceMap_inverse[ipatch][iface];

                            std::ostringstream doubleToOs;                            
                            doubleToOs << std::scientific << std::setprecision(IODigits);
                            doubleToOs << ca_patchP[ipatch][faceIndex];
                            newStr += removeTrailZero(doubleToOs.str());
                            newStr += "\n";
                        }
                        newStr += ")\n";
                    }
                    else if (fileName == "T" && ca_patchT!= nullptr)
                    {
                        newStr += "<scalar>\n";
                        //int nfaces = ca_patchSize[ipatch];
                        newStr += std::to_string(nfaces);
                        newStr += "\n(\n";

                        for(int iface=0; iface<nfaces; iface++)
                        {
                            int faceIndex = ca_patchFaceToFaceMap_inverse[ipatch][iface];

                            std::ostringstream doubleToOs;                            
                            doubleToOs << std::scientific << std::setprecision(IODigits);
                            doubleToOs << ca_patchT[ipatch][faceIndex];
                            newStr += removeTrailZero(doubleToOs.str());
                            newStr += "\n";
                        }
                        newStr += ")\n";
                    }
                    else if (fileName == "rho" && ca_patchRho!= nullptr)
                    {
                        newStr += "<scalar>\n";
                        //int nfaces = ca_patchSize[ipatch];
                        newStr += std::to_string(nfaces);
                        newStr += "\n(\n";

                        for(int iface=0; iface<nfaces; iface++)
                        {
                            int faceIndex = ca_patchFaceToFaceMap_inverse[ipatch][iface];

                            std::ostringstream doubleToOs;                            
                            doubleToOs << std::scientific << std::setprecision(IODigits);
                            doubleToOs << ca_patchRho[ipatch][faceIndex];
                            newStr += removeTrailZero(doubleToOs.str());
                            newStr += "\n";
                        }
                        newStr += ")\n";
                    }
                    else if (fileName == "phi" && ca_patchPhi!= nullptr)
                    {
                        newStr += "<scalar>\n";
                        //int nfaces = ca_patchSize[ipatch];
                        newStr += std::to_string(nfaces);
                        newStr += "\n(\n";

                        for(int iface=0; iface<nfaces; iface++)
                        {
                            int faceIndex = ca_patchFaceToFaceMap_inverse[ipatch][iface];

                            std::ostringstream doubleToOs;                            
                            doubleToOs << std::scientific << std::setprecision(IODigits);
                            doubleToOs << ca_patchPhi[ipatch][faceIndex];
                            newStr += removeTrailZero(doubleToOs.str());
                            newStr += "\n";
                        }
                        newStr += ")\n";
                    }
                    else if (fileName == "rhoUf" && ca_patchRhoUf!= nullptr)
                    {
                        newStr += "<vector>\n";
                        newStr += std::to_string(nfaces);
                        newStr += "\n(\n";

                        for(int iface=0; iface<nfaces; iface++)
                        {
                            int faceIndex = ca_patchFaceToFaceMap_inverse[ipatch][iface];

                            newStr += "(";
                            for(int jcomp=0; jcomp<nComponents; jcomp++)
                            {
                                int localComp = jcomp + faceIndex*nComponents;
                                
                                std::ostringstream doubleToOs;                            
                                doubleToOs << std::scientific << std::setprecision(IODigits);
                                doubleToOs << ca_patchRhoUf[ipatch][localComp];
                                newStr += removeTrailZero(doubleToOs.str());
                                if (jcomp<nComponents-1)
                                    newStr += " ";
                            }
                            newStr += ")\n";
                        }
                        newStr += ")\n";
                    }
                    // Turbulence data ^^^^^^^^^^^^^^
                    else if (fileName == "alphat" && ca_patchAlphaT!= nullptr)
                    {
                        newStr += "<scalar>\n";
                        //int nfaces = ca_patchSize[ipatch];
                        newStr += std::to_string(nfaces);
                        newStr += "\n(\n";

                        for(int iface=0; iface<nfaces; iface++)
                        {
                            int faceIndex = ca_patchFaceToFaceMap_inverse[ipatch][iface];

                            std::ostringstream doubleToOs;                            
                            doubleToOs << std::scientific << std::setprecision(IODigits);
                            doubleToOs << ca_patchAlphaT[ipatch][faceIndex];
                            newStr += removeTrailZero(doubleToOs.str());
                            newStr += "\n";
                        }
                        newStr += ")\n";
                    }
                    else if (fileName == "epsilon" && ca_patchEpsilon!= nullptr)
                    {
                        newStr += "<scalar>\n";
                        //int nfaces = ca_patchSize[ipatch];
                        newStr += std::to_string(nfaces);
                        newStr += "\n(\n";

                        for(int iface=0; iface<nfaces; iface++)
                        {
                            int faceIndex = ca_patchFaceToFaceMap_inverse[ipatch][iface];

                            std::ostringstream doubleToOs;                            
                            doubleToOs << std::scientific << std::setprecision(IODigits);
                            doubleToOs << ca_patchEpsilon[ipatch][faceIndex];
                            newStr += removeTrailZero(doubleToOs.str());
                            newStr += "\n";
                        }
                        newStr += ")\n";
                    }
                    else if (fileName == "omega" && ca_patchOmega!= nullptr)
                    {
                        newStr += "<scalar>\n";
                        //int nfaces = ca_patchSize[ipatch];
                        newStr += std::to_string(nfaces);
                        newStr += "\n(\n";

                        for(int iface=0; iface<nfaces; iface++)
                        {
                            int faceIndex = ca_patchFaceToFaceMap_inverse[ipatch][iface];

                            std::ostringstream doubleToOs;                            
                            doubleToOs << std::scientific << std::setprecision(IODigits);
                            doubleToOs << ca_patchOmega[ipatch][faceIndex];
                            newStr += removeTrailZero(doubleToOs.str());
                            newStr += "\n";
                        }
                        newStr += ")\n";
                    }
                    else if (fileName == "k" && ca_patchK!= nullptr)
                    {
                        newStr += "<scalar>\n";
                        //int nfaces = ca_patchSize[ipatch];
                        newStr += std::to_string(nfaces);
                        newStr += "\n(\n";

                        for(int iface=0; iface<nfaces; iface++)
                        {
                            int faceIndex = ca_patchFaceToFaceMap_inverse[ipatch][iface];

                            std::ostringstream doubleToOs;                            
                            doubleToOs << std::scientific << std::setprecision(IODigits);
                            doubleToOs << ca_patchK[ipatch][faceIndex];
                            newStr += removeTrailZero(doubleToOs.str());
                            newStr += "\n";
                        }
                        newStr += ")\n";
                    }
                    else if (fileName == "nut" && ca_patchNuT!= nullptr)
                    {
                        newStr += "<scalar>\n";
                        //int nfaces = ca_patchSize[ipatch];
                        newStr += std::to_string(nfaces);
                        newStr += "\n(\n";

                        for(int iface=0; iface<nfaces; iface++)
                        {
                            int faceIndex = ca_patchFaceToFaceMap_inverse[ipatch][iface];

                            std::ostringstream doubleToOs;                            
                            doubleToOs << std::scientific << std::setprecision(IODigits);
                            doubleToOs << ca_patchNuT[ipatch][faceIndex];
                            newStr += removeTrailZero(doubleToOs.str());
                            newStr += "\n";
                        }
                        newStr += ")\n";
                    }
                    else if ((fileName == "pointDisplacement" ||
                             fileName == "pointDisplacementNew")
                             && ca_Disp!= nullptr)
                    {
                        newStr += "<vector>\n";
                        
                        int npoints = *ca_patchPointToPointMap_size[ipatch];
                        newStr += std::to_string(npoints);
                        newStr += "\n(\n";

                        //int localIndex = 0;
                        for(int ipoint=0; ipoint<npoints; ipoint++)
                        {
                            int globalPointID = ca_patchPointToPointMap[ipatch][ipoint];

                            newStr += "(";
                            for(int jcomp=0; jcomp<nComponents; jcomp++)
                            {
                                int localComp = jcomp + globalPointID*nComponents;

                                std::ostringstream doubleToOs;                            
                                doubleToOs << std::scientific << std::setprecision(IODigits);
                                doubleToOs << ca_Disp[localComp];
                                newStr += removeTrailZero(doubleToOs.str());
                                if (jcomp<nComponents-1)
                                    newStr += " ";
                            }
                            newStr += ")\n";
                        }
                        newStr += ")\n";
                    }
                    //-------------------------------
                    
                    content.insert(valueStart, newStr);
                }
                // Files that are the in folder "0"" and
                // not being taken care of yet. Need to revisit this later
                else if (fileName == "cellToRegion")
                {}
                else
                {
                    std::cout << "========== WARNING ===============" << std::endl
                              << "Boundary field for file " 
                              << fileName << " is not found." << std::endl;
                }
            }

            std::string newLocalDir = locaParAddr+std::string(ca_timeName);
            fullDir = rootAddr+"/"+newLocalDir;
            fullAddr = fullDir+"/"+fileName;

            BF::path fullPath = fullDir;
            BF::create_directories(fullPath);

            std::ofstream outpuFile;
            outpuFile.open(fullAddr);
            
            outpuFile << content;
            outpuFile.close();
            
            std::cout << "    " << fullAddr
                      << " created." << std::endl;
        }
    } 
    return 0;
}

int comFoam::createUniformTimeFile(const std::string& rootAddr)
{
    namespace BF = boost::filesystem;

    std::string locaParAddr = "";
    if (ca_nProc>1)
    {
        locaParAddr = "processor"+std::to_string(ca_myRank)+"/";
    }

    std::string timeNameStr = std::string(ca_timeName);
    std::string localDir = timeNameStr+"/uniform";

    std::string
    content  = "/*--------------------------------*- C++ -*----------------------------------*\\\n";
    content += "  =========                 |\n";
    content += "  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox\n";
    content += "   \\\\    /   O peration     | Website:  https://openfoam.org\n";
    content += "    \\\\  /    A nd           | Version:  7\n";
    content += "     \\\\/     M anipulation  |\n";
    content += "\\*---------------------------------------------------------------------------*/\n";
    content += "FoamFile\n{\n";
    content += "    version     2.0;\n";
    content += "    format      ascii;\n";
    content += "    class       dictionary;\n";
    content += "    location    \""+localDir+"\";\n";
    content += "    object      time;\n}\n";
    content += "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n";

    std::ostringstream doubleToOs;
    doubleToOs << std::scientific << std::setprecision(IODigits);
    doubleToOs.str("");
    doubleToOs.clear();
    doubleToOs << *ca_time;
    content += "value      "+doubleToOs.str()+";\n";
    
    content += "name       \""+timeNameStr+"\";\n";
    
    content += "index      "+std::to_string(*ca_timeIndex)+";\n";

    doubleToOs.str("");
    doubleToOs.clear();
    doubleToOs << *ca_deltaT;

    content += "deltaT     "+removeTrailZero(doubleToOs.str())+";\n";

    doubleToOs.str("");
    doubleToOs.clear();
    doubleToOs << *ca_deltaT0;

    content += "deltaT0    "+removeTrailZero(doubleToOs.str())+";\n";
    content += "// ************************************************************************* //";

    localDir = locaParAddr+timeNameStr+"/uniform";
    std::string fullDir  = rootAddr+"/"+localDir;
    std::string fullAddr = fullDir+"/time";

    BF::path fullPath = fullDir;
    BF::create_directories(fullPath);

    std::ofstream outpuFile;
    outpuFile.open(fullAddr);
    
    outpuFile << content;
    outpuFile.close();

    std::cout << "    " << fullAddr
              << " created." << std::endl;

    return 0;
}

int comFoam::createPointsFile(const std::string& rootAddr)
{
    namespace BF = boost::filesystem;

    std::string locaParAddr = "";
    if (ca_nProc>1)
    {
        locaParAddr = "processor"+std::to_string(ca_myRank)+"/";
    }

    std::string localDir = "constant/polyMesh";
    std::string fullDir  = rootAddr+"/"+locaParAddr+localDir;
    std::string fullAddr = fullDir+"/points";

    std::string
    content  = "/*--------------------------------*- C++ -*----------------------------------*\\\n";
    content += "  =========                 |\n";
    content += "  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox\n";
    content += "   \\\\    /   O peration     | Website:  https://openfoam.org\n";
    content += "    \\\\  /    A nd           | Version:  7\n";
    content += "     \\\\/     M anipulation  |\n";
    content += "\\*---------------------------------------------------------------------------*/\n";
    content += "FoamFile\n{\n";
    content += "    version     2.0;\n";
    content += "    format      ascii;\n";
    content += "    class       vectorField;\n";
    content += "    location    \""+localDir+"\";\n";
    content += "    object      points;\n}\n";
    content += "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n";

    content += "\n\n";
    int npoints = *ca_nPoints;
    content += std::to_string(npoints);
    content += "\n(\n";
    
    for(int ipoint=0; ipoint<npoints; ipoint++)
    {
        content += "(";
        for(int jcomp=0; jcomp<nComponents; jcomp++)
        {
            int localComp = jcomp + ipoint*nComponents;
            
            std::ostringstream doubleToOs;
            doubleToOs << std::scientific << std::setprecision(IODigits);
            doubleToOs << ca_Points[localComp];

            content += removeTrailZero(doubleToOs.str());
            if (jcomp<nComponents-1)
                content += " ";
        }
        content += ")\n";
    }
    content += ")\n\n\n";
    
    content += "// ************************************************************************* //";

    BF::path fullPath = fullDir;
    BF::create_directories(fullPath);

    std::ofstream outpuFile;
    outpuFile.open(fullAddr);
    
    outpuFile << content;
    outpuFile.close();

    std::cout << "    " << fullAddr
              << " created." << std::endl;


    if (*ca_isDynamicFvMesh == 1)
    {
        std::string timeNameStr = std::string(ca_timeName);
        //std::string
        localDir = timeNameStr+"/polymesh";

        //std::string
        fullDir  = rootAddr+"/"+locaParAddr+localDir;
        //std::string
        fullAddr = fullDir+"/points";

        std::string locationStr = "location";
        size_t locationStart = findWord(fullAddr, content, locationStr);
        if (locationStart != 0 && locationStart != std::string::npos)
        {
            std::string scStr = ";";
            size_t locationEnd = findChar(fullAddr, content, scStr, locationStart);
            size_t length = locationEnd-locationStart;
            content.erase(locationStart, length);

            std::string newStr  = "location    \"";
            newStr += std::string(localDir);
            newStr += "\"";
            content.insert(locationStart, newStr);
        }

        //BF::path
        fullPath = fullDir;
        BF::create_directories(fullPath);

        //std::ofstream outpuFile;
        outpuFile.open(fullAddr);
        
        outpuFile << content;
        outpuFile.close();

        std::cout << "    " << fullAddr
                  << " created." << std::endl;
    }

    return 0;
}


int comFoam::createOwnerFile(const std::string& rootAddr)
{
    namespace BF = boost::filesystem;

    std::string locaParAddr = "";
    if (ca_nProc>1)
    {
        locaParAddr = "processor"+std::to_string(ca_myRank)+"/";
    }
    
    std::string localDir = "constant/polyMesh";
    std::string fullDir  = rootAddr+"/"+locaParAddr+localDir;
    std::string fullAddr = fullDir+"/owner";

    std::string
    content  = "/*--------------------------------*- C++ -*----------------------------------*\\\n";
    content += "  =========                 |\n";
    content += "  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox\n";
    content += "   \\\\    /   O peration     | Website:  https://openfoam.org\n";
    content += "    \\\\  /    A nd           | Version:  7\n";
    content += "     \\\\/     M anipulation  |\n";
    content += "\\*---------------------------------------------------------------------------*/\n";
    content += "FoamFile\n{\n";
    content += "    version     2.0;\n";
    content += "    format      ascii;\n";
    content += "    class       labelList;\n";
    content += "    note        \"nPoints: ";
    content += std::to_string(*ca_nPoints);
    content += " nCells: ";
    content += std::to_string(*ca_nCells);
    content += " nFaces: ";
    content += std::to_string(*ca_nFaces);
    content += " nInternalFaces: ";

    int procStartIndex{0};
    for (int iproc=0; iproc<ca_myRank; iproc++)
    {
        procStartIndex += ca_nPatches[iproc];
    }
    int index  = procStartIndex + 0; //ipatch;
    int nfacesInternal = ca_patchStart[index];

    content += std::to_string(nfacesInternal);
    content += "\";\n";
    content += "    location    \""+localDir+"\";\n";
    content += "    object      owner;\n}\n";
    content += "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n";

    content += "\n\n";
    int nfaces = *ca_nFaces;
    content += std::to_string(nfaces);
    content += "\n(\n";

    for(int iface=0; iface<nfaces; iface++)
    {
        int faceIndex = ca_faceToFaceMap_inverse[iface];

        content += std::to_string(ca_faceOwner[faceIndex]);
        content += "\n";
    }
    content += ")\n\n\n";
    
    content += "// ************************************************************************* //";

    BF::path fullPath = fullDir;
    BF::create_directories(fullPath);

    std::ofstream outpuFile;
    outpuFile.open(fullAddr);
    
    outpuFile << content;
    outpuFile.close();

    std::cout << "    " << fullAddr
              << " created." << std::endl;

    return 0;
}


int comFoam::createNeighborFile(const std::string& rootAddr)
{
    namespace BF = boost::filesystem;

    std::string locaParAddr = "";
    if (ca_nProc>1)
    {
        locaParAddr = "processor"+std::to_string(ca_myRank)+"/";
    }

    std::string localDir = "constant/polyMesh";
    std::string fullDir  = rootAddr+"/"+locaParAddr+localDir;
    std::string fullAddr = fullDir+"/neighbour";

    std::string
    content  = "/*--------------------------------*- C++ -*----------------------------------*\\\n";
    content += "  =========                 |\n";
    content += "  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox\n";
    content += "   \\\\    /   O peration     | Website:  https://openfoam.org\n";
    content += "    \\\\  /    A nd           | Version:  7\n";
    content += "     \\\\/     M anipulation  |\n";
    content += "\\*---------------------------------------------------------------------------*/\n";
    content += "FoamFile\n{\n";
    content += "    version     2.0;\n";
    content += "    format      ascii;\n";
    content += "    class       labelList;\n";
    content += "    note        \"nPoints: ";
    content += std::to_string(*ca_nPoints);
    content += " nCells: ";
    content += std::to_string(*ca_nCells);
    content += " nFaces: ";
    content += std::to_string(*ca_nFaces);
    content += " nInternalFaces: ";

    int procStartIndex{0};
    for (int iproc=0; iproc<ca_myRank; iproc++)
    {
        procStartIndex += ca_nPatches[iproc];
    }
    int index  = procStartIndex + 0; //ipatch;
    int nfacesInternal = ca_patchStart[index];

    content += std::to_string(nfacesInternal);
    content += "\";\n";
    content += "    location    \""+localDir+"\";\n";
    content += "    object      neighbour;\n}\n";
    content += "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n";

    content += "\n\n";
    content += std::to_string(nfacesInternal);
    content += "\n(\n";

    for(int iface=0; iface<nfacesInternal; iface++)
    {
        int faceIndex = ca_faceToFaceMap_inverse[iface];

        content += std::to_string(ca_faceNeighb[faceIndex]);
        content += "\n";
    }
    content += ")\n\n\n";
    
    content += "// ************************************************************************* //";

    BF::path fullPath = fullDir;
    BF::create_directories(fullPath);

    std::ofstream outpuFile;
    outpuFile.open(fullAddr);
    
    outpuFile << content;
    outpuFile.close();

    std::cout << "    " << fullAddr
              << " created." << std::endl;

    return 0;
}


int comFoam::createFacesFile(const std::string& rootAddr)
{
    namespace BF = boost::filesystem;

    std::string locaParAddr = "";
    if (ca_nProc>1)
    {
        locaParAddr = "processor"+std::to_string(ca_myRank)+"/";
    }
    
    std::string localDir = "constant/polyMesh";
    std::string fullDir  = rootAddr+"/"+locaParAddr+localDir;
    std::string fullAddr = fullDir+"/faces";

    std::string
    content  = "/*--------------------------------*- C++ -*----------------------------------*\\\n";
    content += "  =========                 |\n";
    content += "  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox\n";
    content += "   \\\\    /   O peration     | Website:  https://openfoam.org\n";
    content += "    \\\\  /    A nd           | Version:  7\n";
    content += "     \\\\/     M anipulation  |\n";
    content += "\\*---------------------------------------------------------------------------*/\n";
    content += "FoamFile\n{\n";
    content += "    version     2.0;\n";
    content += "    format      ascii;\n";
    content += "    class       faceList;\n";
    content += "    location    \""+localDir+"\";\n";
    content += "    object      faces;\n}\n";
    content += "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n";

    content += "\n\n";
    int nfaces = *ca_nFaces;
    content += std::to_string(nfaces);
    content += "\n(\n";

    for(int iface=0; iface<nfaces; iface++)
    {
        int faceIndex = ca_faceToFaceMap_inverse[iface];
        int ntypes     = *ca_faceToPointConn_types;
        int typeSelect = -1;
        int faceCountFloor = 0;

        for(int itype=0; itype<ntypes; itype++) 
        {
            if (faceIndex < faceCountFloor+ca_faceToPointConn_size[itype])
            {
                typeSelect = itype;
                break;
            }
            else
            {
                faceCountFloor += ca_faceToPointConn_size[itype];
            }
        }

        if (typeSelect == -1)
        {
            std::cout << "========== WARNING ===============" << std::endl
                      << "Face type not found." << std::endl;
        }
        int localFaceIndex = faceIndex - faceCountFloor;
        int npoints = ca_faceToPointConn_map[typeSelect];

        content += std::to_string(npoints);
        content += "(";
        for(int ipoint=0; ipoint<npoints; ipoint++)
        {
            int index = ipoint+localFaceIndex*npoints;

            content += std::to_string(ca_faceToPointConn[typeSelect][index] -1 );
            if (ipoint<npoints-1)
                content += " ";
        }
        content += ")\n";
    }
    content += ")\n\n\n";
    
    content += "// ************************************************************************* //";

    BF::path fullPath = fullDir;
    BF::create_directories(fullPath);

    std::ofstream outpuFile;
    outpuFile.open(fullAddr);

    outpuFile << content;
    outpuFile.close();

    std::cout << "    " << fullAddr
              << " created." << std::endl;

    return 0;
}

int comFoam::createBoundaryFile
            (
                const std::string& rootAddr,
                const std::vector<fileContainer>& vecFile
            )
{
    namespace BF = boost::filesystem;

    std::string locaParAddr = "";
    if (ca_nProc>1)
    {
        locaParAddr = "processor"+std::to_string(ca_myRank)+"/";
    }

    for(int ifile=0; ifile<*ca_nFiles; ifile++)
    {
        std::string fileName = vecFile[ifile].name;
        std::string localDir = vecFile[ifile].path;
        std::string fullDir = rootAddr+"/"+localDir;
        std::string fullAddr = fullDir+"/"+fileName;

        if (fileName == "boundary")
        {
            std::string content = vecFile[ifile].content;
            for(int ipatch=0; ipatch<ca_nPatches[ca_myRank]; ipatch++)
            {
                std::string patchName = patchNameStr[ipatch];
                size_t patchStart = 0;

                std::string brckOpen = "{";
                std::string brckClose = "}";
                size_t brckStart = std::string::npos;
                size_t brckEnd   = std::string::npos;

                bool found = false;
                int startTmp = 0;
                while (!found)
                {
                    patchStart = findWordOnly(fullAddr, content, patchName, startTmp);
                    brckStart = findChar(fullAddr, content, brckOpen, patchStart);
                    brckEnd   = findChar(fullAddr, content, brckClose, brckStart);

                    if (patchStart != std::string::npos)
                    {
                        if (brckStart!=std::string::npos &&
                            brckEnd!=std::string::npos)
                        {
                            if (patchStart<brckStart && brckStart<brckEnd)
                              found = true;  
                        }
                        else
                        {
                            startTmp = patchStart;
                        }
                    }
                    else
                    {
                        std::cout << "========== WARNING ===============" << std::endl
                                  << "patchName " << patchName
                                  << " not found." << std::endl;
                        exit(-1);
                    }
                }
 
                // Update nFaces of each patch in boundary file
                std::string itemStr = "nFaces";
                std::string scStr = ";";
                size_t itemStart = findWordOnly(fullAddr, content, itemStr, brckStart, brckEnd);
                size_t itemEnd = findChar(fullAddr, content, scStr, itemStart);

                if(itemStart == std::string::npos ||
                   !(brckStart<=itemStart && itemStart<=brckEnd) ||
                   !(brckStart<=itemEnd && itemEnd<=brckEnd)
                   )
                {
                    std::cout << "========== WARNING ===============" << std::endl
                              << "Patch " << patchName << "\'s nFaces limit"
                              << " not found." << std::endl;
                    exit(-1);
                }

                int length = itemEnd-itemStart;
                std::string subType = content.substr(itemStart, length);
                std::stringstream iStr(subType);
                std::string tmpStr;
                while(iStr >> tmpStr)
                {
                }
                int tmpInt = std::stoi(tmpStr);

                int procStartIndex{0};
                for (int iproc=0; iproc<ca_myRank; iproc++)
                {
                    procStartIndex += ca_nPatches[iproc];
                }
                int index  = procStartIndex + ipatch;
                int nfaces = ca_patchSize[index];
                if (tmpInt != nfaces)
                {
                    length = itemEnd-itemStart;
                    content.erase(itemStart, length);

                    tmpStr.clear();
                    tmpStr  = "nFaces          ";
                    tmpStr += std::to_string(nfaces);
               
                    content.insert(itemStart, tmpStr);
                }
                
                // Update startFace of each patch in boundary file
                itemStr = "startFace";
                scStr = ";";
                itemStart = findWordOnly(fullAddr, content, itemStr, brckStart, brckEnd);
                itemEnd = findChar(fullAddr, content, scStr, itemStart);

                if(itemStart == std::string::npos ||
                   !(brckStart<=itemStart && itemStart<=brckEnd) ||
                   !(brckStart<=itemEnd && itemEnd<=brckEnd)
                   )
                {
                    std::cout << "========== WARNING ===============" << std::endl
                              << "Patch " << patchName << "\'s startFace limit"
                              << " not found." << std::endl;
                    exit(-1);
                }

                length = itemEnd-itemStart;
                subType = content.substr(itemStart, length);
                iStr.str("");
                iStr.clear();
                
                iStr << subType;
                tmpStr.clear();
                while(iStr >> tmpStr)
                {
                }
                tmpInt = std::stoi(tmpStr);

                int patchStart_ = ca_patchStart[index];
                if (tmpInt != patchStart_)
                {
                    length = itemEnd-itemStart;
                    content.erase(itemStart, length);

                    tmpStr.clear();
                    tmpStr  = "startFace       ";
                    tmpStr += std::to_string(patchStart_);
               
                    content.insert(itemStart, tmpStr);
                }
            }
            
            BF::path fullPath = fullDir;
            BF::create_directories(fullPath);

            std::ofstream outpuFile;
            outpuFile.open(fullAddr);
            
            outpuFile << content;
            outpuFile.close();

            std::cout << "    " << fullAddr
                      << " created." << std::endl;
        }
    } 
    return 0;
}

int comFoam::createConnectivityFiles
            (
                const std::string& rootAddr,
                const std::vector<fileContainer>& vecFile
            )
{
    namespace BF = boost::filesystem;

    std::string locaParAddr = "";
    if (ca_nProc>1)
    {
        locaParAddr = "processor"+std::to_string(ca_myRank)+"/";
    }

    for(int ifile=0; ifile<*ca_nFiles; ifile++)
    {
        std::string fileName = vecFile[ifile].name;
        std::string localDir = vecFile[ifile].path;
        std::string fullDir = rootAddr+"/"+localDir;
        std::string fullAddr = fullDir+"/"+fileName;

        if (
            fileName == "boundaryProcAddressing" ||
            fileName == "cellProcAddressing" ||
            fileName == "faceProcAddressing" ||
            fileName == "pointProcAddressing"||

            fileName  == "cellZones" || //should be directly generated
            fileName  == "faceZones" || //should be directly generated
            fileName  == "pointZones"   //should be directly generated
           )
        {
            std::string content = vecFile[ifile].content;


            BF::path fullPath = fullDir;
            BF::create_directories(fullPath);

            std::ofstream outpuFile;
            outpuFile.open(fullAddr);
            
            outpuFile << content;
            outpuFile.close();

            std::cout << "    " << fullAddr
                      << " created." << std::endl;
        }
    } 
    return 0;
}

std::string comFoam::createBaseFile
(
    std::string name,
    std::string loc,
    std::string type,
    std::string dim
)
{

    if (type!="Scalar" && type!="Vector")
    {
        std::cout << "Watning: Invalid type \""
                  << type << "\" sent to createBaseFile"
                  << std::endl;
        exit(-1);
    }

    if (loc!="vol" && loc!="surface" && loc!="point")
    {
        std::cout << "Watning: Invalid location \""
                  << loc << "\" sent to createBaseFile"
                  << std::endl;
        exit(-1);
    }


    std::string content  = "/*--------------------------------*- C++ -*----------------------------------*\\\n";
    content += "  =========                 |\n";
    content += "  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox\n";
    content += "   \\\\    /   O peration     | Website:  https://openfoam.org\n";
    content += "    \\\\  /    A nd           | Version:  7\n";
    content += "     \\\\/     M anipulation  |\n";
    content += "\\*---------------------------------------------------------------------------*/\n";
    content += "FoamFile\n{\n";
    content += "    version     2.0;\n";
    content += "    format      ascii;\n";

    content += "    class       " + loc + type + "Field;\n";

    content += "    location    \""+std::string("0")+"\";\n";
    content += "    object      "+name+";\n}\n";
    content += "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";

    content += "dimensions      "+dim+";\n\n";

    content += "internalField   uniform 1;\n\n";
    content += "boundaryField\n";
    content += "{\n";

    for(int ipatch=0; ipatch<ca_nPatches[ca_myRank]; ipatch++)
    {
        std::string patchName = patchNameStr[ipatch];
        std::string patchType = patchTypeStr[ipatch];

        int procStartIndex{0};
        for (int iproc=0; iproc<ca_myRank; iproc++)
        {
            procStartIndex += ca_nPatches[iproc];
        }
        int index  = procStartIndex + ipatch;
        int nfaces = ca_patchSize[index];
        
        content += "    "+patchName+"\n";
        content += "    {\n";
        
        if (patchType == "empty")
        {
            content += "        type            "+patchType+";\n";
        }
        else
        {
            if (
                patchType=="processor"      ||
                patchType=="cyclic"         ||
                patchType=="cyclicAMI"      ||
                patchType=="cyclicSlip"     ||
                patchType=="symmetry"       ||
                patchType=="basicSymmetry"  ||
                patchType=="symmetryPlane"  ||
                patchType=="wedge")
            {
                content += "        type            "+patchType+";\n";
            }
            else
            {
                content += "        type            calculated;\n";
            }

            if (nfaces>0)
            {
                if (type == "Scalar")
                {
                    content += "        value           uniform 1;\n";
                }
                else if (type == "Vector")
                {
                    content += "        value           uniform (1 1 1);\n";
                }
            }
            else
            {
                if (type == "Scalar")
                {
                    content += "        value           nonuniform 0();\n";
                }
                else if (type == "Vector")
                {
                    content += "        value           nonuniform 0();\n";
                }
            }
        }
        content += "    }\n";
    }
    content += "}\n\n";
    content += "// ******************************";
    content += "******************************************* //";
    
    return content;
}

int comFoam::deleteTempFiles(const std::string& addr)
{
    namespace BF = boost::filesystem;

    BF::path rootPath(addr);
    BF::remove_all(rootPath);

    return 0;
}

int comFoam::registerFilesData(const char *name)
{
    std::string volName = name+std::string("VOL");
    
    std::cout << "rocFoam.registerFilesData: "
               << "Registering flow data with name "
               << volName
               << std::endl;

    int paneID = ca_myRank+1;// Use this paneID for file data

    std::cout << "procID = " << ca_myRank
              << ", paneID = " << paneID
              << " ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;


    // Genral file data ^^^^^^^^^^^^^^^^^^^^^^^^
    std::string dataName = volName+std::string(".nFiles");
    COM_new_dataitem( dataName, 'p', COM_INT, 1, "");
    COM_set_size(     dataName, paneID, 1);
    COM_set_array(    dataName, paneID, ca_nFiles);
    std::cout << "  " << dataName << " registered." << std::endl;
    
    dataName = volName+std::string(".fileSize");
    COM_new_dataitem( dataName, 'p', COM_INT, 1, "");
    COM_set_size(     dataName, paneID, *ca_nFiles);
    COM_set_array(    dataName, paneID, ca_fileSize);
    //std::cout << dataName << " registered." << std::endl;

    // Register file paths
    for (int ifile=0; ifile<*ca_nFiles; ifile++)
    {
        std::string charToStr = std::string(ca_filePath[ifile]);
        int charSize = charToStr.length()+1;

        dataName = volName+".filePath"+std::to_string(ifile);
        COM_new_dataitem( dataName, 'p', COM_CHAR, 1, "");
        COM_set_size( dataName, paneID, charSize);
        COM_set_array(dataName, paneID, ca_filePath[ifile]);
        //std::cout << dataName << " registered." << std::endl;
    }

    // Register file names
    for (int ifile=0; ifile<*ca_nFiles; ifile++)
    {
        std::string charToStr = std::string(ca_fileName[ifile]);
        int charSize = charToStr.length()+1;

        dataName = volName+".fileName"+std::to_string(ifile);
        COM_new_dataitem( dataName, 'p', COM_CHAR, 1, "");
        COM_set_size( dataName, paneID, charSize);
        COM_set_array(dataName, paneID, ca_fileName[ifile]);
        //std::cout << "  File "
        //          << ca_filePath[ifile]<<"/"<<ca_fileName[ifile]
        //          << " registered."
        //          << std::endl;
    }

    // Register file contents
    for (int ifile=0; ifile<*ca_nFiles; ifile++)
    {
        std::string charToStr = std::string(ca_fileContent[ifile]);
        int charSize = charToStr.length()+1;

        dataName = volName+".fileContent"+std::to_string(ifile);
        COM_new_dataitem( dataName, 'p', COM_CHAR, 1, "");
        COM_set_size( dataName, paneID, charSize);
        COM_set_array(dataName, paneID, ca_fileContent[ifile]);
        //std::cout << dataName << " registered." << std::endl;
    }

    for (int ifile=0; ifile<*ca_nFiles; ifile++)
    {
        std::cout << "  procID = " << ca_myRank
                  << ", " << ca_filePath[ifile]
                  << "/" << ca_fileName[ifile]
                  << " registered." << std::endl;
    }

    std::cout << "----------------------------------------------------"
              << std::endl;

    COM_window_init_done(volName);

    return 0;
}

int comFoam::reconstFilesData(const char *name)
{
    std::string volName = name+std::string("VOL");
    std::cout << "rocFoam.reconstructInitFiles, procID = "
              << ca_myRank
              << ", Retreiving file data form window "
              << volName << "."
              << std::endl;

    int paneID = ca_myRank+1;// Use this paneID for files

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
        }
    }
    std::cout << "  Number of items = " << dataItemNames.size()
              << std::endl << std::endl;

    // Genral file data ^^^^^^^^^^^^^^^^^^^^^^^^
    std::string dataName = std::string("nFiles");
    nameExists(dataItemNames, dataName);
    std::string regName = volName+std::string(".")+dataName;
    COM_get_array(regName.c_str(), paneID, &ca_nFiles);
    std::cout << "  " << dataName.c_str() << " = " << *ca_nFiles << std::endl;

    dataName = std::string("fileSize");
    nameExists(dataItemNames, dataName);
    regName = volName+std::string(".")+dataName;
    int nComp;
    COM_get_array(regName.c_str(), paneID, &ca_fileSize);
    COM_get_size(regName.c_str(), paneID, &nComp);

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
        COM_get_array(regName.c_str(), paneID, &ca_fileName[ifile]);
        COM_get_size(regName.c_str(), paneID, &nComp);

        dataName = std::string("filePath")+std::to_string(ifile);
        nameExists(dataItemNames, dataName);
        regName = volName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_filePath[ifile]);
        COM_get_size(regName.c_str(), paneID, &nComp);

        dataName = std::string("fileContent")+std::to_string(ifile);
        nameExists(dataItemNames, dataName);
        regName = volName+std::string(".")+dataName;
        COM_get_array(regName.c_str(), paneID, &ca_fileContent[ifile]);
        COM_get_size(regName.c_str(), paneID, &nComp);

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
        std::cout << "    file[" << ifile << "], "
                  << " path = " << ca_filePath[ifile] <<"/"<< ca_fileName[ifile] << ", "
                  << " size = " << ca_fileSize[ifile] 
                  << " retreived." << std::endl;
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

    createFilesData();

    return 0;
}


int comFoam::readRecursive
(
    const std::string& locaParAddr,
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

                //if (fileName == "pointDisplacementNew")
                //{
                //    return 0;
                //}

                if (fileShouldBeRead(locaParAddr, localAddr, fileName))
                {
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
                        tmpFile.content = std::string(content);
                        tmpFile.content.resize(size);
                        
                        vecFile.push_back(tmpFile);

                        //std::cout << vecFile[fileCount].name << std::endl;
                        //std::cout << vecFile[fileCount].path << std::endl;
                        //std::cout << vecFile[fileCount].size << std::endl;
                        //std::cout << vecFile[fileCount].content << std::endl;

                        if (content != NULL)
                        {
                            delete [] content;
                            content = NULL;
                        }

                        fileCount++;

                        /*
                        if (tmpFile.name == "pointDisplacement")
                        {
                            fileContainer tmpFile_;
                            tmpFile_.name = fileName+"New";
                            tmpFile_.path = localAddr;
                            tmpFile_.size = size;
                            tmpFile_.content = tmpFile.content;

                            std::ofstream newFile;
                            newFile.open(fullAddr+"New",
                                         std::ofstream::trunc);
                            newFile << tmpFile_.content;
                            newFile.close();

                            vecFile.push_back(tmpFile_);

                            fileCount++;
                        }
                        */
                    }
                }
            }
            else if (BF::is_directory(fullPath))
            {
                for (auto x : BF::directory_iterator(fullPath))
                {
                    BF::path newPath=x.path();
                    std::string neWaddr = BF::canonical(newPath).string();
                    readRecursive(locaParAddr, neWaddr, vecFile, fileCount);
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


bool comFoam::fileShouldBeRead
(
    const std::string& locaParAddr,
    const std::string& localAddr,
    const std::string& fileName
)
{
    // Only keeping these specific files and folders
    bool addFile = false;
    if (locaParAddr == "") 
    {
        if (localAddr == "system" ||
        localAddr == "constant") addFile = true;

        if (localAddr == "0") 
        {
            addFile = true;

            if (fileName  == "pointDisplacement" ||
                fileName  == "pointDisplacementNew")
            {
                if (ca_Disp == nullptr)
                    addFile = false;
            }
        }
        

        if (fileName  == "boundary" ||
            fileName  == "boundaryProcAddressing" ||
            fileName  == "cellProcAddressing" ||
            fileName  == "faceProcAddressing" ||
            fileName  == "pointProcAddressing" ||
            
            fileName  == "cellZones" || //should be directly generated
            fileName  == "faceZones" || //should be directly generated
            fileName  == "pointZones"   //should be directly generated
            ) addFile = true;

    }
    else
    {
        if (Pstream::master())
        {
            if (localAddr == "system" ||
            localAddr == "constant") addFile = true;
        }

        if (localAddr == locaParAddr+"/system" ||
            localAddr == locaParAddr+"/constant") addFile = true;
        if (localAddr == locaParAddr+"/0")
        {
            addFile = true;

            if (fileName  == "pointDisplacement" ||
                fileName  == "pointDisplacementNew")
            {
                if (ca_Disp == nullptr)
                    addFile = false;
            }
        }

        if (localAddr == locaParAddr+"/constant/polyMesh") 
        {
            if (fileName  == "boundary" ||
                fileName  == "boundaryProcAddressing" ||
                fileName  == "cellProcAddressing" ||
                fileName  == "faceProcAddressing" ||
                fileName  == "pointProcAddressing" ||
                
                fileName  == "cellZones" || //should be directly generated
                fileName  == "faceZones" || //should be directly generated
                fileName  == "pointZones"   //should be directly generated
                ) addFile = true;
        }
    }
    
    return addFile;
}


int comFoam::deleteFilesData()
{
    if (ca_fileContent != nullptr)
    {
        for (int ifile=0; ifile<*ca_nFiles; ifile++)
        {
            if (ca_fileContent[ifile] != nullptr)
            {
                delete [] ca_fileContent[ifile];
                ca_fileContent[ifile] = nullptr;
            }
        }
        delete [] ca_fileContent;
        ca_fileContent = nullptr;
    }

    if (ca_filePath != nullptr)
    {
        for (int ifile=0; ifile<*ca_nFiles; ifile++)
        {
            if (ca_filePath[ifile] != nullptr)
            {
                delete [] ca_filePath[ifile];
                ca_filePath[ifile] = nullptr;
            }
        }
        delete [] ca_filePath;
        ca_filePath = nullptr;
    }

    if (ca_fileName != nullptr)
    {
        for (int ifile=0; ifile<*ca_nFiles; ifile++)
        {
            if (ca_fileName[ifile] != nullptr)
            {
                delete [] ca_fileName[ifile];
                ca_fileName[ifile] = nullptr;
            }
        }
        delete [] ca_fileName;
        ca_fileName = nullptr;
    }

    if (ca_fileSize != nullptr)
    {
        delete [] ca_fileSize;
        ca_fileSize = nullptr;
    }

    if (ca_nFiles != nullptr)
    {
        delete ca_nFiles;
        ca_nFiles = nullptr;
    }

    return 0;
}

int comFoam::findGlobalIndex(int* arr, const int& size,  const int& elem)
{
    int index = -1;

    auto itr = std::find(arr, arr+size, elem);
    if (itr != arr+size)
    {
        index = std::distance(arr, itr);
    }
    else
    {
        std::cout << "========== WARNING ===============" << std::endl
                  << "Element is not found." << std::endl;
    }
    return index;
}



size_t comFoam::findChar
(
    const std::string& fullAddr,
    const std::string& content,
    const std::string& exp,
    size_t startPos,
    size_t endPos
)
{
    size_t expStart = content.find(exp, startPos);
    if (expStart==std::string::npos || expStart>endPos)
    {
        std::cout << "Warning: Cannot find \""
                  << exp+"\" keyword"
                  << " in the range "+std::to_string(startPos)
                  << " and "+std::to_string(endPos)
                  << " in file "
                  << fullAddr << std::endl;
    }

    return expStart;
}



size_t comFoam::findWord
(
    const std::string& fullAddr,
    const std::string& content,
    const std::string& exp,
    size_t startPos,
    size_t endPos
)
{
    size_t firstLoc=0;
    int count=0;
    size_t expStart=0;

    while (expStart!=std::string::npos && expStart<endPos)
    {
        expStart = min( content.find(exp, startPos), endPos);

        if (expStart != std::string::npos)
        {
            size_t startTmp  = std::max(static_cast<size_t>(0),
                                        expStart-exp.length());
            size_t endTmp    = std::min(expStart+2*exp.length(),
                                        std::string::npos);
            size_t lengthTmp = endTmp-startTmp+1;

            std::stringstream iStr(content.substr(startTmp, lengthTmp));
            std::string lookUp;

            while(iStr >> lookUp)
            {
                if (lookUp == exp)
                {
                    count++;
                    if (count == 1) firstLoc = expStart;
                    break;
                }
            }
        }
        startPos = expStart+exp.length();
        if (count == 1) break;
    }

    if (count == 0)
    {
        std::cout << "Warning: Cannot find \""
                  << exp+"\" keyword in file "
                 << fullAddr << std::endl;
         //exit(-1);
    }
    else if (count > 1)
    {
        std::cout << "Warning: Found more than one \""
                  << exp+"\" keyword in file "
                    << fullAddr << std::endl;
        //exit(-1);
    }

    return firstLoc;
}


size_t comFoam::findWordOnly
(
    const std::string& fullAddr,
    const std::string& content,
    const std::string& exp,
    size_t startPos,
    size_t endPos
)
{
    size_t firstLoc=0;
    int count=0;
    size_t expStart=0;

    while (expStart!=std::string::npos && expStart<endPos)
    {
        expStart = min( content.find(exp, startPos), endPos);

        if (expStart != std::string::npos)
        {
            size_t startTmp  = std::max(static_cast<size_t>(0),
                                        expStart-exp.length());
            size_t endTmp    = std::min(expStart+2*exp.length(),
                                        std::string::npos);
            size_t lengthTmp = endTmp-startTmp+1;

            std::stringstream iStr(content.substr(startTmp, lengthTmp));
            std::string lookUp;

            while(iStr >> lookUp)
            {
                if (lookUp == exp)
                {
                    count++;
                    if (count == 1) firstLoc = expStart;
                    break;
                }
            }
        }
        startPos = expStart+exp.length();
    }

    if (count == 0)
    {
        std::cout << "Warning: Cannot find \""
                  << exp+"\" keyword in file "
                 << fullAddr << std::endl;
         exit(-1);
    }
    else if (count > 1)
    {
        std::cout << "Warning: Found more than one \""
                  << exp+"\" keyword in file "
                    << fullAddr << std::endl;
        exit(-1);
    }

    return firstLoc;
}


std::string comFoam::removeTrailZero(std::string in)
{
    size_t decLoc = in.find(".");
    size_t powLoc = in.find("e");
    //size_t endLoc = in.size();
    
    if(decLoc != std::string::npos)
    {
        //decLoc++;
        int length = powLoc-decLoc;
        std::string decimal = in.substr(decLoc, length);

        // Remove trailing zeroes
        decimal = decimal.substr(0, decimal.find_last_not_of('0')+1);

        // If the decimal point is now the last character, remove that as well
        if(decimal.find('.') == 0 && decimal.find('.') == decimal.size()-1)
        {
            decimal = decimal.substr(0, decimal.size()-1);
        }

        in.erase(decLoc, length);
        in.insert(decLoc, decimal);
    }    

    return in;
}


/*
std::string comFoam::removeTrailZero(std::string in)
{
    size_t decLoc = in.find(".");
    size_t powLoc = in.find("e");
    size_t endLoc = in.size();
    
    std::string decimal = "";

std::cout << "in = " << in << ", decLoc = " << decLoc 
          << ", powLoc = " << powLoc << ", endLoc = " << endLoc << std::endl;

    
    if(powLoc != std::string::npos)
    {
        powLoc++;
        int length = powLoc-endLoc;
        std::string powStr = in.substr(powLoc, length);
        
        int powInt = std::stoi(powStr);
        
        if (powInt != 0)
        {
            decimal = in.substr(0, powLoc+1) + std::to_string(powInt);
        }
        else
        {
            decimal = in.substr(0, powLoc);
        }
        
        powLoc++;

std::cout << "powStr = " << powStr << ", powInt = " << powInt << ", decimal = " << decimal << std::endl;

    }    
    else
    {
        powLoc = endLoc+1;
        decimal = in;
    }


    
    if(decLoc != std::string::npos)
    {
        decLoc++;
        int length = powLoc-decLoc;
        decimal = decimal.substr(decLoc, length);
        // Remove trailing zeroes
        decimal = decimal.substr(0, decimal.find_last_not_of('0')+1);

        // If the decimal point is now the last character, remove that as well
        if(decimal.find('.') == decimal.size()-1)
        {
            decimal = decimal.substr(0, decimal.size()-1);
        }

        decimal.erase(decLoc, length);
        decimal.insert(decLoc, decimal);
    }    

    return decimal;
}
*/



