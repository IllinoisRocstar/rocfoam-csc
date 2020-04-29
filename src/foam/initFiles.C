int comFoam::readInitFiles(const std::string& rootAddr)
{
    std::string fullAddr=rootAddr;
    std::string locaParAddr = "";
    
    int myRank = 0;
    if (Pstream::parRun())
    {
        myRank = Pstream::myProcNo();
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

    std::cout << std::endl;
    
    createSysConstFiles(rootAddr, vecFile);
    createFieldFiles(rootAddr, vecFile);
    createUniformTimeFile(rootAddr);
    createPointsFile(rootAddr);
    createOwnerFile(rootAddr);
    createNeighborFile(rootAddr);
    createFacesFile(rootAddr);
    createBoundaryFile(rootAddr, vecFile);
    createConnectivityFiles(rootAddr, vecFile);

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
    for(int ifile=0; ifile<*ca_nFiles; ifile++)
    {
        std::string locaParAddr = "";
        if (ca_nProc>1)
        {
            locaParAddr = "processor"+std::to_string(ca_myRank)+"/";
        }

        std::string fileName = vecFile[ifile].name;
        std::string localDir = vecFile[ifile].path;
        std::string fullDir = rootAddr+"/"+localDir;
        std::string fullAddr = fullDir+"/"+fileName;

        //if (localDir=="0")
        if (localDir==locaParAddr+"0")
        {
            // Modify the internalField
            std::string content = vecFile[ifile].content;

            std::string locationStr = "location";
            size_t locationStart = content.find(locationStr);
            if (locationStart == std::string::npos)
            {
                std::cout << "Warning: Cannot find location in file "
                         << fileName << std::endl;
            }
            else
            {
                std::string scStr = ";";
                size_t locationEnd = content.find(scStr, locationStart);
                size_t length = locationEnd-locationStart;
                content.erase(locationStart, length);
                std::string newStr  = "location    \"";
                newStr += std::string(ca_timeName);
                newStr += "\"";
                content.insert(locationStart, newStr);
            }


            std::string intFieldStr = "internalField";
            size_t intFieldStart = content.find(intFieldStr);
            if (intFieldStart == std::string::npos)
            {
                std::cout << "Warning: Cannot find internalField in file"
                         << fullAddr << std::endl;
            }
            std::string scStr = ";";
            size_t intFieldEnd = content.find(scStr, intFieldStart);
            size_t length = intFieldEnd-intFieldStart;
            content.erase(intFieldStart, length);
            
            std::string newStr = "internalField   nonuniform List";
            if (fileName == "U")
            {
                newStr += "<vector>\n";
                int ncells = *ca_nCells;
                newStr += std::to_string(ncells);
                newStr += "\n(\n";
                
                for(int icell=0; icell<ncells; icell++)
                {
                    int cellIndex = findGlobalIndex
                                    (
                                        ca_cellToCellMap,
                                        ncells,
                                        icell
                                    );

                    newStr += "(";
                    for(int jcomp=0; jcomp<nComponents; jcomp++)
                    {
                        int localComp = jcomp + cellIndex*nComponents;
                        
                        std::ostringstream doubleToOs;                            
                        //doubleToOs << std::fixed;
                        doubleToOs << std::setprecision(IODigits);
                        doubleToOs << ca_Vel[localComp];

                        newStr += doubleToOs.str();
                        if (jcomp<nComponents-1)
                            newStr += " ";
                    }
                    newStr += ")\n";
                }
                newStr += ")\n";
                
                content.insert(intFieldStart, newStr);
            }
            else if (fileName == "p")
            {
                newStr += "<scalar>\n";
                int ncells = *ca_nCells;
                newStr += std::to_string(ncells);
                newStr += "\n(\n";
                
                for(int icell=0; icell<ncells; icell++)
                {
                    int cellIndex = findGlobalIndex
                                    (
                                        ca_cellToCellMap,
                                        ncells,
                                        icell
                                    );

                    std::ostringstream doubleToOs;                            
                    doubleToOs << std::setprecision(IODigits);
                    doubleToOs << ca_P[cellIndex];

                    newStr += doubleToOs.str();
                    newStr += "\n";
                }
                newStr += ")\n";
                
                content.insert(intFieldStart, newStr);
            }
            else if (fileName == "T")
            {
                newStr += "<scalar>\n";
                int ncells = *ca_nCells;
                newStr += std::to_string(ncells);
                newStr += "\n(\n";
                
                for(int icell=0; icell<ncells; icell++)
                {
                    int cellIndex = findGlobalIndex
                                    (
                                        ca_cellToCellMap,
                                        ncells,
                                        icell
                                    );

                    std::ostringstream doubleToOs;                            
                    doubleToOs << std::setprecision(IODigits);
                    doubleToOs << ca_T[cellIndex];

                    newStr += doubleToOs.str();
                    newStr += "\n";
                }
                newStr += ")\n";
                
                content.insert(intFieldStart, newStr);
            }
            else if (fileName == "rho")
            {
                newStr += "<scalar>\n";
                int ncells = *ca_nCells;
                newStr += std::to_string(ncells);
                newStr += "\n(\n";
                
                for(int icell=0; icell<ncells; icell++)
                {
                    int cellIndex = findGlobalIndex
                                    (
                                        ca_cellToCellMap,
                                        ncells,
                                        icell
                                    );

                    std::ostringstream doubleToOs;                            
                    doubleToOs << std::setprecision(IODigits);
                    doubleToOs << ca_Rho[cellIndex];

                    newStr += doubleToOs.str();
                    newStr += "\n";
                }
                newStr += ")\n";
                
                content.insert(intFieldStart, newStr);
            }
            else
            {
                std::cout << "========== WARNING ===============" << std::endl
                          << "Solution field for file " 
                          << fileName << "is not found." << std::endl;
            }
            
            for(int ipatch=0; ipatch<*ca_nPatches; ipatch++)
            {
                int nfaces = *ca_patchSize[ipatch];
                if (nfaces<=0)
                {
                    continue;
                }
                
                std::string patchName = patchNameStr[ipatch];

                size_t patchStart = content.find(patchName);
                if (patchStart == std::string::npos)
                {
                    std::cout << "Warning: patchName " << patchName
                              << " was not found" << std::endl;
                }
                
                // Assume only 1 patch of this name
                std::string typeStr = "type";
                scStr = ";";
                size_t typeStart = content.find(typeStr, patchStart);
                size_t typeEnd = content.find(scStr, typeStart);

                length = typeEnd-typeStart;
                std::string subType = content.substr(typeStart, length);
                std::stringstream iStr(subType);

                std::string tmpStr;
                while(iStr >> tmpStr)
                {
                }
                if (tmpStr=="fixedValue"    ||
                    tmpStr=="symmetryPlane" ||
                    tmpStr=="slip"          ||
                    tmpStr=="zeroGradient"  ||
                    tmpStr=="empty"
                   ) continue;

                
                std::string valueStr = "value";
                size_t valueStart = content.find(valueStr, typeStart);
                if (valueStart != std::string::npos)
                {
                    size_t valueEnd = content.find(scStr, valueStart);

                    length = valueEnd-valueStart;
                    content.erase(valueStart, length);

                    newStr.clear();
                    newStr = "value           nonuniform List";
                    if (fileName == "U")
                    {
                        newStr += "<vector>\n";
                        newStr += std::to_string(nfaces);
                        newStr += "\n(\n";

                        for(int iface=0; iface<nfaces; iface++)
                        {
                            int faceIndex = findGlobalIndex
                                            (
                                                ca_patchFaceToFaceMap[ipatch],
                                                nfaces,
                                                iface
                                            );

                            newStr += "(";
                            for(int jcomp=0; jcomp<nComponents; jcomp++)
                            {
                                int localComp = jcomp + faceIndex*nComponents;
                                
                                std::ostringstream doubleToOs;                            
                                //doubleToOs << std::fixed;
                                doubleToOs << std::setprecision(IODigits);
                                doubleToOs << ca_patchVel[ipatch][localComp];

                                newStr += doubleToOs.str();
                                if (jcomp<nComponents-1)
                                    newStr += " ";
                            }
                            newStr += ")\n";
                        }
                        newStr += ")\n";
                    }
                    else if (fileName == "p")
                    {
                        newStr += "<scalar>\n";
                        int nfaces = *ca_patchSize[ipatch];
                        newStr += std::to_string(nfaces);
                        newStr += "\n(\n";

                        for(int iface=0; iface<nfaces; iface++)
                        {
                            int faceIndex = findGlobalIndex
                                            (
                                                ca_patchFaceToFaceMap[ipatch],
                                                nfaces,
                                                iface
                                            );

                            std::ostringstream doubleToOs;                            
                            doubleToOs << std::setprecision(IODigits);
                            doubleToOs << ca_patchP[ipatch][faceIndex];

                            newStr += doubleToOs.str();
                            newStr += "\n";
                        }
                        newStr += ")\n";
                    }
                    else if (fileName == "T")
                    {
                        newStr += "<scalar>\n";
                        int nfaces = *ca_patchSize[ipatch];
                        newStr += std::to_string(nfaces);
                        newStr += "\n(\n";

                        for(int iface=0; iface<nfaces; iface++)
                        {
                            int faceIndex = findGlobalIndex
                                            (
                                                ca_patchFaceToFaceMap[ipatch],
                                                nfaces,
                                                iface
                                            );

                            std::ostringstream doubleToOs;                            
                            doubleToOs << std::setprecision(IODigits);
                            doubleToOs << ca_patchT[ipatch][faceIndex];

                            newStr += doubleToOs.str();
                            newStr += "\n";
                        }
                        newStr += ")\n";
                    }
                    else if (fileName == "rho")
                    {
                        newStr += "<scalar>\n";
                        int nfaces = *ca_patchSize[ipatch];
                        newStr += std::to_string(nfaces);
                        newStr += "\n(\n";

                        for(int iface=0; iface<nfaces; iface++)
                        {
                            int faceIndex = findGlobalIndex
                                            (
                                                ca_patchFaceToFaceMap[ipatch],
                                                nfaces,
                                                iface
                                            );

                            std::ostringstream doubleToOs;                            
                            doubleToOs << std::setprecision(IODigits);
                            doubleToOs << ca_patchRho[ipatch][faceIndex];

                            newStr += doubleToOs.str();
                            newStr += "\n";
                        }
                        newStr += ")\n";
                    }
                    else
                    {
                        std::cout << "========== WARNING ===============" << std::endl
                                  << "Boundary field for file " 
                                  << fileName << "is not found." << std::endl;
                    }
                    content.insert(valueStart, newStr);
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

    std::ostringstream doubleToOs;
    doubleToOs << std::setprecision(IODigits);

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

    doubleToOs.str("");
    doubleToOs.clear();
    doubleToOs << *ca_time;

    content += "value      "+doubleToOs.str()+";\n";
    content += "name       \""+timeNameStr+"\";\n";
    content += "index      "+std::to_string(*ca_timeIndex)+";\n";

    doubleToOs.str("");
    doubleToOs.clear();
    doubleToOs << *ca_deltaT;

    content += "deltaT     "+doubleToOs.str()+";\n";

    doubleToOs.str("");
    doubleToOs.clear();
    doubleToOs << *ca_deltaT0;

    content += "deltaT0    "+doubleToOs.str()+";\n";
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

    std::ostringstream doubleToOs;
    doubleToOs << std::setprecision(IODigits);

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
            //doubleToOs << std::fixed;
            doubleToOs << std::setprecision(IODigits);
            doubleToOs << ca_Points[localComp];

            content += doubleToOs.str();
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

    std::ostringstream doubleToOs;
    doubleToOs << std::setprecision(IODigits);

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
    content += std::to_string(*ca_patchStart[0]);
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
        int faceIndex = findGlobalIndex
                        (
                            ca_faceToFaceMap,
                            nfaces,
                            iface
                        );

        std::ostringstream doubleToOs;                            
        doubleToOs << std::setprecision(IODigits);
        doubleToOs << ca_faceOwner[faceIndex];

        content += doubleToOs.str();
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

    std::ostringstream doubleToOs;
    doubleToOs << std::setprecision(IODigits);

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
    content += std::to_string(*ca_patchStart[0]);
    content += "\";\n";
    content += "    location    \""+localDir+"\";\n";
    content += "    object      neighbour;\n}\n";
    content += "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n";

    content += "\n\n";
    int nfaces = *ca_patchStart[0];
    content += std::to_string(nfaces);
    content += "\n(\n";

    for(int iface=0; iface<nfaces; iface++)
    {
        int faceIndex = findGlobalIndex
                        (
                            ca_faceToFaceMap,
                            nfaces,
                            iface
                        );

        std::ostringstream doubleToOs;                            
        doubleToOs << std::setprecision(IODigits);
        doubleToOs << ca_faceNeighb[faceIndex];

        content += doubleToOs.str();
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

    std::ostringstream doubleToOs;
    doubleToOs << std::setprecision(IODigits);

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
        int faceIndex = findGlobalIndex
                        (
                            ca_faceToFaceMap,
                            nfaces,
                            iface
                        );

        int ntypes = *ca_faceToPointConn_types;
        int typeSelect = -1;
        int faceCountFloor = 0;
        for(int itype=0; itype<ntypes; itype++) 
        {
            if (faceIndex< faceCountFloor+ca_faceToPointConn_size[itype])
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
            
            content += std::to_string(ca_faceToPointConn[typeSelect][index]);
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
            for(int ipatch=0; ipatch<*ca_nPatches; ipatch++)
            {
                std::string patchName = patchNameStr[ipatch];

                size_t patchStart = content.find(patchName);
                if (patchStart == std::string::npos)
                {
                    std::cout << "========== WARNING ===============" << std::endl
                              << "patchName " << patchName
                              << " not found." << std::endl;
                }
                
                // Assume only 1 patch of this name
                std::string startStr = "{";
                std::string endtStr = "}";
                size_t typeStart = content.find(startStr, patchStart);
                size_t typeEnd = content.find(endtStr, patchStart);

                // Update nFaces of each patch in boundary file
                std::string itemStr = "nFaces";
                std::string scStr = ";";
                size_t itemStart = content.find(itemStr, typeStart);
                size_t itemEnd = content.find(scStr, itemStart);

                if(patchStart == std::string::npos ||
                   !(typeStart<=itemStart && itemStart<=typeEnd) ||
                   !(typeStart<=itemEnd && itemEnd<=typeEnd)
                   )
                {
                    std::cout << "========== WARNING ===============" << std::endl
                              << "Patch " << patchName << "\'s nFaces limit"
                              << " not found." << std::endl;
                }
                int length = itemEnd-itemStart;
                std::string subType = content.substr(itemStart, length);
                std::stringstream iStr(subType);
                std::string tmpStr;
                while(iStr >> tmpStr)
                {
                }
                int tmpInt = std::stoi(tmpStr);
                if (tmpInt != *ca_patchSize[ipatch])
                {
                    length = itemEnd-itemStart;
                    content.erase(itemStart, length);

                    tmpStr.clear();
                    tmpStr  = "nFaces          ";
                    tmpStr += std::to_string(*ca_patchSize[ipatch]);
               
                    content.insert(itemStart, tmpStr);
                }
                
                // Update startFace of each patch in boundary file
                itemStr = "startFace";
                scStr = ";";
                itemStart = content.find(itemStr, typeStart);
                itemEnd = content.find(scStr, itemStart);

                if(patchStart == std::string::npos ||
                   !(typeStart<=itemStart && itemStart<=typeEnd) ||
                   !(typeStart<=itemEnd && itemEnd<=typeEnd)
                   )
                {
                    std::cout << "========== WARNING ===============" << std::endl
                              << "Patch " << patchName << "\'s startFace limit"
                              << " not found." << std::endl;
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

                if (tmpInt != *ca_patchStart[ipatch])
                {
                    length = itemEnd-itemStart;
                    content.erase(itemStart, length);

                    tmpStr.clear();
                    tmpStr  = "startFace       ";
                    tmpStr += std::to_string(*ca_patchStart[ipatch]);
               
                    content.insert(itemStart, tmpStr);
                }
            }
            
            /*
            std::string localDir = "constant/polyMesh";
            std::string fullDir  = rootAddr+"/"+localDir;
            std::string fullAddr = fullDir+"/faces";
            */

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
            fileName == "pointProcAddressing"
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

    int paneID = Pstream::myProcNo()+1;// Use this paneID for volume connectivity

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
    for(int ifile=0; ifile<*ca_nFiles; ifile++)
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
    for(int ifile=0; ifile<*ca_nFiles; ifile++)
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
    for(int ifile=0; ifile<*ca_nFiles; ifile++)
    {
        std::string charToStr = std::string(ca_fileContent[ifile]);
        int charSize = charToStr.length()+1;

        dataName = volName+".fileContent"+std::to_string(ifile);
        COM_new_dataitem( dataName, 'p', COM_CHAR, 1, "");
        COM_set_size( dataName, paneID, charSize);
        COM_set_array(dataName, paneID, ca_fileContent[ifile]);
        //std::cout << dataName << " registered." << std::endl;
    }

    for(int ifile=0; ifile<*ca_nFiles; ifile++)
    {
        std::cout << "  procID = " << Pstream::myProcNo()
                  << ", " << ca_filePath[ifile]
                  << "/" << ca_fileName[ifile]
                  << " registered." << std::endl;
    }

    std::cout << "----------------------------------------------------"
              << std::endl;

    return 0;
}

int comFoam::reconstCaInitFiles(const char *name, const std::string& rootAddr)
{
    std::string volName = name+std::string("VOL");
    std::cout << "rocFoam.reconstructInitFiles, proID = "
              << Pstream::myProcNo()
              << ", Retreiving file data form window "
              << volName << "."
              << std::endl;

    int paneID = Pstream::myProcNo()+1;// Use this paneID for files

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
        std::cout << "    file[" << ifile << "], size = "
                  << ca_fileSize[ifile] << ", "
                  << ca_filePath[ifile] <<"/"<< ca_fileName[ifile] 
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

    createInitFiles(rootAddr);

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
                        tmpFile.content = string(content);
                        tmpFile.content.resize(size);
                        
                        vecFile.push_back(tmpFile);

                        //std::cout << vecFile[fileCount].name << std::endl;
                        //std::cout << vecFile[fileCount].path << std::endl;
                        //std::cout << vecFile[fileCount].size << std::endl;
                        //std::cout << vecFile[fileCount].content << std::endl;

                        fileCount++;
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
        localAddr == "constant" ||
        localAddr == "0") addFile = true;

        if (fileName  == "boundary" ||
            fileName  == "boundaryProcAddressing" ||
            fileName  == "cellProcAddressing" ||
            fileName  == "faceProcAddressing" ||
            fileName  == "pointProcAddressing"
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
            localAddr == locaParAddr+"/constant" ||
            localAddr == locaParAddr+"/0") addFile = true;
        
        if (localAddr == locaParAddr+"/constant/polyMesh") 
        {
            if (fileName  == "boundary" ||
                fileName  == "boundaryProcAddressing" ||
                fileName  == "cellProcAddressing" ||
                fileName  == "faceProcAddressing" ||
                fileName  == "pointProcAddressing"
                ) addFile = true;
        }
    }
    
    return addFile;
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
