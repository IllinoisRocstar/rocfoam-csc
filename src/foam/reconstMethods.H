void comFoam::copyWindow(const char *name1, const char *name2)
{
    std::string winName1 = name1; //+string("VOL");
    std::string winName2 = name2; //+string("VOL");
    Info << endl
         << "rocFoam.copyWindow: window names "
         << winName1 << " & "
         << winName2 << endl;


    // Get window information ^^^^^^^^^^^^^^^^^^^
    int nItemps;
    std::string datNames;
    COM_get_dataitems(winName1.c_str(), &nItemps, datNames);

    Info << "rocFoam.copyWindow: nItemps in "
         << "winName1 = "
         << nItemps
         << endl;
    // ------------------------------------------

    //  Copy pane nodes & connectivities ^^^^^^^^
    int nPanes;
    int* paneList;
    COM_get_panes(winName1.c_str(), &nPanes, &paneList);
    Info << "rocFoam.copyWindow: Number of Panes in "
         << "winName1 = "
         << nPanes << endl;

    for(int ipane=0; ipane<nPanes; ipane++)
    {
        COM_clone_dataitem((winName2+".mesh").c_str(),
                           (winName1+".mesh").c_str(),
                           ipane);
    }
    // ------------------------------------------


    //  Window data ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    std::istringstream issDataNames(datNames);
    std::vector<std::string> vecDataNames;
    for (int item=0; item<nItemps; ++item)
    {
        std::string dataName="";

        issDataNames >> dataName;
        vecDataNames.push_back(dataName);

        std::string dataName1 = winName1+string(".")+dataName;
        std::string dataName2 = winName2+string(".")+dataName;

        COM_clone_dataitem(dataName2.c_str(),
                           dataName1.c_str());

        //COM_clone_dataitem(winName2.c_str(), winName1.c_str(), 0, dataName.c_str(), 0);
        //COM_clone_dataitem(winName2, winName1, 0, dataName);

        Info << "rocFoam.copyWindow: dataName["
             << item << "] = " << dataName
             << ", copied from " << winName1
             << " to " << winName2 << endl;
    }
    Info << endl;
    // ------------------------------------------

    COM_window_init_done(winName2); 
}
