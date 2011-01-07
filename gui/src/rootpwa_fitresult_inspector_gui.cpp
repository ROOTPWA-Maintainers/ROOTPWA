//============================================================================
// Name        : rootpwa_fitresult_inspector_gui.cpp
// Author      : Promme@web.de, jasinski@kph.uni-mainz.de
// Version     :
// Copyright   : gnu
// Description : fitresult inspector spinoff from gui
//				 to use it as a standalone program
//============================================================================

#include <iostream>
#include "TApplication.h"
#include "TCanvas.h"
#include "TrpwaPlotAmpsFrame.h"
#include "TrpwaCommonTools.h"
#include <sstream>
#include <string>
#include <vector>

using namespace std;
using namespace TrpwaCommonTools;

int main(int argc, char *argv[]) {
	cout << " __________________________________________" << endl;
	cout << endl;
	cout << "                rootpwa GUI                " << endl;
	cout << "           Fit result inspector            " << endl;
	cout << "                 vers.0.1                  " << endl;
	cout << endl;
	cout << " __________________________________________" << endl;

	vector<string> fit_result_paths, fit_titles, fit_descriptions;

	if (argc < 2) {
		cout << " usage : rootpwa_fitresult_inspector_gui <paths/to/directories/containing/rootfiles/with/fits>" << endl;
		cout << " note  : All root files in the same path will be loaded as one and the same fit result " << endl;
		cout << "         Several fit results in different paths can be compared to each other by providing the paths " << endl;
		return 0;
	} else {
		for (int ipath = 2; ipath <= argc; ipath++){
			string path = argv[ipath-1];
			if (DirExists(path)){
				cout << " including path to fit result #" << ipath-1 << ": " << endl;
				fit_result_paths.push_back(path);
				cout << " " << path << endl;
				stringstream title;
				title << ipath-1 << ": " << RemovePathExtension(path);
				fit_titles.push_back(title.str());
				stringstream description;
				description << " fit number " << ipath-1;
				fit_descriptions.push_back(description.str());
			} else {
				cout << " Error: not a valid directory: " << endl;
				cout << " " << path << endl;
				cout << " skipping! " << endl;
			}
		}
	}

	TApplication app("App", &argc, argv);
	// load the fit result inspector frame
	TrpwaPlotAmpsFrame* fitresultinspector = new TrpwaPlotAmpsFrame(fit_result_paths, fit_titles, fit_descriptions);
	if (!fitresultinspector) cout << " this will be not executed " << endl; // to prevent compiler warnings
	//app.Run();
	return 0;
}
