//============================================================================
// Name        : rootpwa_gui.cpp
// Author      : Promme@web.de, jasinski@kph.uni-mainz.de
// Version     :
// Copyright   : gnu
// Description : gui for rootpwa usage based on root gui packages
//============================================================================

#include <iostream>
#include "TApplication.h"
#include "TCanvas.h"
#include "TrpwaMainFrame.h"

using namespace std;

int main(int argc, char *argv[]) {
	cout << " __________________________________________" << endl;
	cout << endl;
	cout << "                rootpwa GUI                " << endl;
	cout << "                 vers.0.11                 " << endl;
	cout << endl;
	cout << " __________________________________________" << endl;
	TApplication app("App", &argc, argv);
	// load the main frame and let the show begin
	TrpwaMainFrame::Instance();
	app.Run();
	return 0;
}
