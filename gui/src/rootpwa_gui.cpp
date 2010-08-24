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

using namespace std;

int main(int argc, char *argv[]) {
	TApplication app("App", &argc, argv);
	TCanvas canvas;
	canvas.Draw();
	cout << "" << endl; // prints 
	app.Run();
	return 0;
}
