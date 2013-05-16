///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
//
//    This file is part of rootpwa
//
//    rootpwa is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    rootpwa is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with rootpwa.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
// File and Version Information:
// $Rev:: 862                         $: revision of last commit
// $Author:: schmeing                 $: author of last commit
// $Date:: 2012-07-06 13:54:31 +0200 #$: date of last commit
//
// Description:
//      Main file for converter from CompassPWA files to root tree
//
//
// Author List:
//      Stephan Schmeing          TUM            (original author)
//
//
//-------------------------------------------------------------------------

#include <iostream>

#include "TFile.h"
#include "TTree.h"

#include "reportingUtils.hpp"

#include "CompassPwaFileLoader.h"

using namespace std;
using namespace rpwa;

int main(int argc, char *argv[]){ //[1]name of particleDataTable file, [2]name of root output file and [3] to [x] at least one destination of input files
//	CompassPwaFileLoader::SetDebug(true);
//	CompassPwaFileFitResults::SetDebug(true);
	if( argc > 3){
		TFile *ResultFile = new TFile( argv[2],"RECREATE","PWA fit results");
	
		if( ResultFile ){
			CompassPwaFileLoader FileLoader( argv[1] );

			printInfo << "Begin loading data from CompassPWA files\n";
			for( int i = 3; i < argc; ++i ){
				FileLoader.ReadFiles( argv[i] );
			}
			printInfo << "Done loading data from CompassPWA files\n";

			printInfo << "Merging CompassPWA data into root tree\n";
			TTree *data = FileLoader.Merge();
			if( data ){
				printInfo << "Done merging CompassPWA data into root tree\n";

				if( ResultFile->Write() ){
					printSucc << "Root Tree successfully written to " << argv[2] << '\n';
				}
				else{
					printErr << "File could not be written\n";
				}
			}
			else{
				printErr << "Merging failed\n";
			}

			ResultFile->Close();
			delete ResultFile;
		}
		else{
			printErr << "Could not create new TFile\n";
		}
	}
	else{
		printErr << "Not enough arguments (name of particleDataTable file, name of root output file and at least one destination of input files is needed)\n";
	}

	return 0;
}
