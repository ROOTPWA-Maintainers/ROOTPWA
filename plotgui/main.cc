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
// $Rev::                             $: revision of last commit
// $Author::                          $: author of last commit
// $Date::                            $: date of last commit
//
// Description:
//      Main file for pwa gui
//
//
// Author List:
//      Stephan Schmeing          TUM            (original author)
//
//
//-------------------------------------------------------------------------

#include <iostream>

#include "FortranPwaFileObject.h"

using  namespace std;

int main(int argc, char *argv[]){

	//create instance of class Option
	FortranPwaFileObject FileObject;

//	FileObject.ReadFromFile( "/nfs/hicran/project/compass/analysis/sschmeing/PWA/work/integrals_new/PWAPhaseSpaceIntegrals_2500_2510_0000_1000.txt" );
//	FileObject.ReadFromFile( "/nfs/hicran/project/compass/analysis/sschmeing/PWA/work/integrals_new/PWANormIntegralsAcc_1670_1680_0100_1000.txt" );
	FileObject.ReadFromFile( "/nfs/hicran/project/compass/analysis/sschmeing/PWA/work/fits/fit_2008_W37_acc_53waves_new/PWAfitresults_fit_2008_W37_acc_53waves_new_1040_1060_0100_1000_0_fit3.txt" );
/*	if( 1 < argc ){
		FileObject.ReadFromFile( argv[1] );
	}
	if( 1 == argc ){

	}
*/
	return 0;
}

