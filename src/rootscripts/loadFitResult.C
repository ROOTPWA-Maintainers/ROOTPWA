///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2009 Sebastian Neubert
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
//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Adds trees in file list specified by name pattern to chain.
//      The chain is created in case it does not exist.
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------


#include <iostream>
#include <string>

#include "TFileCollection.h"
#include "TChain.h"
#include "THashList.h"

#include "reportingUtils.hpp"


using namespace std;
using namespace rpwa;


TChain*
loadFitResult(const string& fileNamePattern,
              TChain*       chain    = 0,
              const string& treeName = "pwa")
{
	using rpwa::cout;
	// use TFileCollection to expand file name pattern into file list,
	// because TChain::Add() does support wildcards only for the root
	// files themselves (not in directories)
	printInfo << "constructing chain for '" << fileNamePattern << "' ..." << endl;
	TFileCollection fileList("fitresults", "fitresults");
	fileList.Add(fileNamePattern.c_str());
	cout << "    file list contains " << fileList.GetNFiles() << " files." << endl;
	if (!chain)
		chain = new TChain(treeName.c_str(), treeName.c_str());
	if (!chain->AddFileInfoList(fileList.GetList()))
		printWarn << "chain has problems reading file list." << endl;
	cout << "    chain '" << chain->GetName() << "' contains " << chain->GetEntries()
	     << " entries." << endl;
	return chain;
}
