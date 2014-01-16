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
//
// Description:
//      ROOT logon script that loads libraries required for plotting
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------


{
	//gSystem->Load("libGX11TTF.so");

	// load ROOTPWA libraries
	gSystem->Load("libMinuit2.so");
	gSystem->Load("libRootPwa.so");
	gSystem->AddIncludePath("-I$ROOTPWA/src");
	gSystem->AddIncludePath("-I$ROOTPWA/utilities");

	// basic tree routines
	gROOT->ProcessLine(".L loadFitResult.C+");
	//gROOT->ProcessLine(".L convertTFitResultTree.C+");

	// basic plotting routines
	gROOT->ProcessLine(".L plotIntensity.C+");
	//gROOT->ProcessLine(".L plotPhase.C+");
	//gROOT->ProcessLine(".L plotCoherence.C+");

	// summary plots
	gROOT->ProcessLine(".L plotAllIntensities.C+");
	//gROOT->ProcessLine(".L plotSpinTotals.C+");
	//gROOT->ProcessLine(".L plot4.C+");

	//gROOT->ProcessLine(".L loadFit.C+");
	//gROOT->ProcessLine(".L plotwaves.C+");

	// deprecated routines
	// gROOT->ProcessLine(".L loadFit.C+");
	// gROOT->ProcessLine(".L plotwaves.C+");


	gStyle->SetPadColor(0);
	gStyle->SetCanvasColor(0);

	// define grey-scale palette
	UInt_t   nmbColorEndPoints         = 2;
	Double_t red   [nmbColorEndPoints] = {1, 0  };
	Double_t green [nmbColorEndPoints] = {1, 0  };
	Double_t blue  [nmbColorEndPoints] = {1, 0.5};
	Double_t length[nmbColorEndPoints] = {0, 1  };
	Int_t    nmbColors                 = 50;
	TColor::CreateGradientColorTable(nmbColorEndPoints, length, red, green, blue, nmbColors);


	const Int_t NCont = 25;
	gStyle->SetNumberContours(NCont);

	// const UInt_t Number = 3;
	//Double_t Red[Number]    = {1.0, 0.2, 0.0};
	//Double_t Green[Number]  = {1.0, 0.2, 0.0};
	//Double_t Blue[Number]   = {1.0, 0.2, 0.0};
	//Double_t Length[Number] = {0.0, 0.8, 1.0};
	//TColor::CreateGradientColorTable(Number, Length, Red, Green, Blue, NCont);

	//const UInt_t NRGBs = 5;
	//Double_t stops[NRGBs] = {0.00, 0.34, 0.61, 0.84, 1.00};
	//Double_t red[NRGBs]   = {0.00, 0.00, 0.87, 1.00, 0.51};
	//Double_t green[NRGBs] = {0.00, 0.81, 1.00, 0.20, 0.00};
	//Double_t blue[NRGBs]  = {0.51, 1.00, 0.12, 0.00, 0.00};
	//TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);

}
