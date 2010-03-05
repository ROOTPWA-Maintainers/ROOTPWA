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
//      Loads basic libraries and functions needed for analyzing fit
//      results
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------


{
  gSystem->Load("libGX11TTF.so");
  gSystem->Load("librootpwa.so");

  gSystem->AddIncludePath("-I$ROOTPWA/src");

  gROOT->ProcessLine(".L loadFitResult.C+");
  gROOT->ProcessLine(".L plotIntensity.C+");
  gROOT->ProcessLine(".L plotAllIntensities.C+");
  gROOT->ProcessLine(".L plotSpinTotals.C+");
  gROOT->ProcessLine(".L plotPhase.C+");
  gROOT->ProcessLine(".L plotCoherence.C+");
  gROOT->ProcessLine(".L plot4.C+");
  gROOT->ProcessLine(".L convertTFitResultTree.C+");

  gROOT->ProcessLine(".L loadFit.C+");
  gROOT->ProcessLine(".L plotwaves.C+");

  gStyle->SetPadColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPalette(1);
}
