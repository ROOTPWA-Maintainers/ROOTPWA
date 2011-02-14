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
//      ROOT logon script that sets up environment
//
//
// Author List:
//      Boris Grube    TUM            (original author)
//
//
//-----------------------------------------------------------


{
  // load ROOTPWA includes libraries
	gSystem->Load("libRootPwaGen.so");
  gSystem->AddIncludePath("-I$ROOTPWA/utilities");
  gSystem->AddIncludePath("-I$ROOTPWA/pwa2000/libpp");

  gStyle->SetPadColor(0);
  gStyle->SetCanvasColor(0);
}
