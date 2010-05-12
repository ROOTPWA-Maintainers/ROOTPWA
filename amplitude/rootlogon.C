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
//    along with rootpwa. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
// File and Version Information:
// $Rev:: 200                         $: revision of last commit
// $Author:: bgrube                   $: author of last commit
// $Date:: 2010-04-09 11:32:36 +0200 #$: date of last commit
//
// Description:
//      ROOT logon macro that loads libraries needed by other ROOT macros
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


{
  gSystem->AddIncludePath("-I$ROOTPWA/pwa2000/libpp");
  gSystem->AddIncludePath("-I$ROOTPWA/src");

  gSystem->Load("libpp.so");
  gSystem->Load("libRootPwa.so");
  gSystem->Load("libRootPwaAmp.so");
}
