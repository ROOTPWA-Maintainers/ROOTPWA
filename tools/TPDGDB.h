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
// Description:
//      New PDG data base
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

#ifndef TPDGDB_HH
#define TPDGDB_HH

// Base Class Headers ----------------
#include "TObject.h"

// Collaborating Class Headers -------

// Collaborating Class Declarations --
class TString;
class TTree;
class TPDGEntry;

class TPDGDB : public TObject {
public:

  // Constructors/Destructors ---------
  TPDGDB():_tree(NULL){}
  virtual ~TPDGDB(){}

  // Operators
  

  // Accessors -----------------------
  TTree* getTree() { return _tree;}
  void Draw(char* com, char* sel, char* opt, int n, int s);
  

  // Modifiers -----------------------


  // Operations ----------------------
  unsigned int read(const TString& filename, int num=0);
  

private:

  // Private Data Members ------------
  TTree* _tree; // data

  // Private Methods -----------------

public:
    ClassDef(TPDGDB,1)

};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
