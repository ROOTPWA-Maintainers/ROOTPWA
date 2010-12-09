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
//      Data class reading extended PDG table
//      PDG table files available from 
//      http://pdg.lbl.gov/2008/html/computer_read.html
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

#ifndef TPDGENTRY_HH
#define TPDGENTRY_HH

// Base Class Headers ----------------
#include "TObject.h"

// Collaborating Class Headers -------
#include <ostream> // remove if you do not need streaming op
#include <istream>
#include "TString.h"

// Collaborating Class Declarations --


class TPDGEntry : public TObject {
public:

  enum PDG_AntiparticleFlag { B , F , blank };
  enum PDG_BaryonRank { NoBaryon, Poor , Fair , Likely , Certain };

  // Codes in pdg docu:         F             S                   D         R
  enum PDG_Status { NotDefined, FurtherMeson, NeedsConfirmation , Omitted , Established };


  // Constructors/Destructors ---------
  TPDGEntry(){}
  virtual ~TPDGEntry(){};

  // Operators
  friend bool operator== (const TPDGEntry& lhs, const TPDGEntry& rhs);
  friend std::istream& operator>> (std::istream& s, TPDGEntry& me);

  // Accessors -----------------------
  double mass() const {return _mass;}
  double mass_er() const {return 0.5*(_mass_en+_mass_ep);}
  void mass_errors(double& lower, double& upper) const 
  {upper=_mass_ep;lower=_mass_en;}
  double width() const {return _width;}
  double width_er() const {return 0.5*(_width_en+_width_ep);}
  void width_errors(double& lower, double& upper) const 
  {upper=_width_ep;lower=_width_en;}
  double I() const {return _I;}
  double J() const {return _J;}
  int G() const {return _G;}
  int P() const {return _P;}
  int C() const {return _C;}
  PDG_AntiparticleFlag aflag() const {return _aflag;}
  int pdgID() const {return _pdgID;}
  double q() const {return _q;}
  PDG_BaryonRank R() const {return _R;} 
  PDG_Status status() const {return _status;}
  TString name() const {return _name;}
  TString quarks() const {return _quark_content;} 
  bool isLightMeson() const;
  bool isExotic() const;

  // Modifiers -----------------------


  // Operations ----------------------
  virtual void Print(const Option_t* = 0) const;

private:

  // Private Data Members ------------
  double _mass;      // mass
  double _mass_en;   // negative error for mass
  double _mass_ep;   // positive error for mass
  double _width;
  double _width_en; 
  double _width_ep;
  double _I;            // Isospin
  int _G;               // G-Parity
  double _J;            // Total Spin
  int _P;               // Space Parity
  int _C;               // Charge Conjugation Parity
  PDG_AntiparticleFlag _aflag; // see pdg docu for details
  int _pdgID;           // Monte Carlo ID
  double _q;            // electric charge 
  PDG_BaryonRank _R;    // 4: Existance Certain ... 1: Evidence is poor
  PDG_Status _status;   // particle status
  TString _name;
  TString _quark_content; 

  // Private Methods -----------------
  // for printing
  TString Istr() const;
  TString Gstr() const;
  TString Jstr() const;
  TString Pstr() const;
  TString Cstr() const;
  TString Statstr() const;


public:
  ClassDef(TPDGEntry,1)

};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
