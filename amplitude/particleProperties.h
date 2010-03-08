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
// $Rev::                             $: revision of last commit
// $Author::                          $: author of last commit
// $Date::                            $: date of last commit
//
// Description:
//      container class for particle properties
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef PARTICLEPROPERTIES_H
#define PARTICLEPROPERTIES_H


#include <string>
#include <iostream>
#include <sstream>


namespace rpwa {


  class particleProperties {

  public:
		
    particleProperties();
    particleProperties(const particleProperties& partProp);
    virtual ~particleProperties();

    particleProperties& operator =  (const particleProperties& partProp);
    friend bool         operator == (const particleProperties& lhsProp,
				     const particleProperties& rhsProp);
    friend bool         operator != (const particleProperties& lhsProp,
				     const particleProperties& rhsProp) { return !(lhsProp == rhsProp); }

    std::string name()      const { return _name;      }  ///< returns particle name
    double      mass()      const { return _mass;      }  ///< returns particle mass
    double      width()     const { return _width;     }  ///< returns particle width
    int         baryonNmb() const { return _baryonNmb; }  ///< returns particle's baryon number
    int         I()         const { return _I;         }  ///< returns particle's isospin
    int         S()         const { return _S;         }  ///< returns particle's strangeness
    int         C()         const { return _C;         }  ///< returns particle's charm
    int         B()         const { return _B;         }  ///< returns particle's beauty
    int         G()         const { return _G;         }  ///< returns particle's G-parity
    int         J()         const { return _J;         }  ///< returns particle's spin
    int         P()         const { return _P;         }  ///< returns particle's parity
    int         C()         const { return _C;         }  ///< returns particle's C-parity

    bool fillFromDataTable(const std::string& name);

    void setName     (const std::string& name)      { _name      = name;      }  ///< sets particle name
    void setMass     (const double       mass)      { _mass      = mass;      }  ///< sets particle mass
    void setWidth    (const double       width)     { _width     = width;     }  ///< sets particle width
    void setBaryonNmb(const int          baryonNmb) { _baryonNmb = baryonNmb; }  ///< sets particle's baryon number
    void setI        (const int          I)         { _I         = I;         }  ///< sets particle's isospin
    void setS        (const int          S)         { _S         = S;         }  ///< sets particle's strangeness
    void setC        (const int          C)         { _B         = B;         }  ///< sets particle's charm
    void setB        (const int          B)         { _C         = C;         }  ///< sets particle's beauty
    void setG        (const int          G)         { _G         = G;         }  ///< sets particle's G-parity
    void setJ        (const int          J)         { _J         = J;         }  ///< sets particle's spin
    void setP        (const int          P)         { _P         = P;         }  ///< sets particle's parity
    void setC        (const int          C)         { _C         = C;         }  ///< sets particle's C-parity

    void print(std::ostream& out) const;  ///< prints particle data in human-readable form
    void dump (std::ostream& out) const;  ///< dumps particle properties into one text line as in data file
    friend std::ostream& operator << (std::ostream&             out,
				      const particleProperties& partProp);

    bool read(std::istringstream& line);  ///< reads whitespace separated properties from single line
    friend std::istream& operator << (std::istream&       in,
				      particleProperties& partProp);

    static bool debug() const { return _debug; }
    static void setDebug(const bool debug = true) { _debug = debug; }


  private:

    std::string _name;       ///< full PDG name
    double      _mass;       ///< mass [GeV/c^]
    double      _width;      ///< total width [GeV/c^2]
    int         _baryonNmb;  ///< baryon number
    int         _I;          ///< isospin
    int         _S;          ///< strangeness
    int         _C;          ///< charm
    int         _B;          ///< beauty
    int         _G;          ///< G-parity (0 = undefined)
    int         _J;          ///< spin
    int         _P;          ///< parity (0 = undefined)
    int         _C;          ///< C-parity (0 = undefined)

    static bool _debug;  ///< if set to true, debug messages are printed

  };


} // namespace rpwa
	

#endif  // PARTICLEPROPERTIES_H
