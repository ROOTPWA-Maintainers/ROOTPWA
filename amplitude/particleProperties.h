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
//      !NOTE! all potentially half-integer quantum numbers like I and
//             J are in units of hbar / 2 so that they are always integer
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
#include <map>
#include <vector>
#include <iostream>
#include <sstream>


namespace rpwa {


  class particleProperties {

  public:
		
    particleProperties();
    particleProperties(const particleProperties& partProp);
    particleProperties(const std::string&        partName,
		       const int                 isospin,
		       const int                 G,
		       const int                 J,
		       const int                 P,
		       const int                 C);
    
    virtual ~particleProperties();

    virtual particleProperties& operator =(const particleProperties& partProp);
    friend bool operator ==(const particleProperties& lhsProp,
			    const particleProperties& rhsProp);
    friend bool operator !=(const particleProperties& lhsProp,
			    const particleProperties& rhsProp) { return !(lhsProp == rhsProp); }

    friend bool operator ==(const particleProperties& lhsProp,
			    const std::pair<particleProperties, std::string>& rhsProp);
    friend bool operator !=(const particleProperties& lhsProp,
			    const std::pair<particleProperties, std::string>& rhsProp) {return !(lhsProp==rhsProp);}
    

    
    virtual std::string name() const { return _name; }  ///< returns particle name of the corresponding data table entry
    double      mass()        const { return _mass;        }  ///< returns particle mass
    double      width()       const { return _width;       }  ///< returns particle width
    int         baryonNmb()   const { return _baryonNmb;   }  ///< returns particle's baryon number
    int         isospin()     const { return _isospin;     }  ///< returns particle's isospin * 2 (!!!)
    int         strangeness() const { return _strangeness; }  ///< returns particle's strangeness
    int         charm()       const { return _charm;       }  ///< returns particle's charm
    int         beauty()      const { return _beauty;      }  ///< returns particle's beauty
    int         G()           const { return _G;           }  ///< returns particle's G-parity
    int         J()           const { return _J;           }  ///< returns particle's spin * 2 (!!!)
    int         P()           const { return _P;           }  ///< returns particle's parity
    int         C()           const { return _C;           }  ///< returns particle's C-parity

    bool fillFromDataTable(const std::string& name);

    void setName       (const std::string& name)        { _name        = name;        }  ///< sets particle name
    void setMass       (const double       mass)        { _mass        = mass;        }  ///< sets particle mass
    void setWidth      (const double       width)       { _width       = width;       }  ///< sets particle width
    void setBaryonNmb  (const int          baryonNmb)   { _baryonNmb   = baryonNmb;   }  ///< sets particle's baryon number
    void setIsospin    (const int          isospin)     { _isospin     = isospin;     }  ///< sets particle's isospin * 2 (!!!)
    void setStrangeness(const int          strangeness) { _strangeness = strangeness; }  ///< sets particle's strangeness
    void setCharm      (const int          charm)       { _charm       = charm;       }  ///< sets particle's charm
    void setBeauty     (const int          beauty)      { _beauty      = beauty;      }  ///< sets particle's beauty
    void setG          (const int          G)           { _G           = G;           }  ///< sets particle's G-parity
    void setJ          (const int          J)           { _J           = J;           }  ///< sets particle's spin * 2 (!!!)
    void setP          (const int          P)           { _P           = P;           }  ///< sets particle's parity
    void setC          (const int          C)           { _C           = C;           }  ///< sets particle's C-parity

    virtual std::ostream& print(std::ostream& out) const;  ///< prints particle data in human-readable form
    virtual std::ostream& dump (std::ostream& out) const;  ///< dumps particle properties into one text line as in data file

    bool read(std::istringstream& line);  ///< reads whitespace separated properties from single line

    static std::string chargeFromName(const std::string& name,
				      int&               charge);  ///< checks whether last characters of particle name contain charge information and sets charge accordingly; returns name with charge stripped off
    static std::string stripChargeFromName(const std::string& name);  ///< checks whether last characters of particle name contain charge information and returns name with charge stripped off

    static bool debug() { return _debug; }                             ///< returns debug flag
    static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


  private:

    std::string _name;         ///< full PDG name
    double      _mass;         ///< mass [GeV/c^]
    double      _width;        ///< total width [GeV/c^2]
    int         _baryonNmb;    ///< baryon number
    int         _isospin;      ///< isospin * 2 (!!!)
    int         _strangeness;  ///< strangeness
    int         _charm;        ///< charm
    int         _beauty;       ///< beauty
    int         _G;            ///< G-parity (0 = undefined)
    int         _J;            ///< spin * 2 (!!!)
    int         _P;            ///< parity (0 = undefined)
    int         _C;            ///< C-parity (0 = undefined)

    static bool _debug;  ///< if set to true, debug messages are printed

  };


  inline
  std::ostream&
  operator <<(std::ostream&             out,
	      const particleProperties& partProp) { return partProp.print(out); }


  inline
  std::istream&
  operator >>(std::istream&       in,
	      particleProperties& partProp)
  {
    std::string line;
    if (getline(in, line)) {
      // skip comments and empty lines
      while ((line == "") || (line[0] == '#')) {
       	if (partProp.debug())
       	  printInfo << "ignoring line '" << line << "'" << std::endl;
	if (!getline(in, line)) {
	  printWarn << "could not find valid particle entry before end of file" << std::endl;
	  return in;
	}
      }
      std::istringstream lineStream(line);
      partProp.read(lineStream);
    }
    return in;
  }


} // namespace rpwa
	

#endif  // PARTICLEPROPERTIES_H
