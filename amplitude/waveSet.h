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
//      handles wave sets for fitting
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef WAVESET_H
#define WAVESET_H


#include <string>
#include <vector>

#include "TObject.h"

#ifndef __CINT__
#include "isobarAmplitude.h"
#endif
#include "waveName.h"


class TTree;
namespace rpwa {
	class waveDescription;
}


namespace rpwa {


	class waveSet : public TObject {

	public:

		waveSet();
		virtual ~waveSet();

		void clear();

		waveSet& operator =(const waveSet& set);

		unsigned int                                  nmbWaves          () const { return _decayAmpTreeNames.size(); }  ///< returns number of decay amplitudes in wave set
		const std::vector<std::string>&               decayAmpTreeNames () const { return _decayAmpTreeNames;        }  ///< returns array of decay amplitude tree names
		const std::vector<std::pair<double, double> > decayAmpMassRanges() const { return _decayAmpMassRanges;       }  ///< returns array with mass ranges in which decay amplitude should be used [MeV/c^2]
		const std::vector<std::string>&               decayAmpFileNames () const { return _decayAmpFileNames;        }  ///< returns array of decay amplitude file names
		const std::vector<TTree*>&                    decayAmpTrees     () const { return _decayAmpTrees;            }  ///< returns array of decay amplitude trees
		const std::vector<rpwa::waveDescription*>&    waveDescs         () const { return _waveDescs;                }  ///< returns array of wave descriptions
#ifndef __CINT__
		const std::vector<rpwa::isobarAmplitudePtr>&  decayAmps         () const { return _decayAmps;                }  ///< returns array of decay amplitudes
#endif
		const std::vector<rpwa::waveName>&            waveNames         () const { return _waveNames;                }  ///< returns array of wave names


		void setDecayAmpFileNames(const std::vector<std::string>& ampFileNames);  ///< sets list of amplitude file names


		bool parseWaveSetFile(const std::string& waveSetFileName);  ///< constructs wave set from libconfig file
		bool getDecayAmpTrees();    ///< opens list of ROOT files and creates an array of decay amplitude trees ordered as in wave set definition; assumes that in the given set of files there are no two trees with the same name
		bool getWaveDescs();        ///< reads wave descriptions from decay amplitude trees
		bool constructDecayAmps();  ///< constructs decay amplitude objects from wave descriptions
		bool constructWaveNames();  ///< constructs wave names from decay amplitude objects


		std::ostream& print(std::ostream& out) const;  ///< prints wave set parameters in human-readable form

		static bool debug() { return _debug; }                             ///< returns debug flag
		static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


	private:

		std::vector<std::string>                _decayAmpTreeNames;   ///< array of tree names with decay amplitude values
		std::vector<std::pair<double, double> > _decayAmpMassRanges;  ///< array with mass ranges in which decay amplitude should be used [MeV/c^2]
		std::vector<std::string>                _decayAmpFileNames;   ///< array with decay amplitude file names
		std::vector<TTree*>                     _decayAmpTrees;       //! ///< array with decay amplitude trees
		std::vector<waveDescription*>           _waveDescs;           ///< array with decay amplitude wave descriptions
#ifndef __CINT__
		std::vector<isobarAmplitudePtr>         _decayAmps;           ///< array of decay amplitudes
#endif
		std::vector<waveName>                   _waveNames;           ///< array of decay amplitude names

		static bool _debug;  ///< if set to true, debug messages are printed


#ifdef USE_STD_COMPLEX_TREE_LEAFS
		ClassDef(waveSet,1)
#endif

	};


	inline
	std::ostream&
	operator <<(std::ostream&  out,
	            const waveSet& set)
	{
		return set.print(out);
	}


}  // namespace rpwa


#endif  // WAVESET_H
