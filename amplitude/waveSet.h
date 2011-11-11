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

		bool parseWaveSetFile(const std::string& waveSetFileName);  ///< reads wave set parameters from libconfig file

		bool getDecayAmplitudeTrees(const std::vector<std::string>& ampFileNames);  ///< opens given list of files and creates an array of decay amplitude trees ordered according to wave set file; assumes that in the given set of files there are no two trees with the same name

		unsigned int nmbDecayAmps() const { return decayAmpTreeNames.size(); }  ///< returns number of decay amplitudes in wave set

		std::ostream& print(std::ostream& out) const;  ///< prints wave set parameters in human-readable form

		static bool debug() { return _debug; }                             ///< returns debug flag
		static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


	private:

		std::vector<std::string>                decayAmpTreeNames;   ///< array of tree names with decay amplitude values
		std::vector<std::pair<double, double> > decayAmpMassRanges;  ///< array with mass ranges in which decay amplitude should be used [MeV/c^2]
		std::vector<TTree*>                     decayAmpTrees;       //! ///< array with decay amplitude trees
		std::vector<waveDescription*>           decayAmpWaveDescs;   ///< array with decay amplitude wave descriptions

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
