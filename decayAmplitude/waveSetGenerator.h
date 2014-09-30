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
//
// Description:
//      class that generates set of isobar decay topologies from a
//      template and according to user defined criteria
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef WAVESETGENERATOR_H
#define WAVESETGENERATOR_H


#include <boost/tuple/tuple.hpp>

#include "isobarDecayTopology.h"


namespace rpwa {


	class waveSetGenerator {

	public:

		// some typedefs for convenience
		typedef isobarDecayTopology::nodeDesc    nodeDesc;
		typedef isobarDecayTopology::adjIterator adjIterator;


		waveSetGenerator();
		virtual ~waveSetGenerator();

		bool setWaveSetParameters(const std::string& templateKeyFileName);  ///< constructs template topology and reads wave set parameters from key file

		// wave set parameter accessors
		// !note! isospin and angular momentum quantum numbers are in units of hbar / 2
		void setIsospinRange         (const int maxI = 2,
		                              const int minI = 0) { _isospinRange = boost::tuples::make_tuple(minI, maxI); }
		void setJRange               (const int maxJ = 6,
		                              const int minJ = 0) { _JRange       = boost::tuples::make_tuple(minJ, maxJ); }
		void setLRange               (const int maxL = 6,
		                              const int minL = 0) { _LRange       = boost::tuples::make_tuple(minL, maxL); }
		void setSRange               (const int maxS = 6,
		                              const int minS = 0) { _SRange       = boost::tuples::make_tuple(minS, maxS); }
		void setIsobarBlackList      (const std::vector<std::string>& isobarList)
		{ _isobarBlackList = isobarList; }
		void setIsobarWhiteList      (const std::vector<std::string>& isobarList)
		{ _isobarWhiteList = isobarList; }
		void setAllowSpinExotics     (const bool   flag     ) { _allowSpinExotics      = flag;  }
		void setRequireMinIsobarMass (const bool   flag     ) { _requireMinIsobarMass  = flag;  }
		void setForceDecayCheck      (const bool   flag     ) { _forceDecayCheck       = flag;  }
		void setIsobarMassWindowSigma(const double sigma = 1) { _isobarMassWindowSigma = sigma; }

		std::size_t generateWaveSet();  ///< generates wave set from template topology

		std::vector<isobarDecayTopology>&       waveSet()       { return _waveSet; }  ///< returns wave set
		const std::vector<isobarDecayTopology>& waveSet() const { return _waveSet; }  ///< returns wave set

		bool writeKeyFiles(const std::string& dirName                  = "",
		                   const bool         newKeyFileNameConvention = false);  ///< writes key files for wave set into given directory

		virtual void reset();  ///< resets parameters to default values and clears wave set

		virtual std::ostream& print(std::ostream& out) const;  ///< prints parameters

		static bool debug() { return _debug; }                             ///< returns debug flag
		static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


	private:

		const isobarDecayTopologyPtr createNewDecayTopology(const isobarDecayTopology&  parentDecay,
		                                                    const isobarDecayVertexPtr& parentVertex,
		                                                    const int                   L,
		                                                    const int                   S,
		                                                    const particleProperties&   isobar,
		                                                    const int                   parentCharge);

		std::set<std::size_t> findBoseSymDecays() const;  ///< finds decays wave set that are related by Bose symmetrization and returns indices of superfluous waves

		boost::tuples::tuple<int, int> _isospinRange;           ///< range of allowed isobar isospins
		boost::tuples::tuple<int, int> _JRange;                 ///< range of allowed isobar spins
		boost::tuples::tuple<int, int> _spinProjRange;          ///< range of allowed isobar spin projections
		int                            _reflectivity;           ///< if _useReflectivity is set, 0 means that both reflectivities are generated
		bool                           _useReflectivity;        ///< en/disables generation of waves with reflectivity
		bool                           _allowSpinExotics;       ///< flag that allows/forbids spin-exotics to be generated
		boost::tuples::tuple<int, int> _LRange;                 ///< range of allowed orbital angular momenta in isobar decays
		boost::tuples::tuple<int, int> _SRange;                 ///< range of allowed total intrinsic spins in isobar decays
		std::vector<std::string>       _isobarBlackList;        ///< list of particles not to be used as isobars
		std::vector<std::string>       _isobarWhiteList;        ///< list of particles to be used as isobars
		bool                           _requireMinIsobarMass;   ///< flag that en/disables cut on isobar mass
		double                         _isobarMassWindowSigma;  ///< defines width of isobar mass window in units of full widths of parent and daughter resonances
		bool                           _forceDecayCheck;        ///< enables strict decay checking. Particles without defined decays will be discarded

		isobarDecayTopologyPtr _templateTopo;  ///< template topology

		std::vector<isobarDecayTopology> _waveSet;  ///< generated wave set

		static bool _debug;  ///< if set to true, debug messages are printed

	};


	inline
	std::ostream&
	operator <<(std::ostream&           out,
	            const waveSetGenerator& gen)
	{
		return gen.print(out);
	}


} // namespace rpwa


#endif  // WAVESETGENERATOR_H
