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
//      base class that encapsulates naming scheme for partial waves
//
//      in the general case the intensity is parameterized by
//
//      I = sum_j sum_k sum_l |sum_i T_i^{j l} A_i^{j k}|^2
//
//      where
//      T_i^{j l} - transition amplitude
//      A_i^{j k} - decay amplitude
//      i         - coherent sum index common for T and A (e.g. I^G J^PC M [isobar decay chain])
//      j         - incoherent sum index common for T and A (e.g. reflectivity)
//      k         - incoherent sum index for A only (e.g. FS-particle helicity)
//      l         - incoherent sum index for T only (e.g. rank)
//
//      so both T and A are defined by 3 sets of quantum numbers
//      the total amplitude name thus consists of 3 substrings
//
//      the spin-density matrix is given by
//
//      rho^j_{i, i'} = sum_l T_i^{j l} T*_{i'}^{j l}
//
//      so partial waves are characterized by i and j only
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef WAVENAME_H
#define WAVENAME_H


#include <string>

#include "TObject.h"

#include "waveDescription.h"
#ifndef __CINT__
#include "isobarAmplitude.h"
#endif


namespace rpwa {


	class waveName : public TObject {

	public:

		waveName();
		waveName(const waveDescription&    waveDesc);
#ifndef __CINT__
		waveName(const isobarAmplitudePtr& amp);
#endif
		virtual ~waveName();

		void clear();

		waveName& operator =(const waveName& name);

		friend bool operator ==(const waveName& lhsName,
		                        const waveName& rhsName);
		friend bool operator !=(const waveName& lhsName,
		                        const waveName& rhsName) { return not(lhsName == rhsName); }

		std::string fullName   () const { return cohQnLabel() + "," + incohQnLabel(); }  ///< returns full wave name
		std::string operator ()() const { return fullName();                          }  ///< returns full wave name

		std::string cohQnLabel  () const { return _cohQnLabel;   }  ///< returns quantum numbers that are summed coherently for T and A
		std::string incohQnLabel() const { return _incohQnLabel; }  ///< returns quantum numbers that are summed incoherently for T and A

		bool setQnLabel(const waveDescription&    waveDesc);  ///< sets strings for quantum numbers that are common for T and A from wave description
#ifndef __CINT__
		bool setQnLabel(const isobarAmplitudePtr& amp     );  ///< sets strings for quantum numbers that are common for T and A from decay amplitude
#endif

		static bool debug() { return _debug; }                             ///< returns debug flag
		static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


	protected:

		std::string _cohQnLabel;    ///< quantum numbers that are summed coherently for T and A
		std::string _incohQnLabel;  ///< quantum numbers that are summed incoherently for T and A


	private:

#ifndef __CINT__
		static std::string decayChain(const isobarDecayTopology& topo,
		                              const isobarDecayVertex&   currentVertex);
#endif

		static bool _debug;  ///< if set to true, debug messages are printed


		ClassDef(waveName,1)

	};


	inline
	bool
	operator ==(const waveName& lhsName,
	            const waveName& rhsName)
	{
		return (    (lhsName.cohQnLabel  () == rhsName.cohQnLabel  ())
		        and (lhsName.incohQnLabel() == rhsName.incohQnLabel()));
	}


	inline
	std::ostream&
	operator <<(std::ostream&   out,
	            const waveName& name)
	{
		return out << name();
	}


}  // namespace rpwa


#endif  // WAVENAME_H
