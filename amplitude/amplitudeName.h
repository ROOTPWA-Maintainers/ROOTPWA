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
//      base class that encapsulates naming scheme for production and
//      decay amplitudes
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
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef AMPLITUDENAME_H
#define AMPLITUDENAME_H


#include <string>

#include "TObject.h"

#ifndef __CINT__
#include "isobarAmplitude.h"
#endif


namespace rpwa {


	class amplitudeName : public TObject {

	public:

		amplitudeName();
#ifndef __CINT__
		amplitudeName(const isobarAmplitudePtr& amp,
		              const std::string&        incohQnLabel = "");
#endif
		virtual ~amplitudeName();

		void clear();

		amplitudeName& operator =(const amplitudeName& ampName);

		friend bool operator ==(const amplitudeName& lhsName,
		                        const amplitudeName& rhsName);
		friend bool operator !=(const amplitudeName& lhsName,
		                        const amplitudeName& rhsName) { return not(lhsName == rhsName); }

		std::string fullName   () const;                        ///< returns full amplitude name
		std::string operator ()() const { return fullName(); }  ///< returns full amplitude name

		std::string commonCohQnLabel  () const { return _commonCohQnLabel;   }  ///< returns quantum numbers that are summed coherently for T and A
		std::string commonIncohQnLabel() const { return _commonIncohQnLabel; }  ///< returns quantum numbers that are summed incoherently for T and A
		std::string incohQnLabel      () const { return _incohQnLabel;       }  ///< returns quantum numbers that are summed incoherently for either T or A

#ifndef __CINT__
		void setCommonQn(const isobarAmplitudePtr& amp);  ///< sets strings for quantum numbers that are common for T and A from isobar decay amplitude
#endif
		void setIncohQn (const std::string& incohQnLabel) { _incohQnLabel = incohQnLabel; }  ///< sets string for quantum numbers that are summed incoherently for either T or A

		static bool debug() { return _debug; }                             ///< returns debug flag
		static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


	protected:

		std::string _commonCohQnLabel;    ///< quantum numbers that are summed coherently for T and A
		std::string _commonIncohQnLabel;  ///< quantum numbers that are summed incoherently for T and A
		std::string _incohQnLabel;        ///< quantum numbers that are summed incoherently for either T or A


	private:

		static std::string decayChain(const isobarDecayTopology& topo,
		                              const isobarDecayVertex&   currentVertex);

		static bool _debug;  ///< if set to true, debug messages are printed


		ClassDef(amplitudeName,1)

	};


	inline
	bool
	operator ==(const amplitudeName& lhsName,
	            const amplitudeName& rhsName)
	{
		return (    (lhsName.commonCohQnLabel  () == rhsName.commonCohQnLabel  ())
		        and (lhsName.commonIncohQnLabel() == rhsName.commonIncohQnLabel())
		        and (lhsName.incohQnLabel      () == rhsName.incohQnLabel      ()));
	}


}  // namespace rpwa


#endif  // AMPLITUDENAME_H
