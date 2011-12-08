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

#include "waveName.h"


namespace rpwa {


	class amplitudeName : public waveName {

	public:

		amplitudeName();
		amplitudeName(const waveName&        name,
		              const std::string&     incohAmpQnLabel = "");
		amplitudeName(const waveDescription& waveDesc,
		              const std::string&     incohAmpQnLabel = "");
#ifndef __CINT__
		amplitudeName(const isobarAmplitudePtr& amp,
		              const std::string&        incohAmpQnLabel = "");
#endif
		virtual ~amplitudeName();

		void clear();

		amplitudeName& operator =(const amplitudeName& ampName);

		friend bool operator ==(const amplitudeName& lhsName,
		                        const amplitudeName& rhsName);
		friend bool operator !=(const amplitudeName& lhsName,
		                        const amplitudeName& rhsName) { return not(lhsName == rhsName); }

		std::string waveName()    const { return waveName::fullName();                 }  ///< returns full wave name
		std::string fullName()    const { return waveName() + "," + incohAmpQnLabel(); }  ///< returns full amplitude name
		std::string operator ()() const { return fullName();                           }  ///< returns full wave name

		std::string incohAmpQnLabel() const { return _incohAmpQnLabel; }  ///< returns quantum numbers that are summed incoherently for either T or A

		void setIncohAmpQnLabel(const std::string& incohAmpQnLabel) { _incohAmpQnLabel = incohAmpQnLabel; }  ///< sets string for quantum numbers that are summed incoherently for either T or A

		static bool debug() { return _debug; }                             ///< returns debug flag
		static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


	protected:

		std::string _incohAmpQnLabel;  ///< quantum numbers that are summed incoherently for either T or A


	private:

		static bool _debug;  ///< if set to true, debug messages are printed


		ClassDef(amplitudeName,1)

	};


	inline
	bool
	operator ==(const amplitudeName& lhsName,
	            const amplitudeName& rhsName)
	{
		return (    (lhsName.waveName       () == rhsName.waveName       ())
		        and (lhsName.incohAmpQnLabel() == rhsName.incohAmpQnLabel()));
	}


	inline
	std::ostream&
	operator <<(std::ostream&        out,
	            const amplitudeName& name)
	{
		return out << name();
	}


}  // namespace rpwa


#endif  // AMPLITUDENAME_H
