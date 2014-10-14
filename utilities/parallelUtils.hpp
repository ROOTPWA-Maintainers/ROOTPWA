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
//      functions for simple often used parallel computations
//
//
// Author List:
//      Armin Gensler          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef PARALLELUTILS_HPP
#define PARALLELUTILS_HPP

#include <complex>
#include <iostream>
#include <vector>
#include <boost/date_time/posix_time/posix_time.hpp>

#include "Typedefs.hpp"
#include "reportingUtils.hpp"
#include "reportingUtilsRoot.hpp"


namespace rpwa {

	inline
	void
	parallelLorentzRotationInvert(std::vector<LorentzRotation>& lzRot)
	{
		// !! EVENT PARALLEL LOOP
		boost::posix_time::ptime timeBefore = boost::posix_time::microsec_clock::local_time();
		#pragma omp parallel for
		for(unsigned int i = 0; i < lzRot.size(); ++i)
			lzRot[i].Invert();
		boost::posix_time::ptime timeAfter = boost::posix_time::microsec_clock::local_time();
		uint64_t timeDiff = (timeAfter - timeBefore).total_milliseconds();
		std::cout << "EPL: parallelLorentzRotationInvert timediff = " << timeDiff << std::endl;
	}

}  // namespace rpwa


#endif  // PARALLELUTILS_HPP
