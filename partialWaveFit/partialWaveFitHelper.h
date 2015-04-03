///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2015 Sebastian Uhl (TUM)
//
//    This file is part of ROOTPWA
//
//    ROOTPWA is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    ROOTPWA is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with ROOTPWA. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
//
// Description:
//      helper functions related to the partial-wave fit
//
//-------------------------------------------------------------------------


#include <iostream>


namespace rpwa {


	class fitResult;


	namespace partialWaveFitHelper {


		void extractWaveList(const rpwa::fitResult& fitResult, std::ostream& waveList);


	}  // namespace partialWaveFitHelper


}  // namespace rpwa
