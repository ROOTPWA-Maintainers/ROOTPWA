///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2016 Sebastian Uhl (TUM)
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
//    along with ROOTPWA.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
//
// Description:
//      forward declaration of classes use throughout the resonance fit
//      * including shared_ptr definitions
//
//-------------------------------------------------------------------------


#ifndef MASSDEPFITFORWARD_HH
#define MASSDEPFITFORWARD_HH

namespace rpwa {

	namespace massDepFit {

		class component;
		typedef std::shared_ptr<component> componentPtr;
		typedef std::shared_ptr<const component> componentConstPtr;

		class fsmd;
		typedef std::shared_ptr<fsmd> fsmdPtr;
		typedef std::shared_ptr<const fsmd> fsmdConstPtr;

	} // end namespace massDepFit

} // end namespace rpwa

#endif // MASSDEPFITFORWARD_HH
