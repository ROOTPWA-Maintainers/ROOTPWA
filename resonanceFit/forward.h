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


#ifndef RESONANCEFIT_FORWARD_HH
#define RESONANCEFIT_FORWARD_HH

namespace rpwa {

	namespace resonanceFit {

		class component;
		typedef std::shared_ptr<component> componentPtr;
		typedef std::shared_ptr<const component> componentConstPtr;

		class data;
		// 'data' should not be changed after construction, so this is not needed:
		// typedef std::shared_ptr<data> dataPtr;
		typedef std::shared_ptr<const data> dataConstPtr;

		class fsmd;
		typedef std::shared_ptr<fsmd> fsmdPtr;
		typedef std::shared_ptr<const fsmd> fsmdConstPtr;

		class function;
		typedef std::shared_ptr<function> functionPtr;
		typedef std::shared_ptr<const function> functionConstPtr;

		class information;
		// 'information' should not be changed after construction, so this is not needed:
		// typedef std::shared_ptr<information> informationPtr;
		typedef std::shared_ptr<const information> informationConstPtr;

		class model;
		// 'model' should not be changed after construction, so this is not needed:
		// typedef std::shared_ptr<model> modelPtr;
		typedef std::shared_ptr<const model> modelConstPtr;

	} // end namespace resonanceFit

} // end namespace rpwa

#endif // RESONANCEFIT_FORWARD_HH
