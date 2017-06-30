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

#include <memory>

namespace rpwa {

	namespace resonanceFit {

		class baseData;
#ifdef RESONANCEFIT_FORWARD_HH_FROM_PYTHON
		// 'baseData' should not be changed after construction, so this
		// is not needed from C++, however, Python does have no concept
		// of const-ness
		typedef std::shared_ptr<baseData> baseDataPtr;
#endif
		typedef std::shared_ptr<const baseData> baseDataConstPtr;

		class component;
#ifdef RESONANCEFIT_FORWARD_HH_FROM_PYTHON
		// 'component' should not be changed after construction, so this
		// is not needed from C++, however, Python does have no concept
		// of const-ness
		typedef std::shared_ptr<component> componentPtr;
#endif
		typedef std::shared_ptr<const component> componentConstPtr;

		class data;
#ifdef RESONANCEFIT_FORWARD_HH_FROM_PYTHON
		// 'data' should not be changed after construction, so this
		// is not needed from C++, however, Python does have no concept
		// of const-ness
		typedef std::shared_ptr<data> dataPtr;
#endif
		typedef std::shared_ptr<const data> dataConstPtr;

		class fsmd;
#ifdef RESONANCEFIT_FORWARD_HH_FROM_PYTHON
		// 'fsmd' should not be changed after construction, so this
		// is not needed from C++, however, Python does have no concept
		// of const-ness
		typedef std::shared_ptr<fsmd> fsmdPtr;
#endif
		typedef std::shared_ptr<const fsmd> fsmdConstPtr;

		class function;
#ifdef RESONANCEFIT_FORWARD_HH_FROM_PYTHON
		// 'function' should not be changed after construction, so this
		// is not needed from C++, however, Python does have no concept
		// of const-ness
		typedef std::shared_ptr<function> functionPtr;
#endif
		typedef std::shared_ptr<const function> functionConstPtr;

		class input;
#ifdef RESONANCEFIT_FORWARD_HH_FROM_PYTHON
		// 'input' should not be changed after construction, so this
		// is not needed from C++, however, Python does have no concept
		// of const-ness
		typedef std::shared_ptr<input> inputPtr;
#endif
		typedef std::shared_ptr<const input> inputConstPtr;

		class model;
#ifdef RESONANCEFIT_FORWARD_HH_FROM_PYTHON
		// 'model' should not be changed after construction, so this
		// is not needed from C++, however, Python does have no concept
		// of const-ness
		typedef std::shared_ptr<model> modelPtr;
#endif
		typedef std::shared_ptr<const model> modelConstPtr;

	} // end namespace resonanceFit

} // end namespace rpwa

#endif // RESONANCEFIT_FORWARD_HH
