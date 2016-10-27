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
//      information for a single parameter of the resonance fit:
//      - name
//      - start value
//      - limits
//      - step size
//
//-------------------------------------------------------------------------


#ifndef RESONANCEFIT_PARAMETER_HH
#define RESONANCEFIT_PARAMETER_HH

#include <iostream>
#include <string>

namespace rpwa {

	namespace resonanceFit {

		class parameter {

		public:

			parameter();
			~parameter() {}

			const std::string& name() const { return _name; }
			void setName(const std::string& name) { _name = name; }

			double startValue() const { return _startValue; }
			void setStartValue(const double startValue) { _startValue = startValue; }
			double startError() const { return _startError; }
			void setStartError(const double startError) { _startError = startError; }

			bool fixed() const { return _fixed; }
			void setFixed(const bool fixed) { _fixed = fixed; }

			double limitLower() const { return _limitLower; }
			void setLimitLower(const double limitLower) { _limitLower = limitLower; }
			bool limitedLower() const { return _limitedLower; }
			void setLimitedLower(const bool limitedLower) { _limitedLower = limitedLower; }

			double limitUpper() const { return _limitUpper; }
			void setLimitUpper(const double limitUpper) { _limitUpper = limitUpper; }
			bool limitedUpper() const { return _limitedUpper; }
			void setLimitedUpper(const bool limitedUpper) { _limitedUpper = limitedUpper; }

			double step() const { return _step; }
			void setStep(const double step) { _step = step; }

			std::ostream& print(std::ostream& out = std::cout, const bool newLine = true) const;

		private:

			std::string _name;

			double _startValue;
			double _startError;

			bool _fixed;

			double _limitLower;
			bool _limitedLower;

			double _limitUpper;
			bool _limitedUpper;

			double _step;

		};

		std::ostream& operator<< (std::ostream& out, const rpwa::resonanceFit::parameter& parameter);

	} // end namespace resonanceFit

} // end namespace rpwa


inline
std::ostream&
rpwa::resonanceFit::operator<< (std::ostream& out, const rpwa::resonanceFit::parameter& parameter)
{
	return parameter.print(out, false);
}


#endif // RESONANCEFIT_PARAMETER_HH
