#ifndef GENERATOR_HH_
#define GENERATOR_HH_

#include "generatorParameters.hpp"
#include "generatorPickerFunctions.h"

namespace rpwa {

	class generator {

	  public:

		virtual ~generator() { delete _pickerFunction; };

		virtual void setBeam(const rpwa::Beam& beam) { _beam = beam; };
		virtual void setTarget(const rpwa::Target& target) { _target = target; };
		virtual void setTPrimeAndMassPicker(const rpwa::massAndTPrimePicker& pickerFunction) {
			(*_pickerFunction) = pickerFunction;
		}


	  protected:

		rpwa::Beam _beam;
		rpwa::Target _target;
		rpwa::massAndTPrimePicker* _pickerFunction;

	};

}

#endif
