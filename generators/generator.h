#ifndef GENERATOR_HH_
#define GENERATOR_HH_

#include<iostream>
#include<vector>

#include<boost/shared_ptr.hpp>

#include "generatorParameters.hpp"
#include "generatorPickerFunctions.h"
#include "beamAndVertexGenerator.h"

class TVector3;

namespace rpwa {

	class particle;

	class generator {

	  public:

		generator()
			: _pickerFunction(boost::shared_ptr<rpwa::massAndTPrimePicker>()),
			  _beamAndVertexGenerator(boost::shared_ptr<rpwa::beamAndVertexGenerator>()) { }

		virtual ~generator() { }

		virtual unsigned int event() = 0;

		virtual const rpwa::particle& getGeneratedBeam() const { return _beam.particle; }
		virtual const rpwa::particle& getGeneratedRecoil() const { return _target.recoilParticle; }
		virtual const std::vector<rpwa::particle>& getGeneratedFinalState() const { return _decayProducts; }
		virtual const TVector3& getGeneratedVertex() const { return _vertex; }

		virtual void setBeam(const rpwa::Beam& beam) { _beam = beam; }
		virtual void setTarget(const rpwa::Target& target) { _target = target; }
		virtual void setTPrimeAndMassPicker(const boost::shared_ptr<rpwa::massAndTPrimePicker>& pickerFunction) {
			_pickerFunction = pickerFunction;
		}
		virtual void setPrimaryVertexGenerator(const boost::shared_ptr<rpwa::beamAndVertexGenerator>& beamAndVertexGenerator) {
			_beamAndVertexGenerator = beamAndVertexGenerator;
		}
		virtual void setDecayProducts(const std::vector<rpwa::particle>& particles) {
			_decayProducts = particles;
		}
		virtual void addDecayProduct(const rpwa::particle& particle) {
			_decayProducts.push_back(particle);
		}

		static std::ostream& convertEventToAscii(std::ostream& out,
		                                         const rpwa::particle& beam,
		                                         const std::vector<rpwa::particle>& finalState);

		static std::ostream& convertEventToComgeant(std::ostream& out,
		                                            const rpwa::particle& beam,
		                                            const rpwa::particle& recoil,
		                                            const TVector3& vertex,
		                                            const std::vector<rpwa::particle>& finalState,
		                                            bool writeBinary = false);


	  protected:

		rpwa::Beam _beam;
		rpwa::Target _target;
		boost::shared_ptr<rpwa::massAndTPrimePicker> _pickerFunction;
		boost::shared_ptr<rpwa::beamAndVertexGenerator> _beamAndVertexGenerator;
		std::vector<rpwa::particle> _decayProducts;
		TVector3 _vertex;

	};

}

#endif
