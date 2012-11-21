#ifndef GENERATORMANAGER_HH_
#define GENERATORMANAGER_HH_

#include<string>

#include "particleProperties.h"

class TVector3;

namespace rpwa {

	class generatorManager {

	  public:

		generatorManager();
		~generatorManager();

		bool readReactionFile(const std::string& fileName);

		static bool debug() { return _debug; };
		static void setDebug(bool debug = true) { _debug = debug; };

	  private:

		struct Beam {
			rpwa::particleProperties particle;
			double momentum;
			double momentumSigma;
			double DxDz;
			double DxDzSigma;
			double DyDz;
			double DyDzSigma;

			std::ostream& print(std::ostream& out) {
				out << "Beam parameters:" << std::endl;
				out << "    Particle name ..... " << particle.name() << std::endl;
				out << "    Momentum .......... " << momentum << std::endl;
				out << "    Momentum sigma .... " << momentumSigma << std::endl;
				out << "    DxDz .............. " << DxDz << std::endl;
				out << "    DxDz sigma ........ " << DxDzSigma << std::endl;
				out << "    DyDz .............. " << DyDz << std::endl;
				out << "    DyDz sigma ........ " << DyDzSigma << std::endl;
				return out;
			}

		};

		struct Target {
			rpwa::particleProperties particle;
			TVector3 position;
			double length;
			double radius;

			std::ostream& print(std::ostream& out) {
				out << "Target parameters:" << std::endl;
				out << "    Particle name ..... " << particle.name() << std::endl;
				out << "    Target position ... (" << position.X() << ", "
				                                   << position.Y() << ", "
				                                   << position.Z() << ")" << std::endl;
				out << "    Length ............ " << length << std::endl;
				out << "    Radius ............ " << radius << std::endl;
				return out;
			}
		};

		struct FinalState {
			double minimumMass;
			double maximumMass;
			double tPrimeSlope;
			double minimumTPrime;
			std::vector<rpwa::particleProperties> particles;

			std::ostream& print(std::ostream& out) {
				out << "Final state parameters:" << std::endl;
				out << "    Minimum Mass ...... " << minimumMass << std::endl;
				out << "    Maximum Mass ...... " << maximumMass << std::endl;
				out << "    t' Slope .......... " << tPrimeSlope << std::endl;
				out << "    Minimum t' ........ " << minimumTPrime << std::endl;
				out << "    " << particles.size() << " Final State Particles:" << std::endl;
				for(unsigned int i = 0; i < particles.size(); ++i) {
					out << "        - " << particles[i].name() << std::endl;
				}
				return out;
			}
		};

		Beam* _beam;
		Target* _target;
		FinalState* _finalState;

		static bool _debug;

	};

}

#endif

