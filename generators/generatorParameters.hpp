#ifndef GENERATORPARAMETERS_HH_
#define GENERATORPARAMETERS_HH_


#include "particle.h"


namespace rpwa {

	struct Beam {
		rpwa::particle particle;
		double momentum;
		double momentumSigma;
		double DxDz;
		double DxDzSigma;
		double DyDz;
		double DyDzSigma;

		std::ostream& print(std::ostream& out) const
		{
			out << "Beam parameters:" << std::endl;
			out << "    Particle name .......... " << particle.name() << std::endl;
			out << "    Momentum ............... " << momentum << std::endl;
			out << "    Momentum sigma ......... " << momentumSigma << std::endl;
			out << "    DxDz ................... " << DxDz << std::endl;
			out << "    DxDz sigma ............. " << DxDzSigma << std::endl;
			out << "    DyDz ................... " << DyDz << std::endl;
			out << "    DyDz sigma ............. " << DyDzSigma << std::endl;
			return out;
		}

	};


	inline std::ostream& operator<< (std::ostream& out, const Beam& beam)
	{
		return beam.print(out);
	}


	struct Target {
		rpwa::particle targetParticle;
		rpwa::particle recoilParticle;
		TVector3 position;
		double length;
		double radius;
		double interactionLength;

		std::ostream& print(std::ostream& out) const
		{
			out << "Target parameters:" << std::endl;
			out << "    Target particle name ... " << targetParticle.name() << std::endl;
			out << "    Recoil particle name ... " << recoilParticle.name() << std::endl;
			out << "    Target position ........ (" << position.X() << ", "
			                                        << position.Y() << ", "
			                                        << position.Z() << ")" << std::endl;
			out << "    Length ................. " << length << std::endl;
			out << "    Radius ................. " << radius << std::endl;
			out << "    Interaction length ..... " << interactionLength << std::endl;
			return out;
		}
	};


	inline std::ostream& operator<< (std::ostream& out, const Target& target)
	{
		return target.print(out);
	}


	struct FinalState {
		std::vector<rpwa::particle> particles;

		std::ostream& print(std::ostream& out) const
		{
			out << "Final state parameters:" << std::endl;
			out << "    " << particles.size() << " Final State Particles:" << std::endl;
			for(unsigned int i = 0; i < particles.size(); ++i) {
				out << "        - " << particles[i].name() << std::endl;
			}
			return out;
		}
	};


	inline std::ostream& operator<< (std::ostream& out, const FinalState& finalState)
	{
		return finalState.print(out);
	}


}

#endif
