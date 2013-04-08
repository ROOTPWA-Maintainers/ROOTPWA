/** @addtogroup generators
 * @{
 */

#ifndef TPRIMARYVERTEXGEN_HH
#define TPRIMARYVERTEXGEN_HH

#include <string>

#include <TVector3.h>
#include <TLorentzVector.h>

class TFile;
class TTree;


namespace rpwa {

	class generator;
	struct Target;

	class beamAndVertexGenerator {

	  public:

		beamAndVertexGenerator(std::string rootFileName,
		                       double massBeamParticle,
		                       Target target);

		~beamAndVertexGenerator();

		bool check();

		bool event(const rpwa::generator& generator);

		const TVector3& getVertex() { return _vertex; }

		const TLorentzVector& getBeam() { return _beam; }

		std::ostream& print(std::ostream& out);

	  private:

		std::string _rootFileName;
		TFile* _rootFile;
		TTree* _beamTree;
		TVector3 _vertex;
		TLorentzVector _beam;

		// pairs with [value, sigma]
		std::pair<double, double> _vertexX;
		std::pair<double, double> _vertexY;
		std::pair<double, double> _beamMomentumX;
		std::pair<double, double> _beamMomentumY;
		std::pair<double, double> _beamMomentumZ;

		double _massBeamParticle;
		bool _sigmasPresent;

	};

}

#endif
/* @} **/
