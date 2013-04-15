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

		beamAndVertexGenerator();

		virtual ~beamAndVertexGenerator();

		virtual bool loadBeamFile(const std::string& beamFileName);

		virtual bool check();

		virtual bool event(const rpwa::Target& target, const rpwa::Beam& beam);

		const TVector3& getVertex() { return _vertex; }

		const TLorentzVector& getBeam() { return _beam; }

		void setSigmaScalingFactor(const double& scalingFactor) { _sigmaScalingFactor = scalingFactor; }

		std::ostream& print(std::ostream& out);

	  protected:

		virtual double getVertexZ(const rpwa::Target& target) const;

	  private:

		bool _simpleSimulation;

		std::string _beamFileName;
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

		bool _sigmasPresent;
		double _sigmaScalingFactor;

	};

}

#endif
/* @} **/
