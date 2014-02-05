/** @addtogroup generators
 * @{
 */

#ifndef TPRIMARYVERTEXGEN_HH
#define TPRIMARYVERTEXGEN_HH

#include <string>

#include <boost/shared_ptr.hpp>

#include <TVector3.h>
#include <TLorentzVector.h>

class TFile;
class TTree;


namespace rpwa {

	class generator;
	struct Target;
	struct Beam;

	class beamAndVertexGenerator;
	typedef boost::shared_ptr<beamAndVertexGenerator> beamAndVertexGeneratorPtr;

	class beamAndVertexGenerator {

	  public:

		beamAndVertexGenerator();

		virtual ~beamAndVertexGenerator();

		virtual bool loadBeamFile(const std::string& beamFileName);
		virtual void setBeamfileSequentialReading(bool sequentialReading = true) { _readBeamfileSequentially = sequentialReading; }
		virtual void randomizeBeamfileStartingPosition();

		virtual bool check() const;

		virtual bool event(const rpwa::Target& target, const rpwa::Beam& beam);

		virtual const TVector3& getVertex() const { return _vertex; }

		virtual const TLorentzVector& getBeam() const { return _beam; }

		void setSigmaScalingFactor(const double& scalingFactor) { _sigmaScalingFactor = scalingFactor; }

		std::ostream& print(std::ostream& out) const;

	  protected:

		virtual double getVertexZ(const rpwa::Target& target) const;

		std::string _beamFileName;

		TVector3 _vertex;
		TLorentzVector _beam;

		bool _readBeamfileSequentially;
		long _currentBeamfileEntry;

	  private:

		bool _simpleSimulation;

		TFile* _rootFile;
		TTree* _beamTree;

		// pairs with [value, sigma]
		std::pair<double, double> _vertexX;
		std::pair<double, double> _vertexY;
		std::pair<double, double> _beamMomentumX;
		std::pair<double, double> _beamMomentumY;
		std::pair<double, double> _beamMomentumZ;

		bool _sigmasPresent;
		double _sigmaScalingFactor;

	};

	inline std::ostream& operator<< (std::ostream& out, const beamAndVertexGenerator& beamAndVertGen)
	{
		return beamAndVertGen.print(out);
	}

}

#endif
/* @} **/
