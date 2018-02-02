/** @addtogroup generators
 * @{
 */

#ifndef TPRIMARYVERTEXGEN_HH
#define TPRIMARYVERTEXGEN_HH

#include <string>

#include <boost/shared_ptr.hpp>

#include <TVector3.h>
#include <TLorentzVector.h>
#include <TMatrixDSym.h>

class TFile;
class TTree;


namespace rpwa {

	class generator;
	struct Target;
	struct Beam;

	class beamAndVertexGenerator;
	typedef boost::shared_ptr<beamAndVertexGenerator> beamAndVertexGeneratorPtr;

	class beamAndVertexGenerator {
		enum simulationMode {simpleSimulation, simulationFromMomenta, simulationFromInclinations};
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

		virtual bool loadBeamFileWithMomenta();
		virtual bool loadBeamFileWithInclinations();
		virtual bool eventSimple(const rpwa::Target& target, const rpwa::Beam& beam);
		virtual bool eventFromMomenta(const rpwa::Target& target, const rpwa::Beam& beam);
		virtual bool eventFromInclinations(const rpwa::Target& target, const rpwa::Beam& beam);
		/**
		 * Load next event from beam file into the member variables of this object
		 * @return true if loading was successful
		 */
		virtual bool loadNextEventFromBeamfile();

	  private:

		simulationMode _simulationMode;

		TFile* _rootFile;
		TTree* _beamTree;

		double _vertexX;
		double _vertexY;
		double _vertexZ;
		double _beamMomentumX;
		double _beamMomentumY;
		double _beamMomentumZ;
		double _beamdXdZ;
		double _beamdYdZ;
		double _beamMomentum;
		double _vertexXSigma;
		double _vertexYSigma;
		double _beamMomentumXSigma;
		double _beamMomentumYSigma;
		double _beamMomentumZSigma;
		double _beamMomentumSigma;
		TMatrixDSym* _vertexCovarianceMatrix;

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
