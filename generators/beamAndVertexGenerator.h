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

		/**
		 * Take the z position for the beam events from data
		 * instead of drawing them randomly according to the target definition
		 */
		void setTakeZpositionFromData(const bool takeZpositionFromData = true);

		/**
		 * Set the options to correct the beam momentum for the resolution
		 */
		void setMomentumResolutionCorrection(const double resolution, const double momentumPDFsigma, const double momentumPDFmean);

		std::ostream& print(std::ostream& out) const;

	  protected:

		virtual double getVertexZ(const rpwa::Target& target) const;

		std::string _beamFileName;

		TVector3 _vertex;
		TLorentzVector _beam;

		bool _readBeamfileSequentially;
		long _currentBeamfileEntry;

		/**
		 * Load beam file with momentum given as x,y,z coordinates
		 * @return true if loading was successful
		 */
		virtual bool loadBeamFileWithMomenta();
		/**
		 * Load beam file with momentum given as magnitude and x,y inclination
		 * @return true if loading was successful
		 */
		virtual bool loadBeamFileWithInclinations();
		virtual bool eventSimple(const rpwa::Target& target, const rpwa::Beam& beam);
		virtual bool eventFromMomenta(const rpwa::Target& target, const rpwa::Beam& beam);
		virtual bool eventFromMomentaCorrectMomentumResolution(const rpwa::Target& target, const rpwa::Beam& beam);
		virtual bool eventFromInclinations(const rpwa::Target& target, const rpwa::Beam& beam);
		/**
		 * Load next event from beam file into the member variables of this object
		 * @return true if loading was successful
		 */
		virtual bool loadNextEventFromBeamfile();

	  private:

		enum simulationMode {simpleSimulation, simulationFromMomenta, simulationFromInclinations, simulationFromMomentaCorrectMomentumResolution};

		simulationMode _simulationMode;

		TFile* _rootFile;
		TTree* _beamTree;

		// input data from beam file
		double _vertexX;                      // vertex x-position
		double _vertexY;                      // vertex y-position
		double _vertexZ;                      // vertex z-position
		// beam momentum as x,y,z coordinates
		double _beamMomentumX;                // beam-particle momentum in x-direction
		double _beamMomentumY;                // beam-particle momentum in y-direction
		double _beamMomentumZ;                // beam-particle momentum in z-direction
		// beam momentum as magnitude and x,y inclination
		double _beamdXdZ;                     // beam inclination in x-direction
		double _beamdYdZ;                     // beam inclination in y-direction
		double _beamMomentum;                 // absolute beam momentum
		// Gaussian uncertainties independent for each variable for momenta from momenta
		double _vertexXSigma;                 // uncertainty in vertex x-position
		double _vertexYSigma;                 // uncertainty in vertex y-position
		double _beamMomentumXSigma;           // uncertainty in beam-particle momentum in x-direction
		double _beamMomentumYSigma;           // uncertainty in beam-particle momentum in y-direction
		double _beamMomentumZSigma;           // uncertainty in beam-particle momentum in z-direction
		// correlated uncertainties for momenta from inclinations
		double _beamMomentumSigma;            // uncertainty in beam-particle momentum
		TMatrixDSym* _vertexCovarianceMatrix; // full covariance matrix of vertex parameters (x,y,z,dXdZ, dYdZ)

		bool _sigmasPresent;         // uncertainty is present in the beam file
		double _sigmaScalingFactor;  // scaling factor for all uncertainties
		bool _takeZpositionFromData; // take the z-position of the vertex from the beam file instad of generating one
		double _momentumResolution;  // assumption about momentum resolution for momentum resolution correction of the beam file
		double _momentumPDFSigma;    // assumption about std.-deviation of momentum distribution for momentum resolution correction of the beam file
		double _momentumPDFMean;     // assumption about mean of momentum distribution for momentum resolution correction of the beam file

	};

	inline std::ostream& operator<< (std::ostream& out, const beamAndVertexGenerator& beamAndVertGen)
	{
		return beamAndVertGen.print(out);
	}

}

#endif
/* @} **/
