//-----------------------------------------------------------
//
// Description:
//    (BW) Component of mass dependent fit
//      
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

#ifndef MASSDEPFITMODEL_HH
#define MASSDEPFITMODEL_HH

#include <complex>
#include <vector>

class TF1;

namespace rpwa {

	namespace massDepFit {

		class component;

		class model {

		public:

			model();
			~model() {}

			void add(rpwa::massDepFit::component* comp);

			TF1* getFsmdFunction() const { return _fsmdFunction; }
			void setFsmdFunction(TF1* fsmdFunction);

			bool init(const std::vector<std::string>& waveNames,
			          const std::vector<double>& massBinCenters,
			          const std::string& anchorWaveName,
			          const std::string& anchorComponentName);





    size_t n() const {return _comp.size();}
    unsigned int numPar() const {return _numpar;}
    
    void setPar(const double* par); // set parameters
    void getPar(double* par) const;       // return parameters 
    unsigned int nFreeFsmdPar() const {return _fsmdFreeParameters.size();}
    double getFreeFsmdPar(unsigned int i) const;
    void getFreeFsmdLimits(unsigned int i, double& lower, double& upper) const;


    const rpwa::massDepFit::component* operator[](size_t i) const {return _comp[i];}
    const std::vector<std::pair<size_t, size_t> >&
      getCompChannel(size_t idx) const { return _compChannel[idx]; }


			std::complex<double> productionAmplitude(const size_t idxWave,
			                                         const double mass,
			                                         const size_t idxMass = std::numeric_limits<size_t>::max()) const;
			double intensity(const size_t idxWave,
			                 const double mass,
			                 const size_t idxMass = std::numeric_limits<size_t>::max()) const;
			double phaseAbsolute(const size_t idxWave,
			                     const double mass,
			                     const size_t idxMass = std::numeric_limits<size_t>::max()) const;
			std::complex<double> spinDensityMatrix(const size_t idxWave,
			                                       const size_t jdxWave,
			                                       const double mass,
			                                       const size_t idxMass = std::numeric_limits<size_t>::max()) const;
			double phase(const size_t idxWave,
			             const size_t jdxWave,
			             const double mass,
			             const size_t idxMass = std::numeric_limits<size_t>::max()) const;

			double calcFsmd(const double mass,
			                const size_t idxMass = std::numeric_limits<size_t>::max()) const;

			std::ostream& print(std::ostream& out) const;

		private:

			bool initMapping(const std::string& anchorWaveName,
			                 const std::string& anchorComponentName);
			bool initFsmd(const std::vector<double>& massBinCenters);

			std::vector<std::string> _waveNames;

    std::vector<rpwa::massDepFit::component*> _comp;
    unsigned int _numpar;

			size_t _idxAnchorWave;
			size_t _idxAnchorComponent;
			size_t _idxAnchorChannel;

			TF1* _fsmdFunction;
			bool _fsmdFixed;
			std::vector<unsigned int> _fsmdFreeParameters;
			std::vector<double> _fsmdValues;

    // mapping for wave -> which components with which channel
    // wavelist in same order as given by wavelist
    std::vector<std::vector<std::pair<size_t, size_t> > > _compChannel;    

		};

	} // end namespace massDepFit

} // end namespace rpwa

inline std::ostream& operator<< (std::ostream& out, const rpwa::massDepFit::model& fitModel) {
	return fitModel.print(out);
}

#endif // MASSDEPFITMODEL_HH
