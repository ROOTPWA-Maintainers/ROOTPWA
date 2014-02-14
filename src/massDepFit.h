#ifndef MASSDEPFIT_HH
#define MASSDEPFIT_HH

#include <boost/numeric/ublas/matrix.hpp>

namespace rpwa {


	typedef boost::numeric::ublas::matrix<std::complex<double> > spinDensityMatrixType;
	typedef boost::numeric::ublas::matrix<double> complexCovarianceMatrixType;
	typedef boost::numeric::ublas::matrix<complexCovarianceMatrixType> spinDensityCovarianceMatrixType;


	class massDepFit {

	public:

		massDepFit();
		~massDepFit() {};

// FIXME: make private
		bool prepareMassLimits();

// FIXME: make private
		bool readConfigInput(const libconfig::Setting* configInput);

// FIXME: make private
		bool readInFile(const std::string& valTreeName   = "pwa",
		                const std::string& valBranchName = "fitResult_v2");

// FIXME: get rid
		const std::string& getInFileName() const { return _inFileName; }
// FIXME: get rid
		const std::vector<std::string>& getSysFileNames() const { return _sysFileNames; }
// FIXME: get rid
		bool getSysPlotting() const { return _sysPlotting; }
// FIXME: get rid
		const std::vector<double>& getMassBinCenters() const { return _massBinCenters; }
// FIXME: get rid
		const std::vector<std::string>& getWaveNames() const { return _waveNames; }
// FIXME: get rid
		const std::map<std::string, size_t>& getWaveIndices() const { return _waveIndices; }
// FIXME: get rid
		const std::vector<std::pair<double, double> >& getWaveMassLimits() const { return _waveMassLimits; }
// FIXME: get rid
		const std::vector<std::vector<std::pair<size_t, size_t> > >& getWavePairMassBinLimits() const { return _wavePairMassBinLimits; }
// FIXME: get rid
		const std::vector<spinDensityMatrixType>& getInSpinDensityMatrices() const { return _inSpinDensityMatrices; }
// FIXME: get rid
		const std::vector<spinDensityCovarianceMatrixType>& getInSpinDensityCovarianceMatrices() const { return _inSpinDensityCovarianceMatrices; }

		static void setDebug(bool debug) { _debug = debug; }

	private:

		bool readConfigInputFitResults(const libconfig::Setting* configInputFitResults);
		bool readConfigInputWaves(const libconfig::Setting* configInputWaves);
		bool readConfigInputSystematics(const libconfig::Setting* configInputSystematics);

		bool checkFitResultMassBins(TTree* tree,
		                            rpwa::fitResult* fit,
		                            std::vector<Long64_t>& mapping) const;
		bool readFitResultMassBins(TTree* tree,
		                           rpwa::fitResult* fit);
		bool readFitResultMatrices(TTree* tree,
		                           rpwa::fitResult* fit,
		                           const std::vector<Long64_t>& mapping,
		                           std::vector<spinDensityMatrixType>& spinDensityMatrices,
		                           std::vector<spinDensityCovarianceMatrixType>& spinDensityCovarianceMatrices) const;

		std::string _inFileName;

		bool _sysPlotting;
		std::vector<std::string> _sysFileNames;

		double _massMax;
		double _massMin;
		double _massStep;
		std::vector<double> _massBinCenters;

		std::vector<std::string> _waveNames;
		std::map<std::string, size_t> _waveIndices;
		std::vector<std::pair<double, double> > _waveMassLimits;
		std::vector<std::pair<size_t, size_t> > _waveMassBinLimits;

		std::vector<std::vector<std::pair<size_t, size_t> > > _wavePairMassBinLimits;

		std::vector<spinDensityMatrixType> _inSpinDensityMatrices;
		std::vector<spinDensityCovarianceMatrixType> _inSpinDensityCovarianceMatrices;

		static bool _debug;

	};


} // end namespace rpwa

#endif
