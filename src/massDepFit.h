#ifndef MASSDEPFIT_HH
#define MASSDEPFIT_HH

#include <boost/multi_array.hpp>

namespace rpwa {


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
		const boost::multi_array<std::pair<size_t, size_t>, 2>& getWavePairMassBinLimits() const { return _wavePairMassBinLimits; }
// FIXME: get rid
		const boost::multi_array<std::complex<double>, 3>& getInSpinDensityMatrices() const { return _inSpinDensityMatrices; }
// FIXME: get rid
		const boost::multi_array<double, 5>& getInSpinDensityCovarianceMatrices() const { return _inSpinDensityCovarianceMatrices; }
// FIXME: get rid
		const std::vector<std::vector<double> >& getInPhaseSpaceIntegrals() const { return _inPhaseSpaceIntegrals; }

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
		                           boost::multi_array<std::complex<double>, 3>& spinDensityMatrices,
		                           boost::multi_array<double, 5>& spinDensityCovarianceMatrices) const;
		bool readFitResultIntegrals(TTree* tree,
		                            rpwa::fitResult* fit,
		                            const std::vector<Long64_t>& mapping,
		                            std::vector<std::vector<double> >& phaseSpaceIntegrals) const;

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

		boost::multi_array<std::pair<size_t, size_t>, 2> _wavePairMassBinLimits;

		boost::multi_array<std::complex<double>, 3> _inSpinDensityMatrices;
		boost::multi_array<double, 5> _inSpinDensityCovarianceMatrices;
		std::vector<std::vector<double> > _inPhaseSpaceIntegrals;

		static bool _debug;

	};


} // end namespace rpwa

#endif
