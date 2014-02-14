#ifndef MASSDEPFIT_HH
#define MASSDEPFIT_HH

namespace rpwa {

	class massDepFit {

	public:

		massDepFit();
		~massDepFit() {};

		bool prepareMassLimits();

		bool readConfigInput(const libconfig::Setting* configInput);

		bool readInFile(const std::string& valTreeName   = "pwa",
		                const std::string& valBranchName = "fitResult_v2");

		const std::string& getInFileName() const { return _inFileName; }
		const std::vector<std::string>& getWaveNames() const { return _waveNames; }
		const std::vector<std::pair<double, double> >& getWaveMassLimits() const { return _waveMassLimits; }
		const std::vector<std::string>& getSysFileNames() const { return _sysFileNames; }
		bool getSysPlotting() const { return _sysPlotting; }

		static void setDebug(bool debug) { _debug = debug; }

	private:

		bool readConfigInputFitResults(const libconfig::Setting* configInputFitResults);
		bool readConfigInputWaves(const libconfig::Setting* configInputWaves);
		bool readConfigInputSystematics(const libconfig::Setting* configInputSystematics);

		bool checkFitResultMassBins(TTree* tree, rpwa::fitResult* fit, std::vector<size_t>& mapping) const;
		bool readFitResultMassBins(TTree* tree, rpwa::fitResult* fit);

		std::string _inFileName;

		double _massMax;
		double _massMin;
		double _massStep;
		std::vector<double> _massBinCenters;

		std::vector<std::string> _waveNames;
		std::vector<std::pair<double, double> > _waveMassLimits;
		std::vector<std::pair<size_t, size_t> > _waveMassBinLimits;

		std::vector<std::vector<std::pair<size_t, size_t> > > _wavePairMassBinLimits;

		bool _sysPlotting;
		std::vector<std::string> _sysFileNames;

		static bool _debug;

	};

} // end namespace

#endif
