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

#ifndef MASSDEPFITCOMPONENTS_HH
#define MASSDEPFITCOMPONENTS_HH

#include <map>

#include <boost/multi_array.hpp>

#include <Math/Interpolator.h>

namespace libconfig {
	class Setting;
}
namespace ROOT {
	namespace Math {
		class Minimizer;
	}
}

namespace rpwa {

	namespace massDepFit {

		class channel {

		public:

			channel(const std::string& waveName,
			        const size_t nrBins,
			        const std::vector<std::complex<double> >& couplings,
			        const std::vector<double>& massBinCenters,
			        const boost::multi_array<double, 2>& phaseSpace);
			channel(const rpwa::massDepFit::channel& ch);
			~channel();

			const std::string& getWaveName() const { return _waveName; }

			bool isAnchor() const { return _anchor; }
			void setAnchor(const bool anchor) { _anchor = anchor; }

			size_t getNrBins() const { return _nrBins; }

			std::complex<double> getCoupling(const size_t idxBin) const { return _couplings[idxBin]; }
			const std::vector<std::complex<double> >& getCouplings() const { return _couplings; }
			void setCouplings(const std::vector<std::complex<double> >& couplings) { _couplings = couplings; }

			double getPhaseSpace(const size_t idxBin,
			                     const double mass,
			                     const size_t idxMass = std::numeric_limits<size_t>::max()) const;

			std::complex<double> getCouplingPhaseSpace(const size_t idxBin,
			                                           const double mass,
			                                           const size_t idxMass = std::numeric_limits<size_t>::max()) const;

		private:

			std::string _waveName;

			bool _anchor;

			size_t _nrBins;

			std::vector<std::complex<double> > _couplings;

			std::vector<double> _massBinCenters;
			boost::multi_array<double, 2> _phaseSpace;
			std::vector<ROOT::Math::Interpolator*> _interpolator;

		};

		class component {

		public:

			component(const std::string& name,
			          const size_t nrParameters);
			virtual ~component() {};

			const std::string& getName() const { return _name; }

			virtual bool init(const libconfig::Setting* configComponent,
			                  const size_t nrBins,
			                  const std::vector<double>& massBinCenters,
			                  const std::map<std::string, size_t>& waveIndices,
			                  const boost::multi_array<double, 3>& phaseSpaceIntegrals,
			                  const bool debug);

			virtual bool update(const libconfig::Setting* configComponent,
			                    const ROOT::Math::Minimizer* minimizer,
			                    const bool debug) const;

			const size_t getNrChannels() const { return _channels.size(); }
			const std::vector<channel>& getChannels() const { return _channels; }
			const channel& getChannel(const size_t i) const { return _channels[i]; }
			void setChannelAnchor(const size_t i, const bool anchor) { _channels[i].setAnchor(anchor); }
			const std::string& getChannelWaveName(const size_t i) const { return _channels[i].getWaveName(); }

			size_t getCouplings(double* par) const;
			size_t setCouplings(const double* par);

			size_t getNrParameters() const { return _nrParameters; }
			virtual size_t getParameters(double* par) const;
			virtual size_t setParameters(const double* par);

			virtual double getParameter(const size_t idx) const { return _parameters[idx]; }
			virtual bool getParameterFixed(const size_t idx) const { return _parametersFixed[idx]; }
			virtual double getParameterLimitLower(const size_t idx) const { return _parametersLimitLower[idx]; }
			virtual bool getParameterLimitedLower(const size_t idx) const { return _parametersLimitedLower[idx]; }
			virtual double getParameterLimitUpper(const size_t idx) const { return _parametersLimitUpper[idx]; }
			virtual bool getParameterLimitedUpper(const size_t idx) const { return _parametersLimitedUpper[idx]; }
			virtual const std::string& getParameterName(const size_t idx) const { return _parametersName[idx]; }
			virtual double getParameterStep(const size_t idx) const { return _parametersStep[idx]; }

			virtual std::complex<double> val(const size_t idxBin,
			                                 const double m) const = 0;

			virtual std::ostream& print(std::ostream& out) const;

		private:

			const std::string _name;

			std::vector<channel> _channels;

			const size_t _nrParameters;

		protected:

			std::vector<double> _parameters;
			std::vector<bool> _parametersFixed;
			std::vector<double> _parametersLimitLower;
			std::vector<bool> _parametersLimitedLower;
			std::vector<double> _parametersLimitUpper;
			std::vector<bool> _parametersLimitedUpper;
			std::vector<std::string> _parametersName;
			std::vector<double> _parametersStep;

		};

		class fixedWidthBreitWigner : public component {

		public:

			fixedWidthBreitWigner(const std::string& name);

			virtual bool init(const libconfig::Setting* configComponent,
			                  const size_t nrBins,
			                  const std::vector<double>& massBinCenters,
			                  const std::map<std::string, size_t>& waveIndices,
			                  const boost::multi_array<double, 3>& phaseSpaceIntegrals,
			                  const bool debug);

			virtual std::complex<double> val(const size_t idxBin,
			                                 const double m) const;

			virtual std::ostream& print(std::ostream& out) const;

		};

		class relativisticBreitWigner : public component {

		public:

			relativisticBreitWigner(const std::string& name);

			virtual bool init(const libconfig::Setting* configComponent,
			                  const size_t nrBins,
			                  const std::vector<double>& massBinCenters,
			                  const std::map<std::string, size_t>& waveIndices,
			                  const boost::multi_array<double, 3>& phaseSpaceIntegrals,
			                  const bool debug);

			virtual std::complex<double> val(const size_t idxBin,
			                                 const double m) const;

			virtual std::ostream& print(std::ostream& out) const;

		private:

			std::vector<double> _ratio;
			std::vector<int> _l;
			std::vector<double> _m1;
			std::vector<double> _m2;

		};

		class parameterizationA1Bowler : public component {

		public:

			parameterizationA1Bowler(const std::string& name);

			virtual bool init(const libconfig::Setting* configComponent,
			                  const size_t nrBins,
			                  const std::vector<double>& massBinCenters,
			                  const std::map<std::string, size_t>& waveIndices,
			                  const boost::multi_array<double, 3>& phaseSpaceIntegrals,
			                  const bool debug);

			virtual std::complex<double> val(const size_t idxBin,
			                                 const double m) const;

			virtual std::ostream& print(std::ostream& out) const;

		};

		class exponentialBackground : public component {

		public:

			exponentialBackground(const std::string& name);

			virtual bool init(const libconfig::Setting* configComponent,
			                  const size_t nrBins,
			                  const std::vector<double>& massBinCenters,
			                  const std::map<std::string, size_t>& waveIndices,
			                  const boost::multi_array<double, 3>& phaseSpaceIntegrals,
			                  const bool debug);

			virtual std::complex<double> val(const size_t idxBin,
			                                 const double m) const;

			virtual std::ostream& print(std::ostream& out) const;

		private:

			double _m1;
			double _m2;

		};

		class tPrimeDependentBackground : public component {

		public:

			tPrimeDependentBackground(const std::string& name);

			virtual bool setTPrimeMeans(const std::vector<double> tPrimeMeans);

			virtual bool init(const libconfig::Setting* configComponent,
			                  const size_t nrBins,
			                  const std::vector<double>& massBinCenters,
			                  const std::map<std::string, size_t>& waveIndices,
			                  const boost::multi_array<double, 3>& phaseSpaceIntegrals,
			                  const bool debug);

			virtual std::complex<double> val(const size_t idxBin,
			                                 const double m) const;

			virtual std::ostream& print(std::ostream& out) const;

		private:

			std::vector<double> _tPrimeMeans;
			double _m1;
			double _m2;

		};

	} // end namespace massDepFit

} // end namespace rpwa


inline
double
rpwa::massDepFit::channel::getPhaseSpace(const size_t idxBin,
                                         const double mass,
                                         const size_t idxMass) const
{
	if(idxMass != std::numeric_limits<size_t>::max()) {
		return _phaseSpace[idxBin][idxMass];
	}

	return _interpolator[idxBin]->Eval(mass);
}


inline
std::complex<double>
rpwa::massDepFit::channel::getCouplingPhaseSpace(const size_t idxBin,
                                                 const double mass,
                                                 const size_t idxMass) const
{
	return getCoupling(idxBin) * getPhaseSpace(idxBin, mass, idxMass);
}


inline
std::ostream&
operator<< (std::ostream& out, const rpwa::massDepFit::component& component)
{
	return component.print(out);
}


#endif // MASSDEPFITCOMPONENTS_HH
