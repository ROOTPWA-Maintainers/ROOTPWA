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

#include <iostream>
#include <map>

#include <boost/multi_array.hpp>

#include <Math/Interpolator.h>

#include "massDepFitParameters.h"

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
			        const std::vector<double>& massBinCenters,
			        const boost::multi_array<double, 2>& phaseSpace);
			channel(const rpwa::massDepFit::channel& ch);
			~channel();

			const std::string& getWaveName() const { return _waveName; }

			bool isAnchor() const { return _anchor; }
			void setAnchor(const bool anchor) { _anchor = anchor; }

			size_t getNrBins() const { return _nrBins; }

			double getPhaseSpace(const size_t idxBin,
			                     const double mass,
			                     const size_t idxMass = std::numeric_limits<size_t>::max()) const;

		private:

			std::string _waveName;

			bool _anchor;

			size_t _nrBins;

			std::vector<double> _massBinCenters;
			boost::multi_array<double, 2> _phaseSpace;
			std::vector<ROOT::Math::Interpolator*> _interpolator;

		};

		class component {

		public:

			component(const size_t id,
			          const std::string& name,
			          const size_t nrParameters);
			virtual ~component() {};

			size_t getId() const { return _id; }
			const std::string& getName() const { return _name; }

			virtual bool init(const libconfig::Setting* configComponent,
			                  rpwa::massDepFit::parameters& fitParameters,
			                  const size_t nrBins,
			                  const std::vector<double>& massBinCenters,
			                  const std::map<std::string, size_t>& waveIndices,
			                  const boost::multi_array<double, 3>& phaseSpaceIntegrals,
			                  const bool useBranchings,
			                  const bool debug);

			virtual bool update(const libconfig::Setting* configComponent,
			                    const rpwa::massDepFit::parameters& fitParameters,
			                    const ROOT::Math::Minimizer* minimizer,
			                    const bool useBranchings,
			                    const bool debug) const;

			size_t getNrChannels() const { return _channels.size(); }
			const std::vector<channel>& getChannels() const { return _channels; }
			const channel& getChannel(const size_t i) const { return _channels[i]; }
			void setChannelAnchor(const size_t i, const bool anchor) { _channels[i].setAnchor(anchor); }
			const std::string& getChannelWaveName(const size_t i) const { return _channels[i].getWaveName(); }
			size_t getChannelIdxCoupling(const size_t i) const { return _channelsCoupling[i]; }
			size_t getChannelIdxBranching(const size_t i) const { return _channelsBranching[i]; }

			size_t getNrCouplings() const { return _nrCouplings; }
			size_t importCouplings(const double* par,
			                       rpwa::massDepFit::parameters& fitParameters);

			size_t getNrBranchings() const { return _nrBranchings; }
			size_t importBranchings(const double* par,
			                        rpwa::massDepFit::parameters& fitParameters);

			size_t getNrParameters() const { return _nrParameters; }
			virtual size_t importParameters(const double* par,
			                                rpwa::massDepFit::parameters& fitParameters);

			virtual bool getParameterFixed(const size_t idxParameter) const { return _parametersFixed[idxParameter]; }
			virtual double getParameterLimitLower(const size_t idxParameter) const { return _parametersLimitLower[idxParameter]; }
			virtual bool getParameterLimitedLower(const size_t idxParameter) const { return _parametersLimitedLower[idxParameter]; }
			virtual double getParameterLimitUpper(const size_t idxParameter) const { return _parametersLimitUpper[idxParameter]; }
			virtual bool getParameterLimitedUpper(const size_t idxParameter) const { return _parametersLimitedUpper[idxParameter]; }
			virtual const std::string& getParameterName(const size_t idxParameter) const { return _parametersName[idxParameter]; }
			virtual double getParameterStep(const size_t idxParameter) const { return _parametersStep[idxParameter]; }

			virtual std::complex<double> val(const rpwa::massDepFit::parameters& fitParameters,
			                                 const size_t idxBin,
			                                 const double m) const = 0;

			std::complex<double> getCouplingPhaseSpace(const rpwa::massDepFit::parameters& fitParameters,
			                                           const size_t idxChannel,
			                                           const size_t idxBin,
			                                           const double mass,
			                                           const size_t idxMass = std::numeric_limits<size_t>::max()) const;

			virtual std::ostream& print(std::ostream& out = std::cout) const;

		private:

			const size_t _id;
			const std::string _name;

			std::vector<channel> _channels;
			std::vector<size_t> _channelsCoupling;
			std::vector<size_t> _channelsFromCoupling;
			std::vector<size_t> _channelsBranching;

			const size_t _nrParameters;
			size_t _nrCouplings;
			size_t _nrBranchings;

		protected:

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

			fixedWidthBreitWigner(const size_t id,
			                      const std::string& name);

			virtual bool init(const libconfig::Setting* configComponent,
			                  rpwa::massDepFit::parameters& fitParameters,
			                  const size_t nrBins,
			                  const std::vector<double>& massBinCenters,
			                  const std::map<std::string, size_t>& waveIndices,
			                  const boost::multi_array<double, 3>& phaseSpaceIntegrals,
			                  const bool useBranchings,
			                  const bool debug);

			virtual std::complex<double> val(const rpwa::massDepFit::parameters& fitParameters,
			                                 const size_t idxBin,
			                                 const double m) const;

			virtual std::ostream& print(std::ostream& out = std::cout) const;

		};

		class dynamicWidthBreitWigner : public component {

		public:

			dynamicWidthBreitWigner(const size_t id,
			                        const std::string& name);

			virtual bool init(const libconfig::Setting* configComponent,
			                  rpwa::massDepFit::parameters& fitParameters,
			                  const size_t nrBins,
			                  const std::vector<double>& massBinCenters,
			                  const std::map<std::string, size_t>& waveIndices,
			                  const boost::multi_array<double, 3>& phaseSpaceIntegrals,
			                  const bool useBranchings,
			                  const bool debug);

			virtual std::complex<double> val(const rpwa::massDepFit::parameters& fitParameters,
			                                 const size_t idxBin,
			                                 const double m) const;

			virtual std::ostream& print(std::ostream& out = std::cout) const;

		private:

			std::vector<double> _ratio;
			std::vector<int> _l;
			std::vector<double> _m1;
			std::vector<double> _m2;

		};

		class parameterizationA1Bowler : public component {

		public:

			parameterizationA1Bowler(const size_t id,
			                         const std::string& name);
			virtual ~parameterizationA1Bowler();

			virtual bool init(const libconfig::Setting* configComponent,
			                  rpwa::massDepFit::parameters& fitParameters,
			                  const size_t nrBins,
			                  const std::vector<double>& massBinCenters,
			                  const std::map<std::string, size_t>& waveIndices,
			                  const boost::multi_array<double, 3>& phaseSpaceIntegrals,
			                  const bool useBranchings,
			                  const bool debug);

			virtual std::complex<double> val(const rpwa::massDepFit::parameters& fitParameters,
			                                 const size_t idxBin,
			                                 const double m) const;

			virtual std::ostream& print(std::ostream& out = std::cout) const;

		private:

			ROOT::Math::Interpolator* _interpolator;

		};

		class exponentialBackground : public component {

		public:

			exponentialBackground(const size_t id,
			                      const std::string& name);

			virtual bool init(const libconfig::Setting* configComponent,
			                  rpwa::massDepFit::parameters& fitParameters,
			                  const size_t nrBins,
			                  const std::vector<double>& massBinCenters,
			                  const std::map<std::string, size_t>& waveIndices,
			                  const boost::multi_array<double, 3>& phaseSpaceIntegrals,
			                  const bool useBranchings,
			                  const bool debug);

			virtual std::complex<double> val(const rpwa::massDepFit::parameters& fitParameters,
			                                 const size_t idxBin,
			                                 const double m) const;

			virtual std::ostream& print(std::ostream& out = std::cout) const;

		private:

			double _m1;
			double _m2;

		};

		class tPrimeDependentBackground : public component {

		public:

			tPrimeDependentBackground(const size_t id,
			                          const std::string& name);

			virtual bool setTPrimeMeans(const std::vector<double> tPrimeMeans);

			virtual bool init(const libconfig::Setting* configComponent,
			                  rpwa::massDepFit::parameters& fitParameters,
			                  const size_t nrBins,
			                  const std::vector<double>& massBinCenters,
			                  const std::map<std::string, size_t>& waveIndices,
			                  const boost::multi_array<double, 3>& phaseSpaceIntegrals,
			                  const bool useBranchings,
			                  const bool debug);

			virtual std::complex<double> val(const rpwa::massDepFit::parameters& fitParameters,
			                                 const size_t idxBin,
			                                 const double m) const;

			virtual std::ostream& print(std::ostream& out = std::cout) const;

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
rpwa::massDepFit::component::getCouplingPhaseSpace(const rpwa::massDepFit::parameters& fitParameters,
                                                   const size_t idxChannel,
                                                   const size_t idxBin,
                                                   const double mass,
                                                   const size_t idxMass) const
{
	const std::complex<double> coupling = fitParameters.getCoupling(_id, _channelsCoupling[idxChannel], idxBin);
	const std::complex<double> branching = fitParameters.getBranching(_id, _channelsBranching[idxChannel]);

	const channel& channelPhaseSpace = _channels[idxChannel];
	const double phaseSpace = channelPhaseSpace.getPhaseSpace(idxBin, mass, idxMass);

	return coupling * branching * phaseSpace;
}


inline
std::ostream&
operator<< (std::ostream& out, const rpwa::massDepFit::component& component)
{
	return component.print(out);
}


#endif // MASSDEPFITCOMPONENTS_HH
