///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010-2012 Sebastian Neubert (TUM)
//    Copyright 2014,2015 Sebastian Uhl (TUM)
//
//    This file is part of ROOTPWA
//
//    ROOTPWA is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    ROOTPWA is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with ROOTPWA.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
//
// Description:
//      components and channels for the resonance fit
//      - channels are used to keep track of the coupling of a component to
//        one particular decay channel in the partial-wave fit
//      - several types of components have been implemented
//        * Breit-Wigner functions with various types of (dynamic) widths
//        * two background parameterizations
//
//-------------------------------------------------------------------------


#ifndef RESONANCEFIT_COMPONENTS_HH
#define RESONANCEFIT_COMPONENTS_HH

#include <iostream>
#include <map>

#include <boost/multi_array.hpp>

#include <Math/Interpolator.h>

#include "cache.h"
#include "parameters.h"

namespace YAML {
	class Emitter;
	class Node;
}

namespace rpwa {

	namespace resonanceFit {

		class channel {

		public:

			channel(const size_t waveIdx,
			        const std::string& waveName,
			        const std::vector<size_t>& bins,
			        const std::vector<size_t>& nrMassBins,
			        const boost::multi_array<double, 2>& massBinCenters,
			        const boost::multi_array<double, 2>& phaseSpace);

			size_t getWaveIdx() const { return _waveIdx; }
			const std::string& getWaveName() const { return _waveName; }

			bool isAnchor(const size_t idxBin) const { return _anchor[idxBin]; }
			void setAnchor(const size_t idxBin) { _anchor[idxBin] = true; }

			const std::vector<size_t>& getBins() const { return _bins; }

			double getPhaseSpace(const size_t idxBin,
			                     const double mass,
			                     const size_t idxMass = std::numeric_limits<size_t>::max()) const;

		private:

			const size_t _waveIdx;
			const std::string _waveName;

			std::vector<bool> _anchor;

			std::vector<size_t> _bins;

			boost::multi_array<double, 2> _phaseSpace;
			std::vector<std::shared_ptr<const ROOT::Math::Interpolator> > _interpolator;

		};

		class component {

		public:

			component(const size_t id,
			          const std::string& name,
			          const std::string& type,
			          const size_t nrParameters);
			virtual ~component() {}

			size_t getId() const { return _id; }
			const std::string& getName() const { return _name; }
			const std::string& getType() const { return _type; }

			virtual bool init(const YAML::Node& configComponent,
			                  rpwa::resonanceFit::parameters& fitParameters,
			                  rpwa::resonanceFit::parameters& fitParametersError,
			                  const std::vector<size_t>& nrMassBins,
			                  const boost::multi_array<double, 2>& massBinCenters,
			                  const std::map<std::string, size_t>& waveIndices,
			                  const std::map<std::string, std::vector<size_t> >& waveBins,
			                  const boost::multi_array<double, 3>& phaseSpaceIntegrals,
			                  const bool useBranchings,
			                  const bool debug);
			virtual bool readDecayChannel(const YAML::Node& decayChannel,
			                              const size_t idxDecayChannel,
			                              const bool debug);

			virtual bool write(YAML::Emitter& yamlOutput,
			                   const rpwa::resonanceFit::parameters& fitParameters,
			                   const rpwa::resonanceFit::parameters& fitParametersError,
			                   const bool debug) const;
			virtual bool writeDecayChannel(YAML::Emitter& yamlOutput,
			                               const size_t idxDecayChannel,
			                               const bool debug) const;

			size_t getNrChannels() const { return _channels.size(); }
			const std::vector<channel>& getChannels() const { return _channels; }
			const channel& getChannel(const size_t i) const { return _channels[i]; }
			const channel& getChannelFromBranchingIdx(const size_t i) const { return _channels[_channelsFromBranching[i]]; }
			const channel& getChannelFromCouplingIdx(const size_t i) const { return _channels[_channelsFromCoupling[i]]; }
			void setChannelAnchor(const std::vector<size_t>& channels);

			size_t getNrCouplings() const { return _nrCouplings; }
			size_t importCouplings(const double* par,
			                       rpwa::resonanceFit::parameters& fitParameters,
			                       rpwa::resonanceFit::cache& cache);

			size_t getNrBranchings() const { return _nrBranchings; }
			size_t importBranchings(const double* par,
			                        rpwa::resonanceFit::parameters& fitParameters,
			                        rpwa::resonanceFit::cache& cache);

			size_t getNrParameters() const { return _nrParameters; }
			size_t importParameters(const double* par,
			                        rpwa::resonanceFit::parameters& fitParameters,
			                        rpwa::resonanceFit::cache& cache);

			bool isBranchingFixed(const size_t idxBranching) const { return _branchingsFixed[idxBranching]; }

			double getParameterStart(const size_t idxParameter) const { return _parametersStart[idxParameter]; }
			double getParameterError(const size_t idxParameter) const { return _parametersError[idxParameter]; }

			bool getParameterFixed(const size_t idxParameter) const { return _parametersFixed[idxParameter]; }
			double getParameterLimitLower(const size_t idxParameter) const { return _parametersLimitLower[idxParameter]; }
			bool getParameterLimitedLower(const size_t idxParameter) const { return _parametersLimitedLower[idxParameter]; }
			double getParameterLimitUpper(const size_t idxParameter) const { return _parametersLimitUpper[idxParameter]; }
			bool getParameterLimitedUpper(const size_t idxParameter) const { return _parametersLimitedUpper[idxParameter]; }
			const std::string& getParameterName(const size_t idxParameter) const { return _parametersName[idxParameter]; }
			double getParameterStep(const size_t idxParameter) const { return _parametersStep[idxParameter]; }

			std::complex<double> val(const rpwa::resonanceFit::parameters& fitParameters,
			                         rpwa::resonanceFit::cache& cache,
			                         const size_t idxBin,
			                         const double m,
			                         const size_t idxMass = std::numeric_limits<size_t>::max()) const;

			std::complex<double> getCouplingPhaseSpace(const rpwa::resonanceFit::parameters& fitParameters,
			                                           rpwa::resonanceFit::cache& cache,
			                                           const size_t idxChannel,
			                                           const size_t idxBin,
			                                           const double mass,
			                                           const size_t idxMass = std::numeric_limits<size_t>::max()) const;

			virtual std::ostream& print(std::ostream& out = std::cout) const;

		private:

			virtual std::complex<double> val(const rpwa::resonanceFit::parameters& fitParameters,
			                                 const size_t idxBin,
			                                 const double m) const = 0;

			const size_t _id;
			const std::string _name;
			const std::string _type;

			std::vector<channel> _channels;
			std::vector<size_t> _channelsCoupling;
			std::vector<size_t> _channelsFromCoupling;
			std::vector<size_t> _channelsBranching;
			std::vector<size_t> _channelsFromBranching;

			const size_t _nrParameters;
			size_t _nrCouplings;
			size_t _nrBranchings;

			std::vector<bool> _branchingsFixed;

		protected:

			std::vector<std::vector<size_t> > _binsEqualValues;

			std::vector<double> _parametersStart;
			std::vector<double> _parametersError;

			std::vector<bool> _parametersFixed;
			std::vector<double> _parametersLimitLower;
			std::vector<bool> _parametersLimitedLower;
			std::vector<double> _parametersLimitUpper;
			std::vector<bool> _parametersLimitedUpper;
			std::vector<std::string> _parametersName;
			std::vector<double> _parametersStep;

		};

		std::ostream& operator<< (std::ostream& out, const rpwa::resonanceFit::component& component);

		class fixedWidthBreitWigner : public component {

		public:

			fixedWidthBreitWigner(const size_t id,
			                      const std::string& name);

			virtual bool init(const YAML::Node& configComponent,
			                  rpwa::resonanceFit::parameters& fitParameters,
			                  rpwa::resonanceFit::parameters& fitParametersError,
			                  const std::vector<size_t>& nrMassBins,
			                  const boost::multi_array<double, 2>& massBinCenters,
			                  const std::map<std::string, size_t>& waveIndices,
			                  const std::map<std::string, std::vector<size_t> >& waveBins,
			                  const boost::multi_array<double, 3>& phaseSpaceIntegrals,
			                  const bool useBranchings,
			                  const bool debug);

			virtual bool write(YAML::Emitter& yamlOutput,
			                   const rpwa::resonanceFit::parameters& fitParameters,
			                   const rpwa::resonanceFit::parameters& fitParametersError,
			                   const bool debug) const;

			virtual std::ostream& print(std::ostream& out = std::cout) const;

		private:

			virtual std::complex<double> val(const rpwa::resonanceFit::parameters& fitParameters,
			                                 const size_t idxBin,
			                                 const double m) const;

		};

		class dynamicWidthBreitWigner : public component {

		public:

			dynamicWidthBreitWigner(const size_t id,
			                        const std::string& name);

			virtual bool init(const YAML::Node& configComponent,
			                  rpwa::resonanceFit::parameters& fitParameters,
			                  rpwa::resonanceFit::parameters& fitParametersError,
			                  const std::vector<size_t>& nrMassBins,
			                  const boost::multi_array<double, 2>& massBinCenters,
			                  const std::map<std::string, size_t>& waveIndices,
			                  const std::map<std::string, std::vector<size_t> >& waveBins,
			                  const boost::multi_array<double, 3>& phaseSpaceIntegrals,
			                  const bool useBranchings,
			                  const bool debug);
			virtual bool readDecayChannel(const YAML::Node& decayChannel,
			                              const size_t idxDecayChannel,
			                              const bool debug);

			virtual bool write(YAML::Emitter& yamlOutput,
			                   const rpwa::resonanceFit::parameters& fitParameters,
			                   const rpwa::resonanceFit::parameters& fitParametersError,
			                   const bool debug) const;
			virtual bool writeDecayChannel(YAML::Emitter& yamlOutput,
			                               const size_t idxDecayChannel,
			                               const bool debug) const;

			virtual std::ostream& print(std::ostream& out = std::cout) const;

		private:

			virtual std::complex<double> val(const rpwa::resonanceFit::parameters& fitParameters,
			                                 const size_t idxBin,
			                                 const double m) const;

			std::vector<double> _ratio;
			std::vector<int> _l;
			std::vector<double> _m1;
			std::vector<double> _m2;

		};

		class integralWidthBreitWigner : public component {

		public:

			integralWidthBreitWigner(const size_t id,
			                         const std::string& name);

			virtual bool init(const YAML::Node& configComponent,
			                  rpwa::resonanceFit::parameters& fitParameters,
			                  rpwa::resonanceFit::parameters& fitParametersError,
			                  const std::vector<size_t>& nrMassBins,
			                  const boost::multi_array<double, 2>& massBinCenters,
			                  const std::map<std::string, size_t>& waveIndices,
			                  const std::map<std::string, std::vector<size_t> >& waveBins,
			                  const boost::multi_array<double, 3>& phaseSpaceIntegrals,
			                  const bool useBranchings,
			                  const bool debug);
			virtual bool readDecayChannel(const YAML::Node& decayChannel,
			                              const size_t idxDecayChannel,
			                              const bool debug);

			virtual bool write(YAML::Emitter& yamlOutput,
			                   const rpwa::resonanceFit::parameters& fitParameters,
			                   const rpwa::resonanceFit::parameters& fitParametersError,
			                   const bool debug) const;
			virtual bool writeDecayChannel(YAML::Emitter& yamlOutput,
			                               const size_t idxDecayChannel,
			                               const bool debug) const;

			virtual std::ostream& print(std::ostream& out = std::cout) const;

		private:

			virtual std::complex<double> val(const rpwa::resonanceFit::parameters& fitParameters,
			                                 const size_t idxBin,
			                                 const double m) const;

			std::vector<std::vector<double> > _masses;
			std::vector<std::vector<double> > _values;
			std::vector<std::shared_ptr<const ROOT::Math::Interpolator> > _interpolator;
			std::vector<double> _ratio;

		};

		class constantBackground : public component {

		public:

			constantBackground(const size_t id,
			                   const std::string& name);

			virtual bool init(const YAML::Node& configComponent,
			                  rpwa::resonanceFit::parameters& fitParameters,
			                  rpwa::resonanceFit::parameters& fitParametersError,
			                  const std::vector<size_t>& nrMassBins,
			                  const boost::multi_array<double, 2>& massBinCenters,
			                  const std::map<std::string, size_t>& waveIndices,
			                  const std::map<std::string, std::vector<size_t> >& waveBins,
			                  const boost::multi_array<double, 3>& phaseSpaceIntegrals,
			                  const bool useBranchings,
			                  const bool debug);

			virtual bool write(YAML::Emitter& yamlOutput,
			                   const rpwa::resonanceFit::parameters& fitParameters,
			                   const rpwa::resonanceFit::parameters& fitParametersError,
			                   const bool debug) const;

			virtual std::ostream& print(std::ostream& out = std::cout) const;

		private:

			virtual std::complex<double> val(const rpwa::resonanceFit::parameters& fitParameters,
			                                 const size_t idxBin,
			                                 const double m) const;

		};

		class exponentialBackground : public component {

		public:

			exponentialBackground(const size_t id,
			                      const std::string& name);

			virtual bool init(const YAML::Node& configComponent,
			                  rpwa::resonanceFit::parameters& fitParameters,
			                  rpwa::resonanceFit::parameters& fitParametersError,
			                  const std::vector<size_t>& nrMassBins,
			                  const boost::multi_array<double, 2>& massBinCenters,
			                  const std::map<std::string, size_t>& waveIndices,
			                  const std::map<std::string, std::vector<size_t> >& waveBins,
			                  const boost::multi_array<double, 3>& phaseSpaceIntegrals,
			                  const bool useBranchings,
			                  const bool debug);

			virtual bool write(YAML::Emitter& yamlOutput,
			                   const rpwa::resonanceFit::parameters& fitParameters,
			                   const rpwa::resonanceFit::parameters& fitParametersError,
			                   const bool debug) const;

			virtual std::ostream& print(std::ostream& out = std::cout) const;

		private:

			virtual std::complex<double> val(const rpwa::resonanceFit::parameters& fitParameters,
			                                 const size_t idxBin,
			                                 const double m) const;

			int _l;
			double _m1;
			double _m2;
			double _exponent;

			double _norm;

		};

		class tPrimeDependentBackground : public component {

		public:

			tPrimeDependentBackground(const size_t id,
			                          const std::string& name);

			bool setTPrimeMeans(const std::vector<double> tPrimeMeans);

			virtual bool init(const YAML::Node& configComponent,
			                  rpwa::resonanceFit::parameters& fitParameters,
			                  rpwa::resonanceFit::parameters& fitParametersError,
			                  const std::vector<size_t>& nrMassBins,
			                  const boost::multi_array<double, 2>& massBinCenters,
			                  const std::map<std::string, size_t>& waveIndices,
			                  const std::map<std::string, std::vector<size_t> >& waveBins,
			                  const boost::multi_array<double, 3>& phaseSpaceIntegrals,
			                  const bool useBranchings,
			                  const bool debug);

			virtual bool write(YAML::Emitter& yamlOutput,
			                   const rpwa::resonanceFit::parameters& fitParameters,
			                   const rpwa::resonanceFit::parameters& fitParametersError,
			                   const bool debug) const;

			virtual std::ostream& print(std::ostream& out = std::cout) const;

		private:

			virtual std::complex<double> val(const rpwa::resonanceFit::parameters& fitParameters,
			                                 const size_t idxBin,
			                                 const double m) const;

			std::vector<double> _tPrimeMeans;
			int _l;
			double _m1;
			double _m2;
			double _exponent;

			double _norm;

		};

		class exponentialBackgroundIntegral : public component {

		public:

			exponentialBackgroundIntegral(const size_t id,
			                              const std::string& name);

			virtual bool init(const YAML::Node& configComponent,
			                  rpwa::resonanceFit::parameters& fitParameters,
			                  rpwa::resonanceFit::parameters& fitParametersError,
			                  const std::vector<size_t>& nrMassBins,
			                  const boost::multi_array<double, 2>& massBinCenters,
			                  const std::map<std::string, size_t>& waveIndices,
			                  const std::map<std::string, std::vector<size_t> >& waveBins,
			                  const boost::multi_array<double, 3>& phaseSpaceIntegrals,
			                  const bool useBranchings,
			                  const bool debug);

			virtual bool write(YAML::Emitter& yamlOutput,
			                   const rpwa::resonanceFit::parameters& fitParameters,
			                   const rpwa::resonanceFit::parameters& fitParametersError,
			                   const bool debug) const;

			virtual std::ostream& print(std::ostream& out = std::cout) const;

		private:

			virtual std::complex<double> val(const rpwa::resonanceFit::parameters& fitParameters,
			                                 const size_t idxBin,
			                                 const double m) const;

			std::vector<double> _masses;
			std::vector<double> _values;
			std::shared_ptr<const ROOT::Math::Interpolator> _interpolator;
			double _exponent;

			double _norm;

		};

		class tPrimeDependentBackgroundIntegral : public component {

		public:

			tPrimeDependentBackgroundIntegral(const size_t id,
			                                  const std::string& name);

			bool setTPrimeMeans(const std::vector<double> tPrimeMeans);

			virtual bool init(const YAML::Node& configComponent,
			                  rpwa::resonanceFit::parameters& fitParameters,
			                  rpwa::resonanceFit::parameters& fitParametersError,
			                  const std::vector<size_t>& nrMassBins,
			                  const boost::multi_array<double, 2>& massBinCenters,
			                  const std::map<std::string, size_t>& waveIndices,
			                  const std::map<std::string, std::vector<size_t> >& waveBins,
			                  const boost::multi_array<double, 3>& phaseSpaceIntegrals,
			                  const bool useBranchings,
			                  const bool debug);

			virtual bool write(YAML::Emitter& yamlOutput,
			                   const rpwa::resonanceFit::parameters& fitParameters,
			                   const rpwa::resonanceFit::parameters& fitParametersError,
			                   const bool debug) const;

			virtual std::ostream& print(std::ostream& out = std::cout) const;

		private:

			virtual std::complex<double> val(const rpwa::resonanceFit::parameters& fitParameters,
			                                 const size_t idxBin,
			                                 const double m) const;

			std::vector<double> _tPrimeMeans;
			std::vector<double> _masses;
			std::vector<double> _values;
			std::shared_ptr<const ROOT::Math::Interpolator> _interpolator;
			double _exponent;

			double _norm;

		};

	} // end namespace resonanceFit

} // end namespace rpwa


inline
double
rpwa::resonanceFit::channel::getPhaseSpace(const size_t idxBin,
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
rpwa::resonanceFit::component::getCouplingPhaseSpace(const rpwa::resonanceFit::parameters& fitParameters,
                                                     rpwa::resonanceFit::cache& cache,
                                                     const size_t idxChannel,
                                                     const size_t idxBin,
                                                     const double mass,
                                                     const size_t idxMass) const
{
	if (idxMass != std::numeric_limits<size_t>::max()) {
		const std::complex<double> couplingPhaseSpace = cache.getCoupling(_id, idxChannel, idxBin, idxMass);
		if (couplingPhaseSpace != 0.) {
			return couplingPhaseSpace;
		}
	}

	const std::complex<double> coupling = fitParameters.getCoupling(_id, _channelsCoupling[idxChannel], idxBin);
	const std::complex<double> branching = fitParameters.getBranching(_id, _channelsBranching[idxChannel]);

	const channel& channelPhaseSpace = _channels[idxChannel];
	const double phaseSpace = channelPhaseSpace.getPhaseSpace(idxBin, mass, idxMass);

	const std::complex<double> couplingPhaseSpace = coupling * branching * phaseSpace;

	if (idxMass != std::numeric_limits<size_t>::max()) {
		cache.setCoupling(_id, idxChannel, idxBin, idxMass, couplingPhaseSpace);
	}

	return couplingPhaseSpace;
}


inline
std::ostream&
rpwa::resonanceFit::operator<< (std::ostream& out, const rpwa::resonanceFit::component& component)
{
	return component.print(out);
}


#endif // RESONANCEFIT_COMPONENTS_HH
