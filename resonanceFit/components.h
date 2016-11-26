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
#include "parameter.h"
#include "parameters.h"

namespace rpwa {

	namespace resonanceFit {

		class component {

		public:

			class channel {

			public:

				channel(const size_t waveIdx,
				        const std::string& waveName,
				        const std::vector<size_t>& bins,
				        const std::vector<size_t>& nrMassBins,
				        const boost::multi_array<double, 2>& massBinCenters,
				        const boost::multi_array<double, 2>& phaseSpaceIntegrals);

				size_t getWaveIdx() const { return _waveIdx; }
				const std::string& getWaveName() const { return _waveName; }

				bool isAnchor(const size_t idxBin) const { return _anchors[idxBin]; }
				void setAnchor(const size_t idxBin) { _anchors[idxBin] = true; }

				const std::vector<size_t>& getBins() const { return _bins; }

				double getPhaseSpaceIntegral(const size_t idxBin,
				                             const double mass,
				                             const size_t idxMass = std::numeric_limits<size_t>::max()) const;

			private:

				const size_t _waveIdx;
				const std::string _waveName;

				std::vector<bool> _anchors;

				const std::vector<size_t> _bins;

				const boost::multi_array<double, 2> _phaseSpaceIntegrals;
				std::vector<std::shared_ptr<const ROOT::Math::Interpolator> > _interpolators;

			};

			component(const size_t id,
			          const std::string& name,
			          const std::vector<rpwa::resonanceFit::parameter>& parameters,
			          const std::vector<rpwa::resonanceFit::component::channel>& decayChannels,
			          const std::vector<size_t>& nrMassBins,
			          const boost::multi_array<double, 2>& massBinCenters,
			          const bool useBranchings);
			virtual ~component() {}

			virtual std::string getType() const = 0;

			size_t getId() const { return _id; }
			const std::string& getName() const { return _name; }

			size_t getNrChannels() const { return _channels.size(); }
			const channel& getChannel(const size_t idxDecayChannel) const { return _channels[idxDecayChannel]; }
			const channel& getChannelFromCouplingIdx(const size_t idxCoupling) const { return _channels[_channelsFromCoupling[idxCoupling]]; }
			const channel& getChannelFromBranchingIdx(const size_t idxBranching) const { return _channels[_channelsFromBranching[idxBranching]]; }
			void setChannelAnchor(const std::vector<size_t>& channels);

			virtual size_t getTotalNrChannels() const { return _channels.size(); }

			size_t mapChannelToCoupling(const size_t idxDecayChannel) const { return _channelsCoupling[idxDecayChannel]; }
			size_t mapChannelToBranching(const size_t idxDecayChannel) const { return _channelsBranching[idxDecayChannel]; }

			size_t mapCouplingToMasterChannel(const size_t idxCoupling) const { return _channelsFromCoupling[idxCoupling]; }
			size_t mapBranchingToMasterChannel(const size_t idxBranching) const { return _channelsFromBranching[idxBranching]; }

			size_t getNrCouplings() const { return _nrCouplings; }
			size_t importCouplings(const double* par,
			                       rpwa::resonanceFit::parameters& fitParameters,
			                       rpwa::resonanceFit::cache& cache) const;

			size_t getNrBranchings() const { return _nrBranchings; }
			size_t importBranchings(const double* par,
			                        rpwa::resonanceFit::parameters& fitParameters,
			                        rpwa::resonanceFit::cache& cache) const;

			size_t getNrParameters() const { return _parameters.size(); }
			size_t importParameters(const double* par,
			                        rpwa::resonanceFit::parameters& fitParameters,
			                        rpwa::resonanceFit::cache& cache) const;

			bool isBranchingFixed(const size_t idxBranching) const { return _branchingsFixed[idxBranching]; }

			const rpwa::resonanceFit::parameter& getParameter(const size_t idxParameter) const { return _parameters[idxParameter]; }

			std::complex<double> val(const rpwa::resonanceFit::parameters& fitParameters,
			                         rpwa::resonanceFit::cache& cache,
			                         const size_t idxBin,
			                         const double mass,
			                         const size_t idxMass = std::numeric_limits<size_t>::max()) const;

			std::complex<double> getCouplingPhaseSpace(const rpwa::resonanceFit::parameters& fitParameters,
			                                           rpwa::resonanceFit::cache& cache,
			                                           const size_t idxChannel,
			                                           const size_t idxBin,
			                                           const double mass,
			                                           const size_t idxMass = std::numeric_limits<size_t>::max()) const;

			virtual std::ostream& print(std::ostream& out = std::cout, const bool newLine = true) const;

		private:

			virtual std::complex<double> val(const rpwa::resonanceFit::parameters& fitParameters,
			                                 const size_t idxBin,
			                                 const double mass) const = 0;

			const size_t _id;
			const std::string _name;

			std::vector<channel> _channels;
			std::vector<size_t> _channelsCoupling;
			std::vector<size_t> _channelsBranching;

			size_t _nrCouplings;
			std::vector<size_t> _channelsFromCoupling;

			size_t _nrBranchings;
			std::vector<size_t> _channelsFromBranching;
			std::vector<bool> _branchingsFixed;

			const std::vector<rpwa::resonanceFit::parameter> _parameters;

		protected:

			std::vector<std::vector<size_t> > _binsEqualValues;

		};

		std::ostream& operator<< (std::ostream& out, const rpwa::resonanceFit::component& component);

		class fixedWidthBreitWigner : public component {

		public:

			static std::vector<rpwa::resonanceFit::parameter> getDefaultParameters();

			fixedWidthBreitWigner(const size_t id,
			                      const std::string& name,
			                      const std::vector<rpwa::resonanceFit::parameter>& parameters,
			                      const std::vector<rpwa::resonanceFit::component::channel>& decayChannels,
			                      const std::vector<size_t>& nrMassBins,
			                      const boost::multi_array<double, 2>& massBinCenters,
			                      const bool useBranchings);

			virtual std::string getType() const { return "fixedWidthBreitWigner"; }

		private:

			virtual std::complex<double> val(const rpwa::resonanceFit::parameters& fitParameters,
			                                 const size_t idxBin,
			                                 const double mass) const;

		};

		class dynamicWidthBreitWigner : public component {

		public:

			static std::vector<rpwa::resonanceFit::parameter> getDefaultParameters();

			dynamicWidthBreitWigner(const size_t id,
			                        const std::string& name,
			                        const std::vector<rpwa::resonanceFit::parameter>& parameters,
			                        const std::vector<rpwa::resonanceFit::component::channel>& decayChannels,
			                        const std::vector<size_t>& nrMassBins,
			                        const boost::multi_array<double, 2>& massBinCenters,
			                        const bool useBranchings,
			                        const std::vector<double>& branchingRatio,
			                        const std::vector<int>& relAngularMom,
			                        const std::vector<double>& mIsobar1,
			                        const std::vector<double>& mIsobar2);

			virtual std::string getType() const { return "dynamicWidthBreitWigner"; }

			virtual size_t getTotalNrChannels() const { return _ratio.size(); }

			const std::vector<double>& branchingRatio() const { return _ratio; }
			const std::vector<int>& relAngularMom() const { return _l; }
			const std::vector<double>& mIsobar1() const { return _m1; }
			const std::vector<double>& mIsobar2() const { return _m2; }

			virtual std::ostream& print(std::ostream& out = std::cout, const bool newLine = true) const;

		private:

			virtual std::complex<double> val(const rpwa::resonanceFit::parameters& fitParameters,
			                                 const size_t idxBin,
			                                 const double mass) const;

			std::vector<double> _ratio;
			const std::vector<int> _l;
			const std::vector<double> _m1;
			const std::vector<double> _m2;

		};

		class integralWidthBreitWigner : public component {

		public:

			static std::vector<rpwa::resonanceFit::parameter> getDefaultParameters();

			integralWidthBreitWigner(const size_t id,
			                         const std::string& name,
			                         const std::vector<rpwa::resonanceFit::parameter>& parameters,
			                         const std::vector<rpwa::resonanceFit::component::channel>& decayChannels,
			                         const std::vector<size_t>& nrMassBins,
			                         const boost::multi_array<double, 2>& massBinCenters,
			                         const bool useBranchings,
			                         const std::vector<double>& branchingRatio,
			                         const std::vector<std::vector<double> >& masses,
			                         const std::vector<std::vector<double> >& values);

			virtual std::string getType() const { return "integralWidthBreitWigner"; }

			virtual size_t getTotalNrChannels() const { return _ratio.size(); }

			const std::vector<double>& branchingRatio() const { return _ratio; }
			const std::vector<std::vector<double> >& masses() const { return _masses; }
			const std::vector<std::vector<double> >& values() const { return _values; }

		private:

			virtual std::complex<double> val(const rpwa::resonanceFit::parameters& fitParameters,
			                                 const size_t idxBin,
			                                 const double mass) const;

			std::vector<double> _ratio;
			const std::vector<std::vector<double> > _masses;
			const std::vector<std::vector<double> > _values;
			std::vector<std::shared_ptr<const ROOT::Math::Interpolator> > _interpolator;

		};

		class constantBackground : public component {

		public:

			static std::vector<rpwa::resonanceFit::parameter> getDefaultParameters();

			constantBackground(const size_t id,
			                   const std::string& name,
			                   const std::vector<rpwa::resonanceFit::parameter>& parameters,
			                   const std::vector<rpwa::resonanceFit::component::channel>& decayChannels,
			                   const std::vector<size_t>& nrMassBins,
			                   const boost::multi_array<double, 2>& massBinCenters,
			                   const bool useBranchings);

			virtual std::string getType() const { return "constantBackground"; }

		private:

			virtual std::complex<double> val(const rpwa::resonanceFit::parameters& fitParameters,
			                                 const size_t idxBin,
			                                 const double mass) const;

		};

		class exponentialBackground : public component {

		public:

			static std::vector<rpwa::resonanceFit::parameter> getDefaultParameters();

			exponentialBackground(const size_t id,
			                      const std::string& name,
			                      const std::vector<rpwa::resonanceFit::parameter>& parameters,
			                      const std::vector<rpwa::resonanceFit::component::channel>& decayChannels,
			                      const std::vector<size_t>& nrMassBins,
			                      const boost::multi_array<double, 2>& massBinCenters,
			                      const bool useBranchings,
			                      const int relAngularMom,
			                      const double mIsobar1,
			                      const double mIsobar2,
			                      const double exponent);

			virtual std::string getType() const { return "exponentialBackground"; }

			int relAngularMom() const { return _l; }
			double mIsobar1() const { return _m1; }
			double mIsobar2() const { return _m2; }
			double exponent() const { return _exponent; }

			virtual std::ostream& print(std::ostream& out = std::cout, const bool newLine = true) const;

		private:

			virtual std::complex<double> val(const rpwa::resonanceFit::parameters& fitParameters,
			                                 const size_t idxBin,
			                                 const double mass) const;

			const int _l;
			const double _m1;
			const double _m2;
			const double _exponent;

			double _norm;

		};

		class tPrimeDependentBackground : public component {

		public:

			static std::vector<rpwa::resonanceFit::parameter> getDefaultParameters();

			tPrimeDependentBackground(const size_t id,
			                          const std::string& name,
			                          const std::vector<rpwa::resonanceFit::parameter>& parameters,
			                          const std::vector<rpwa::resonanceFit::component::channel>& decayChannels,
			                          const std::vector<size_t>& nrMassBins,
			                          const boost::multi_array<double, 2>& massBinCenters,
			                          const bool useBranchings,
			                          const std::vector<double>& tPrimeMeans,
			                          const int relAngularMom,
			                          const double mIsobar1,
			                          const double mIsobar2,
			                          const double exponent);

			virtual std::string getType() const { return "tPrimeDependentBackground"; }

			int relAngularMom() const { return _l; }
			double mIsobar1() const { return _m1; }
			double mIsobar2() const { return _m2; }
			double exponent() const { return _exponent; }

			virtual std::ostream& print(std::ostream& out = std::cout, const bool newLine = true) const;

		private:

			virtual std::complex<double> val(const rpwa::resonanceFit::parameters& fitParameters,
			                                 const size_t idxBin,
			                                 const double mass) const;

			const std::vector<double> _tPrimeMeans;

			const int _l;
			const double _m1;
			const double _m2;
			const double _exponent;

			double _norm;

		};

		class exponentialBackgroundIntegral : public component {

		public:

			static std::vector<rpwa::resonanceFit::parameter> getDefaultParameters();

			exponentialBackgroundIntegral(const size_t id,
			                              const std::string& name,
			                              const std::vector<rpwa::resonanceFit::parameter>& parameters,
			                              const std::vector<rpwa::resonanceFit::component::channel>& decayChannels,
			                              const std::vector<size_t>& nrMassBins,
			                              const boost::multi_array<double, 2>& massBinCenters,
			                              const bool useBranchings,
			                              const std::vector<double>& masses,
			                              const std::vector<double>& values,
			                              const double exponent);

			virtual std::string getType() const { return "exponentialBackgroundIntegral"; }

			const std::vector<double>& masses() const { return _masses; }
			const std::vector<double>& values() const { return _values; }
			double exponent() const { return _exponent; }

			virtual std::ostream& print(std::ostream& out = std::cout, const bool newLine = true) const;

		private:

			virtual std::complex<double> val(const rpwa::resonanceFit::parameters& fitParameters,
			                                 const size_t idxBin,
			                                 const double mass) const;

			const std::vector<double> _masses;
			const std::vector<double> _values;
			std::shared_ptr<const ROOT::Math::Interpolator> _interpolator;
			const double _exponent;

			double _norm;

		};

		class tPrimeDependentBackgroundIntegral : public component {

		public:

			static std::vector<rpwa::resonanceFit::parameter> getDefaultParameters();

			tPrimeDependentBackgroundIntegral(const size_t id,
			                                  const std::string& name,
			                                  const std::vector<rpwa::resonanceFit::parameter>& parameters,
			                                  const std::vector<rpwa::resonanceFit::component::channel>& decayChannels,
			                                  const std::vector<size_t>& nrMassBins,
			                                  const boost::multi_array<double, 2>& massBinCenters,
			                                  const bool useBranchings,
			                                  const std::vector<double>& tPrimeMeans,
			                                  const std::vector<double>& masses,
			                                  const std::vector<double>& values,
			                                  const double exponent);

			virtual std::string getType() const { return "tPrimeDependentBackgroundIntegral"; }

			const std::vector<double>& masses() const { return _masses; }
			const std::vector<double>& values() const { return _values; }
			double exponent() const { return _exponent; }

			virtual std::ostream& print(std::ostream& out = std::cout, const bool newLine = true) const;

		private:

			virtual std::complex<double> val(const rpwa::resonanceFit::parameters& fitParameters,
			                                 const size_t idxBin,
			                                 const double mass) const;

			const std::vector<double> _tPrimeMeans;

			const std::vector<double> _masses;
			const std::vector<double> _values;
			std::shared_ptr<const ROOT::Math::Interpolator> _interpolator;
			const double _exponent;

			double _norm;

		};

	} // end namespace resonanceFit

} // end namespace rpwa


inline
double
rpwa::resonanceFit::component::channel::getPhaseSpaceIntegral(const size_t idxBin,
                                                              const double mass,
                                                              const size_t idxMass) const
{
	if(idxMass != std::numeric_limits<size_t>::max()) {
		return _phaseSpaceIntegrals[idxBin][idxMass];
	}

	return _interpolators[idxBin]->Eval(mass);
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
	if(idxMass != std::numeric_limits<size_t>::max()) {
		const std::complex<double> couplingPhaseSpace = cache.getCoupling(_id, idxChannel, idxBin, idxMass);
		if(couplingPhaseSpace != 0.) {
			return couplingPhaseSpace;
		}
	}

	const std::complex<double> coupling = fitParameters.getCoupling(_id, _channelsCoupling[idxChannel], idxBin);
	const std::complex<double> branching = fitParameters.getBranching(_id, _channelsBranching[idxChannel]);

	const channel& channelPhaseSpaceIntegral = _channels[idxChannel];
	const double phaseSpaceIntegral = channelPhaseSpaceIntegral.getPhaseSpaceIntegral(idxBin, mass, idxMass);

	const std::complex<double> couplingPhaseSpace = coupling * branching * phaseSpaceIntegral;

	if(idxMass != std::numeric_limits<size_t>::max()) {
		cache.setCoupling(_id, idxChannel, idxBin, idxMass, couplingPhaseSpace);
	}

	return couplingPhaseSpace;
}


inline
std::ostream&
rpwa::resonanceFit::operator<< (std::ostream& out, const rpwa::resonanceFit::component& component)
{
	return component.print(out, false);
}


#endif // RESONANCEFIT_COMPONENTS_HH
