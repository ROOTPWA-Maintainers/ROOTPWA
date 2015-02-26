///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2015 Sebastian Uhl (TUM)
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
//      implementation of cache results during the resonance fit
//
//-------------------------------------------------------------------------


#include "massDepFitCache.h"


rpwa::massDepFit::cache::cache(const size_t maxWaves,
                               const size_t maxComponents,
                               const size_t maxChannels,
                               const size_t maxBins,
                               const size_t maxMassBins)
	: _couplings(boost::extents[maxComponents][maxChannels][maxBins][maxMassBins]),
	  _components(boost::extents[maxComponents][maxBins][maxMassBins]),
	  _prodAmps(boost::extents[maxWaves][maxBins][maxMassBins])
{
}


void
rpwa::massDepFit::cache::setCoupling(const size_t idxComponent,
                                     const size_t idxChannel,
                                     const size_t idxBin,
                                     const size_t idxMassBin,
                                     const std::complex<double> coupling)
{
	if (idxBin == std::numeric_limits<size_t>::max() && idxMassBin == std::numeric_limits<size_t>::max()) {
		for (size_t idx=0; idx<*(_couplings.shape()+2) ; ++idx) {
			for (size_t jdx=0; jdx<*(_couplings.shape()+3) ; ++jdx) {
				_couplings[idxComponent][idxChannel][idx][jdx] = coupling;
			}
		}
	} else if (idxBin == std::numeric_limits<size_t>::max()) {
		for (size_t idx=0; idx<*(_couplings.shape()+2) ; ++idx) {
			_couplings[idxComponent][idxChannel][idx][idxMassBin] = coupling;
		}
	} else if (idxMassBin == std::numeric_limits<size_t>::max()) {
		for (size_t idx=0; idx<*(_couplings.shape()+3) ; ++idx) {
			_couplings[idxComponent][idxChannel][idxBin][idx] = coupling;
		}
	} else {
		_couplings[idxComponent][idxChannel][idxBin][idxMassBin] = coupling;
	}
}


void
rpwa::massDepFit::cache::setComponent(const size_t idxComponent,
                                      const size_t idxBin,
                                      const size_t idxMassBin,
                                      const std::complex<double> component)
{
	if (idxBin == std::numeric_limits<size_t>::max() && idxMassBin == std::numeric_limits<size_t>::max()) {
		for (size_t idx=0; idx<*(_components.shape()+1) ; ++idx) {
			for (size_t jdx=0; jdx<*(_components.shape()+2) ; ++jdx) {
				_components[idxComponent][idx][jdx] = component;
			}
		}
	} else if (idxBin == std::numeric_limits<size_t>::max()) {
		for (size_t idx=0; idx<*(_components.shape()+1) ; ++idx) {
			_components[idxComponent][idx][idxMassBin] = component;
		}
	} else if (idxMassBin == std::numeric_limits<size_t>::max()) {
		for (size_t idx=0; idx<*(_components.shape()+2) ; ++idx) {
			_components[idxComponent][idxBin][idx] = component;
		}
	} else {
		_components[idxComponent][idxBin][idxMassBin] = component;
	}
}


void
rpwa::massDepFit::cache::setProdAmp(const size_t idxWave,
                                    const size_t idxBin,
                                    const size_t idxMassBin,
                                    const std::complex<double> prodAmp)
{
	if (idxBin == std::numeric_limits<size_t>::max() && idxMassBin == std::numeric_limits<size_t>::max()) {
		for (size_t idx=0; idx<*(_prodAmps.shape()+1) ; ++idx) {
			for (size_t jdx=0; jdx<*(_prodAmps.shape()+2) ; ++jdx) {
				_prodAmps[idxWave][idx][jdx] = prodAmp;
			}
		}
	} else if (idxBin == std::numeric_limits<size_t>::max()) {
		for (size_t idx=0; idx<*(_prodAmps.shape()+1) ; ++idx) {
			_prodAmps[idxWave][idx][idxMassBin] = prodAmp;
		}
	} else if (idxMassBin == std::numeric_limits<size_t>::max()) {
		for (size_t idx=0; idx<*(_prodAmps.shape()+2) ; ++idx) {
			_prodAmps[idxWave][idxBin][idx] = prodAmp;
		}
	} else {
		_prodAmps[idxWave][idxBin][idxMassBin] = prodAmp;
	}
}


std::ostream&
rpwa::massDepFit::cache::print(std::ostream& out) const
{
	out << "cache for couplings: "
	    << *(_couplings.shape()) << " components, "
	    << *(_couplings.shape()+1) << " channels, "
	    << *(_couplings.shape()+2) << " bins, "
	    << *(_couplings.shape()+3) << " mass bins"
	    << std::endl;
	out << "cache for component values: "
	    << *(_components.shape()) << " components, "
	    << *(_components.shape()+1) << " bins, "
	    << *(_components.shape()+2) << " mass bins"
	    << std::endl;
	out << "cache for production amplitudes: "
	    << *(_prodAmps.shape()) << " waves, "
	    << *(_prodAmps.shape()+1) << " bins, "
	    << *(_prodAmps.shape()+2) << " mass bins"
	    << std::endl;

	return out;
}
