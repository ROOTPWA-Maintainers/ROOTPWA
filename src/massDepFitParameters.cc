#include "massDepFitParameters.h"


#include "reportingUtils.hpp"


rpwa::massDepFit::parameters::parameters()
	: _fixed(false)
{
}


rpwa::massDepFit::parameters::parameters(const size_t maxComponents,
                                         const size_t maxChannels,
                                         const size_t maxParameters,
                                         const size_t maxBins)
	: _fixed(true),
	  _branchings(boost::extents[maxComponents][maxChannels]),
	  _couplings(boost::extents[maxComponents][maxChannels][maxBins]),
	  _parameters(boost::extents[maxComponents][maxParameters])
{
}


void
rpwa::massDepFit::parameters::resize(const size_t maxComponents,
                                     const size_t maxChannels,
                                     const size_t maxParameters,
                                     const size_t maxBins)
{
	if (_fixed) {
		printErr << "cannot resize 'parameters' object that has been initialized with a size." << std::endl;
		throw;
	}

	std::vector<size_t> branchingsSize(_branchings.shape(), _branchings.shape()+_branchings.num_dimensions());
	std::vector<size_t> couplingsSize(_couplings.shape(), _couplings.shape()+_couplings.num_dimensions());
	std::vector<size_t> parametersSize(_parameters.shape(), _parameters.shape()+_parameters.num_dimensions());

	// check that all three arrays know about the same number of components
	// (strictly speaking _parameters should have one more than the other
	// two, but for simplicity's sake that's the way it is)
	assert(branchingsSize[0] == couplingsSize[0]);
	assert(branchingsSize[0] == parametersSize[0]);
	assert(couplingsSize[0] == parametersSize[0]);

	// check that the same number of channels is present
	assert(branchingsSize[1] == couplingsSize[1]);

	const size_t newComponents = std::max(maxComponents, couplingsSize[0]);
	const size_t newChannels   = std::max(maxChannels,   couplingsSize[1]);
	const size_t newParameters = std::max(maxParameters, parametersSize[1]);
	const size_t newBins       = std::max(maxBins,       couplingsSize[2]);

	_branchings.resize(boost::extents[newComponents][newChannels]);
	_couplings.resize(boost::extents[newComponents][newChannels][newBins]);
	_parameters.resize(boost::extents[newComponents][newParameters]);
}


std::ostream&
rpwa::massDepFit::parameters::print(std::ostream& out) const
{
	out << "branchings: " << *(_branchings.shape()) << " components, " << *(_branchings.shape()+1) << " channels" << std::endl;
	out << "couplings: " << *(_couplings.shape()) << " components, " << *(_couplings.shape()+1) << " channels, " << *(_couplings.shape()+2) << " bins" << std::endl;
	out << "parameters: " << *(_parameters.shape()) << " components, " << *(_parameters.shape()+1) << " parameters" << std::endl;

	return out;
}
