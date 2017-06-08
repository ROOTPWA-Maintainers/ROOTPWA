#include "stlContainers_py.h"

#include "boost/python/suite/indexing/vector_indexing_suite.hpp"
#include "boost/python/suite/indexing/map_indexing_suite.hpp"

#include<amplitudeMetadata.h>
#include<interactionVertex.h>
#include<isobarDecayVertex.h>
#include<particle.h>
#include<particleProperties.h>
#include<pwaLikelihood.h>
#include<waveDescription.h>

namespace bp = boost::python;

void rpwa::py::exportStlContainers() {

	// std::pair<std::string, rpwa::particleProperties>
	typedef std::pair< const std::string, rpwa::particleProperties > stdpair_int_particleProperties;
	bp::class_< stdpair_int_particleProperties >("__stdpair_int_particleProperties")
		.add_property("first", &stdpair_int_particleProperties::first)
		.def_readwrite("second", &stdpair_int_particleProperties::second);

	// std::vector<rpwa::particleProperties>
	bp::class_<std::vector<rpwa::particleProperties> >("__vector_particleProperties")
		.def(bp::vector_indexing_suite<std::vector<rpwa::particleProperties> >());

	// std::vector<rpwa::particle>
	bp::class_<std::vector<rpwa::particle> >("__vector_particle")
		.def(bp::vector_indexing_suite<std::vector<rpwa::particle>, true>());

	// std::vector<rpwa::particlePtr>
	bp::class_<std::vector<rpwa::particlePtr> >("__vector_particlePtr")
		.def(bp::vector_indexing_suite<std::vector<rpwa::particlePtr>, true>());

	// std::vector<rpwa::interactionVertexPtr>
	bp::class_<std::vector<rpwa::interactionVertexPtr> >("__vector_interactionVertexPtr")
		.def(bp::vector_indexing_suite<std::vector<rpwa::interactionVertexPtr>, true>());

	// std::vector<rpwa::isobarDecayVertexPtr>
	bp::class_<std::vector<rpwa::isobarDecayVertexPtr> >("__vector_isobarDecayVertexPtr")
		.def(bp::vector_indexing_suite<std::vector<rpwa::isobarDecayVertexPtr>, true>());

	// std::vector<std::string>
	bp::class_<std::vector<std::string> >("__vector_string")
		.def(bp::vector_indexing_suite<std::vector<std::string> >());

	// std::vector<unsigned int>
	bp::class_<std::vector<unsigned int> >("__vector_unsigned_int")
		.def(bp::vector_indexing_suite<std::vector<unsigned int> >());

	// std::vector<double>
	bp::class_<std::vector<double> >("__vector_double")
		.def(bp::vector_indexing_suite<std::vector<double> >());

	// std::vector<std::complex<double> >
	bp::class_<std::vector<std::complex<double> > >("__vector_complex_double")
		.def(bp::vector_indexing_suite<std::vector<std::complex<double> > >());

	// std::vector<rpwa::amplitudeMetadata*>
	bp::class_<std::vector<rpwa::amplitudeMetadata*> >("__vector_amplitudeMetadata")
		.def(bp::vector_indexing_suite<std::vector<rpwa::amplitudeMetadata*> >());

	// std::vector<rpwa::waveDescriptionPtr>
	bp::class_<std::vector<rpwa::waveDescriptionPtr> >("__vector_waveDescription")
		.def(bp::vector_indexing_suite<std::vector<rpwa::waveDescriptionPtr>, true>());

	// std::vector<rpwa::pwaLikelihood<double>::fitParameter>
	bp::class_<std::vector<rpwa::pwaLikelihood<std::complex<double> >::fitParameter> >("__vector_pwaLikelihood_fitParameter")
		.def(bp::vector_indexing_suite<std::vector<rpwa::pwaLikelihood<std::complex<double> >::fitParameter>, true>());

}


rpwa::multibinBoundariesType
rpwa::py::convertMultibinBoundariesFromPy(const boost::python::dict& pyMultibinBoundaries) {
	rpwa::multibinBoundariesType multibinBoundaries;
	boost::python::list keys = pyMultibinBoundaries.keys();
	for (int i = 0; i < boost::python::len(keys); ++i) {
		rpwa::boundaryType element;
		if (not rpwa::py::convertBPObjectToPair<double, double>(pyMultibinBoundaries[keys[i]], element)) {
			PyErr_SetString(PyExc_TypeError, "Got invalid pair for multibin boundaries");
			boost::python::throw_error_already_set();
		}
		boost::python::extract<std::string> getString(keys[i]);
		if (not getString.check()) {
			PyErr_SetString(PyExc_TypeError, "Got invalid key for multibin boundaries");
			boost::python::throw_error_already_set();
		}
		multibinBoundaries[getString()] = element;
	}
	return multibinBoundaries;
}


boost::python::dict
rpwa::py::convertMultibinBoundariesToPy(const rpwa::multibinBoundariesType& multibinBoundaries) {
	boost::python::dict pyMultibinBoundaries;
	for (const auto& variable_boundaries : multibinBoundaries) {
		const boost::python::str variablePy(variable_boundaries.first);
		pyMultibinBoundaries[variablePy] = boost::python::make_tuple(variable_boundaries.second.first, variable_boundaries.second.second);
	}
	return pyMultibinBoundaries;
}

boost::python::object
rpwa::py::convertMultibinBoundariesToPyMultibin(const rpwa::multibinBoundariesType& multibinBoundaries) {
			static bp::object mod = bp::import("pyRootPwa.utils");

			const bp::dict pyMultibinBoundaries = rpwa::py::convertMultibinBoundariesToPy(multibinBoundaries);

			return mod.attr("multiBin")(pyMultibinBoundaries);
}
