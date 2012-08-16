
#include <Python.h>

#include <boost/python/class.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>

#include <isobarHelicityAmplitude.h>

BOOST_PYTHON_MODULE(libRootPwaPy)
{
using namespace boost::python;
//class_<rpwa::isobarAmplitude, boost::noncopyable>("isobarAmplitude", no_init);
class_<rpwa::isobarHelicityAmplitude>("isobarAmplitude")
	.def("name", &rpwa::isobarHelicityAmplitude::name);
}


