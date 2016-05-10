#include"getMassShapes.h"
#include"physUtils.hpp"
#include"decayGraph.hpp"

namespace {
	std::complex<double> getSingleMassShape(const rpwa::isobarDecayVertexPtr& vertex,           // vertex
	                                        const double                      mass,             // mass
	                                        const bool                        useBarrierFactors) {

		rpwa::particlePtr& parent    = vertex->parent();
		rpwa::particlePtr& daughter1 = vertex->daughter1();
		rpwa::particlePtr& daughter2 = vertex->daughter2();
		parent->setLzVec(TLorentzVector(0.,0.,0.,mass));
		double m1 = daughter1->mass();
		double m2 = daughter2->mass();
		double lambda = mass*mass*mass*mass + m1*m1*m1*m1 + m2*m2*m2*m2 - 2*(mass*mass*m1*m1 + mass*mass*m2*m2 + m1*m1*m2*m2);
		if (lambda < 0.) {
			printWarn << "kinematically not allowed point: mParent = " << mass << " mDaughter1 = " << m1 << "mDaughter2 = " << m2 << std::endl;
			return std::complex<double>(0.,0.);
		}
		double q = std::pow(lambda, .5)/(2*mass);

		daughter1->setLzVec(TLorentzVector(0.,0., q,std::pow(m1*m1+q*q, .5)));
		daughter2->setLzVec(TLorentzVector(0.,0.,-q,std::pow(m2*m2+q*q, .5)));

		std::complex<double> amp = vertex->massDepAmplitude();

		if (useBarrierFactors) {
			const int L = vertex->L();
			amp *= rpwa::barrierFactor(L, q, false);
		}
		return amp;
	}


	std::vector<std::complex<double> > getMassShapesRecursive(const rpwa::isobarDecayVertexPtr& vertex,            // vertex
	                                                          rpwa::isobarDecayTopologyPtr&     topo,
	                                                          const double                      mass,              // masss
	                                                          const bool                        useBarrierFactors) {

		std::vector<std::complex<double> > shapes(1,getSingleMassShape(vertex, mass, useBarrierFactors));
		rpwa::particlePtr daughter1 = vertex->daughter1();
		rpwa::isobarDecayVertexPtr daughter1Vertex = boost::dynamic_pointer_cast<rpwa::isobarDecayVertex>(topo->toVertex(daughter1));
		if (daughter1Vertex) {
			std::vector<std::complex<double> > daughter1Amps = getMassShapesRecursive(daughter1Vertex, topo, mass, useBarrierFactors);
			for (size_t i = 0; i < daughter1Amps.size(); ++i) {
				shapes.push_back(daughter1Amps[i]);
			}
		}
		rpwa::particlePtr daughter2 = vertex->daughter2();
		rpwa::isobarDecayVertexPtr daughter2Vertex = boost::dynamic_pointer_cast<rpwa::isobarDecayVertex>(topo->toVertex(daughter2));
		if (daughter2Vertex) {
			std::vector<std::complex<double> > daughter2Amps = getMassShapesRecursive(daughter2Vertex, topo, mass, useBarrierFactors);
			for (size_t i = 0; i < daughter2Amps.size(); ++i) {
				shapes.push_back(daughter2Amps[i]);
			}
		}
		return shapes;
	}
}

std::vector<std::complex<double> > rpwa::getMassShapes(isobarDecayTopologyPtr &topo,
	                                               const double            mass,
	                                               const bool              useBarrierFactors) {

	return ::getMassShapesRecursive(topo->XIsobarDecayVertex(), topo, mass, useBarrierFactors);
}
