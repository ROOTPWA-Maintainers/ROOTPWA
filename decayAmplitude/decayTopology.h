///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
//
//    This file is part of rootpwa
//
//    rootpwa is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    rootpwa is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with rootpwa. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
//
// Description:
//      container class that holds all external information for
//      amplitude calculation
//      internally the decay process is represented as a graph
//      the graph is constraint to contain exactly one production
//      vertex and at least one interaction vertex; in addtion for
//      each final-state particle a corresponding final-state vertex
//      is created for internal use
//
//      "final-state" particles are the measured decay daughters;
//      additional final-state particles that belong to the production
//      process are handled by the production vertex
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef DECAYTOPOLOGY_H
#define DECAYTOPOLOGY_H


#include <vector>
#include <map>

#include "TLorentzRotation.h"

#include "particle.h"
#include "productionVertex.h"
#include "interactionVertex.h"
#include "fsVertex.h"
#include "decayGraph.hpp"


class TClonesArray;
class TVector3;


namespace rpwa {


	class decayTopology;
	typedef boost::shared_ptr<decayTopology> decayTopologyPtr;
	typedef decayGraph<interactionVertex, particle> decayTopologyGraphType;


	class decayTopology : public decayTopologyGraphType {

	public:

		decayTopology();
		decayTopology(const productionVertexPtr&               productionVertex,
		              const std::vector<interactionVertexPtr>& decayVertices,
		              const std::vector<particlePtr>&          fsParticles,
		              const bool                               performTopologyCheck = true);
		decayTopology(const decayTopology&                     topo);
		decayTopology(const decayTopologyGraphType&            graph);
		virtual ~decayTopology();

		decayTopology& operator =(const decayTopology&          topo);
		decayTopology& operator =(const decayTopologyGraphType& graph);
		decayTopologyPtr clone(const bool cloneFsParticles    = false,
		                       const bool cloneProdKinematics = false) const  ///< creates deep copy of decay topology; must not be virtual
		{ return decayTopologyPtr(doClone(cloneFsParticles, cloneProdKinematics)); }

		virtual void clear();  ///< deletes all information

		unsigned int nmbDecayVertices() const { return _decayVertices.size(); }  ///< returns number of decay vertices
		unsigned int nmbFsParticles  () const { return _fsParticles.size();   }  ///< returns number of final-state particles
		std::map<std::string, unsigned int> nmbIndistFsParticles() const;  ///< returns multiplicities of indistinguishable final-state particles

		int fsParticlesIntrinsicParity() const;  ///< returns intrinsic parity of final-state particles
		int spaceInvEigenValue()         const;  ///< returns eigenvalue of decay under space inversion
		int reflectionEigenValue()       const;  ///< returns eigenvalue of decay under reflection through production plane

		const std::vector<particlePtr>&          fsParticles  () const { return _fsParticles;   }  ///< returns final-state particles ordered depth-first
		const std::vector<interactionVertexPtr>& decayVertices() const { return _decayVertices; }  ///< returns decay vertices ordered depth-first

		const particlePtr&          XParticle       () const { return productionVertex()->XParticle(); }  ///< returns X particle
		const productionVertexPtr&  productionVertex() const { return _prodVertex;                     }  ///< returns production vertex
		const interactionVertexPtr& XDecayVertex    () const { return _decayVertices[0];               }  ///< returns X-decay vertex

		void transformFsParticles(const TLorentzRotation& L);  ///< applies Lorentz-transformation to all final-state particles

		bool isVertex          (const interactionVertexPtr& vert) const;  ///< returns whether given vertex is a vertex in this topology
		bool isParticle        (const particlePtr&          part) const;  ///< returns whether given particle is a particle in this topology
		bool isProductionVertex(const interactionVertexPtr& vert) const;  ///< returns whether given vertex is the production vertex
		bool isDecayVertex     (const interactionVertexPtr& vert) const;  ///< returns whether given vertex is one of the interaction vertices
		int  decayVertexIndex  (const interactionVertexPtr& vert) const;  ///< returns index of given vertex in decay-vertex array; -1 means vertex is not a decay vertex
		bool isFsVertex        (const interactionVertexPtr& vert) const;  ///< returns whether given vertex is one of the final-state vertices
		bool isFsParticle      (const particlePtr&          part) const;  ///< returns whether given particle is one of the final-state particles
		int  fsParticlesIndex  (const particlePtr&          part) const;  ///< returns index of given particle in final-state particle array; -1 means particle is not a final-state particle

		bool checkTopology   () const;                  ///< returns whether decay has the correct topology
		bool checkConsistency() const { return true; }  ///< checks consistency of information in vertices

		decayTopology subDecay(const nodeDesc& startNd,
		                       const bool      linkToMotherTopo = false);  ///< returns sub-decay tree that starts at given node
		decayTopology subDecay(const interactionVertexPtr& startVert,
		                       const bool                  linkToParentTopo = false)
		{ return subDecay(node(startVert), linkToParentTopo); }  ///< returns sub-decay tree that starts at given vertex

		void addDecay(const decayTopology& topo);  ///< copies all vertices and particles into this topology

		void setProductionVertex(const productionVertexPtr& productionVertex);  ///< (re)defines production vertex

		bool initKinematicsData(const std::vector<std::string>& prodKinParticles,
		                        const std::vector<std::string>& decayKinParticles);  ///< initializes input data

		bool initKinematicsData(const TClonesArray& prodKinParticles,
		                        const TClonesArray& decayKinParticles);  ///< initializes input data

		bool readKinematicsData(const std::vector<TVector3>& prodKinMomenta,
		                        const std::vector<TVector3>& decayKinMomenta);  ///< reads production and decay kinematics data and sets respective 4-momenta

		bool readKinematicsData(const TClonesArray& prodKinMomenta,
		                        const TClonesArray& decayKinMomenta);    ///< reads production and decay kinematics data and sets respective 4-momenta

		void fillKinematicsDataCache();  ///< copies kinematics data into cache; needed for Bose symmetrization

		bool revertMomenta();  ///< resets momenta to the values of last event read
		bool revertMomenta(const std::vector<unsigned int>& fsPartPermMap);  ///< resets momenta to the values of last event read, but reordering them according to index map

		void saveDecayToVertices(const decayTopologyPtr decay);


		virtual std::ostream& print(std::ostream& out) const;  ///< prints decay topology in human-readable form
		virtual std::ostream& printProdKinParticles (std::ostream& out) const;  ///< prints production kinematics data in human-readable form
		virtual std::ostream& printDecayKinParticles(std::ostream& out) const;  ///< prints decay kinematics data in human-readable form

		static bool debug() { return _debug; }                             ///< returns debug flag
		static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


	protected:

		virtual decayTopology* doClone(const bool cloneFsParticles,
		                               const bool cloneProdKinematics) const;  ///< helper function to use covariant return types with smart pointers; needed for public clone()

		decayTopology& constructDecay(const productionVertexPtr&               productionVertex,
		                              const std::vector<interactionVertexPtr>& decayVertices,
		                              const std::vector<particlePtr>&          fsParticles,
		                              const bool                               performTopologyCheck = true);  ///< constructs the decay graph based on production vertex, intermediate vertices, and final-state particles

		void buildInternalData();  ///< (re)builds internal data structure of vertex and particle pointers

		virtual interactionVertexPtr cloneNode(const nodeDesc& nd,
		                                       const bool      cloneInParticles  = false,
		                                       const bool      cloneOutParticles = false);
		virtual particlePtr          cloneEdge(const edgeDesc& ed);

	private:

		productionVertexPtr               _prodVertex;     ///< pointer to production vertex
		std::vector<interactionVertexPtr> _decayVertices;  ///< array of decay vertices; ordered depth-first
		std::vector<particlePtr>          _fsParticles;    ///< array of final-state particles; ordered depth-first

		std::map<unsigned int, unsigned int> _fsDataPartIndexMap;  ///< final-state particle indices in input data array
		std::vector<TVector3>                _fsDataPartMomCache;  ///< caches final-state momenta of last event read from input data; allows to "reset" kinematics for multiple passes over the same data

		static bool _debug;  ///< if set to true, debug messages are printed

	};


	inline
	decayTopologyPtr
	createDecayTopology(const productionVertexPtr&               productionVertex,
	                    const std::vector<interactionVertexPtr>& decayVertices,
	                    const std::vector<particlePtr>&          fsParticles,
	                    const bool                               performTopologyCheck = true)
	{
		decayTopologyPtr topo(new decayTopology(productionVertex, decayVertices,
		                                        fsParticles, performTopologyCheck));
		return topo;
	}


	inline
	decayTopologyPtr
	createDecayTopology(const decayTopology& topo)
	{
		decayTopologyPtr topoCopy(new decayTopology(topo));
		return topoCopy;
	}


	inline
	decayTopologyPtr
	createDecayTopology(const decayTopologyGraphType& graph)
	{
		decayTopologyPtr topoCopy(new decayTopology(graph));
		return topoCopy;
	}


	inline
	std::ostream&
	operator <<(std::ostream&        out,
	            const decayTopology& topo)
	{
		return topo.print(out);
	}


}  // namespace rpwa


#endif  // DECAYTOPOLOGY_H
