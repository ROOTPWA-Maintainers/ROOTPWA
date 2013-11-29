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
//      amplitude calculation of isobar decay
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef ISOBARDECAYTOPOLOGY_H
#define ISOBARDECAYTOPOLOGY_H

#include <boost/tuple/tuple.hpp>

#include "isobarDecayVertex.h"
#include "decayTopology.h"


namespace rpwa {

	// amplitude symmetrization info
	struct symTermMap {
		symTermMap(const std::complex<double>&      f,
		           const std::vector<unsigned int>& m)
			: factor       (f),
			  fsPartPermMap(m)
		{ }
		std::complex<double>      factor;         ///< factor to be applied to symmetrization term
		std::vector<unsigned int> fsPartPermMap;  ///< final-state-particle permutation map
	};


	class isobarDecayTopology;
	typedef boost::shared_ptr<isobarDecayTopology> isobarDecayTopologyPtr;


	class isobarDecayTopology : public decayTopology {

	public:

		isobarDecayTopology();
		isobarDecayTopology(const productionVertexPtr&               productionVertex,
		                    const std::vector<isobarDecayVertexPtr>& isobarDecayVertices,
		                    const std::vector<particlePtr>&          fsParticles,
		                    const bool                               performTopologyCheck = true);
		isobarDecayTopology(const productionVertexPtr&               productionVertex,
		                    const std::vector<interactionVertexPtr>& isobarDecayVertices,
		                    const std::vector<particlePtr>&          fsParticles,
		                    const bool                               performTopologyCheck = true);
		isobarDecayTopology(const isobarDecayTopology&               topo);
		isobarDecayTopology(const decayTopology&                     topo);
		virtual ~isobarDecayTopology();

		isobarDecayTopology& operator =(const isobarDecayTopology& topo);
		isobarDecayTopology& operator =(const decayTopology&       topo);
		isobarDecayTopologyPtr clone(const bool cloneFsParticles    = false,
		                             const bool cloneProdKinematics = false) const  ///< creates deep copy of isobar decay topology; must not be virtual
		{ return isobarDecayTopologyPtr(doClone(cloneFsParticles, cloneProdKinematics)); }
		virtual void clear();  ///< deletes all information

		const std::vector<isobarDecayVertexPtr>& isobarDecayVertices() const { return _isobarVertices;    }  ///< returns all isobar decay vertices ordered by depth-first; first vertex is X-decay vertex
		const isobarDecayVertexPtr&              XIsobarDecayVertex () const { return _isobarVertices[0]; }  ///< returns X-decay vertex

		bool checkTopology   () const;  ///< returns whether decay has the correct topology
		bool checkConsistency() const;  ///< checks conservation rules on all vertices

		isobarDecayTopology subDecay(const nodeDesc& startNd,
		                             const bool      linkToParentTopo = false);  ///< returns sub-decay tree that starts at given node
		isobarDecayTopology subDecay(const isobarDecayVertexPtr& startVert,
		                             const bool                  linkToParentTopo = false)
		{ return subDecay(node(startVert), linkToParentTopo); }  ///< returns sub-decay tree that starts at given vertex
		isobarDecayTopology subDecayConsistent(const isobarDecayVertexPtr& startVertex); ///< returns sub-decay tree that starts at given vertex, but adds a nonInteractionVertex as production vertex to have a consistent topology

		void addDecay(const isobarDecayTopology& topo);  ///< returns sub-decay tree that starts at given vertex

		static isobarDecayTopology joinDaughterDecays
		(const isobarDecayVertexPtr&             parentVertex,
		 const std::vector<isobarDecayTopology>& daughterDecays);  ///< joins daughter decay graphs and connects them to a common parent vertex
		static isobarDecayTopology joinDaughterDecays
		(const isobarDecayVertexPtr& parentVertex,
		 const isobarDecayTopology&  daughter1Decay,
		 const isobarDecayTopology&  daughter2Decay);  ///< joins daughter decay graphs and connects them to a common parent vertex

		const TLorentzVector& calcIsobarLzVec();  ///< (re)calculates Lorentz-vectors of all isobars in the decay from final-state particles and returns Lorentz-vector of X-system

		void calcIsobarCharges   (bool quiet = false);  ///< sets isobar charges as defined by final-state particles
		void calcIsobarBaryonNmbs();                    ///< sets isobar baryon numbers as defined by final-state particles

		virtual std::ostream& print(std::ostream& out) const;  ///< prints decay topology in human-readable form

		virtual std::ostream& writeGraphViz(std::ostream&      out);          ///< writes graph in GraphViz DOT format
		virtual bool          writeGraphViz(const std::string& outFileName);  ///< writes graph in GraphViz DOT format

		static bool debug() { return _debug; }                             ///< returns debug flag
		static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag

		double getIsospinClebschGordanProduct(isobarDecayVertexPtr vertex = isobarDecayVertexPtr()) const;  ///< returns product of isospin Clebsch-Gordans for all two-body decays in the topology

		std::vector<symTermMap> getIsospinSymmetrization();  ///< returns all final-state permutations needed for isospin symmetrization
		std::vector<symTermMap> getBoseSymmetrization() const;     ///< returns all final-state permutations needed for Bose symmetrization

		bool isobarIsAffectedByPermutation(const isobarDecayVertexPtr& vertex, const std::vector<unsigned int>& permutation) const; ///< returns true if the isobar decaying to the vertex is changed by the permutation
		bool daughtersAreAffectedByPermutation(const isobarDecayVertexPtr& vertex, const std::vector<unsigned int>& permutation) const; ///< returns true if the daughters the vertex decays into are changed by the permutation
		std::vector<unsigned int> getFsPartIndicesConnectedToVertex(const isobarDecayVertexPtr& vertex) const; ///< returns the indices of the final state particles which are 'below' the given vertex

		std::vector<unsigned int> findIsobarBoseSymVertices() const;  ///< returns indices of all isobar vertices that have isobar daughters that decay into the same final state


	protected:

		isobarDecayTopology& constructDecay(const productionVertexPtr&               productionVertex,
		                                    const std::vector<isobarDecayVertexPtr>& isobarDecayVertices,
		                                    const std::vector<particlePtr>&          fsParticles,
		                                    const bool                               performTopologyCheck = true);  ///< constructs the decay graph based on final-state particles and vertices
		isobarDecayTopology& constructDecay(const productionVertexPtr&               productionVertex,
		                                    const std::vector<interactionVertexPtr>& isobarDecayVertices,
		                                    const std::vector<particlePtr>&          fsParticles,
		                                    const bool                               performTopologyCheck = true);  ///< constructs the decay graph based on final-state particles and vertices


	private:

		virtual isobarDecayTopology* doClone(const bool cloneFsParticles,
		                                     const bool cloneProdKinematics) const;  ///< helper function to use covariant return types with smart pointers; needed for public clone()

		void buildIsobarVertexArray();  ///< (re)builds array of isobar decay vertices

		virtual void genBoseSymTermMaps
		(const std::map<std::string, std::vector<unsigned int> >&     origFsPartIndices,
		 const std::map<std::string, std::vector<unsigned int> >&     newFsPartIndices,
		 std::map<std::string, std::vector<unsigned int> >::iterator& newFsPartIndicesEntry,
		 std::vector<symTermMap>&                                     symTermMaps) const;  ///< recursive function that generates all permutation maps of indistinguishable final state particles

		std::vector<isobarDecayVertexPtr> _isobarVertices;  ///< array of isobar-decay vertices excluding production vertex; ordered depth-first; this is a copy of the respective array in decayTopology

		static bool _debug;  ///< if set to true, debug messages are printed

	};


	inline
	isobarDecayTopologyPtr
	createIsobarDecayTopology(const productionVertexPtr&               productionVertex,
	                          const std::vector<isobarDecayVertexPtr>& isobarDecayVertices,
	                          const std::vector<particlePtr>&          fsParticles,
	                          const bool                               performTopologyCheck = true)
	{
		isobarDecayTopologyPtr topo(new isobarDecayTopology(productionVertex, isobarDecayVertices,
		                                                    fsParticles, performTopologyCheck));
		return topo;
	}


	inline
	isobarDecayTopologyPtr
	createIsobarDecayTopology(const productionVertexPtr&               productionVertex,
	                          const std::vector<interactionVertexPtr>& isobarDecayVertices,
	                          const std::vector<particlePtr>&          fsParticles,
	                          const bool                               performTopologyCheck = true)
	{
		isobarDecayTopologyPtr topo(new isobarDecayTopology(productionVertex, isobarDecayVertices,
		                                                    fsParticles, performTopologyCheck));
		return topo;
	}


	inline
	isobarDecayTopologyPtr
	createIsobarDecayTopology(const isobarDecayTopology& topo)
	{
		isobarDecayTopologyPtr topoCopy(new isobarDecayTopology(topo));
		return topoCopy;
	}


	inline
	isobarDecayTopologyPtr
	createIsobarDecayTopology(const decayTopology& topo)
	{
		isobarDecayTopologyPtr topoCopy(new isobarDecayTopology(topo));
		return topoCopy;
	}


	inline
	std::ostream&
	operator <<(std::ostream&              out,
	            const isobarDecayTopology& topo)
	{
		return topo.print(out);
	}


} // namespace rpwa


#endif  // ISOBARDECAYTOPOLOGY_H
