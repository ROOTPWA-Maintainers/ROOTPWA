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
// File and Version Information:
// $Rev::                             $: revision of last commit
// $Author::                          $: author of last commit
// $Date::                            $: date of last commit
//
// Description:
//      container class that holds all external information for
//      amplitude calculation
//      internally the decay process is represented as a graph
//      the graph is constraint to contain exactly one production
//      vertex and at least one interaction vertex; in addtion for
//      each final state particle a corresponding final state vertex
//      is created for internal use
//
//      "final state" particles are the measured decay daughters;
//      additional final state particles that belong to the production
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

#include "particle.h"
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
    decayTopology(const interactionVertexPtr&              productionVertex,
		  const std::vector<interactionVertexPtr>& interactionVertices,
		  const std::vector<particlePtr>&          fsParticles);
    decayTopology(const decayTopology&                     topo);
    decayTopology(const decayTopologyGraphType&            graph);
    virtual ~decayTopology();

    virtual decayTopology& operator =(const decayTopology&          topo);
    virtual decayTopology& operator =(const decayTopologyGraphType& graph);
    decayTopologyPtr clone(const bool cloneFsParticles    = false,
			   const bool cloneProdKinematics = false) const  ///< creates deep copy of decay topology; must not be virtual
    { return decayTopologyPtr(doClone(cloneFsParticles, cloneProdKinematics)); }

    virtual void clear();  ///< deletes all information
    
    unsigned int nmbInteractionVertices() const { return _intVertices.size(); }  ///< returns number of interaction vertices
    unsigned int nmbFsParticles        () const { return _fsParticles.size(); }  ///< returns number of final state particles
    std::map<std::string, unsigned int> nmbIndistFsParticles() const;  ///< returns multiplicities of indistinguishable final state particles

    const std::vector<particlePtr>&          fsParticles        () const { return _fsParticles; }  ///< returns final state particles ordered depth-first
    const std::vector<interactionVertexPtr>& interactionVertices() const { return _intVertices; }  ///< returns interaction vertices (excluding production vertex) ordered depth-first

    const interactionVertexPtr& productionVertex  () const { return _prodVertex;     }  ///< returns production vertex
    const interactionVertexPtr& xInteractionVertex() const { return _intVertices[0]; }  ///< returns X-decay vertex

    bool isProductionVertex (const interactionVertexPtr& vert) const { return (vert == _prodVertex); }  ///< returns whether given vertex is the production vertex
    bool isInteractionVertex(const interactionVertexPtr& vert) const;  ///< returns whether given vertex is one of the interaction vertices
    bool isFsVertex         (const interactionVertexPtr& vert) const;  ///< returns whether given vertex is one of the final state vertices
    bool isFsParticle       (const particlePtr&          part) const;  ///< returns whether given particle is one of the final state particles
    int  fsParticleIndex    (const particlePtr&          part) const;  ///< returns index in final state particle array; -1 means particle is not a final state particle

    bool checkTopology   () const;                  ///< returns whether decay has the correct topology
    bool checkConsistency() const { return true; }  ///< checks consistency of information in vertices

    decayTopology subDecay(const nodeDesc& startNd,
			   const bool      linkToMotherTopo = false);  ///< returns sub-decay tree that starts at given vertex

    void addDecay(const decayTopology& topo);  ///< copies all vertices and particles into this topology

    void setProductionVertex(const interactionVertexPtr& productionVertex);  ///< (re)defines production vertex

    bool readData(const TClonesArray& prodKinParticles,
		  const TClonesArray& prodKinMomenta,
		  const TClonesArray& decayKinParticles,
		  const TClonesArray& decayKinMomenta);  ///< reads production and decay kinematics data and sets respective 4-momenta

    bool revertMomenta();  ///< resets momenta to the values of last event read
    bool revertMomenta(const std::vector<unsigned int>& indexMap);  ///< resets momenta to the values of last event read, but reordering them according to index map


    virtual std::ostream& print(std::ostream& out) const;  ///< prints decay topology in human-readable form
    virtual std::ostream& printProdKinParticles (std::ostream& out) const;  ///< prints production kinematics data in human-readable form
    virtual std::ostream& printDecayKinParticles(std::ostream& out) const;  ///< prints decay kinematics data in human-readable form

    static bool debug() { return _debug; }                             ///< returns debug flag
    static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag


  protected:

    virtual decayTopology* doClone(const bool cloneFsParticles,
				   const bool cloneProdKinematics) const;  ///< helper function to use covariant return types with smart pointers; needed for public clone()

    decayTopology& constructDecay(const interactionVertexPtr&              productionVertex,
				  const std::vector<interactionVertexPtr>& interactionVertices,
				  const std::vector<particlePtr>&          fsParticles);  ///< constructs the decay graph based on, production vertex, intermediate vertices, and final state particles

    void buildInternalData();  ///< (re)builds internal data structure of vertex and particle pointers

    virtual interactionVertexPtr cloneNode(const nodeDesc& nd,
					   const bool      cloneInParticles  = false,
					   const bool      cloneOutParticles = false);
    virtual particlePtr          cloneEdge(const edgeDesc& ed);

  private:

    interactionVertexPtr              _prodVertex;      ///< pointer to production vertex
    std::vector<interactionVertexPtr> _intVertices;     ///< array of interaction vertices excluding production vertex; ordered depth-first
    std::vector<particlePtr>          _fsParticles;     ///< array of final state particles; ordered depth-first
    std::vector<TVector3>             _fsPartMomCache;  ///< caches final state momenta of last event read from input data; allows to "reset" kinematics for multiple passes over the same data
    
    static bool _debug;  ///< if set to true, debug messages are printed
    
  };
  

  inline
  decayTopologyPtr
  createDecayTopology(const interactionVertexPtr&              productionVertex,
		      const std::vector<interactionVertexPtr>& interactionVertices,
		      const std::vector<particlePtr>&          fsParticles)
  {
    decayTopologyPtr t(new decayTopology(productionVertex, interactionVertices, fsParticles));
    return t;
  }


  inline
  decayTopologyPtr
  createDecayTopology(const decayTopology& topo)
  {
    decayTopologyPtr t(new decayTopology(topo));
    return t;
  }


  inline
  decayTopologyPtr
  createDecayTopology(const decayTopologyGraphType& graph)
  {
    decayTopologyPtr t(new decayTopology(graph));
    return t;
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
