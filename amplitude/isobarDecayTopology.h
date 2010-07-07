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


#include <boost/shared_ptr.hpp>

#include "isobarDecayVertex.h"
#include "decayTopology.h"


namespace rpwa {  


  class isobarDecayTopology;
	typedef boost::shared_ptr<isobarDecayTopology> isobarDecayTopologyPtr;
	

  class isobarDecayTopology : public decayTopology {
  
  public:
      
    isobarDecayTopology();
    isobarDecayTopology(const interactionVertexPtr&              productionVertex,
                        const std::vector<isobarDecayVertexPtr>& isobarDecayVertices,
                        const std::vector<particlePtr>&          fsParticles);
    isobarDecayTopology(const interactionVertexPtr&              productionVertex,
                        const std::vector<interactionVertexPtr>& isobarDecayVertices,
                        const std::vector<particlePtr>&          fsParticles);
    isobarDecayTopology(const isobarDecayTopology&               topo);
    isobarDecayTopology(const decayTopology&                     topo);
    virtual ~isobarDecayTopology();

    virtual isobarDecayTopology& operator =(const isobarDecayTopology& topo);
    virtual isobarDecayTopology& operator =(const decayTopology&       topo);
    isobarDecayTopologyPtr clone(const bool cloneFsParticles    = false,
                                 const bool cloneProdKinematics = false) const  ///< creates deep copy of isobar decay topology; must not be virtual
	  { return isobarDecayTopologyPtr(doClone(cloneFsParticles, cloneProdKinematics)); }
	  virtual void clear();  ///< deletes all information

	  const std::vector<isobarDecayVertexPtr>& isobarDecayVertices() const { return _isobarVertices;    }  ///< returns all isobar decay vertices ordered by depth-first; first vertex is X-decay vertex
	  const isobarDecayVertexPtr               XIsobarDecayVertex () const { return _isobarVertices[0]; }  ///< returns X-decay vertex

    bool checkTopology   () const;  ///< returns whether decay has the correct topology
    bool checkConsistency() const;  ///< checks conservation rules on all vertices

    isobarDecayTopology subDecay(const nodeDesc& startNd,
                                 const bool      linkToMotherTopo = false);  ///< returns sub-decay tree that starts at given vertex

    void addDecay(const isobarDecayTopology& topo);  ///< returns sub-decay tree that starts at given vertex

    static isobarDecayTopology joinDaughterDecays(const isobarDecayVertexPtr&             motherVertex,
                                                  const std::vector<isobarDecayTopology>& daughterDecays);  ///< joins daughter decay graphs and connects them to a common mother vertex
	  static isobarDecayTopology joinDaughterDecays(const isobarDecayVertexPtr& motherVertex,
	                                                const isobarDecayTopology&  daughter1Decay,
	                                                const isobarDecayTopology&  daughter2Decay);  ///< joins daughter decay graphs and connects them to a common mother vertex

    std::vector<isobarDecayTopology> possibleDecays(const int  minI           = 0,
                                                    const int  maxI           = 2,
                                                    const int  minJ           = 0,
                                                    const int  maxJ           = 8,
                                                    const int  minL           = 0,
                                                    const int  maxL           = 6,
                                                    const int  minS           = 0,
                                                    const int  maxS           = 6,
                                                    const bool allowJpcExotic = false);  ///< constructs set of all possible decays given the final state particles and the constraints on I, J, L, and S
	  
	  const TLorentzVector& calcIsobarLzVec();  ///< (re)calculates Lorentz-vectors of all isobars in the decay from final state particles and returns Lorentz-vector of X-system
	  
	  void calcIsobarCharges();  ///< sets isobar charges as defined by final state particles
	  
	  virtual std::ostream& print(std::ostream& out) const;  ///< prints decay topology in human-readable form
	  
	  virtual std::ostream& writeGraphViz(std::ostream&      out);          ///< writes graph in GraphViz DOT format
	  virtual bool          writeGraphViz(const std::string& outFileName);  ///< writes graph in GraphViz DOT format
	  
	  static bool debug() { return _debug; }                             ///< returns debug flag
	  static void setDebug(const bool debug = true) { _debug = debug; }  ///< sets debug flag
	  

  private:

	  virtual isobarDecayTopology* doClone(const bool cloneFsParticles,
	                                       const bool cloneProdKinematics) const;  ///< helper function to use covariant return types with smart pointers; needed for public clone()

    isobarDecayTopology& constructDecay(const interactionVertexPtr&              productionVertex,
                                        const std::vector<isobarDecayVertexPtr>& isobarDecayVertices,
                                        const std::vector<particlePtr>&          fsParticles);  ///< constructs the decay graph based on final state particles and vertices
	  isobarDecayTopology& constructDecay(const interactionVertexPtr&              productionVertex,
	                                      const std::vector<interactionVertexPtr>& isobarDecayVertices,
	                                      const std::vector<particlePtr>&          fsParticles);  ///< constructs the decay graph based on final state particles and vertices
	  
    void buildIsobarVertexArray();  ///< (re)builds array of isobar decay vertices
	  
	  std::vector<isobarDecayVertexPtr> _isobarVertices;  ///< array of isobar-decay vertices excluding production vertex; ordered depth-first; this is a copy of the respective array in decayTopology
	  
    static bool _debug;  ///< if set to true, debug messages are printed

  };
	

  inline
  isobarDecayTopologyPtr
  createIsobarDecayTopology(const interactionVertexPtr&              productionVertex,
                            const std::vector<isobarDecayVertexPtr>& isobarDecayVertices,
                            const std::vector<particlePtr>&          fsParticles)
  {
    isobarDecayTopologyPtr t(new isobarDecayTopology(productionVertex,
                                                     isobarDecayVertices, fsParticles));
    return t;
  }


  inline
  isobarDecayTopologyPtr
  createIsobarDecayTopology(const interactionVertexPtr&              productionVertex,
                            const std::vector<interactionVertexPtr>& isobarDecayVertices,
                            const std::vector<particlePtr>&          fsParticles)
  {
	  isobarDecayTopologyPtr t(new isobarDecayTopology(productionVertex,
	                                                   isobarDecayVertices, fsParticles));
	  return t;
  }
	
	
	inline
	isobarDecayTopologyPtr
	createIsobarDecayTopology(const isobarDecayTopology& topo)
	{
		isobarDecayTopologyPtr t(new isobarDecayTopology(topo));
		return t;
	}
	
	
	inline
	isobarDecayTopologyPtr
	createIsobarDecayTopology(const decayTopology& topo)
	{
		isobarDecayTopologyPtr t(new isobarDecayTopology(topo));
		return t;
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
