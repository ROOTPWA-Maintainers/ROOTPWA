#ifndef NONINTERACTIONVERTEX_H
#define NONINTERACTIONVERTEX_H

#include "productionVertex.h"

namespace rpwa {

	class nonInteractionVertex;
	typedef boost::shared_ptr<nonInteractionVertex> nonInteractionVertexPtr;

	class nonInteractionVertex : public productionVertex {

	public:

		nonInteractionVertex(const particlePtr& XParticle);
		virtual ~nonInteractionVertex() { }

		virtual const std::vector<TLorentzVector>& referenceLzVec() const { return XParticle()->lzVecs(); }
		inline const particlePtr&                  XParticle     () const { return outParticles()[0];     }

		virtual bool initKinematicsData(const TClonesArray& prodKinPartNames);  ///< initializes input data
		virtual bool readKinematicsData(const std::vector<std::vector<TVector3> >& prodKinMomenta);    ///< add input data event

		virtual bool revertMomenta();  ///< resets momenta to the values of last event read

		virtual void setXFlavorQN() { }  ///< does not do anything, as incoming and outgoing particles are the same

		virtual bool addInParticle (const particlePtr&);  ///< disabled; all incoming particles have to be specified at construction
		virtual bool addOutParticle(const particlePtr&);  ///< disabled; all outgoing particles have to be specified at construction

		virtual inline std::vector<particlePtr>&       inParticles()       { return outParticles(); }  ///< returns array of outgoing particles
		virtual inline const std::vector<particlePtr>& inParticles() const { return outParticles(); }  ///< returns array of outgoing particles

		virtual std::ostream& print(std::ostream& out) const;  ///< prints vertex parameters in human-readable form

		virtual std::string name() const { return "nonInteractionVertex"; }  ///< returns label used in graph visualization, reporting, and key file

		static bool debug() { return _debug; }
		static void setDebug(bool debug = true) { _debug = debug; }

	private:

		std::vector<TVector3> _XParticleCache;

		static bool _debug;

	};

	inline
	nonInteractionVertexPtr
	createNonInteractionVertex(const particlePtr& XParticle)
	{
		nonInteractionVertexPtr vert(new nonInteractionVertex(XParticle));
		return vert;
	}


	inline
	std::ostream&
	operator <<(std::ostream&                out,
	            const nonInteractionVertex& vert)
	{
		return vert.print(out);
	}

}

#endif
