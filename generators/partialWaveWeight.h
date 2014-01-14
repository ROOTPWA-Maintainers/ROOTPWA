/** @addtogroup generators
 * @{
 */

#ifndef TPWWEIGHT_HH
#define TPWWEIGHT_HH

#include <string>
#include <vector>
#include <map>
#include <complex>
#include <event.h>
#include <Tgamp.h>
#include "productionAmp.h"

#include <integral.h>


/** @brief Calculates the partial wave weight for a given event
 *  @author Sebastian Neubert TUM (original author)
 *
 */

namespace rpwa {

	class partialWaveWeight {

	  public:

		partialWaveWeight(): m_hasInt(false){}
		~partialWaveWeight(){}

		/** @brief Add a partial wave to model
		 *
		 *  keyfilename = name of gamp keyfile describing this wave
		 *  amp         = production amplitude for this wave
		 *  vectori     = index of production vector
		 *                (different indices are added incoherently)
		 */
		void addWave(const std::string& keyfilename,
		             rpwa::productionAmp* amp,
		             const std::complex<double>& branching,
		             unsigned int vectori=0);


		/** @brief Load normalization integrals for a mass bin
		 *
		 *  several integrals can be loaded to cover a larger range in mass
		 */
		void loadIntegrals(const std::string& filename, double mass);

		/** @brief Calculate weight for an event
		*/
		double weight(event& e);

		/** @brief get production amplitude
		 *
		 *  iv: production vector
		 *  iw: wave
		 */
		std::complex<double> prodAmp(unsigned int iv,
		                             unsigned int iw,
		                             event& e);


	  private:

		// Private Data Members ------------
		std::vector<std::vector<std::string> > m_waves;

		// Production amplitudes
		std::vector<std::vector<std::complex<double> > > m_branchings;
		std::vector<std::vector<rpwa::productionAmp*> > m_amps;

		std::vector<Tgamp> m_gamp;
		std::map<std::string,double> m_relphase;
		std::map<double,integral> m_normInt;

		bool m_hasInt;

			// Private Methods -----------------

	};

}

#endif
/* @} **/

