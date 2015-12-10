#ifndef AMPINTEGRALMATRIXMETADATA_H
#define AMPINTEGRALMATRIXMETADATA_H

#include"ampIntegralMatrix.h"
#include<TObject.h>

namespace rpwa {
	class ampIntegralMatrixMetadata : public TObject {
			friend class integralFileWriter;
		public:
			~ampIntegralMatrixMetadata();

			bool add(const ampIntegralMatrixMetadata& metaData);

			size_t nmbWaves();
			std::string getWaveName(size_t iWave);
			std::string getKeyFileContent(size_t iWave);

			const std::string& contentHash() const { return _contentHash; }
			std::string recalculateHash(const bool& printProgress = false) const;

			bool addAmplitudeMetadata(const rpwa::amplitudeMetadata& amplMeta);
			bool addEventMetadata(const rpwa::eventMetadata& evtMeta, size_t indexMin, size_t indexMax); 

			bool Write();
#if defined(__CINT__) || defined(__CLING__) || defined(G__DICTIONARY)
		// root needs a public default constructor
		public:
#else
		private:
#endif
			ampIntegralMatrixMetadata();
		private:
			std::string                          _contentHash;
			mutable ampIntegralMatrix*           _ampIntegralMatrix;

			std::vector<std::string>             _amplitudeHashes;
			std::vector<std::string>             _keyFileContents;
			
			typedef std::vector<std::pair<rpwa::eventMetadata, std::vector<std::pair<size_t, size_t> > > > evtMetasEventList;
			evtMetasEventList                    _evtMetas;
	};
};
#endif//AMPINTEGRALMATRIXMETADATA_H
