#ifndef AMPINTEGRALMATRIXMETADATA_H
#define AMPINTEGRALMATRIXMETADATA_H

#include"eventMetadata.h"
#include"ampIntegralMatrix.h"
#include<TObject.h>
#include<map>

namespace rpwa {
	typedef std::map<std::string, std::pair<double, double> > binningMapType;
	class ampIntegralMatrixMetadata : public TObject {
			friend class ampIntegralFileWriter;

		public:
			ampIntegralMatrixMetadata();
			~ampIntegralMatrixMetadata();

	//		bool add(const ampIntegralMatrixMetadata& metaData);

			const std::vector<std::string>& getKeyFileContents() const {return _keyFileContents;};
			const std::vector<std::string>& getAmplitudeHashes() const {return _amplitudeHashes;};

			const std::string& contentHash()       const {return _contentHash;};
			const std::string& rootpwaGitHash()    const {return _rootpwaGitHash;};
			const std::string& objectBaseName()    const {return _objectBaseName;};
			const std::map<std::string, std::pair<double, double> >& binningMap()     const {return _binningMap;};
			const std::vector<std::pair<rpwa::eventMetadata, std::vector<std::pair<size_t, size_t> > > > evtMetas() const {return _evtMetas;};

			std::string recalculateHash() const;

			std::ostream& print(std::ostream& out) const;

			static const ampIntegralMatrixMetadata* readIntegralFile(TFile* inputFile,
			                                                         const std::string& objectBaseName,
			                                                         const bool& quiet = false);

			ampIntegralMatrix* getAmpIntegralMatrix() const { return _ampIntegralMatrix;};

			Int_t Write(const char* name = 0, Int_t option = 0, Int_t bufsize = 0) {return ((const ampIntegralMatrixMetadata*)this)->Write(name, option, bufsize);};
			Int_t Write(const char* name = 0, Int_t option = 0, Int_t bufsize = 0) const;

			bool mergeIntegralMatrix(const ampIntegralMatrixMetadata& second);

			bool setAmpIntegralMatrix(ampIntegralMatrix* matrix);
			bool setHash();
			bool setGitHash(const std::string &gitHash);
			bool setObjectBaseName(const std::string &baseName);
			bool addKeyFileContent(const std::string &content);
			bool setBinningMap(const std::map<std::string, std::pair<double, double> > &binningMapIn);
			bool mergeBinningMap(const std::map<std::string, std::pair<double, double> > &binnignMapIn);
			bool addAmplitudeHash(const std::string &hash);
			bool addEventMetadata(const rpwa::eventMetadata& evtMeta, size_t eventMin, size_t eventMax);

			bool check() const;
			bool hasAmplitudeHash(const std::string& hash) const;
			bool hasKeyFileContent(const std::string& content) const;
			bool writeToFile(TFile* outputFile);

			bool setAllZeroHash();
		private:
			std::string                          _contentHash;
			std::string                          _rootpwaGitHash;
			std::string                          _objectBaseName;
			std::string                          _allZeroHash;
			mutable ampIntegralMatrix*           _ampIntegralMatrix; //!

			std::vector<std::string>             _amplitudeHashes;
			std::vector<std::string>             _keyFileContents;

			std::map<std::string, std::pair<double, double> >                       _binningMap;
			std::vector<std::pair<rpwa::eventMetadata, std::vector<std::pair<size_t, size_t> > > >                _evtMetas;


			static std::pair<std::string, std::string> getObjectNames(const std::string& objectBaseName);
			ClassDef(ampIntegralMatrixMetadata, 1);
	};// class ampIntegralMatrixMetadata

	inline
	std::ostream&
	operator <<(std::ostream&                       out,
	            const ampIntegralMatrixMetadata&    metadata) {
		return metadata.print(out);
	};  // class ampIntegralMatrixMetadata 
};  // namespace rpwa
#endif//AMPINTEGRALMATRIXMETADATA_H  
