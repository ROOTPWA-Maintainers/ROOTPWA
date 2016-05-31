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

			const std::vector<std::string>& getKeyFileContents() const { return _keyFileContents; }
			const std::vector<std::string>& getAmplitudeHashes() const { return _amplitudeHashes; }

			const std::string& contentHash()       const { return _contentHash; }
			const std::string& rootpwaGitHash()    const { return _rootpwaGitHash; }
			const std::map<std::string, std::pair<double, double> >& binningMap() const { return _binningMap; }
			const std::vector<rpwa::eventMetadata>& evtMetas() const { return _evtMetas; }

			std::string recalculateHash() const;

			std::ostream& print(std::ostream& out) const;

			static const ampIntegralMatrixMetadata* readIntegralFile(TFile* inputFile,
			                                                         const bool& quiet = false);

			ampIntegralMatrix* getAmpIntegralMatrix() const { return _ampIntegralMatrix; }

			Int_t Write(const char* name = 0, Int_t option = 0, Int_t bufsize = 0) {
				return ((const ampIntegralMatrixMetadata*)this)->Write(name, option, bufsize);
			}
			Int_t Write(const char* name = 0, Int_t option = 0, Int_t bufsize = 0) const;

			bool mergeIntegralMatrix(const ampIntegralMatrixMetadata& second);

			bool setAmpIntegralMatrix(ampIntegralMatrix* matrix);
			bool setHash();
			void setGitHash(const std::string& gitHash) { _rootpwaGitHash = gitHash; }
			void setBinningMap(const std::map<std::string, std::pair<double, double> >& binningMap) { _binningMap = binningMap; }
			bool addKeyFileContent(const std::string& content);
			bool mergeBinningMap(const std::map<std::string, std::pair<double, double> >& binnignMapIn);
			bool addAmplitudeHash(const std::string& hash);
			bool addEventMetadata(const rpwa::eventMetadata& evtMeta);

			bool check() const;
			bool hasAmplitudeHash(const std::string& hash) const;
			bool hasKeyFileContent(const std::string& content) const;
			bool writeToFile(TFile* outputFile);

			bool setAllZeroHash();

			static const std::string objectNameInFile;
		private:
			std::string                          _contentHash;
			std::string                          _rootpwaGitHash;
			std::string                          _allZeroHash;
			mutable ampIntegralMatrix*           _ampIntegralMatrix; //!

			std::vector<std::string>             _amplitudeHashes;
			std::vector<std::string>             _keyFileContents;

			std::map<std::string, std::pair<double, double> > _binningMap;
			std::vector<rpwa::eventMetadata> _evtMetas;

			ClassDef(ampIntegralMatrixMetadata, 2);
	};  // class ampIntegralMatrixMetadata

	inline
	std::ostream&
	operator<< (std::ostream&                       out,
	            const ampIntegralMatrixMetadata&    metadata)
	{
		return metadata.print(out);
	}

}  // namespace rpwa
#endif  //AMPINTEGRALMATRIXMETADATA_H
