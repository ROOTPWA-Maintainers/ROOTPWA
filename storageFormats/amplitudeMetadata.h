#ifndef AMPLITUDEMETADATA_H
#define AMPLITUDEMETADATA_H

#include <TObject.h>
#include <TTree.h>

#include "eventMetadata.h"


namespace rpwa {

	class amplitudeMetadata : public TObject {
		friend class amplitudeFileWriter;

	  public:

		~amplitudeMetadata();

		const std::string& contentHash() const { return _contentHash; }
		const std::vector<rpwa::eventMetadata>& eventMetadata() const { return _eventMetadata; }
		const std::string& keyfileContent() const { return _keyfileContent; }
		const std::string& rootpwaGitHash() const { return _rootpwaGitHash; }
		const std::string& objectBaseName() const { return _objectBaseName; }

		std::string recalculateHash(const bool& printProgress = false) const;

		std::ostream& print(std::ostream& out) const;

		static const amplitudeMetadata* readAmplitudeFile(TFile* inputFile,
		                                                  const std::string& objectBaseName,
		                                                  const bool& quiet = false);

		Int_t Write(const char* name = 0, Int_t option = 0, Int_t bufsize = 0) { return ((const amplitudeMetadata*)this)->Write(name, option, bufsize); }
		Int_t Write(const char* name = 0, Int_t option = 0, Int_t bufsize = 0) const;

#if defined(__CINT__) || defined(__CLING__) || defined(G__DICTIONARY)
	// root needs a public default constructor
	  public:
#else
	  private:
#endif
		amplitudeMetadata();

	  private:

		void setContentHash(const std::string& contentHash) { _contentHash = contentHash; }
		void setEventMetadata(const std::vector<rpwa::eventMetadata>& eventMetadata) { _eventMetadata = eventMetadata; }
		void setKeyfileContent(const std::string& keyfileContent) { _keyfileContent = keyfileContent; }
		void setRootpwaGitHash(const std::string& rootpwaGitHash) { _rootpwaGitHash = rootpwaGitHash; }
		void setObjectBaseName(const std::string& objectBaseName) { _objectBaseName = objectBaseName; }

		static std::pair<std::string, std::string> getObjectNames(const std::string& objectBaseName);
		std::pair<std::string, std::string> getObjectNames() const { return amplitudeMetadata::getObjectNames(objectBaseName()); }

		std::string _contentHash;
		std::vector<rpwa::eventMetadata> _eventMetadata;
		std::string _keyfileContent;
		std::string _rootpwaGitHash;
		std::string _objectBaseName;

		mutable TTree* _amplitudeTree; //!

		ClassDef(amplitudeMetadata, 1);

	};

	inline
	std::ostream&
	operator <<(std::ostream&               out,
	            const amplitudeMetadata&    metadata)
	{
		return metadata.print(out);
	}


} // namespace rpwa

#endif
