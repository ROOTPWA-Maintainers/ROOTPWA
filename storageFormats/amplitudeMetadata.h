
#ifndef AMPLITUDEMETADATA_H
#define AMPLITUDEMETADATA_H

#include <vector>
#include <complex>

#include <TObject.h>

#include "eventMetadata.h"

class TFile;
class TTree;


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

		TTree* amplitudeTree() const { return _amplitudeTree; } // changing this tree is not allowed (it should be const, but then you can't read it...)

		Int_t Write(const char* name = 0, Int_t option = 0, Int_t bufsize = 0) { return ((const amplitudeMetadata*)this)->Write(name, option, bufsize); }
		Int_t Write(const char* name = 0, Int_t option = 0, Int_t bufsize = 0) const;

		static const std::string amplitudeLeafName;

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

	}; // class amplitudeMetadata


	inline
	std::ostream&
	operator <<(std::ostream&               out,
	            const amplitudeMetadata&    metadata)
	{
		return metadata.print(out);
	}

std::vector<std::vector<std::complex<double>>>
loadAmplitudes(const std::vector<std::string>& ampFilenames,
               const std::vector<std::string>& waveNames,
               const std::string&              eventFilename,
               const multibinBoundariesType&   otfBin,
               unsigned long                   maxNmbEvents = 0);

} // namespace rpwa

#endif
