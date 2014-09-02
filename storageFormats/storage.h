
#ifndef RPWA_STORAGE_H
#define RPWA_STORAGE_H

#include <TObject.h>
#include <TTree.h>


namespace rpwa {

	template<typename T>
	class storage : public TObject
	{
		friend class eventFileWriter;

	  public:

		storage();
		virtual ~storage();

		const T&      metadata() const { return _metadata; }
		      T&      metadata()       { return _metadata; }
		      TTree* data()       { return &_data; }
		const TTree* data() const { return &_data; }

		virtual std::string hash(const bool& printProgress = false) = 0;

		virtual Long64_t Merge(TCollection* list, Option_t* option = "") = 0;

	  protected:

		void setMetadata(const T& metadata) { _metadata = metadata; }

		T _metadata;
		TTree _data;

		ClassDef(storage, 1);

	};

} // namespace rpwa

#endif
