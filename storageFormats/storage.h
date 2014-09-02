
#ifndef RPWA_STORAGE_H
#define RPWA_STORAGE_H

#include <TObject.h>
#include <TTree.h>


namespace rpwa {

	template<typename T>
	class storage : public TObject
	{

	  public:

		const T&     metadata() { return _metadata; }
		const TTree* data()     { return _data; }

		virtual Long64_t Merge(TCollection* list, Option_t* option = "") = 0;

	  private:

		T _metadata;
		TTree* _data;


	};

} // namespace rpwa

#endif
