
#include "eventMetadata.h"
#include "storage.h"

template<typename T>
rpwa::storage<T>::storage()
	: _data() { }

template<typename T>
rpwa::storage<T>::~storage() { }

template class rpwa::storage<rpwa::eventMetadata>;
