#ifndef RPWA_PRIMENUMBERS_HHH
#define RPWA_PRIMENUMBERS_HHH

#include <cstddef>
#include <string>
#include <vector>

class TFile;
class TTree;


namespace rpwa {

	class primeNumbers {


	  public:

		typedef unsigned long long entryType;

		static primeNumbers& instance();

		const entryType& primeNumber(const size_t& index);

		bool readCacheFile(const std::string& fileName);
		~primeNumbers();

		static bool isPrimeNumber(const entryType& number);

		static const std::string TREE_NAME;
		static const std::string BRANCH_NAME;

	  private:

		primeNumbers();

		TFile* _cacheFile;
		TTree* _cacheTree;
		const entryType _treeEntry;
		std::vector<entryType> _primeNumberCache;

		static primeNumbers* _instance;
		const static size_t _blockSize;
		const static size_t _emergencyCacheSize;

	};

}

#endif
