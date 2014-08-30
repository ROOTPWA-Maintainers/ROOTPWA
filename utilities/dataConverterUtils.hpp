
#ifndef DATACONVERTERUTILS_HPP
#define DATACONVERTERUTILS_HPP

#include <map>

#include <TObject.h>

namespace rpwa {

	class dataMetadata : public TObject {

	  public:

		const std::string& getUserString() { return _userString; }
		void setUserString(const std::string& userString) { _userString = userString; }

		const std::string& getContentHash() { return _contentHash; }

		void setBinningVariableLables(std::vector<std::string> labels)
		{
			for(unsigned int i = 0; i < labels.size(); ++i) {
				_binningMap[labels[i]] = std::pair<double, double>(0., 0.);
			}
		}

		void setBinningVariableRange(const std::string& label, const std::pair<double, double>& range)
		{
			_binningMap[label] = range;
		}

		void setBinningMap(const std::map<std::string, std::pair<double, double> > binningMap)
		{
			_binningMap = binningMap;
		}

	  private:

		std::string _userString;
		std::string _contentHash;

		std::map<std::string, std::pair<double, double> > _binningMap;

	};

}

#endif
