#ifndef MULTIBINTYPES_H
#define MULTIBINTYPES_H

#include <map>

namespace rpwa {

	typedef std::pair<double, double> boundaryType;
	typedef std::map<std::string, boundaryType> multibinBoundariesType;
	typedef std::map<std::string, double> multibinCenterType;

} // namespace rpwa

#endif
