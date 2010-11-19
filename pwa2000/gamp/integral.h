#ifndef INTEGRAL_H
#define INTEGRAL_H


#include <complex>
#include <iostream>
#include <string>
#include <list>
#include <map>

#include "matrix.h"


class integral {

public:
        
	integral();
	integral(char**          files);
	integral(const integral& ni);
	virtual ~integral();

	integral& operator = (const integral& ni);

	integral& files(char**                        files);
	integral& files(const std::list<std::string>& files);
	std::list<std::string> files()       const;
	char**                 files_c_str() const;
	void weightfile(const std::string& fileName) { _weightFileName = fileName; }

	integral& integrate();
	integral& renormalize(const int n);
	integral& max(const int m);
	integral& events(const int n);
	int       nevents() const { return _nevents; }

	std::complex<double>& el(const std::string& iName,
	                         const std::string& jName)
	{ return (_sum.el(_index[iName],_index[jName])); }
	std::complex<double> val(const std::string& iName,
	                         const std::string& jName);

	integral get(char**                        files);
	integral get(const std::list<std::string>& files);

	int index(const std::string& s) { return _index[s]; }
	int index(const char*        s) { return _index[s]; }

	matrix<std::complex<double> > mat();

	const integral& print       (std::ostream& os = std::cout) const;
	const integral& print_events(std::ostream& os = std::cout) const;

	integral& scan(std::istream& is = std::cin);

private:

	matrix<std::complex<double> > _sum;
	std::map<std::string, int>    _index;
	int                           _nwaves;
	int                           _nevents;
	int                           _maxEvents;
	std::string                   _weightFileName;

};


#endif  // INTEGRAL_H
