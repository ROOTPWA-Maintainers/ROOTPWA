#include <cstdlib>

#include "sumAccumulators.hpp"
#include "integral.h"


using namespace std;
using namespace boost::accumulators;

    
integral::integral()
	: _nwaves   (0),
	  _nevents  (0),
	  _maxEvents(0)
{
}


integral::integral(char** fileList)
	: _nwaves   (0),
	  _nevents  (0),
	  _maxEvents(0)
{
	files(fileList);
	_sum = matrix<complex<double> >(_nwaves, _nwaves);    
}


integral::integral(const integral& ni)
{
	_nwaves    = ni._nwaves;
	_nevents   = ni._nevents;
	_maxEvents = ni._maxEvents;
	_index     = ni._index;
	_sum       = ni._sum;
}


integral::~integral()
{
}


integral&
integral::operator = (const integral& ni)
{
	_nwaves         = ni._nwaves;
	_nevents        = ni._nevents;
	_maxEvents      = ni._maxEvents;
	_index          = ni._index;
	_sum            = ni._sum;
	_weightFileName = ni._weightFileName;
	return *this;
}


integral&
integral::files(char** fileList)
{
	list<string> fList;
	while (*fileList) {
		fList.push_back(*fileList);
		++fileList;
	}
	files(fList);
	return *this;
}


integral&
integral::files(const list<string>& fileList)
{
	list<string>::const_iterator file = fileList.begin();
	while (file != fileList.end()) {
		_index[*file] = _nwaves++;
		file++;
	}
	_sum = matrix<complex<double> >(_nwaves, _nwaves);    
	return *this;
}


list<string>
integral::files() const
{
	list<string> fileList;
	map<string, int>::const_iterator i = _index.begin();
	while (i != _index.end() ) {
		fileList.push_back(i->first);
		++i;
	}
	return fileList;
}


char** integral::files_c_str() const
{
	char** fileList = (char**)malloc((_nwaves + 1) * sizeof(char*));
	int    index = 0;
	map<string, int>::const_iterator i = _index.begin();
	while (i != _index.end()) {
		const string fileName = i->first;
		fileList[index] = (char*)malloc((fileName.size() + 1) * sizeof(char));
		strcpy(fileList[index], fileName.c_str());
		++i;
		++index;
	}
	fileList[index] = NULL;
	return fileList;
}


integral&
integral::integrate()
{
	if (_nwaves == 0)
		throw "no waves";

	ifstream weightFile;
	bool     hasWeight = false;
	if (_weightFileName.size() != 0) {
		weightFile.open(_weightFileName.c_str());
		if (!weightFile) { 
			cerr << "error: cannot open " << _weightFileName << endl;
			throw "file not found";
		}
		hasWeight = true;
	}

	ifstream* ampfile = new ifstream [_nwaves];
	for (map<string, int>::const_iterator i = _index.begin(); i != _index.end(); ++i) {
		const string fileName  = i->first;
		const int    fileIndex = i->second;
		ampfile[fileIndex].open((fileName).c_str());
		if(!ampfile[fileIndex]) {
			cerr << "error: cannot open " << fileName << endl;
			throw "file not found";
		}
	}

	complex<double>* amps  = new complex<double>[_nwaves];
	int              nRead = 0;
	int              eof   = 0;
	accumulator_set<double,          stats<tag::sum(compensated)> > weightInt;
	accumulator_set<complex<double>, stats<tag::sum(compensated)> > sums[_nwaves][_nwaves];
	while (!eof && ((_maxEvents) ? nRead < _maxEvents : true)) {
		double w = 1;
		if (hasWeight)
			weightFile >> w;
		const double weight = 1. / w; // we have to de-weight the events!!!
		weightInt(weight);

		for (map<string, int>::const_iterator i = _index.begin(); i != _index.end(); ++i) {
			const int index = i->second;
			ampfile[index].read((char*)&amps[index], sizeof(complex<double>));
			if ((eof = ampfile[index].eof()))
				break;
		}
		if (eof)
			break;
		++_nevents;
		++nRead;

		if (!(nRead % 100))
			cerr << nRead << "\r" << flush;

		for (map<string, int>::const_iterator i = _index.begin(); i != _index.end(); ++i) {
			const int indexI = i->second;
			for (map<string, int>::const_iterator j = _index.begin(); j != _index.end(); ++j) {
				const int indexJ = j->second;
				complex<double> val = amps[indexI] * conj(amps[indexJ]);
				if (hasWeight)
					val *= weight;
				sums[indexI][indexJ](val);
			}
		}
	}
	for (int i = 0; i < _nwaves; ++i)
		for (int j = 0; j < _nwaves; ++j)
			_sum.el(i, j) = sum(sums[i][j]);

	if (hasWeight) {
		// renormalize to importance sampling weight integral:
		const double weightNorm = sum(weightInt) / (double)nRead;
		cerr << "Weight integral= " << weightNorm << endl;
		_sum *= 1. / weightNorm;
	}

	delete [] ampfile;
	delete [] amps;
	return *this;
}


integral&
integral::renormalize(const int n)
{
	_sum     = ((complex<double>) ((double) n / (double)_nevents)) * _sum;
	_nevents = n;
	return *this;
}


integral&
integral::max(const int m)
{
	_maxEvents = m; 
	return *this; 
}


integral& integral::events(const int n)
{
	_nevents = n; 
	return *this; 
}


complex<double>
integral::val(const string& iName,
              const string& jName)
{
	if (_index.find(iName) == _index.end()) {
		cerr << "error: " << iName << " not in integral" << endl;
		throw "bad wave access";
	}
	if (_index.find(jName) == _index.end()) {
		cerr << "error: " << jName << " not in integral" << endl;
		throw "bad wave access";
	}
	return el(iName, jName) / ((double)_nevents);
}


integral
integral::get(char** fileList)
{
	list<string> fList;
	while (*fileList) {
		fList.push_back(*fileList);
		++fileList;
	}
	return get(fList);
}


integral
integral::get(const list<string>& fileList)
{
	// need to check that all requested files are in list
	for (list<string>::const_iterator i = fileList.begin(); i != fileList.end(); ++i)
		if( _index.find(*i) == _index.end() ) {
			cerr << "error: " << *i << " not in integral" << endl;
			throw "bad wave access";
		}
	integral ret;
	ret.files(fileList);
	ret.events(_nevents);
	for (list<string>::const_iterator i = fileList.begin(); i != fileList.end(); ++i)
		for (list<string>::const_iterator j = fileList.begin(); j != fileList.end(); ++j)
			ret.el(*i, *j) = el(*i, *j);
	return ret;
}


matrix<complex<double> >
integral::mat()
{
	return ((complex<double>)(1.0 / ((double)_nevents))) * _sum;
}


const integral&
integral::print(ostream& os) const
{
	os << _nwaves << endl
	   << _nevents << endl
	   << _sum
	   << _index.size() << endl;
	map<string, int>::const_iterator i = _index.begin();
	while (i != _index.end()) {
		os << i->first << " " << i->second << endl;
		++i;
	}
	return *this;
}


const integral&
integral::print_events(ostream& os) const
{
	os << _nevents << endl;
	return *this;
}


integral&
integral::scan(istream& is)
{
	int    indexSize = 0, index = 0;
	string name;
	is >> _nwaves >> _nevents >> _sum >> indexSize;
	while (indexSize--) {
		is >> name >> index;
		_index[name] = index;
	}
	return *this;
}
