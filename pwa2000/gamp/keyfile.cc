#include "event.h"
#include "keyfile.h"


using namespace std;


extern particleDataTable table;
extern FILE*             keyin;
extern char*             fname;
event                    e;

extern "C++" int keylex();
extern "C++" int keyparse(std::complex<double>&);


keyfile::keyfile()
  : _filename("")
{
}


keyfile::keyfile(const string& filename)
  : _filename(filename)
{
  fname = (char*) malloc(((filename.length()) + 1) * sizeof(char));
  strcpy(fname, filename.c_str());
}


keyfile::keyfile(const keyfile& kf) 
{
  _atomicWaves = kf._atomicWaves;
  _filename    = kf._filename;
  _file        = kf._file;
}


keyfile::~keyfile()
{
}


keyfile
keyfile::operator = (const keyfile& kf)
{
  _atomicWaves = kf._atomicWaves;
  _filename    = kf._filename;
  _file        = kf._file;
  return *this;
}


keyfile
keyfile::open(const string& filename)
{
  fname = (char*) malloc(((filename.length()) + 1) * sizeof(char));
  strcpy(fname, filename.c_str());
  _file = fopen(fname, "r");
  return *this;
}


keyfile
keyfile::run() const
{
  keyin = _file;
  complex<double> result;
  keyparse(result);
  return *this;
}


keyfile
keyfile::run(const event& ev) const
{
  e     = ev; 
  keyin = _file;
  complex<double> result;
  keyparse(result);
  return *this;
}


keyfile
keyfile::run(const event& ev,
	     std::complex<double>& result) const
{
  e     = ev; 
  keyin = _file;
  keyparse(result);
  return *this;
}


wave
keyfile::addWave(const wave& w)
{
  this->_atomicWaves.push_back(w);
  return w;
}


wave
keyfile::operator[](const int index) const
{
  int i = 0;
  list<wave>::const_iterator w = _atomicWaves.begin();
  while (w != _atomicWaves.end()) {
    if (i == index)
      return *w;
    ++w;
    ++i;
  }
  return *w;
}
