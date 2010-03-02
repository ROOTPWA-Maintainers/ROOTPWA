#include "event.h"
#include "keyfile.h"


using namespace std;


extern FILE* gKeyInFile;
extern char* gKeyInFileName;
event        gEvent;
bool         gSuppressKeyParseOutput = false;

extern "C++" int keyparse(complex<double>& result);


keyfile::keyfile()
  : _filename("")
{
}


keyfile::keyfile(const string& filename)
  : _filename(filename)
{
  gKeyInFileName = (char*) malloc(((filename.length()) + 1) * sizeof(char));
  strcpy(gKeyInFileName, filename.c_str());
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
  gKeyInFileName = (char*) malloc(((filename.length()) + 1) * sizeof(char));
  strcpy(gKeyInFileName, filename.c_str());
  _file = fopen(gKeyInFileName, "r");
  return *this;
}


keyfile
keyfile::run() const
{
  gKeyInFile = _file;
  complex<double> result;
  keyparse(result);
  return *this;
}


keyfile
keyfile::run(const event& ev) const
{
  gEvent     = ev; 
  gKeyInFile = _file;
  complex<double> result;
  keyparse(result);
  return *this;
}


keyfile
keyfile::run(const event&     ev,
	     complex<double>& result,
	     const bool       suppressOutput) const
{
  gEvent     = ev; 
  gKeyInFile = _file;
  const bool flag = gSuppressKeyParseOutput;
  gSuppressKeyParseOutput = suppressOutput;
  keyparse(result);
  gSuppressKeyParseOutput = flag;
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
