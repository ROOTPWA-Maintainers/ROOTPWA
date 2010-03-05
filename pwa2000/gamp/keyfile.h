#ifndef KEYFILE_H
#define KEYFILE_H


#include <string>
#include <list>
#include <cstdio>
#include <complex>

#include "wave.h"


class keyfile {

 public:

  keyfile();
  keyfile(const std::string& filename);
  keyfile(const keyfile&     kf);
  virtual ~keyfile();

  keyfile operator = (const keyfile&);
		
  keyfile open(const std::string& filename);
  void    rewind() { fseek(_file, 0L, SEEK_SET); }
  void    close()  { fclose(_file);              }
  FILE*   file()   { return _file;               }

  keyfile run()                             const;
  keyfile run(const event& ev)              const;
  keyfile run(const event&          ev,
	      std::complex<double>& result,
	      const bool            suppressOutput = false) const;
		
  wave addWave(const wave& w);
  wave operator [] (const int index) const;
  int  nWaves() const { return _atomicWaves.size(); }

 private:

  std::list<wave> _atomicWaves;
  std::string     _filename;
  FILE*           _file;
};


keyfile readKey (keyfile);


#endif  // KEYFILE_H
