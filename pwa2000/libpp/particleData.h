#ifndef __PARTICLEDATA_H_
#define __PARTICLEDATA_H_


#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>

#include <pputil.h>


#define getsign(x) (((x)>0)?"+":"-")

	
class particleData {

private:

  static int _particle_data_count;
  static int _particle_data_debug;
  std::string _name;
  double _mass;
  double _width;
  int _isospin;
  int _gparity;
  int _spin;
  int _parity;
  int _cparity;

  void _init(const std::string&, double, double, int, int, int, int, int);

public:
		
  particleData();
  particleData(const particleData&);
  particleData(const std::string& n,double m, double w, int i, int g, int j, int p, int c);
  ~particleData();
  int OK() const;

  particleData& setName(const std::string&);

  particleData& setMass(double m)
  {
    this->_mass = m; 
    return(*this);
  }

  particleData& setWidth(double w)
  {
    this->_width = w; 
    return(*this);
  }

  particleData& setI(int i)
  {
    this->_isospin = i;
    return(*this);
  }

  particleData& setG(int g)
  {
    this->_gparity = g;
    return(*this);
  }

  particleData& setJ(int j)
  {
    this->_spin = j;
    return(*this);
  }

  particleData& setP(int p)
  {
    this->_parity = p;
    return(*this);
  }

  particleData& setC(int c)
  {
    this->_cparity = c;
    return(*this);
  }

  std::string Name() const;

  double Mass()  const { return(this->_mass);    }
  double Width() const { return(this->_width);   }
  int    I()     const { return(this->_isospin); }
  int    G()     const { return(this->_gparity); }
  int    J()     const { return(this->_spin);    }
  int    P()     const { return(this->_parity);  }
  int    C()     const { return(this->_cparity); }

  void print() const;
  void dump() const;
	
  particleData& operator = (const particleData& p);
	
  void debug(int d = 1) { _particle_data_debug = d; }
	
  friend std::ostream & operator << (std::ostream& out, particleData& p);
};


class tableEntry {
	
private:

  particleData particle;
  tableEntry*  nextparticle;
	
  void _init(const particleData& p, tableEntry* n);
	
public:

  tableEntry(const particleData& p, tableEntry* n);
  tableEntry(const tableEntry& te);
  ~tableEntry();
	
  tableEntry& operator = (const tableEntry& te);
	
  tableEntry* next() const;
  particleData Particle() const;
  void print() const;
  void dump() const;
};


class particleDataTable {
	
private:

  tableEntry* head;
	
public:

  particleDataTable(tableEntry* p = NULL);
	
  void initialize();
  void initialize(const char* PDTfile);
  void insert(const particleData& p);
	
  particleData get(const std::string& _name) const;
  double mass(const std::string& _name) const;
  double width(const std::string& _name) const;
  int ListLen() const;
  char** List() const;
  void print() const;
  void dump() const;
};
	

#endif  // __PARTICLEDATA_H_
