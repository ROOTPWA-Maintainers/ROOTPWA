#line 18 "particleData.nw"
#ifndef __PARTICLEDATA_H_
#define __PARTICLEDATA_H_

#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <pputil.h>
#define getsign(x) (((x)>0)?"+":"-")

	
#line 35 "particleData.nw"
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

		void _init(std::string, double, double, int, int, int, int, int);

	public:
		
#line 63 "particleData.nw"
	particleData();
	particleData(const particleData&);
	particleData(std::string n,double m, double w, int i, int g, int j, int p, int c);
	~particleData();
	int OK() const;

#line 76 "particleData.nw"
	particleData& setName(std::string);

	particleData& setMass(double m) {
		this->_mass = m; 
		return(*this);
	}

	particleData& setWidth(double w) {
		this->_width = w; 
		return(*this);
	}

	particleData& setI(int i) {
		this->_isospin = i;
		return(*this);
	}

	particleData& setG(int g) {
		this->_gparity = g;
		return(*this);
	}

	particleData& setJ(int j) {
		this->_spin = j;
		return(*this);
	}

	particleData& setP(int p) {
		this->_parity = p;
		return(*this);
	}

	particleData& setC(int c) {
		this->_cparity = c;
		return(*this);
	}

	std::string Name() const;

	double Mass() const{
		return(this->_mass);
	}

	double Width() const {
		return(this->_width);
	}

	int I() const {
		return(this->_isospin);
	}

	int G() const {
		return(this->_gparity);
	}

	int J() const {
		return(this->_spin);
	}

	int P() const {
		return(this->_parity);
	}

	int C() const {
		return(this->_cparity);
	}

	void print() const;
	void dump() const;
	
	particleData& operator=(const particleData&);
	
	void debug(int d = 1) {
		_particle_data_debug = d;
	}
	
	friend std::ostream & operator << (std::ostream &, particleData &);


#line 53 "particleData.nw"
};

#line 29 "particleData.nw"
	
#line 275 "particleData.nw"
	class tableEntry {
	
	private:
		particleData particle;
		tableEntry  *nextparticle;
	
		void _init(const particleData& p, tableEntry *n);
	
	public:
		tableEntry(particleData p, tableEntry *n);
		tableEntry(const tableEntry& te);
		~tableEntry();
	
		tableEntry& operator=(const tableEntry& te);
	
		tableEntry* next() const;
		particleData Particle() const;
		void print() const;
		void dump() const;
	};

#line 30 "particleData.nw"
	
#line 370 "particleData.nw"
	class particleDataTable {
	
	private:
		tableEntry *head;
	
	public:
		particleDataTable( tableEntry *p = NULL );
	
		void initialize();
		void initialize(char* PDTfile);
		void insert(particleData p);
	
		particleData get(std::string _name) const;
		double mass(std::string _name) const;
		double width(std::string _name) const;
		int ListLen() const;
		char** List() const;
		void print() const;
		void dump() const;
	};
	
#line 31 "particleData.nw"
#endif


