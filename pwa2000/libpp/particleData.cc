#line 161 "particleData.nw"
#include <particleData.h>
	
#line 169 "particleData.nw"
	using std::string;
	using std::ifstream;
	using std::cout;
	using std::cerr;
	using std::endl;


#line 163 "particleData.nw"
	int particleData::_particle_data_debug = 0;
	int particleData::_particle_data_count = 0;
	particleDataTable PDGtable;
	
#line 180 "particleData.nw"
	void particleData::_init(string nm,double m,double w,int i,int g,int j,int p,int c) {
		_name = nm;
		_mass = m;
		_width = w;
		_isospin = i;
		_gparity = g;
		_spin = j;
		_parity = p;
		_cparity = c;
	}

#line 193 "particleData.nw"
	particleData::particleData() {
		_particle_data_count++;
		if (_particle_data_debug) { 
			cerr << "in particleData(" << this << ")::particleData()\t_particle_data_count = " << _particle_data_count << endl;
		}
		_init("",0.0,0.0,0,0,0,0,0);
	}

	particleData::particleData(const particleData& p) {
		_particle_data_count++;
		if (_particle_data_debug){ 
			cerr << "in particleData(" << this << ")::particleData(const particleData& p=" << p.Name() << ")\t_particle_data_count = " << _particle_data_count << endl;
		}
		_init(p._name,p._mass,p._width,p._isospin,p._gparity,p._spin,p._parity,p._cparity);
	}

	particleData::particleData(string n,double m, double w, int i, int g, int j, int p, int c) {
		_particle_data_count++;
		if (_particle_data_debug) { 
			cerr << "in particleData(" << this << ")::particleData(string ,double, double, int, int, int, int, int)\t_particle_data_count = " << _particle_data_count << endl;
		}
		_init(n, m, w, i, g, j, p, c);
	}

	particleData& particleData::operator=(const particleData& p) {
		if (_particle_data_debug) { 
			cerr << "in particleData(" << this << ")::operator=(const particleData&(" << &p << ")=" << p.Name() << ")\t_particle_data_count = " << _particle_data_count << endl;
		}
		_init(p._name,p._mass,p._width,p._isospin,p._gparity,p._spin,p._parity,p._cparity);
		return(*this);
	}

#line 227 "particleData.nw"
	particleData::~particleData() {
		_particle_data_count--;
		if (_particle_data_debug) { 
			cerr << "in particleData(" << this << ")=" << this->Name() << ")::~particleData()\t_particle_data_count = " << _particle_data_count << endl;
		}
	}

#line 236 "particleData.nw"
	string particleData::Name() const {
		return(this->_name);
	}

#line 242 "particleData.nw"
	particleData& particleData::setName(string nm) {
		this->_name = nm;
		return(*this);
	}

	void particleData::print() const {
		ptab();
		cout << this->_name << ":\tmass=" << this->_mass << "\twidth=" << this->_width;
		cout << "\t" << this->_isospin << getsign(this->_gparity);
		cout << "(" << this->_spin << getsign(this->_parity) << getsign(this->_cparity) << ")";
		cout << endl;
	}

	void particleData::dump() const {
		cout << this->_name << "\t" \
			<< this->_mass << "\t" \
			<< this->_width << "\t" \
			<< this->_isospin << "\t" \
			<< this->_gparity << "\t" \
			<< this->_spin << "\t" \
			<< this->_parity << "\t" 
			<< this->_cparity << endl;
	}



#line 298 "particleData.nw"
	
#line 302 "particleData.nw"
	void tableEntry::_init(const particleData& p, tableEntry *n) {
		particle = p; 
		nextparticle = n;
	}

#line 309 "particleData.nw"
	tableEntry::tableEntry(particleData p, tableEntry *n) {
		_init(p,n);
	}
	
	tableEntry::tableEntry(const tableEntry& te) {
		_init(te.particle,te.nextparticle);
	}

	tableEntry::~tableEntry() {
		;
	}
	
	tableEntry& tableEntry::operator=(const tableEntry& te) {
		_init(te.particle,te.nextparticle);
		return(*this);
	}

#line 328 "particleData.nw"
	particleData tableEntry::Particle() const {
		return (particle);
	}

#line 334 "particleData.nw"
	tableEntry* tableEntry::next() const {
		return nextparticle;
	}

	void tableEntry::print() const {
		particle.print();
	}

	void tableEntry::dump() const {
		particle.dump();
	}



#line 393 "particleData.nw"
	
#line 399 "particleData.nw"
	particleDataTable::particleDataTable( tableEntry *p ) {
		head = p;
	}

#line 405 "particleData.nw"
	void particleDataTable::initialize() {
		insert( particleData("e",0.00051,0.0,1,0,1,+1,0));
		insert( particleData("gamma",0.0,0,0,0,2,-1,-1));
		insert( particleData("pi",0.13957018,0,2,-1,0,-1,1));
		insert(particleData("pi0",0.1349766,0,2,-1,0,-1,1));
		insert(particleData("eta",0.54730,0.00000118,0,1,0,-1,1));
		insert(particleData("sigma",0.800,0.800,0,1,0,1,1));
		insert(particleData("rho(770)",0.7693,0.1502,2,1,2,-1,-1));
		insert(particleData("omega(782)",0.78257,0.00844,0,-1,2,-1,-1));
		insert(particleData("eta'(958)",0.95778,0.000202,0,1,0,-1,1));
		insert(particleData("f0(980)",0.980,0.070,0,1,0,1,1));
		insert(particleData("a0(980)",0.9931,0.071,2,-1,0,1,1));
		insert(particleData("phi(1020)",1.019417,0.004468,0,-1,2,-1,-1));
		insert(particleData("h1(1170)",1.170,0.360,0,-1,2,1,-1));
		insert(particleData("b1(1235)",1.229,0.142,2,1,2,1,-1));
		insert(particleData("a1(1269)",1.230,0.425,2,-1,2,1,1));
		insert(particleData("f2(1270)",1.2754,0.1851,0,1,4,1,1));
		insert(particleData("f1(1285)",1.2819,0.024,0,1,2,1,1));
		insert(particleData("eta(1295)",1.297,0.053,0,1,0,-1,1));
		insert(particleData("pi(1300)",1.300,0.400,2,-1,0,-1,1));
		insert(particleData("a2(1320)",1.318,0.107,2,-1,4,1,1));
		insert(particleData("f0(1370)",1.350,0.350,0,1,0,1,1));
		insert(particleData("f1(1420)",1.4263,0.0555,0,1,2,1,1));
		insert(particleData("omega(1420)",1.419,0.174,0,-1,2,-1,-1));
		insert(particleData("eta(1440)",1.420,0.060,0,1,0,-1,1));
		insert(particleData("etaL(1405)",1.405,0.056,0,1,0,-1,1));
		insert(particleData("etaH(1475)",1.475,0.081,0,1,0,-1,1));
		insert(particleData("a0(1450)",1.474,0.265,2,-1,0,1,1));
		insert(particleData("rho(1450)",1.465,0.310,2,1,2,-1,-1));
		insert(particleData("f0(1500)",1.500,0.112,0,1,0,1,1));
		insert(particleData("f1(1510)",1.510,0.073,0,1,2,1,1));
		insert(particleData("f2'(1525)",1.525,0.076,0,1,4,1,1));
		insert(particleData("omega(1650)",1.649,0.220,0,-1,2,-1,-1));
		insert(particleData("omega3(1670)",1.667,0.168,0,-1,6,-1,-1));
		insert(particleData("pi2(1670)",1.670,0.259,2,-1,4,-1,1));
		insert(particleData("phi(1680)",1.680,0.150,0,-1,2,-1,-1));
		insert(particleData("rho3(1690)",1.691,0.161,2,1,6,-1,-1));
		insert(particleData("rho(1700)",1.700,0.240,2,1,2,-1,-1));
		insert(particleData("f0(1700)",1.715,0.125,0,1,0,1,1));
		insert(particleData("pi(1800)",1.801,0.210,2,-1,0,-1,1));
		insert(particleData("phi3(1850)",1.854,0.087,0,-1,6,-1,-1));
		insert(particleData("f2(2010)",2.011,0.202,0,1,4,1,1));
		insert(particleData("a4(2040)",2.014,0.361,2,-1,8,1,1));
		insert(particleData("f4(2050)",2.034,0.222,0,1,8,1,1));
		insert(particleData("f2(2300)",2.297,0.149,0,1,4,1,1));
		insert(particleData("f2(2340)",2.339,0.319,0,1,4,1,1));

		insert(particleData("K",0.493677,0.0,1,0,0,-1,0));
		insert(particleData("K0",0.497672,0.0,1,0,0,-1,0));
		insert(particleData("Kstar(892)",0.89166,0.0508,1,0,2,-1,0));
		insert(particleData("Kstar(892)0",0.8961,0.0507,1,0,2,-1,0));
		insert(particleData("K1(1270)",1.273,0.090,1,0,2,1,0));
		insert(particleData("K1(1400)",1.402,0.174,1,0,2,1,0));
		insert(particleData("Kstar(1410)",1.414,0.232,1,0,2,-1,0));
		insert(particleData("Kstar0(1430)",1.412,0.294,1,0,0,+1,0));
		insert(particleData("Kstar2(1430)",1.4256,0.0985,1,0,4,+1,0));
		insert(particleData("Kstar2(1430)0",1.4324,0.109,1,0,4,+1,0));
		insert(particleData("Kstar(1680)",1.717,0.322,1,0,2,-1,0));
		insert(particleData("K2(1770)",1.773,0.186,1,0,4,-1,0));
		insert(particleData("Kstar3(1780)",1.776,0.159,1,0,6,-1,0));
		insert(particleData("K2(1820)",1.816,0.276,1,0,4,-1,0));
		insert(particleData("Kstar4(2045)",2.045,0.198,1,0,8,+1,0));

		insert(particleData("p",0.938272,0.0,1,0,1,+1,0));
		insert(particleData("pbar",0.938272,0.0,1,0,1,+1,0));
		insert(particleData("n",0.93956533,0.0,1,0,1,+1,0));
		insert(particleData("d",1.875612762,0.0,0,1,0,1,1));
		insert(particleData("N(1440)",1.440,0.350,1,0,1,+1,0));
		insert(particleData("N(1520)",1.520,0.120,1,0,3,-1,0));
		insert(particleData("N(1535)",1.535,0.150,1,0,1,-1,0));
		insert(particleData("N(1650)",1.650,0.150,1,0,1,-1,0));
		insert(particleData("N(1675)",1.675,0.150,1,0,5,-1,0));
		insert(particleData("N(1680)",1.680,0.130,1,0,5,+1,0));
		insert(particleData("N(1700)",1.700,0.100,1,0,3,-1,0));
		insert(particleData("N(1710)",1.710,0.100,1,0,1,+1,0));
		insert(particleData("N(1720)",1.720,0.150,1,0,3,+1,0));
		insert(particleData("N(2190)",2.190,0.450,1,0,7,-1,0));
		insert(particleData("N(2220)",2.220,0.400,1,0,9,+1,0));
		insert(particleData("N(2250)",2.250,0.400,1,0,9,-1,0));
		insert(particleData("N(2600)",2.600,0.650,1,0,11,-1,0));


		insert(particleData("Delta(1232)",1.232,0.120,3,0,3,+1,0));
		insert(particleData("Delta(1600)",1.600,0.350,3,0,3,+1,0));
		insert(particleData("Delta(1620)",1.600,0.150,3,0,1,-1,0));
		insert(particleData("Delta(1700)",1.600,0.300,3,0,3,-1,0));
		insert(particleData("Delta(1905)",1.905,0.350,3,0,5,+1,0));
		insert(particleData("Delta(1910)",1.910,0.250,3,0,1,+1,0));
		insert(particleData("Delta(1920)",1.920,0.200,3,0,3,+1,0));
		insert(particleData("Delta(1930)",1.930,0.350,3,0,5,-1,0));
		insert(particleData("Delta(1950)",1.950,0.300,3,0,7,+1,0));
		insert(particleData("Delta(2420)",2.420,0.400,3,0,11,+1,0));


		insert(particleData("lambda",1.115684,0.0,0,0,1,+1,0));
	}

#line 504 "particleData.nw"
	void particleDataTable::initialize(char* PDTfile) {
		string name;
		double mass, width;
		int isospin, gparity, spin, parity,cparity;
		ifstream ifs(PDTfile);
	
		while ( !(ifs >> name).eof() ) {
			ifs >> mass;
			ifs >> width;
			ifs >> isospin;
			ifs >> gparity;
			ifs >> spin;
			ifs >> parity;
			ifs >> cparity;
			insert(particleData(name,mass,width,isospin,gparity,spin,parity,cparity));
		}
	
	}

#line 525 "particleData.nw"
	void particleDataTable::insert(particleData p) {
		head = new tableEntry( p, head );
	}

#line 531 "particleData.nw"
	void particleDataTable::print() const {
		tableEntry *te = this->head;
		while(te!=NULL) {
			te->print();
			te = te->next();
		}
	}

#line 543 "particleData.nw"
	void particleDataTable::dump() const {
		tableEntry *te = this->head;
		while(te!=NULL) {
			te->dump();
			te = te->next();
		}
	}

#line 553 "particleData.nw"
	int particleDataTable::ListLen() const {
		int len=0;
		tableEntry *te = this->head;
		while(te!=NULL) {
			len++;
			te = te->next();
		}
		return(len);
	}
	
	char** particleDataTable::List() const {
		int particleno = 0;
		int listlen = ListLen();
		particleData p;
		char** list;
		list = (char**) malloc(listlen*sizeof(char*));
		tableEntry *te = this->head;
		while(te!=NULL) {
			p = te->Particle();
			list[particleno] = (char*) malloc( ((p.Name().length())+1)*sizeof(char));
			strcpy(list[particleno],p.Name().c_str());
//			p.Name().to_char_type(list[particleno]);
			particleno++;
			te = te->next();
		}
		return(list);
	}
 
#line 583 "particleData.nw"
	particleData particleDataTable::get(string name) const {
		tableEntry *te = this->head;
		while(te!=NULL) {
			if( (te->Particle()).Name() == name )
				return te->Particle();
			te = te->next();
		}
		return particleData(); 
		// null particle
	}

#line 596 "particleData.nw"
	double particleDataTable::mass(string name) const {
		return this->get(name).Mass();
	}

	double particleDataTable::width(string name) const {
		return this->get(name).Width();
	}

