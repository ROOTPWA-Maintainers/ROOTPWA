#line 382 "Vec.nw"
	
#line 390 "Vec.nw"
#include <Vec.h>

#line 383 "Vec.nw"
	
#line 393 "Vec.nw"
	using std::ostream;
	using std::istream;
	using std::endl;
	using std::cerr;

#line 384 "Vec.nw"
	
	
#line 403 "Vec.nw"
	void threeVec::_init(double x, double y, double z) {
		this->_x = x;
		this->_y = y;
		this->_z = z;
	}
#line 411 "Vec.nw"
	threeVec::threeVec() {
		this->_init(0.0, 0.0, 0.0);
	}

	threeVec::threeVec(double x, double y, double z) {
		this->_init(x, y, z);
	}

	threeVec::threeVec(const threeVec& V) {
		this->_init(V._x, V._y, V._z);
	}

#line 426 "Vec.nw"
	threeVec::~threeVec() {
		;
	}

#line 433 "Vec.nw"
	threeVec threeVec::operator+(const threeVec& V) const {
		return threeVec(this->_x + V._x,
			this->_y + V._y,
			this->_z + V._z);
	}

	threeVec threeVec::operator-(const threeVec& V) const {
		return threeVec(this->_x - V._x,
			this->_y - V._y,
			this->_z - V._z);
	}

	threeVec threeVec::operator-() const {
		return threeVec( -this->_x, -this->_y, -this->_z);
	}

	double threeVec::operator*(const threeVec& V) const {
		return (this->_x*V._x + this->_y*V._y + this->_z*V._z);
	}

	threeVec threeVec::operator/(const threeVec& V) const {
		return threeVec(this->_y*V._z - this->_z*V._y,
			this->_z*V._x - this->_x*V._z,
			this->_x*V._y - this->_y*V._x);
	}

#line 462 "Vec.nw"
	threeVec operator*(double a, const threeVec& V) {
		return threeVec(a*V._x, a*V._y, a*V._z);
	}

	threeVec operator*(const threeVec& V, double a) {
		return ( a*V );
	}

#line 473 "Vec.nw"
	threeVec& threeVec::operator=(const threeVec& V) {
		this->_x = V._x;
		this->_y = V._y;
		this->_z = V._z;
		return *this;
	}

	threeVec& threeVec::operator+=(const threeVec& V) {
		this->_x += V._x;
		this->_y += V._y;
		this->_z += V._z;
		return *this;
	}

	threeVec& threeVec::operator-=(const threeVec& V) {
		this->_x -= V._x;
		this->_y -= V._y;
		this->_z -= V._z;
		return *this;
	}

	threeVec& threeVec::operator*=(double a) {
		this->_x *= a;
		this->_y *= a;
		this->_z *= a;
		return *this;
	}

#line 504 "Vec.nw"
	const threeVec& threeVec::print(ostream& os) const {
		os << this->_x << "\t" << this->_y << "\t" << this->_z << endl;
		return *this;
	}

	threeVec& threeVec::scan(istream& is) {
		is >> this->_x;
		is >> this->_y;
		is >> this->_z;
		return *this;
	}

#line 519 "Vec.nw"
	ostream& operator<<(ostream& os, const threeVec& V) {
		V.print(os);
		return os;
	}

	istream& operator>>(istream& is, threeVec& V) {
		V.scan(is);
		return is;
	}

#line 538 "Vec.nw"
	threeVec threeVec::write(ostream& os) const {
		os.write((char*) &(this->_x), sizeof(this->_x) );
		os.write((char*) &(this->_y), sizeof(this->_y) );
		os.write((char*) &(this->_z), sizeof(this->_z) );
		os.flush();
		return *this;
	}
	
	threeVec threeVec::read(istream& is) {
		is.read((char*) &(this->_x), sizeof(this->_x) );
		is.read((char*) &(this->_y), sizeof(this->_y) );
		is.read((char*) &(this->_z), sizeof(this->_z) );
		return *this;
	}

#line 556 "Vec.nw"
	double& threeVec::operator[](int index) {
		assert(index >=0 && index <=2);
		switch(index) {
			case 0:
				return _x;
				break;
			case 1:
				return _y;
				break;
			case 2:
				return _z;
				break;
			default:
				cerr << "threeVec: error: index " << index << " out of bounds";
				cerr << endl;
				break;
		}
		return _x;
	}


#line 578 "Vec.nw"
	threeVec threeVec::set(double x,double y,double z) {
		this->_x = x;
		this->_y = y;
		this->_z = z;
		return *this;
	}

#line 588 "Vec.nw"
	int threeVec::operator==(const threeVec& V) const {
		return (this->_x==V._x
			&& this->_y==V._y
			&& this->_z==V._z);
	}

	int threeVec::operator!=(const threeVec& V) const {
		return ( !(*this==V) );
	}

	int threeVec::operator<(const threeVec& V) const {
		return ( this->lenSq() < V.lenSq() );
	}

	int threeVec::operator>(const threeVec& V) const {
		return ( this->lenSq() > V.lenSq() );
	}

	int threeVec::operator>=(const threeVec& V) const {
		return ( this->lenSq() >= V.lenSq() );
	}

	int threeVec::operator<=(const threeVec& V) const {
		return ( this->lenSq() <= V.lenSq() );
	}

#line 617 "Vec.nw"
	double threeVec::x() const {
		return this->_x;
	}

	double threeVec::y() const {
		return this->_y;
	}

	double threeVec::z() const {
		return this->_z;
	}

	double threeVec::r() const {
		return pow( this->lenSq(),0.5 );
	}

	double threeVec::theta() const {
		return acos( this->cosTheta() );
	}

	double threeVec::cosTheta() const {
		return ( this->_z/this->r() );
	}

	double threeVec::phi() const {
		return atan2( this->_y, this->_x );
	}

	double threeVec::len() const {
		return pow(this->lenSq(),0.5);
	}

	double threeVec::lenSq() const {
		return (this->_x*this->_x
			+ this->_y*this->_y
			+ this->_z*this->_z);
	}

	double threeVec::operator~() const {
		return pow(this->lenSq(),0.5);
	}

#line 662 "Vec.nw"
	threeVec& threeVec::x(double x) {
		this->_x = x;
		return *this;
	}
	
	threeVec& threeVec::y(double y) {
		this->_y = y;
		return *this;
	}

	threeVec& threeVec::z(double z) {
		this->_z = z;
		return *this;
	}

	threeVec& threeVec::cartesian(double x,double y,double z) {
		this->_init(x, y, z);
		return *this;
	}

	threeVec& threeVec::polar(double r,double theta,double phi) {
		this->_x = r*sin(theta)*cos(phi);
		this->_y = r*sin(theta)*sin(phi);
		this->_z = r*cos(theta);
		return *this;
	}

#line 386 "Vec.nw"
	
#line 694 "Vec.nw"
	void fourVec::_init(double t, threeVec V) {
		this->_t = t;
		this->_V = V;
	}
#line 701 "Vec.nw"
	fourVec::fourVec() {
		this->_init(0.0,threeVec(0.0, 0.0, 0.0));
	}

	fourVec::fourVec(double t, threeVec V) {
		this->_init(t,V);
	}

	fourVec::fourVec(const fourVec& v) {
		this->_init(v._t,v._V);
	}

#line 716 "Vec.nw"
	fourVec::~fourVec() {
		;
	}

#line 723 "Vec.nw"
	fourVec fourVec::operator+(const fourVec& v) const {
		return fourVec( this->_t + v._t,
			this->_V + v._V );
	}

	fourVec fourVec::operator-(const fourVec& v) const {
		return fourVec( this->_t - v._t,
			this->_V - v._V );
	}

	fourVec fourVec::operator-() const {
		return fourVec( -this->_t, -this->_V );
	}

	double fourVec::operator*(const fourVec& v) const {
		return ( this->_t*v._t - this->_V*v._V );
	}

	threeVec fourVec::operator/(const fourVec& v) const {
		return (this->_V/v._V );
	}

#line 748 "Vec.nw"
	fourVec operator*(double a, const fourVec& v) {
		return fourVec(a*v._t, a*v._V);
	}

	fourVec operator*(const fourVec& v, double a) {
		return fourVec(a*v._t, a*v._V);
	}

#line 759 "Vec.nw"
	fourVec& fourVec::operator=(const fourVec& v) {
		this->_t = v._t;
		this->_V = v._V;
		return *this;
	}

	fourVec& fourVec::operator+=(const fourVec& v) {
		this->_t += v._t;
		this->_V += v._V;
		return *this;
	}

	fourVec& fourVec::operator-=(const fourVec& v) {
		this->_t -= v._t;
		this->_V -= v._V;
		return *this;
	}

	fourVec& fourVec::operator*=(double a) {
		this->_t *= a;
		this->_V *= a;
		return *this;
	}

#line 786 "Vec.nw"
	const fourVec& fourVec::print(ostream& os) const {
		os << this->_t << "\t";
		this->_V.print(os);
		return *this;
	}

	fourVec& fourVec::scan(istream& is) {
		is >> this->_t;
		this->_V.scan(is);
		return *this;
	}

#line 801 "Vec.nw"
	ostream& operator<<(ostream& os, const fourVec& v) {
		v.print(os);
		return os;
	}

	istream& operator>>(istream& is, fourVec& v) {
		v.scan(is);
		return is;
	}

#line 820 "Vec.nw"
	fourVec fourVec::write(ostream& os) const {
		os.write((char*) &(this->_t), sizeof(this->_t) );
		this->_V.write(os);
		os.flush();
		return *this;
	}
	
	fourVec fourVec::read(istream& is) {
		is.read((char*) &(this->_t), sizeof(this->_t) );
		this->_V.read(is);
		return *this;
	}

#line 836 "Vec.nw"
	double& fourVec::operator[](int index) {
		assert(index >=0 && index <=3);
		switch(index) {
			case 0:
				return _t;
				break;
			case 1:
			case 2:
			case 3:
				return _V[index-1];
				break;
			default:
				cerr << "fourVec: warning: index " << index << " out of bounds";
				cerr << endl;
				break;
		}
		return _t;
	}

#line 857 "Vec.nw"
	fourVec fourVec::set(double t,double x,double y,double z) {
		this->set(t,threeVec(x,y,z));
		return *this;
	}

#line 863 "Vec.nw"
	fourVec fourVec::set(double t,threeVec V) {
		this->_t = t;
		this->_V = V;
		return *this;
	}

#line 872 "Vec.nw"
	int fourVec::operator==(const fourVec& v) const {
		return (this->_t==v._t && (this->_V == v._V) );
	}

	int fourVec::operator!=(const fourVec& v) const {
		return ( !(*this==v) );
	}

	int fourVec::operator<(const fourVec& v) const {
		return ( this->lenSq() < v.lenSq() );
	}

	int fourVec::operator>(const fourVec& v) const {
		return ( this->lenSq() > v.lenSq() );
	}

	int fourVec::operator>=(const fourVec& v) const {
		return ( this->lenSq() >= v.lenSq() );
	}

	int fourVec::operator<=(const fourVec& v) const {
		return ( this->lenSq() <= v.lenSq() );
	}

#line 899 "Vec.nw"
	threeVec fourVec::V() const {
		return this->_V;
	}

	double fourVec::t() const {
		return this->_t;
	}

	double fourVec::x() const {
		return this->_V.x();
	}

	double fourVec::y() const {
		return this->_V.y();
	}

	double fourVec::z() const {
		return this->_V.z();
	}

#line 922 "Vec.nw"
	double fourVec::r() const {
		return pow( this->_V.lenSq(),0.5 );
	}

	double fourVec::theta() const {
		return acos( this->_V.cosTheta() );
	}

	double fourVec::cosTheta() const {
		return ( this->_V.cosTheta() );
	}

	double fourVec::phi() const {
		return ( this->_V.phi() );
	}

	double fourVec::len() const {
		return pow(this->lenSq(),0.5);
	}

	double fourVec::lenSq() const {
		return ( this->_t*this->_t - this->_V.lenSq() );
	}

	double fourVec::operator~() const {
		return pow(this->lenSq(),0.5);
	}

#line 953 "Vec.nw"
	fourVec& fourVec::V(threeVec V) {
		this->_V = V;
		return *this;
	}
	
	fourVec& fourVec::x(double x) {
		this->_V.x(x);
		return *this;
	}

	fourVec& fourVec::y(double y) {
		this->_V.y(y);
		return *this;
	}

	fourVec& fourVec::z(double z) {
		this->_V.z(z);
		return *this;
	}

	fourVec& fourVec::t(double t) {
		this->_t = t;
		return *this;
	}

	fourVec& fourVec::cartesian(double x,double y,double z) {
		this->_V.cartesian(x, y, z);
		return *this;
	}

	fourVec& fourVec::polar(double r,double theta,double phi) {
		this->_V.polar(r, theta, phi);
		return *this;
	}

#line 991 "Vec.nw"
	fourVec fourVec::mass(double m) {
		this->_t = pow( this->_V.lenSq()+m*m, 0.5 );
		return *this;
	}


