#include "Vec.h"


using namespace std;


void threeVec::_init(double x, double y, double z) {
  this->_x = x;
  this->_y = y;
  this->_z = z;
}


threeVec::threeVec() {
  this->_init(0.0, 0.0, 0.0);
}


threeVec::threeVec(double x, double y, double z) {
  this->_init(x, y, z);
}


threeVec::threeVec(const threeVec& V) {
  this->_init(V._x, V._y, V._z);
}


threeVec::~threeVec() {
  ;
}


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


threeVec operator*(double a, const threeVec& V) {
  return threeVec(a*V._x, a*V._y, a*V._z);
}


threeVec operator*(const threeVec& V, double a) {
  return ( a*V );
}


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


ostream& operator<<(ostream& os, const threeVec& V) {
  V.print(os);
  return os;
}


istream& operator>>(istream& is, threeVec& V) {
  V.scan(is);
  return is;
}


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


double&
threeVec::operator [] (const int index)
{
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


const double&
threeVec::operator [] (const int index) const
{
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


threeVec threeVec::set(double x,double y,double z) {
  this->_x = x;
  this->_y = y;
  this->_z = z;
  return *this;
}


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


void fourVec::_init(double t, threeVec V) {
  this->_t = t;
  this->_V = V;
}


fourVec::fourVec() {
  this->_init(0.0,threeVec(0.0, 0.0, 0.0));
}


fourVec::fourVec(double t, threeVec V) {
  this->_init(t,V);
}


fourVec::fourVec(const fourVec& v) {
  this->_init(v._t,v._V);
}


fourVec::~fourVec() {
  ;
}


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


fourVec operator*(double a, const fourVec& v) {
  return fourVec(a*v._t, a*v._V);
}


fourVec operator*(const fourVec& v, double a) {
  return fourVec(a*v._t, a*v._V);
}


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


ostream& operator<<(ostream& os, const fourVec& v) {
  v.print(os);
  return os;
}


istream& operator>>(istream& is, fourVec& v) {
  v.scan(is);
  return is;
}


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


double&
fourVec::operator [] (const int index)
{
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


const double&
fourVec::operator [] (const int index) const
{
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


fourVec fourVec::set(double t,double x,double y,double z) {
  this->set(t,threeVec(x,y,z));
  return *this;
}


fourVec fourVec::set(double t,threeVec V) {
  this->_t = t;
  this->_V = V;
  return *this;
}


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


fourVec fourVec::mass(double m) {
  this->_t = pow( this->_V.lenSq()+m*m, 0.5 );
  return *this;
}
