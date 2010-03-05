#ifndef __VEC_H_
#define __VEC_H_


#include <fstream>
#include <cmath>
#include <cassert>

#include <pputil.h>
	

class threeVec {

private:
			
  double _x, _y, _z;

  void _init(double, double, double);

public:
			
  threeVec();
  threeVec(double, double, double);
  threeVec(const threeVec&);
  ~threeVec();

  threeVec operator + (const threeVec&) const;
  threeVec operator - (const threeVec&) const;
  threeVec operator - () const;
  double   operator * (const threeVec&) const;
  threeVec operator / (const threeVec&) const;
  friend threeVec operator * (double, const threeVec&);
  friend threeVec operator * (const threeVec&, double);
  threeVec& operator =  (const threeVec&);
  threeVec& operator += (const threeVec&);
  threeVec& operator -= (const threeVec&);
  threeVec& operator *= (double);

  const threeVec& print(std::ostream& = std::cout) const;
  threeVec& scan(std::istream& is = std::cin);

  threeVec write(std::ostream&) const;
  threeVec read(std::istream&);

  double&       operator [] (const int i);
  const double& operator [] (const int i) const;
  double&       el(const int i)       { return operator[](i); }
  const double& el(const int i) const { return operator[](i); }
  threeVec set(double, double, double);

  int operator == (const threeVec&) const;
  int operator != (const threeVec&) const;
  int operator <  (const threeVec&) const;
  int operator >  (const threeVec&) const;
  int operator >= (const threeVec&) const;
  int operator <= (const threeVec&) const;

  double x() const;
  double y() const;
  double z() const;

  double r() const;
  double theta() const;
  double cosTheta() const;
  double phi() const;

  double len() const;
  double lenSq() const;

  threeVec& x(double x);
  threeVec& y(double y);
  threeVec& z(double z);

  threeVec& cartesian(double x, double y, double z);
  threeVec& polar(double r, double theta, double phi);

  double operator~() const;

};


std::ostream& operator << (std::ostream& os, const threeVec& V);
std::istream& operator >> (std::istream& is, threeVec& V);


class fourVec {

private:
			
  double   _t;
  threeVec _V;

  void     _init(double, threeVec);

public:
			
  fourVec();
  fourVec(double, threeVec);
  fourVec(const fourVec&);
  ~fourVec();

  fourVec  operator + (const fourVec&) const;
  fourVec  operator - (const fourVec&) const;
  fourVec  operator - () const;
  double   operator * (const fourVec&) const;
  threeVec operator / (const fourVec&) const;
  friend fourVec operator*(double,const fourVec&);
  friend fourVec operator*(const fourVec&,double);
  fourVec& operator =  (const fourVec&);
  fourVec& operator += (const fourVec&);
  fourVec& operator -= (const fourVec&);
  fourVec& operator *= (double);

  const fourVec& print(std::ostream& = std::cout) const;
  fourVec& scan(std::istream& = std::cin);

  fourVec write(std::ostream&) const;
  fourVec read(std::istream&);

  fourVec mass(double);

  double&       operator [] (const int);
  const double& operator [] (const int) const;
  double&       el(const int i)       { return this->operator[](i); }
  const double& el(const int i) const { return this->operator[](i); }
  fourVec set(double, double, double, double);
  fourVec set(double, threeVec);

  int operator == (const fourVec&) const;
  int operator != (const fourVec&) const;
  int operator <  (const fourVec&) const;
  int operator >  (const fourVec&) const;
  int operator >= (const fourVec&) const;
  int operator <= (const fourVec&) const;

  threeVec V() const;
  double x() const;
  double y() const;
  double z() const;
  double t() const;

  double r() const;
  double theta() const;
  double cosTheta() const;
  double phi() const;

  fourVec& V(threeVec V);
  fourVec& x(double x);
  fourVec& y(double y);
  fourVec& z(double z);
  fourVec& t(double t);

  fourVec& cartesian(double x, double y, double z);
  fourVec& polar(double r, double theta, double phi);

  double len() const;
  double lenSq() const;


  double operator~() const;

};


std::ostream& operator << (std::ostream& os, const fourVec& v);
std::istream& operator >> (std::istream& is, fourVec& v);


#endif  // __VEC_H_
