#line 77 "Vec.nw"
#ifndef __VEC_H_
#define __VEC_H_
#include <fstream>
#include <cmath>
#include <cassert>
#include <pputil.h>
	
#line 370 "Vec.nw"
//	class matrix<double>;
#line 84 "Vec.nw"
	
#line 94 "Vec.nw"
	class threeVec {
		private:
			
#line 106 "Vec.nw"
	double _x, _y, _z;
	void _init(double, double, double);

#line 97 "Vec.nw"
		public:
			
#line 117 "Vec.nw"
	threeVec();
	threeVec(double, double, double);
	threeVec(const threeVec&);
	~threeVec();

#line 125 "Vec.nw"
	threeVec operator+(const threeVec&) const;
	threeVec operator-(const threeVec&) const;
	threeVec operator-() const;
	double operator*(const threeVec&) const;
#line 134 "Vec.nw"
	threeVec operator/(const threeVec&) const;
#line 139 "Vec.nw"
	friend threeVec operator*(double,const threeVec&);
	friend threeVec operator*(const threeVec&,double);
#line 144 "Vec.nw"
	threeVec& operator=(const threeVec&);
	threeVec& operator+=(const threeVec&);
	threeVec& operator-=(const threeVec&);
	threeVec& operator*=(double);

#line 152 "Vec.nw"
	const threeVec& print(std::ostream& = std::cout) const;
	threeVec& scan(std::istream& is = std::cin);

#line 158 "Vec.nw"
	threeVec write(std::ostream&) const;
	threeVec read(std::istream&);


#line 165 "Vec.nw"
	double& operator[](int);
#line 169 "Vec.nw"
	double& el(int i) { return this->operator[](i); }
#line 173 "Vec.nw"
	threeVec set(double,double,double);



#line 180 "Vec.nw"
	int operator==(const threeVec&) const;
	int operator!=(const threeVec&) const;
	int operator<(const threeVec&) const;
	int operator>(const threeVec&) const;
	int operator>=(const threeVec&) const;
	int operator<=(const threeVec&) const;

#line 190 "Vec.nw"
	double x() const;
	double y() const;
	double z() const;

	double r() const;
	double theta() const;
	double cosTheta() const;
	double phi() const;

	double len() const;
	double lenSq() const;

#line 205 "Vec.nw"
	threeVec& x(double x);
	threeVec& y(double y);
	threeVec& z(double z);

	threeVec& cartesian(double x,double y,double z);
	threeVec& polar(double r,double theta,double phi);

#line 215 "Vec.nw"
	double operator~() const;


#line 99 "Vec.nw"
	};

#line 85 "Vec.nw"
	
#line 223 "Vec.nw"
	class fourVec {
		private:
			
#line 235 "Vec.nw"
	double _t;
	threeVec _V;
	void _init(double, threeVec);

#line 226 "Vec.nw"
		public:
			
#line 247 "Vec.nw"
	fourVec();
	fourVec(double, threeVec);
	fourVec(const fourVec&);
	~fourVec();

#line 255 "Vec.nw"
	fourVec operator+(const fourVec&) const;
	fourVec operator-(const fourVec&) const;
	fourVec operator-() const;
	double operator*(const fourVec&) const;
#line 265 "Vec.nw"
	threeVec operator/(const fourVec&) const;
#line 270 "Vec.nw"
	friend fourVec operator*(double,const fourVec&);
	friend fourVec operator*(const fourVec&,double);

#line 277 "Vec.nw"
	fourVec& operator=(const fourVec&);
	fourVec& operator+=(const fourVec&);
	fourVec& operator-=(const fourVec&);
	fourVec& operator*=(double);

#line 285 "Vec.nw"
	const fourVec& print(std::ostream& = std::cout) const;
	fourVec& scan(std::istream& = std::cin);

#line 291 "Vec.nw"
	fourVec write(std::ostream&) const;
	fourVec read(std::istream&);

#line 297 "Vec.nw"
	fourVec mass(double);

#line 305 "Vec.nw"
	double& operator[](int);
#line 309 "Vec.nw"
	double& el(int i) { return this->operator[](i); }
#line 313 "Vec.nw"
	fourVec set(double,double,double,double);
#line 317 "Vec.nw"
	fourVec set(double,threeVec);

#line 322 "Vec.nw"
	int operator==(const fourVec&) const;
	int operator!=(const fourVec&) const;
	int operator<(const fourVec&) const;
	int operator>(const fourVec&) const;
	int operator>=(const fourVec&) const;
	int operator<=(const fourVec&) const;

#line 332 "Vec.nw"
	threeVec V() const;
	double x() const;
	double y() const;
	double z() const;
	double t() const;

	double r() const;
	double theta() const;
	double cosTheta() const;
	double phi() const;

#line 346 "Vec.nw"
	fourVec& V(threeVec V);
	fourVec& x(double x);
	fourVec& y(double y);
	fourVec& z(double z);
	fourVec& t(double t);

	fourVec& cartesian(double x,double y,double z);
	fourVec& polar(double r,double theta,double phi);

#line 358 "Vec.nw"
	double len() const;
	double lenSq() const;


#line 365 "Vec.nw"
	double operator~() const;

#line 372 "Vec.nw"
//	fourVec operator*=(const matrix<double>&);
//		*this = *this * T;
//		return *this;
//	}

#line 228 "Vec.nw"
	};

#line 86 "Vec.nw"
	
#line 532 "Vec.nw"
	std::ostream& operator<<(std::ostream& os, const threeVec& V);
	std::istream& operator>>(std::istream& is, threeVec& V);

#line 814 "Vec.nw"
	std::ostream& operator<<(std::ostream& os, const fourVec& v);
	std::istream& operator>>(std::istream& is, fourVec& v);

#line 87 "Vec.nw"
#endif

