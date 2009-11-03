#line 341 "../integral.nw"
#ifndef INTEGRAL_H
#define INTEGRAL_H

#include <complex>
#include <iostream>
#include <string>
#include <list>
#include <map>
#include <cstdlib>
#include <matrix.h>

using namespace std;

class integral {
    private:
        matrix<complex<double> > _sum;
        map<string, int> _index;
        int _nwaves;
        int _nevents;
        int _maxEvents;
    public:
        
#line 68 "../integral.nw"
    
#line 85 "../integral.nw"
        integral();

#line 94 "../integral.nw"
        integral(char**);

#line 101 "../integral.nw"
        integral(const integral&);

#line 108 "../integral.nw"
        ~integral();

#line 115 "../integral.nw"
        integral& operator=(const integral&);

#line 69 "../integral.nw"
    
#line 134 "../integral.nw"
        integral& integrate();

#line 143 "../integral.nw"
        integral& files(char**);

#line 152 "../integral.nw"
        integral& files(list<string>);

#line 175 "../integral.nw"
        integral& renormalize(int n);

#line 182 "../integral.nw"
        integral& max(int m);

#line 190 "../integral.nw"
        integral& events(int n);

#line 200 "../integral.nw"
        complex<double>& el(string, string);


#line 70 "../integral.nw"
    
#line 216 "../integral.nw"
        int nevents() const;
#line 222 "../integral.nw"
        list<string> files() const;

#line 231 "../integral.nw"
        char** files_c_str() const;

#line 239 "../integral.nw"
        complex<double> val(string, string);

#line 247 "../integral.nw"
        integral get(char** flist);

#line 259 "../integral.nw"
        integral get(list<string> flist);

#line 266 "../integral.nw"
        int index(string s);

#line 273 "../integral.nw"
        int index(char* s);

#line 285 "../integral.nw"
        matrix<complex<double> > mat();


#line 71 "../integral.nw"
    
#line 301 "../integral.nw"
        const integral& print(ostream& os = cout) const;

#line 309 "../integral.nw"
        const integral& print_events(ostream& os = cout) const;

#line 317 "../integral.nw"
        integral& scan(istream& is = cin);



#line 363 "../integral.nw"
};

#endif

