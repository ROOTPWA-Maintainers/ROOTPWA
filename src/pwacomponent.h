//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//    (BW) Component of mass dependent fit
//      
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

#ifndef PWACOMPONENT_HH
#define PWACOMPONENT_HH

// Base Class Headers ----------------


// Collaborating Class Headers -------
#include <ostream>
#include <string>
#include <map>
#include <vector>
#include <complex>

// Collaborating Class Declarations --

namespace rpwa {

  class pwacomponent {
  public:
    pwacomponent(const std::string& name,
		 double m0, double gamma,
		 const std::map<std::string,std::complex<double> >& channels);
    
    virtual ~pwacomponent(){}
    
    std::string name() const {return _name;}
    virtual std::complex<double> val(double m) const ;
    void setPar(double m0, double gamma){_m0=m0;_m02=m0*m0;_gamma=gamma;}
    unsigned int numChannels() const {return _channels.size();}
    unsigned int numPar() const {return numChannels()*2+2;}
    void setCouplings(const double* par);
    void getCouplings(double* par);
    void setLimits(double mmin, double mmax, double gmin, double gmax)
    {_m0min=mmin;_m0max=mmax;_gammamin=gmin;_gammamax=gmax;}
    void setFixed(bool mflag=true, bool gflag=true)
    {_fixm=mflag;_fixgamma=gflag;}
    void getLimits(double& mmin, double& mmax, double& gmin, double& gmax)const
    {mmin=_m0min;mmax=_m0max;gmin=_gammamin;gmax=_gammamax;}
    bool fixM() const {return _fixm;}
    bool fixGamma() const {return _fixgamma;}

    double m0() const {return _m0;}
    double gamma() const {return _gamma;}
    const std::map<std::string,std::complex<double> >& channels()const {return _channels;}

    friend std::ostream& operator<< (std::ostream& o,const rpwa::pwacomponent& c);

  private:
    std::string _name;
    double _m0;
    double _m02;
    double _m0min,_m0max;
    double _gamma;
    double _gammamin,_gammamax;
    bool _fixm;
    bool _fixgamma;
    std::map<std::string,std::complex<double> > _channels;

  };



  class pwacompset {
  public:
    pwacompset():_numpar(0){}
    ~pwacompset(){}

    void add(const pwacomponent& comp){_comp.push_back(comp);_numpar+=comp.numPar();}
    
    unsigned int n() const {return _comp.size();}
    unsigned int numPar() const {return _numpar;}
    
    std::vector<std::string> wavelist()const;

    void setPar(const double* par); // set parameters
    void getPar(double* par);       // return parameters 

    const pwacomponent& operator[](unsigned int i) const {return _comp[i];}

    friend std::ostream& operator<< (std::ostream& o,const rpwa::pwacompset& cs);
    
    double intensity(std::string wave, double m);
    
    
  private:
    std::vector<pwacomponent> _comp;
    unsigned int _numpar;

  };




} //  end namespace
 


#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
