//-----------------------------------------------------------
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
#include <TGraph.h>

// Collaborating Class Declarations --
class TF1;

namespace rpwa {

  class pwachannel {
  public:
    pwachannel() : _C(0,0),_ps(NULL){}
    pwachannel(std::complex<double> coupling,
	       TGraph* phasespace)
      : _C(coupling), _ps(phasespace)
    {
    }

    pwachannel(const rpwa::pwachannel& ch); /// cp ctor

    // accessors
    std::complex<double> C() const {return _C;}
    std::complex<double> CsqrtPS(double m) const {return _C*sqrt(_ps->Eval(m));}

    TGraph* ps() const {return _ps;}
    double ps(double m) const {return _ps->Eval(m);}

    //modifiers
    void setCoupling(std::complex<double> c){_C=c;}


  private:
    std::complex<double> _C;
    TGraph* _ps;
  };


  class pwacomponent {
  public:
    pwacomponent(const std::string& name,
		 double m0, double gamma,
		 const std::map<std::string,pwachannel >& channels);

    virtual ~pwacomponent(){}

    std::string name() const {return _name;}
    virtual std::complex<double> val(double m) const ;
    virtual void setPar(double m0, double gamma){_m0=m0;_m02=m0*m0;_gamma=gamma;}
    unsigned int numChannels() const {return _channels.size();}
    const std::string& getChannelName(unsigned int i) const {return _channelname[i];}
    const pwachannel& getChannel(unsigned int i) const {return *_vchannels[i];}
    unsigned int numPar() const {return numChannels()*2+2;}
    void setCouplings(const double* par);
    void getCouplings(double* par);
    void setLimits(double mmin, double mmax, double gmin, double gmax)
    {_m0min=mmin;_m0max=mmax;_gammamin=gmin;_gammamax=gmax;}
    void setFixed(bool mflag=true, bool gflag=true)
    {_fixm=mflag;_fixgamma=gflag;}
    void setConstWidth(bool flag=true){_constWidth=flag;}
    void getLimits(double& mmin, double& mmax, double& gmin, double& gmax)const
    {mmin=_m0min;mmax=_m0max;gmin=_gammamin;gmax=_gammamax;}
    bool fixM() const {return _fixm;}
    bool fixGamma() const {return _fixgamma;}
    bool constWidth() const {return _constWidth;}

    double m0() const {return _m0;}
    double gamma() const {return _gamma;}
    const std::map<std::string,pwachannel >& channels()const {return _channels;}

    friend std::ostream& operator<< (std::ostream& o,const rpwa::pwacomponent& c);

  protected:
    std::string _name;
    double _m0;
    double _m02;
    double _m0min,_m0max;
    double _gamma;
    double _gammamin,_gammamax;
    bool _fixm;
    bool _fixgamma;
    bool _constWidth;
    std::map<std::string,pwachannel > _channels;
    std::vector<pwachannel*> _vchannels; // vector of channels
    std::vector<std::string> _channelname;

  };


  class pwabkg : public pwacomponent{
  public:
    pwabkg(const std::string& name,
	   double m0, double gamma,
	   const std::map<std::string,pwachannel >& channels)
      : pwacomponent(name,m0,gamma,channels){}

    virtual ~pwabkg(){}

    virtual std::complex<double> val(double m) const ;

    void setIsobars(double m1,double m2){_m1=m1;_m2=m2;}
    void getIsobars(double& m1, double& m2){m1=_m1;m2=_m2;}

  private:
    // isobar masses
    double _m1;
    double _m2;

  };


  class pwacompset {
  public:
    pwacompset():_numpar(0),_funcFsmd(NULL){}
    ~pwacompset(){}

		const std::vector<std::string>& getWaveList() const {
			return _waveList;
		}
		void setWaveList(const std::vector<std::string>& waveList) {
			_waveList = waveList;
		}



    void add(pwacomponent* comp){_comp.push_back(comp);_numpar+=comp->numPar();}
    void setFuncFsmd(TF1* funcFsmd);
		bool doMapping(); // necessary for performance. to be called after all
		                  // components have been added


    unsigned int n() const {return _comp.size();}
    unsigned int numPar() const {return _numpar;}
    
    void setPar(const double* par); // set parameters
    void getPar(double* par);       // return parameters 
    unsigned int nFreeFsmdPar() const {return _freeFsmdPar.size();}
    double getFreeFsmdPar(unsigned int i) const;
    void getFreeFsmdLimits(unsigned int i, double& lower, double& upper) const;


    const pwacomponent* operator[](unsigned int i) const {return _comp[i];}
    std::vector<std::pair<unsigned int,unsigned int> >
      getCompChannel(const std::string& wave) const;


    friend std::ostream& operator<< (std::ostream& o,const rpwa::pwacompset& cs);
    double calcFsmd(double m);
    double intensity(const std::string& wave, double m);
    double phase(const std::string& wave, double m);
    double phase(const std::string& wave1,
		 const std::string& wave2,
		 double m);
    std::complex<double> overlap(const std::string& wave1,
		 const std::string& wave2,
		 double m);
    std::complex<double> overlap(unsigned int wave1,
				 unsigned int wave2,
				 double m,
                                 const size_t idxMass = std::numeric_limits<size_t>::max());

  private:
		std::vector<std::string> _waveList;
    std::vector<pwacomponent*> _comp;
    unsigned int _numpar;
    TF1* _funcFsmd;
    std::vector<unsigned int> _freeFsmdPar; // parameters of phase space to keep floating
    // mapping for wave -> which components with which channel
    // wavelist in same order as given by wavelist
    std::vector<std::vector<std::pair<unsigned int,unsigned int> > > _compChannel;


  };




} //  end namespace



#endif // PWACOMPONENT_HH