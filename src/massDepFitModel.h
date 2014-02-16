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

#ifndef MASSDEPFITMODEL_HH
#define MASSDEPFITMODEL_HH

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

  class pwacomponent;

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



    void add(pwacomponent* comp);
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
 


#endif // MASSDEPFITMODEL_HH
