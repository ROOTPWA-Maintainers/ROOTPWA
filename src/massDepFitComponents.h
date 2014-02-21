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

#ifndef MASSDEPFITCOMPONENTS_HH
#define MASSDEPFITCOMPONENTS_HH

// Base Class Headers ----------------


// Collaborating Class Headers -------
#include <ostream>
#include <string>
#include <map>
#include <vector>
#include <complex>

#include <boost/multi_array.hpp>

#include <TGraph.h>

// Collaborating Class Declarations --
namespace libconfig {
	class Setting;
}
class TF1;

namespace rpwa {


	class pwachannel {

	public:

    pwachannel(const std::string& waveName,
               std::complex<double> coupling,
	       const std::vector<double>& massBinCenters,
	       const std::vector<double>& phaseSpace)
      : _waveName(waveName), _C(coupling), _phaseSpace(phaseSpace)
    {
        _ps = new TGraph(_phaseSpace.size(), &(massBinCenters[0]), &(_phaseSpace[0]));
    }

    pwachannel(const rpwa::pwachannel& ch); /// cp ctor

    // accessors
    std::complex<double> C() const {return _C;}
    std::complex<double> CsqrtPS(double m) const {return _C*sqrt(_ps->Eval(m));}
    
		const std::string& getWaveName() const { return _waveName; }

    TGraph* ps() const {return _ps;}
    double ps(double m) const {return _ps->Eval(m);}

    //modifiers
    void setCoupling(std::complex<double> c){_C=c;}

	private:

		std::string _waveName;
    std::complex<double> _C;
    std::vector<double> _phaseSpace;
    TGraph* _ps;

	};


	class massDepFitComponent {

	public:

		massDepFitComponent(const std::string& name,
		                    const size_t nrParameters);
		virtual ~massDepFitComponent() {};

		const std::string& getName() const { return _name; }

		virtual bool init(const libconfig::Setting* configComponent,
		                  const std::vector<double>& massBinCenters,
		                  const std::map<std::string, size_t>& waveIndices,
		                  const boost::multi_array<double, 2>& phaseSpaceIntegrals,
		                  const bool debug);

		const size_t getNrChannels() const { return _channels.size(); }
		const std::vector<pwachannel>& getChannels() const { return _channels; }
		const pwachannel& getChannel(const size_t i) const { return _channels[i]; } 
		const std::string& getChannelWaveName(const size_t i) const { return _channels[i].getWaveName(); }

		void getCouplings(double* par) const;
		void setCouplings(const double* par);

		size_t getNrParameters() const { return _nrParameters; }
		virtual void getParameters(double* par) const = 0;
		virtual void setParameters(const double* par) = 0;

		virtual double getParameter(const size_t idx) const = 0;
		virtual void setParameter(const size_t idx, const double value) = 0;
		virtual bool getParameterFixed(const size_t idx) const = 0;
		virtual void setParameterFixed(const size_t idx, const bool fixed) = 0;
		virtual std::pair<double, double> getParameterLimits(const size_t idx) const = 0;
		virtual void setParameterLimits(const size_t idx, const std::pair<double, double>& limits) = 0;

		virtual std::complex<double> val(const double m) const = 0;

		virtual std::ostream& print(std::ostream& out) const;

	private:

		const std::string _name;

		std::vector<pwachannel> _channels;

		const size_t _nrParameters;

	};


	class pwacomponent : public massDepFitComponent {

	public:

		pwacomponent(const std::string& name);
    
		virtual bool init(const libconfig::Setting* configComponent,
		                  const std::vector<double>& massBinCenters,
		                  const std::map<std::string, size_t>& waveIndices,
		                  const boost::multi_array<double, 2>& phaseSpaceIntegrals,
		                  const bool debug);

		virtual void getParameters(double* par) const;
		virtual void setParameters(const double* par);

		virtual double getParameter(const size_t idx) const;
		virtual void setParameter(const size_t idx, const double value);
		virtual bool getParameterFixed(const size_t idx) const;
		virtual void setParameterFixed(const size_t idx, const bool fixed);
		virtual std::pair<double, double> getParameterLimits(const size_t idx) const ;
		virtual void setParameterLimits(const size_t idx, const std::pair<double, double>& limits);

		virtual std::complex<double> val(const double m) const;

		std::ostream& print(std::ostream& out) const;

	private:

		double _m0;
		double _m02;
		double _m0min;
		double _m0max;
		bool _fixm;

		double _gamma;
		double _gammamin;
		double _gammamax;
		bool _fixgamma;

		bool _constWidth;

	};


	class pwabkg : public massDepFitComponent {

	public:

		pwabkg(const std::string& name);
    
		virtual bool init(const libconfig::Setting* configComponent,
		                  const std::vector<double>& massBinCenters,
		                  const std::map<std::string, size_t>& waveIndices,
		                  const boost::multi_array<double, 2>& phaseSpaceIntegrals,
		                  const bool debug);

		virtual void getParameters(double* par) const;
		virtual void setParameters(const double* par);

		virtual double getParameter(const size_t idx) const;
		virtual void setParameter(const size_t idx, const double value);
		virtual bool getParameterFixed(const size_t idx) const;
		virtual void setParameterFixed(const size_t idx, const bool fixed);
		virtual std::pair<double, double> getParameterLimits(const size_t idx) const;
		virtual void setParameterLimits(const size_t idx, const std::pair<double, double>& limits);

		virtual std::complex<double> val(const double m) const;

  private:

		double _m0;
		double _m0min;
		double _m0max;
		bool _fixm;

		double _gamma;
		double _gammamin;
		double _gammamax;
		bool _fixgamma;

		double _m1;
		double _m2;
    
	};


} // end namespace rpwa

inline std::ostream& operator<< (std::ostream& out, const rpwa::massDepFitComponent& component) {
	return component.print(out);
}

#endif // MASSDEPFITCOMPONENTS_HH
