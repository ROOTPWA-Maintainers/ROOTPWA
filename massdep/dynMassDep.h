
///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2009 Sebastian Neubert
//
//    This file is part of rootpwa
//
//    rootpwa is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    rootpwa is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with rootpwa.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////

/** @brief Dynamical Mass Dependency (Propagator)
 */

#ifndef DYNMASSDEP_HH
#define DYNMASSDEP_HH

#include <complex>
#include "TGraph.h"
#include "absMassDep.h"
#include "absDecayChannel.h"
#include "mcPhaseSpace.h"


namespace rpwa {
  
  const double gChargedPionMass = 0.13957018;  // charged pion rest mass [GeV/c^2] [PDG08
  
  class dynMassDep : public absMassDep
  {
  public:
    dynMassDep(double M, double width, 
	       unsigned int nparticles, double* masses);
    ~dynMassDep() {if(ps!=NULL)delete ps;ps=NULL;}
    
    cd val(double m){return val(m,_channel);}
    cd val(double m, unsigned int i);
    cd val_static(double m){return val_static(m,_channel);}
    cd val_static(double m, unsigned int i);
    cd val_nodisperse(double m){return val_nodisperse(m,_channel);}
    cd val_nodisperse(double m, unsigned int i);
    cd val_bnl(double m, unsigned int i);
    mcPhaseSpace* phasespace() const {return ps;}
    double get_rho0(unsigned int i)const;
    double get_rho(double m, unsigned int i)const {return ps->eval(m,i);}
    void store_ms(double maxM, unsigned int steps=400);
    double get_ms(double s, unsigned int i) const ;
    double calc_ms(double s, unsigned int i) const ;
    double disperse(double* x, double* par) const;
int l(unsigned int i) const {if(_channels.size()==0)return 0; else return _channels[i]->l();}
    int l()const {return l(_channel);}

    void addDecayChannel(absDecayChannel* ch);
    void setFixedChannel(unsigned int channel){_channel=channel;}

    TGraph* graph_ms(unsigned int i)const {return _graphms[i];}

  private:
    double mS;
    double mM;
    double mWidth;
    mcPhaseSpace* ps;
    std::vector<TGraph*> _graphms;
    double rho0;
    std::vector<absDecayChannel*> _channels;
    unsigned int _channel; // selected channel;
    
  };
  
  
} // end namespace

#endif
