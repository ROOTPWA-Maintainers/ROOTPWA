#include <sstream>

#include "event.h"
#include "wave.h"


using namespace std;


wave::wave()
  : particle(),
    _b(0),
    _t(0)
{
}


wave::wave(const wave& wv)
  : particle(wv)
{
  _m       = wv._m;
  _epsilon = wv._epsilon;
  _beam    = wv._beam;
  _target  = wv._target;
  _t       = wv._t;
  _b       = wv._b;
}


wave::~wave()
{
}


wave&
wave::operator = (const wave& wv)
{
  particle::operator = (wv);
  _m       = wv._m;
  _epsilon = wv._epsilon;
  _beam    = wv._beam;
  _target  = wv._target;
  _t       = wv._t;
  _b       = wv._b;
  return *this;
}


wave&
wave::operator *= (const lorentzTransform& L)
{
  _beam   *= L;
  _target *= L;
  particle::operator *= (L);
  return *this;
}


wave
wave::setM(const int m)
{
  _m = m;
  return *this;
}


wave
wave::setSlope(const double b)
{
  _b = b;
  return *this;
}


wave
wave::setT(const double t)
{
  _t = t;
  return *this;
}


wave
wave::channel(const char* ch)
{
  _channel = ch;
  return *this;
}


wave
wave::fill(const event& e,
	   const int    debug)
{
  if (debug) {
    cout << "Filling wave" << endl
	 << "Setting beam p to:" << endl;
    e.beam().get4P().print();
  }
  _beam = e.beam().get4P();
  if (debug) {
    cout << "Setting target p to:" << endl;
    e.target().get4P().print();
  }
  _target = e.target().get4P();
  pwa2000::decay* d = Decay();
  if (debug)
    cout << "Calling fill for wave: " << endl;
  fourVec p = d->fill(e, debug);
  if (debug) {
    cout << "Setting p of wave to: " << endl;
    p.print();
  }
  set4P(p);
  return *this;	
}


wave&
wave::setupFrames(int debug)
{
  if (debug)
    cout << "put wave into Gottfried-Jackson frame:" << endl;
  if (Stable()) {
    ;
  } else {
    fourVec tempX      = get4P();
    fourVec tempBeam   = _beam;
    fourVec tempTarget = _target;
    fourVec tempChild1 = Decay()->_children.begin()->get4P();

    if (debug) {
      cout << "initially in lab:" << endl
	   << "tempX: ";
      tempX.print();
      cout << "tempBeam: ";
      tempBeam.print();
      cout << "tempTarget: ";
      tempTarget.print();
      cout << "tempChild1: ";
      tempChild1.print();
    }

    threeVec N;
    if (channel() == "t")
      // put normal to production plane along y
      N = tempBeam.V () / tempX.V ();
    else if ((channel() == "s") || (channel() == "expt"))
      // use lab y
      N = threeVec(0.0, 1.0, 0.0);
    lorentzTransform T;
    rotation         R;
    T.set(R.set(N.phi(), N.theta() - M_PI / 2.0, -M_PI / 2.0));
    lorentzTransform L = T;
    tempX      *= T;
    tempBeam   *= T;
    tempTarget *= T;
    tempChild1 *= T;
    if (debug) {
      cout << "put normal to PP along y:" << endl
	   << "tempX: ";
      tempX.print();
      cout << "tempBeam: ";
      tempBeam.print();
      cout << "tempTarget: ";
      tempTarget.print();
      cout << "tempChild1: ";
      tempChild1.print();
    }

    // boost to X rest frame
    T.set (tempX);
    matrix <double> X(4, 4);
    X = T * L;
    // gives error: dereferencing pointer 'X.95' does break strict-aliasing rules
    //L = *((lorentzTransform*) &X);
    // !!! BG workaround
    L = X;
    tempX      *= T;
    tempBeam   *= T;
    tempTarget *= T;
    tempChild1 *= T;
    if (debug) {
      cout << "boost to XRF:" << endl
	   << "tempX: ";
      tempX.print();
      cout << "tempBeam: ";
      tempBeam.print();
      cout << "tempTarget: ";
      tempTarget.print();
      cout << "tempChild1: ";
      tempChild1.print();
    }

    // put beam along z
    // T.set (R.set (tempBeam.V ()));
    T.set(R.set(0.0, signof(tempBeam.x()) * tempBeam.V().theta(), 0.0));
    X = T * L;
    // gives error: dereferencing pointer 'X.95' does break strict-aliasing rules
    //L = *((lorentzTransform*) &X);
    // !!! BG workaround
    L = X;
    tempX      *= T;
    tempBeam   *= T;
    tempTarget *= T;
    tempChild1 *= T;
    if (debug) {
      cout << "put beam along z:" << endl
	   << "tempX: ";
      tempX.print();
      cout << "tempBeam: ";
      tempBeam.print();
      cout << "tempTarget: ";
      tempTarget.print();
      cout << "tempChild1: ";
      tempChild1.print();
    }

    // boost the beam and the target
    _beam   *= L;
    _target *= L;

    // setupFrames of children
    Decay()->setupFrames(L, debug);
  }
  return *this;
}


complex<double>
wave::decayAmp(const int debug)
{
  complex<double> a;
  if (Stable())
    a = complex<double>(1, 0);
  else if (_b != 0.0) {
    if (debug)
      cout << "calculate decay amplitude for expt wave b=" << _b << " t=" << _t << endl;
    pwa2000::decay* d = Decay();
    a = d->expt_amp(_b, _t, debug);
  } else {
    if (debug)
      cout << "calculate decay amplitude for wave J=" << J() << " m=" << _m << endl;
    pwa2000::decay* d = Decay();
    a = d->amp(J(), _m, debug);
  }
  return a;
}


string
wave::sprint(const string& space) const
{
  stringstream s;
  s << "J=" << J() << space << "P=" << P() << space << "M=" << M() << space << "{";
  for (list<particle>::const_iterator i = Decay()->_children.begin();
       i != Decay()->_children.end(); ++i)
    s << space << i->sprint(space);
  s << space << Decay()->L() << space << Decay()->S() << space << "};";
  return s.str();
}


void
wave::print() const
{
  cout << "wave: " << endl;
  cout << "beam: ";
  _beam.print();
  cout << "target: ";
  _target.print();
  if (_b == 0.0) {
    cout << I() << getsign(G());
    cout << "(" << J() << getsign(P()) << getsign(C()) << ")";
    cout << " m = " << _m << " eps = " << _epsilon << endl;
  } else
    cout << "b: " << _b << " t: " << _t;
  cout << "momentum: ";
  get4P().print();
  cout << "decays to:" << endl;
  if(!Stable())
    Decay()->print();
}


void
wave::printFrames() const
{
  cout << "wave: " << endl;
  cout << "beam: ";
  _beam.print();
  cout << "target: ";
  _target.print();
  cout << I() << getsign(G());
  cout << "(" << J() << getsign(P()) << getsign(C()) << ")";
  cout << " m = " << _m << " eps = " << _epsilon << endl;
  cout << "momentum: ";
  get4P().print();
  cout << "decays to:" << endl;
  if(!Stable())
    Decay()->printFrames();
}
