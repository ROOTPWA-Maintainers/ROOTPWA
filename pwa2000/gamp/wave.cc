#include "wave.h"


wave::wave(const wave& wv)
  : particle(wv)
{
  this->_m       = wv._m;
  this->_epsilon = wv._epsilon;
  this->_beam    = wv._beam;
  this->_target  = wv._target;
  this->_t       = wv._t;
  this->_b       = wv._b;
}


wave&
wave::operator = (const wave& wv)
{
  particle::operator = (wv);
  this->_m       = wv._m;
  this->_epsilon = wv._epsilon;
  this->_beam    = wv._beam;
  this->_target  = wv._target;
  this->_t       = wv._t;
  this->_b       = wv._b;
  return(*this);
}


void
wave::printFrames() const
{
  cout << "wave: " << endl;
  cout << "beam: ";
  this->_beam.print();
  cout << "target: ";
  this->_target.print();
  cout << this->I() << getsign(this->G());
  cout << "(" << this->J() << getsign(this->P()) << getsign(this->C()) << ")";
  cout << " m = " << this->_m << " eps = " << this->_epsilon << endl;
  cout << "momentum: ";
  this->get4P().print();
  cout << "decays to:" << endl;
  if(!this->Stable())
    this->Decay()->printFrames();
}


void
wave::print() const
{
  cout << "wave: " << endl;
  cout << "beam: ";
  this->_beam.print();
  cout << "target: ";
  this->_target.print();
  if (this->_b == 0.0) {
    cout << this->I() << getsign(this->G());
    cout << "(" << this->J() << getsign(this->P()) << getsign(this->C()) << ")";
    cout << " m = " << this->_m << " eps = " << this->_epsilon << endl;
  } else
    cout << "b: " << this->_b << " t: " << this->_t;
  cout << "momentum: ";
  this->get4P().print();
  cout << "decays to:" << endl;
  if(!this->Stable())
    this->Decay()->print();
}


wave
wave::fill(const event& e, int debug)
{
  decay* d;
  fourVec p;
  if (debug) {
    cout << "Filling wave" << endl;
    cout << "Setting beam p to:" << endl;
    e.beam().get4P().print();
  }
  this->_beam = e.beam().get4P();
  if (debug) {
    cout << "Setting target p to:" << endl;
    e.target().get4P().print();
  }
  this->_target = e.target().get4P();
  d = this->Decay();
  if (debug)
    cout << "Calling fill for wave: " << endl;
  p = d->fill(e, debug);
  if (debug) {
    cout << "Setting p of wave to: " << endl;
    p.print();
  }
  this->set4P(p);
  return *this;	
}


wave&
wave::setupFrames (int debug)
{
  if (debug)
    cout << "put wave into Gottfried-Jackson frame:" << endl;
  if (Stable()) {
    ;
  } else {
    lorentzTransform L, T;
    matrix <double> X(4, 4);
    rotation R;
    fourVec tempX, tempBeam, tempTarget, tempChild1;
    list<particle>::iterator child1 = this->Decay()->_children.begin();
    threeVec N;

    tempX      = this->get4P ();
    tempBeam   = this->_beam;
    tempTarget = this->_target;
    tempChild1 = (*child1).get4P();

    if (debug) {
      cout << "initially in lab:" << endl;
      cout << "tempX: ";
      tempX.print();
      cout << "tempBeam: ";
      tempBeam.print();
      cout << "tempTarget: ";
      tempTarget.print();
      cout << "tempChild1: ";
      tempChild1.print();
    }

    if (this->channel() == "t")
      // put normal to production plane along y
      N = tempBeam.V () / tempX.V ();
    else if ((this->channel() == "s") || (this->channel() == "expt"))
      // use lab y
      N = threeVec(0.0, 1.0, 0.0);
    T.set(R.set(N.phi(), N.theta() - M_PI / 2.0, -M_PI / 2.0));
    L = T;
    tempX      *= T;
    tempBeam   *= T;
    tempTarget *= T;
    tempChild1 *= T;
    if (debug) {
      cout << "put normal to PP along y:" << endl;
      cout << "tempX: ";
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
      cout << "boost to XRF:" << endl;
      cout << "tempX: ";
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
      cout << "put beam along z:" << endl;
      cout << "tempX: ";
      tempX.print();
      cout << "tempBeam: ";
      tempBeam.print();
      cout << "tempTarget: ";
      tempTarget.print();
      cout << "tempChild1: ";
      tempChild1.print();
    }

    // boost the beam and the target
    this->_beam *= L;
    this->_target *= L;

    // setupFrames of children
    this->Decay()->setupFrames(L, debug);
  }
  return *this;
}


complex<double>
wave::decayAmp(int debug)
{
  complex<double> a;
  if (Stable())
    a = complex<double>(1, 0);
  else if (this->_b != 0.0) {
    if (debug) {
      cout << "calculate decay amplitude for expt wave b=" << this->_b;
      cout << " t=" << this->_t << endl;
    }
    decay* d = this->Decay();
    a = d->expt_amp(this->_b, this->_t,debug);
  } else {
    if (debug) {
      cout << "calculate decay amplitude for wave J=" << this->J();
      cout << " m=" << this->_m << endl;
    }
    decay* d = this->Decay();
    a = d->amp(this->J(), this->_m,debug);
  }
  return a;
}


wave&
wave::operator *= (const lorentzTransform& L)
{
  this->_beam *= L;
  this->_target *= L;
  particle::operator *=(L);
  return *this;
}


string
wave::sprint(string space)
{
  string s;
  s = "J=";
  s += itos(this->J());
  s += space;
  s += "P=";
  s += itos(this->P());
  s += space;
  s += "M=";
  s += itos(this->M());
  s += space;
  s += "{";
  list<particle>::iterator c;
  for (c = this->Decay()->_children.begin(); c != this->Decay()->_children.end(); ++c) {
    s += space;
    s += c->sprint(space);
  }
  s += space;
  s += itos(this->Decay()->L());
  s += space;
  s += itos(this->Decay()->S());
  s += space;
  s += "};";
  return s;
}
