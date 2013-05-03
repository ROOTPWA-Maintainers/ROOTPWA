#include <sstream>

#include "massDep.h"
#include "matrix.h"

#include "event.h"
#include "particle.h"

  
using namespace std;


#define MAXPRECISION(val) setprecision(numeric_limits<double>::digits10 + 1) << val


int particle::_particle_debug = 0;
int decay::_decay_debug       = 0;


particle::particle()
	: particleData(),
	  _lambda     (0),
	  _charge     (0),
	  _p          (fourVec(0, threeVec(0,0,0))),
	  _decay      (NULL),
	  _index      (0),
	  _inRestFrame(0),
    _massDep    (NULL)
{
	if (_particle_debug)
		cout << "in particle(" << this << ")::particle()" << endl;
}


particle::particle(const particle& p)
	: particleData(p)
{
  if (_particle_debug)
	  cout << "in particle(" << this << ")::particle(const particle& part("
	       << &p << ")=" << p.Name() << ")" << endl;
  _lambda = p._lambda;
  _charge = p._charge;
  _p      = p._p;
  if (p._decay)
    _decay = new decay(*(p._decay));
  else
    _decay = NULL;
  _index       = p._index;
  _helicities  = p._helicities;
  _inRestFrame = p._inRestFrame;
  if (p._massDep)
	  _massDep = p._massDep->clone();
  else
	  _massDep = NULL;
}


particle::particle(const particleData& data,
                   const int           charge)
  : particleData(data),
    _lambda     (0),
    _charge     (charge),
    _p          (fourVec(0, threeVec(0,0,0))),
    _decay      (NULL),
    _index      (0),
    _inRestFrame(0),
    _massDep    (NULL)
{
  if (_particle_debug)
    cout << "in particle(" << this << ")::particle(const particleData& data("
         << &data << ")=" << data.Name() << ", int c)" << endl;
}

  
particle::~particle()
{
  if (_particle_debug)
	  cout << "in particle((" << this << ")=" << Name() << ")::~particle()" << endl;
  delete _decay;
  delete _massDep;
}


particle&
particle::operator = (const particle& p)
{
  if (_particle_debug)
    cout << "in particle(" << this << ")::operator=(const particle& particle("
         << &p << ")=" << p.Name() << ")" << endl;
  if (this != &p) {
    particleData::operator = (p);
    _charge      = p._charge;
    _p           = p._p;
    _index       = p._index;
    _lambda      = p._lambda;
    _inRestFrame = p._inRestFrame;
    _helicities  = p._helicities;
    delete _decay;
    if (p._decay)
      _decay = new decay(*(p._decay));
    else
      _decay = NULL;
    delete _massDep;
    if (p._massDep) {
      _massDep = p._massDep->clone();
    } else
      _massDep = NULL;
  }
  return *this;
}


particle&
particle::operator *= (const lorentzTransform& L)
{
  _p *= L;
  if (_decay)
    *_decay *= L;
  return *this;
}
  

particle
operator * (const lorentzTransform& L,
      const particle&         p)
{
  particle part = p;
  part._p *= L;
  return part;
}
  

particle&
particle::setCharge(const int charge)
{
  _charge = charge;
  return *this;
}


particle&
particle::set4P(const fourVec& p4)
{
  _p = p4;
  return *this;
}


particle&
particle::set3P(const threeVec& p3)
{
  _p = fourVec(sqrt(p3.lenSq() + Mass() * Mass()), p3);
  return *this;
}


particle&
particle::Index(const int i)
{
  _index = i;
  return *this;
}

  
particle&
particle::setDecay(const decay& d)
{
  if (_decay)
    delete _decay;
  _decay = new decay(d);
  return *this;
}


particle&
particle::setMassDep(massDep* md)
{
  _massDep = md;
  return *this;
}


fourVec*
particle::get4P(particle* part,
                const int debug)
{
  fourVec* ret = NULL;
  if (   (Name()   == part->Name())
      && (Charge() == part->Charge())  
      && (Index()  == part->Index())) {
    if (debug)
	    cout << "found particle " << part->Name() <<  part->Charge() << "[" << part->Index()
	         << "]" << endl << "returning fourVec:" << ret << endl;
    ret = new fourVec(_p);
  } else if (_decay) {
	  if (debug)
		  cout << "I'm a " << Name()       <<  Charge()       << "[" << Index()       << "]" << endl
		       << "not a " << part->Name() <<  part->Charge() << "[" << part->Index() << "]" << endl
		       << "checking my children..." << endl;
	  ret = _decay->get4P(part, debug);
  }
  return ret;
}


list<int>&
particle::helicities()
{
  if (_helicities.size() == 0)
    for (int lam = -J(); lam <= J(); lam += 2)
      _helicities.push_back(lam);
  return(_helicities);
}

  
particle&
particle::addHelicity(const int lam)
{
  _helicities.push_back(lam);
  return *this;
}

  
int
particle::is(const string& name) const
{
  if (Name() == name)
    return 1;
  else
    return 0;
}


fourVec
particle::setupFrames(const int debug)
{
  fourVec ret;
  if (_decay) {
    if (debug) {
	    cout << "put " << Name() << Charge() << "[" << Index () << "]"
	         << " into helicity frame:" << endl << "momentum: ";
	    _p.print ();
    }
    
    // i should be in my parents rest frame when this is called
    // so this->_p.get3P() should be a breakup momentum
    ret = _p;
    
    fourVec tempP = _p;

    // make y perpendicular to z_old and p
    const threeVec   normal = threeVec(0, 0, 1) / tempP.V();
    lorentzTransform T;
    rotation         R;
    T.set(R.set(normal.phi(), normal.theta() - M_PI / 2, -M_PI / 2));
    lorentzTransform L = T;
    tempP *= T;

    // make z_new parallel to p
    T.set(R.set(0, signof(tempP.x()) * tempP.theta(), 0));
    matrix<double> X(4, 4);
    X = T * L;
    L = lorentzTransform(X);
    tempP *= T;

    // boost into p rest frame
    T.set(tempP);
    X = T * L;
    L = lorentzTransform(X);
    tempP *= T;

    _decay->setupFrames(L, debug);
    _inRestFrame = 1;
  } else {
	  if (debug)
		  cout << "found stable particle " << Name() << Charge() << "[" << Index() << "]" << endl;
	  ret = _p;
  }
  return ret;
}


double
particle::q() const
{
  list<particle>::const_iterator child = _decay->_children.begin();

  if (_inRestFrame)
    return ~(child->get3P());
  else {
    assert(_decay->_children.size() == 2);
    const particle& child1 = *child;
    ++child;
    const particle& child2 = *child;
    const double    Msq    = _p.lenSq();
    const double    m1sq   = child1.get4P().lenSq();
    const double    m2sq   = child2.get4P().lenSq();
    const double    lam    = lambda(Msq, m1sq, m2sq);
    return sqrt(fabs(lam / (4 * Msq)));
  }
}


double
particle::q0() const
{
  list<particle>::const_iterator child = _decay->_children.begin();

  assert(_decay->_children.size() == 2);
  const particle& child1 = *child;
  child++;
  const particle& child2 = *child;
  const double    Msq    = Mass() * Mass();
  //!!! this is wrong! q_0 should be a constant and should not depend
  //!!! on the actual but on the nominal masses of the daughters
  const double    m1sq   = child1.get4P().lenSq();
  const double    m2sq   = child2.get4P().lenSq();
  // const double    m1sq   = child1.Mass() * child1.Mass();
  // const double    m2sq   = child2.Mass() * child2.Mass();
  const double    lam    = lambda(Msq, m1sq, m2sq);
  return sqrt(fabs(lam / (4 * Msq)));  // the fabs is probably wrong
}


string
particle::sprint(const string& space) const
{
  stringstream s;
  s << Name() << chargetos(Charge()) << "[" << _index << "]";
  if (_decay) {
    s << space << "{";
    for (list<particle>::const_iterator i = _decay->_children.begin();
   i != _decay->_children.end(); ++i)
	    s << space << i->sprint(space);
    s << space << _decay->L() << space << _decay->S() << space << "}";
  }
  return s.str();
}


complex<double>
particle::breitWigner() const
{
  const double m0     = Mass();
  const double Gamma0 = Width();
  const double q      = this->q();
  const double q0     = this->q0();
  const double m      = ~(get4P());
  const int    l      = _decay->L();
  const double GammaV = Gamma0 * (m0 / m) * (q / q0) * (pow (F (l, q), 2) / pow (F (l, q0), 2));
  const complex<double> i   = complex < double >(0, 1);
  const complex<double> ret = (m0 * Gamma0) / (m0 * m0 - m * m - i * m0 * GammaV);
  return ret;
}


complex<double>
particle::decayAmp(const int lambda,
       const int debug)
{
  complex<double> a, bw;
  if (Stable())
    a = 1;
  else {
    bw = _massDep->val(*this);
    if (debug) {
      ptab();
      cout << "calculate decay amplitude for " << Name() << Charge() << "[" << Index() << "] "
           << "{bw=" << MAXPRECISION(bw) << "}" << endl;
    }
    a = _decay->amp(J(), lambda, debug) * bw;
  }
  if (debug) {
    ptab();
    cout << Name() << " decay amp = " << a << endl;
  }
  return a;
}


void
particle::print() const
{
  particleData::print();
  ptab();
  cout << "charge: " << _charge << "\tid: " << _index << endl;
  ptab();
  cout << "momentum: ";
  _p.print();
  if (_decay) {
    addtab();
    cout << "mass dependance: ";
    _massDep->print();
    cout << endl;
    _decay->print();
    subtab();
  }
}


void
particle::printFrames() const
{
  particleData::print();
  ptab();
  cout << "charge: " << _charge << "\tid: " << _index << endl;
  ptab();
  cout << "momentum: ";
  _p.print();
  if (_decay) {
    addtab();
    _decay->printFrames();
    subtab();
  }
}




decay::decay()
  : _l(0),
    _s(0)
{
}


void
decay::_init(const list<particle>& children,
             const int             l,
             const int             s,
             const double          mass)
{
  _children = children;
  _l        = l;
  _s        = s;
  _mass     = mass;
}


decay::decay(const decay& d)
{
	if (_decay_debug)
		cout << "in decay(" << this << ")::decay(const decay& d)" << endl;
	_init(d._children, d._l, d._s, d._mass);
}


decay::~decay()
{
  if (_decay_debug)
    cout << "in decay(" << this << ")::~decay()" << endl;
}


decay&
decay::addChild(const particle& p)
{
  _children.push_back(p);
  return *this;
}


decay&
decay::setL(const int l)
{
  _l = l;
  return *this;
}


decay&
decay::setS(const int s)
{
  _s = s;
  return *this;
}


decay&
decay::calculateS()
{
  list<particle>::const_iterator child = _children.begin();
  int spin       = 0;
  int numNonZero = 0;
  while (child != _children.end()) {
    if (child->J() != 0) {
      spin += child->J();
      numNonZero++;
    }
    child++;
  }
  if (numNonZero > 1) {
    cerr << "final state spin is undetermined in decay: " << endl;
    print();
    exit(1);
  }
  _s = spin;
  return *this;
}


decay&
decay::operator = (const decay& d)
{
	if (_decay_debug)
		cout << "in decay(" << this << ")::operator=(const decay& d)" << endl;
	_init(d._children, d._l, d._s, d._mass);
  return *this;
}


decay&
decay::operator *= (const lorentzTransform& L)
{
  for (list<particle>::iterator i = _children.begin(); i != _children.end(); ++i)
    *i *= L;
  return *this;
}


fourVec*
decay::get4P(particle* part,
             const int debug)
{
  if (debug)
    cout << "looking for " << part->Name() << endl;
  list<particle>::iterator child = _children.begin();
  fourVec* p = NULL;
  while (child != _children.end()) {
    if (debug)
      cout << "checking against " << child->Name() << endl;
    p = child->get4P(part, debug);
    if (p != NULL) {
	    if (debug)
		    cout << "found it! " << child->Name() << " == " << part->Name() << endl;
      break;  
    }
    ++child;
  }
  return p;
}


fourVec
decay::fill(const event& e,
            const int    debug)
{
  fourVec p;
  for (list<particle>::iterator i = _children.begin(); i != _children.end(); ++i) {
	  fourVec v;
	  if (i->Stable()) {
      if (debug)
	      cout << "Found stable particle " << i->Name() << endl;
      v = e.getPartPFinal(i->Name(), i->Charge(), i->Index(), debug);
	  } else {
      if (debug)
	      cout << "Found unstable particle " << i->Name() << endl
	           << "Calling fill for " << i->Name() << endl;
      v = i->Decay()->fill(e, debug);
	  }
	  if (debug) {
		  cout << "Setting p of " << i->Name() << " to:"<< endl;
		  v.print();
	  }
	  i->set4P(v);
	  p += v;
  }
  return p;
}


decay&
decay::setupFrames(const lorentzTransform& T,
                   const int               debug)
{
  // boost children into correct frame
  list<particle>::iterator child = _children.begin();
  while (child != _children.end() ) {
    if (debug) {
	    cout << "boosting child (" << child->Name() << ") into correct frame" << endl
	         << "p before: " << endl;
	    child->get4P().print();
    }
    *child *= T;
    if (debug) {
	    cout << "p after: " << endl;
      child->get4P().print();
    }
    ++child;
  }

  // setup decay for children
  fourVec p(0, threeVec(0, 0, 0));
  child = _children.begin();
  while (child != _children.end() ) {
	  p += child->setupFrames(debug);
	  ++child;
  }
  _mass = ~p;

  if (debug) {
	  child = _children.begin();
	  cout << "decay mass: " << _mass << endl
	       << "decay analyzer: " << child->Name() << child->Charge() << "[" << child->Index() << "]"
	       << endl << "decay angles: theta: " << child->get3P().theta() << " "
	       << "phi: " << child->get3P().phi() << endl;
  }

  return *this;
}


complex<double>
decay::expt_amp(const double b,
                const double t,
                const int    debug)
{
  assert(b >= 0);

  list<particle>::iterator child = _children.begin();
  particle&       child1      = *child;
  ++child;
  particle&       child2      = *child;
  const int       s1          = child1.J();
  const int       s2          = child2.J();
  const list<int> helicities1 = child1.helicities();
  const list<int> helicities2 = child2.helicities();
  //const particle& analyzer    = child1;

  addtab();
  if (debug) {
	  ptab();
	  cout << "My children are: " << child1.Name() << " and " << child2.Name() << endl;
	  ptab();
	  cout << "My childrens spins are: " << s1 << " and " << s2 << endl;
    ptab();
    cout << "child1's helicity list is " << helicities1.size() << " elements long" << endl;
    ptab();
    cout << "child1's helicities are: ";
    list<int>::const_iterator i = helicities1.begin();
    while (i != helicities1.end()) {
	    cout << *i << " ";
	    ++i;
    }
    cout << endl;
    ptab();
    cout << "child2's helicity list is " << helicities2.size() << " elements long" << endl;
    ptab();
    cout << "child2's helicities are: ";
    i = helicities2.begin();
    while (i != helicities2.end()) {
      cout << *i << " ";
      ++i;
    }
    cout << endl;
  }

  // calculate my amplitude
  complex<double>           amp     = 0;
  list<int>::const_iterator lambda1 = helicities1.begin();
  while (lambda1 != helicities1.end()) {
    list<int>::const_iterator lambda2 = helicities2.begin();
    while (lambda2 != helicities2.end()) {
	    if (debug) {
		    ptab();
		    cout << "lambda1: " << *lambda1 << " lambda2: " << *lambda2 << endl;
		    cout << "exp(-" << b << "|" << t << "|)" << endl;
	    }
	    complex<double> a = complex<double>(exp(-(b / 2.0) * fabs(t)), 0.0);
	    if (debug) {
		    ptab();
		    cout << "a = " << a << endl;
		    ptab();
		    cout << "calculating amp for " << child1.Name() << endl;
      }
	    a *= child1.decayAmp(*lambda1, debug);
	    if (debug) {
		    ptab();
		    cout << "a *= child1.decayAmp = " << a << endl;
		    ptab();
		    cout << "calculating amp for " << child2.Name() << endl;
	    }
	    a *= child2.decayAmp(*lambda2, debug);
      if (debug) {
	      ptab();
	      cout << "a *= child2.decayAmp = " << a << endl;
	      ptab();
	      cout << "PLUS" << endl;
      }
      amp += a;
      if (debug) {
	      ptab();
	      cout << "amp += a = " << amp << endl;
      }
      ++lambda2;
    }
    ++lambda1;
  }
  subtab();
  return amp;
}


complex<double>
decay::amp(const int j,
           const int m,
           const int debug)
{
	assert(j >= 0);
  assert(m <= j);
  
  list<particle>::iterator child = _children.begin();
  particle&       child1      = *child;
  ++child;
  particle&       child2      = *child;
  const int       s1          = child1.J();
  const int       s2          = child2.J();
  const list<int> helicities1 = child1.helicities();
  const list<int> helicities2 = child2.helicities();
  const particle& analyzer    = child1;

  addtab();
  if (debug) {
	  ptab();
	  cout << "My children are: " << child1.Name() << " and " << child2.Name() << endl;
    ptab();
    cout << "My childrens spins are: " << s1 << " and " << s2 << endl;
    ptab();
    cout << "child1's helicity list is " << helicities1.size() << " elements long" << endl;
    ptab();
    cout << "child1's helicities are: ";
    list<int>::const_iterator i = helicities1.begin();
    while (i != helicities1.end()) {
	    cout << *i << " ";
      ++i;
    }
    cout << endl;
    ptab();
    cout << "child2's helicity list is " << helicities2.size() << " elements long" << endl;
    ptab();
    cout << "child2's helicities are: ";
    i = helicities2.begin();
    while (i != helicities2.end()) {
	    cout << *i << " ";
	    ++i;
    }
    cout << endl;
  }
  
  // calculate my amplitude
  complex<double>           amp     = 0;
  list<int>::const_iterator lambda1 = helicities1.begin();
  while (lambda1 != helicities1.end()) {
	  list<int>::const_iterator lambda2 = helicities2.begin();
	  while (lambda2 != helicities2.end()) {
      const int lambda = *lambda1 - *lambda2;
      if (abs(lambda) <= j) {
	      if (debug) {
		      ptab();
		      cout << "lambda1: " << *lambda1 << " lambda2: " << *lambda2 << endl;
	      }
	      double phi   = 0;
	      double theta = 0;
	      if (_children.size() == 2) {
		      phi   = analyzer.get3P().phi();
		      theta = analyzer.get3P().theta();
	      } else if (_children.size() == 3) {  
		      // omega case, use normal to decay plane
		      const threeVec normal = child1.get3P() / child2.get3P();
		      phi   = normal.phi();
		      theta = normal.theta();
	      }
	      const double          tildeFactor    = tilde(_l);
	      const complex<double> df             = conj(D(phi, theta, 0, j, m, lambda));
	      const double          CGcoefficient1 = clebsch(_l, _s,  j,        0,    lambda, lambda);
	      const double          CGcoefficient2 = clebsch(s1, s2, _s, *lambda1, -*lambda2, lambda);
	      double                barrierFactor  = 1;
	      double                lambdaFactor   = 1;
	      if (_children.size() == 2) {
		      barrierFactor = F(_l, ~(analyzer.get3P()));
		      lambdaFactor  = 1;
	      } else if (_children.size() == 3) {  
		      // omega case
		      // instead of barrier factor use lambda factor
		      particle piZero, piPlus, piMinus;
		      child = _children.begin();
		      while (child != _children.end()) {
			      if (child->is("pi0"))
				      piZero = *child;
			      else if (child->is("pi")) {
				      switch (child->Charge()) {
				      case -1:
					      piMinus = *child;
					      break;
				      case 1:
					      piPlus = *child;
					      break;
				      default:
					      cerr << "bad child for omega: "
					           << child->Name() << endl;
				      }
			      } else {
				      cerr << "bad child for omega: "
				           << child->Name() << endl;
			      }
			      ++child;
		      }
		      if (debug) {
			      ptab();
			      cout << "calculate sqrt(lambda)" << endl;
			      ptab();
			      cout << "P_pi+: " << piPlus.get3P() << endl;
			      ptab();
			      cout << "P_pi-: " << piMinus.get3P() << endl;
			      ptab();
			      cout << "| P_pi+ X P_pi- |: " << (piPlus.get3P() / piMinus.get3P()).len() << endl;
			      ptab();
			      cout << "M(pi^+ pi^- pi^0): " << _mass << endl;
			      ptab();
			      cout << "M(pi^-): " << piMinus.get4P().lenSq() << endl;
			      ptab();
			      cout << "numerator: " << (piPlus.get3P() / piMinus.get3P()).len() << endl;
			      ptab();
			      cout << "denominator: "
			           << sqrt(3.0 / 4.0) * (pow(_mass / 3.0, 2.0) - piMinus.get4P().lenSq()) << endl;
		      }
		      lambdaFactor = (piPlus.get3P() / piMinus.get3P()).len()
			      / (sqrt(3.0 / 4.0) * (pow(_mass / 3.0, 2.0) - piMinus.get4P().lenSq()));
		      barrierFactor = 1;
	      }
	      
	      if (debug) {
		      ptab();
		      cout << "tilde(" << _l << "){=" << MAXPRECISION(tildeFactor) << "}"
		           << "D[" << j << ", " << m << ", " << lambda << "]"
		           << "(" << phi << ", " << theta << ", 0){=" << MAXPRECISION(df) << "}"
		           << "( " << _l << " 0 " << _s << " " << lambda
		           << " | " << j << " " << lambda << " ){=" << MAXPRECISION(CGcoefficient1) << "}"
		           << "( " << s1 << " " << *lambda1 << " " << s2 << " " << -*lambda2
		           << " | " << _s << " " << lambda << " ){=" << MAXPRECISION(CGcoefficient2) << "}"
		           << "F_" << _l << "(" << ~(analyzer.get3P()) << "){=" << MAXPRECISION(barrierFactor)
		           << "}sqrt(lambda)" << "{=" << MAXPRECISION(lambdaFactor) << "}"
		           << "Delta(" << _mass << ")"
		           << endl;
	      }
	      complex<double> a = tildeFactor * df * CGcoefficient1 * CGcoefficient2
		      * barrierFactor * lambdaFactor;
	      if (debug) {
		      ptab();
		      cout << "a = " << MAXPRECISION(a) << endl;
		      ptab();
		      cout << "calculating amp for " << child1.Name() << endl;
	      }
	      a *= child1.decayAmp(*lambda1, debug);
	      if (debug) {
		      ptab();
		      cout << "a *= child1.decayAmp = " << a << endl;
		      ptab();
		      cout << "calculating amp for " << child2.Name() << endl;
	      }
	      a *= child2.decayAmp(*lambda2, debug);
	      if (debug) {
		      ptab();
		      cout << "a *= child2.decayAmp = " << a << endl;
		      ptab();
		      cout << "PLUS" << endl;
	      }
	      amp += a;
	      if (debug) {
		      ptab();
		      cout << "amp += a = " << amp << endl;
	      }
      }
      ++lambda2;
	  }
	  ++lambda1;
  }
  subtab();
  return amp;
}


void
decay::print() const
{
  ptab();
  cout << "children : {" << endl;
  for (list<particle>::const_iterator i = _children.begin(); i != _children.end(); ++i)
	  i->print();
  ptab();
  cout << "L: " << _l << " S: " << _s << endl;
  ptab();
  cout << "}" << endl;
}


void
decay::printFrames() const
{
	ptab();
	cout << "i have " << _childrenInFrames.size() << " children." << endl;
	ptab();
	cout << "children in decay frame : {" << endl;
	for (list<particle>::const_iterator i = _children.begin(); i != _children.end(); ++i)
		i->printFrames();
	ptab();
	cout << "L: " << _l << " S: " << _s << endl;
	ptab();
	cout << "}" << endl;
}
