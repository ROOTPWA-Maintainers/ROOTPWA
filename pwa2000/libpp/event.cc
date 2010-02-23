#include <iomanip>
#include <limits>

#include <event.h>


using namespace std;
	

extern particleDataTable PDGtable;


event::event()
  : _beam     (NULL),
    _target   (NULL),
    _ioversion(1)
{
}


event::event(const event& e)
{
  _final     = e._final;
  _initial   = e._initial;
  _beam      = new particle(*e._beam);
  _target    = new particle(*e._target);
  _ioversion = e._ioversion;
}


event::~event()
{
  if (_beam)
    delete _beam;
  if (_target)
    delete _target;
}


event&
event::operator = (const event& e)
{
  if (_beam)
    delete _beam;
  if (_target)
    delete _target;
  _final     = e._final;
  _initial   = e._initial;
  _beam      = new particle(*e._beam);
  _target    = new particle(*e._target);
  _ioversion = e._ioversion;
  return *this;
}


event
operator * (const lorentzTransform& L,
	    const event&            e)
{
  event r;
  r.beam  (L * e.beam());
  r.target(L * e.target());
  list<particle>::const_iterator p = e._initial.begin();
  while (p != e._initial.end()) {
    r.addinitial(L * (*p));
    ++p;
  }
  p = e._final.begin();
  while (p != e._final.end()) {
    r.addfinal(L * (*p));
    ++p;
  }
  return r;
}


event&
event::addfinal(const particle& p)
{
  _final.push_back(p);
  return *this;
}
		

event&
event::addinitial(const particle& p)
{
  _initial.push_back(p);
  return *this;
}


event&
event::erase()
{
  if (!_initial.empty())
    _initial.erase(_initial.begin(), _initial.end());
  if (!_final.empty())
    _final.erase(_final.begin(), _final.end());
  return *this;
}


event&
event::beam(const particle& p)
{
  _initial.push_front(p);
  if (_beam)
    *_beam = p;
  else
    _beam = new particle(p);
  return *this;
}
	

event&
event::target(const particle& p)
{
  _initial.push_back(p);
  if (_target)
    *_target = p;
  else
    _target = new particle(p);
  return *this;
}


int
event::OK(const double epsilon) const
{
  list<particle>::const_iterator p;
  int     q_initial = 0;
  int     q_final   = 0;
  fourVec p_initial, p_final;
  int     q_conserved, p_conserved;

  p = _initial.begin();
  while (p != _initial.end()) {
    q_initial += p->Charge();
    p_initial += p->get4P();
    ++p;
  }
  p = _final.begin();
  while (p != _final.end()) {
    q_final += p->Charge();
    p_final += p->get4P();
    ++p;
  }
  q_conserved = q_initial == q_final;
  p_conserved = (p_initial - p_final).lenSq() < epsilon;
		
  return(q_conserved && p_conserved);
}


fourVec
event::getPartPFinal(const string& name,
		     const int     charge,
		     const int     index,
		     const int     debug) const
{
  int i = 0;
  if (debug)
    cout << "Looking for " << name << charge << "[" << index << "] in event" << endl;
  list<particle>::const_iterator p = _final.begin();
  while (p != _final.end()) {
    if (debug)
      cout << "checking against " << p->Name() << p->Charge() << endl;
    if ((p->Name() == name) && (p->Charge() == charge)) {
      ++i;
      if (debug)
	cout << "found one" << endl
	     << "checking against index " << i << endl;
      if (i == index) {
	if (debug) {
	  cout << "found the right one, getting 4p" << endl
	       << "4p:" << endl;
	  p->get4P().print();
	}
	return p->get4P();
      }
    }
    ++p;
  }
  throw("PartNotFound");
  return(fourVec(0, threeVec(0, 0, 0)));
}


fourVec
event::getPartPInitial(const string& name,
		       const int     charge,
		       const int     index) const
{
  int i = 1;
  list<particle>::const_iterator p = _initial.begin();
  while (p != _initial.end())
    if ((p->Name() == name) && (i++ == index))
      return p->get4P();
  throw("PartNotFound");
  return(fourVec(0, threeVec(0, 0, 0)));
}


int
event::f_charge() const
{
  int q = 0;
  list<particle>::const_iterator p = _final.begin();
  while (p != _final.end()) {
    q += p->Charge();
    ++p;
  }
  return q;
}


double
event::f_mass() const
{
  list<particle>::const_iterator it = _final.begin();
  fourVec pX;
  fourVec p;
  while (it != _final.end()) {
    pX = it->get4P();
    p += pX;
    ++it;
  }
  // sort by mass
  double m = p.len();
  return m;
}


list<particle>
event::f_mesons() const
{
  list<particle> l;
  list<particle>::const_iterator p = _final.begin();
  while (p != _final.end()) {
    if (p->J()%2 == 0)
      l.push_back(*p);
    ++p;
  }
  return l;
}


list<particle>
event::f_baryons() const
{
  list<particle> l;
  list<particle>::const_iterator p = _final.begin();
  while (p != _final.end()) {
    if (p->J()%2 == 1)
      l.push_back(*p);
    ++p;
  }
  return l;
}


particle
event::f_particle(const string& name,
		  const int     charge,
		  const int     index) const
{
  int i = 0;
  list<particle>::const_iterator p = _final.begin();
  while (p != _final.end()) {
    if ((p->Name() == name) && (p->Charge() == charge)) {
      ++i;
      if (i == index)
	return *p;
    }
    ++p;
  }
  throw("PartNotFound");
}


void
event::set_f_mesons(const list<particle>& l)
{
  _final.clear();
  _final = l;
}


int
event::i_charge() const
{
  int q = 0;
  list<particle>::const_iterator p = _initial.begin();
  while (p != _initial.end()) {
    q += p->Charge();
    ++p;
  }
  return q;
}


list<particle>
event::i_mesons() const
{
  list<particle> l;
  list<particle>::const_iterator p = _initial.begin();
  while (p != _initial.end()) {
    if (p->J()%2 == 0)
      l.push_back(*p);
    ++p;
  }
  return l;
}


list<particle>
event::i_baryons() const
{
  list<particle> l;
  list<particle>::const_iterator p = _initial.begin();
  while (p != _initial.end()) {
    if (p->J()%2 == 1)
      l.push_back(*p);
    ++p;
  }
  return l;
}


particle
event::i_particle(const string& name,
		  const int     charge,
		  const int     index) const
{
  int i = 0;
  list<particle>::const_iterator p = _initial.begin();
  while (p != _initial.end()) {
    if ((p->Name() == name) && (p->Charge() == charge)) {
      ++i;
      if (i == index)
	return *p;
    }
    ++p;
  }
  throw("PartNotFound");
}


threeVec
event::mesonPlane() const
{
  list<particle> i = i_mesons();
  threeVec       A;
  list<particle>::const_iterator p = i.begin();
  while (p != i.end()) {
    A += p->get3P();
    ++p;
  }
  list<particle> f = f_mesons();
  threeVec       C;
  p = f.begin();
  while (p != f.end()) {
    C += p->get3P();
    ++p;
  }
  if ((A < threeVec(1e-4, 0, 0)) || (C < threeVec(1e-4, 0, 0)))
    return threeVec(0, 0, 0);
  threeVec N = A / C;
  N *= (1 / N.len());
  return N;
}


threeVec
event::baryonPlane() const
{
  list<particle> i = i_baryons();
  threeVec       B;
  list<particle>::const_iterator p = i.begin();
  while (p != i.end()) {
    B += p->get3P();
    ++p;
  }
  list<particle> f = f_baryons();
  threeVec       D;
  p = f.begin();
  while (p != f.end()) {
    D += p->get3P();
    ++p;
  }
  if ((B < threeVec(1e-4, 0, 0)) || (D < threeVec(1e-4, 0, 0)))
    return threeVec(0, 0, 0);
  threeVec N  = B / D;
  N *= (1 / N.len());
		
  return N;
}


void
event::print() const
{
  cout << "beam: ";
  _beam->get4P().print();
  cout << "target: ";
  _target->get4P().print();
  cout << "final particles: ";
  cout << endl;
  list<particle>::const_iterator p = _final.begin();
  while (p != _final.end()) {
    p->print();
    ++p;
  }
}


ostream&
operator << (ostream& os, const event& e)
{
  switch (e._ioversion) {
  case 1:
    return e.write1(os);
    break;
  case 2:
    return e.write2(os);
    break;
  default:
    throw("badIOVersion");
  }
}


ostream&
event::write1(ostream& out) const
{
  const unsigned int nmbDigits = numeric_limits<double>::digits10 + 1;
  ostringstream      s;
  s.precision(nmbDigits);
  s.setf(ios_base::scientific, ios_base::floatfield);
  s << _final.size() + 1 << endl;
  fourVec v = _beam->get4P();
  s << name2id(_beam->Name(), _beam->Charge()) << " " 
    << _beam->Charge() << " " 
    << v.x() << " " << v.y() << " " << v.z() << " " 
    << v.t() << endl;
  list<particle>::const_iterator part = _final.begin();
  while (part != _final.end()) {
    v = part->get4P();
    s << name2id(part->Name(), part->Charge()) << " "
      << part->Charge() << " "
      << v.x() << " " << v.y() << " " << v.z() << " "
      << v.t() << endl;
    ++part;
  }
  out << s.str();
  return out;
}


ostream&
event::write2(ostream& out) const
{
  const unsigned int nmbDigits = numeric_limits<double>::digits10 + 1;
  ostringstream      s;
  s.precision(nmbDigits);
  s.setf(ios_base::scientific, ios_base::floatfield);
  fourVec v = _beam->get4P();
  s << "B " << name2id(_beam->Name(), _beam->Charge()) << " " 
    << _beam->Charge() << " " 
    << v.t() << " "
    << v.x() << " " << v.y() << " " << v.z() << " " 
    << endl;
  v = _target->get4P();
  s << "T " << name2id(_target->Name(), _target->Charge()) << " " 
    << _target->Charge() << " " 
    << v.t() << " "
    << v.x() << " " << v.y() << " " << v.z() << " " 
    << endl;
  list<particle>::const_iterator part = _final.begin();
  while (part != _final.end()) {
    v = part->get4P();
    s << "F " << name2id(part->Name(), part->Charge()) << " "
      << part->Charge() << " "
      << v.t() << " "
      << v.x() << " " << v.y() << " " << v.z() << " "
      << endl;
    ++part;
  }
  s << "E" << endl;
  out << s.str();
  return out;
}


istream&
operator >> (istream& is, event& e)
{
  switch (e._ioversion) {
  case 1:
    return e.read1(is);
    break;
  case 2:
    return e.read2(is);
    break;
  default:
    throw("badIOVersion");
  }
}


istream&
event::read1(istream& is)
{
  erase();
  particle Target(PDGtable.get("p"), 1);
  Target.set4P(fourVec(Target.Mass(), threeVec(0, 0, 0)));
  target(Target);
  int nparticles = 0;
  is >> nparticles;
  for (int i = 0; i < nparticles; ++i) {
    int    ptype, q;
    double px, py, pz, t;
    is >> ptype >> q >> px >> py >> pz >> t;
    string name = id2name((Geant_ID)ptype);
    particle part(PDGtable.get(name), q);
    part.setName(name);
    part.set4P(fourVec(t, threeVec(px, py, pz)));
    if (i == 0)
      beam(part);
    else
      addfinal(part);
  }
  return is;
}


istream&
event::read2(istream& is)
{
  erase();
  char Tag = 0;
  while (!(is >> Tag).eof()) {
    int    ptype, q;
    double px, py, pz, t;
    is >> ptype >> q >> t >> px >> py >> pz;
    string name;
    name = id2name((Geant_ID)ptype);
    particle part(PDGtable.get(name), q);
    part.set4P(fourVec(t, threeVec(px, py, pz)));
    switch (Tag) {
    case 'I':
      addinitial(part);
      break;
    case 'F':
      addfinal(part);
      break;
    case 'B':
      beam(part);
      break;
    case 'T':
      target(part);
      break;
    case 'E':
      return is;
    }
  }
  return is;
}
		

event&
event::setIOVersion(int ver)
{
  if ((ver >= 1) && (ver <= 2))
    _ioversion = ver;
  else {
    cerr << "unknown io version " << ver << endl;
    throw ("UnknownIOVersion");
  }
  return *this;
}
