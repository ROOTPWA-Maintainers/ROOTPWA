#include <iomanip>
#include <limits>

#include <event.h>


using namespace std;
	

event::event()
{
  this->_beam      = NULL;
  this->_target    = NULL;
  this->_ioversion = 1;
}


event::~event()
{
  if (this->_beam)
    delete this->_beam;
  if (this->_target)
    delete this->_target;
}


event::event(const event& e)
{
  this->_final     = e._final;
  this->_initial   = e._initial;
  this->_beam      = new particle(*e._beam);
  this->_target    = new particle(*e._target);
  this->_ioversion = e._ioversion;
}


event&
event::operator = (const event& e)
{
  if (this->_beam)
    delete this->_beam;
  if (this->_target)
    delete this->_target;
  this->_final     = e._final;
  this->_initial   = e._initial;
  this->_beam      = new particle(*e._beam);
  this->_target    = new particle(*e._target);
  this->_ioversion = e._ioversion;
  return *this;
}


event&
event::addfinal(const particle& p)
{
  this->_final.push_back(p);
  return *this;
}
		

event&
event::addinitial(const particle& p)
{
  this->_initial.push_back(p);
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
  this->_initial.push_front(p);
  if (this->_beam)
    *(this->_beam) = p;
  else
    this->_beam = new particle(p);
  return *this;
}
	

event&
event::target(const particle& p)
{
  this->_initial.push_back(p);
  if (this->_target)
    *(this->_target) = p;
  else
    this->_target = new particle(p);
  return *this;
}


int
event::OK(double epsilon = 1e-6) const
{
  list<particle>::const_iterator p;
  int q_initial = 0;
  int q_final   = 0;
  fourVec p_initial, p_final;
  int q_conserved, p_conserved;

  p = this->_initial.begin();
  while (p != this->_initial.end()) {
    q_initial += p->Charge();
    p_initial += p->get4P();
    ++p;
  }
  p = this->_final.begin();
  while (p != this->_final.end()) {
    q_final += p->Charge();
    p_final += p->get4P();
    ++p;
  }
  q_conserved = q_initial == q_final;
  p_conserved = (p_initial - p_final).lenSq() < epsilon;
		
  return(q_conserved && p_conserved);
}


particle
event::beam() const
{
  return *(this->_beam);
}
	

particle
event::target() const
{
  return *(this->_target);
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
  list<particle>::const_iterator p = this->_final.begin();
  while (p != this->_final.end()) {
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
  list<particle>::const_iterator p = this->_initial.begin();
  while (p != this->_initial.end())
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


void
event::set_f_mesons(const list<particle>& l)
{
  _final.clear();
  _final = l;
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


list<particle>
event::f_particles() const
{
  return _final;
}


particle
event::f_particle(const string& name,
		  const int     charge,
		  const int     index) const
{
  int i = 0;
  list<particle>::const_iterator p = this->_final.begin();
  while (p != this->_final.end()) {
    if ((p->Name() == name) && (p->Charge() == charge)) {
      ++i;
      if (i == index)
	return *p;
    }
    ++p;
  }
  throw("PartNotFound");
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


list<particle>
event::i_particles() const
{
  return _initial;
}


particle
event::i_particle(const string& name,
		  const int     charge,
		  const int     index) const
{
  int i = 0;
  list<particle>::const_iterator p = this->_initial.begin();
  while (p != this->_initial.end()) {
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
  threeVec A, C, N;
  list<particle> i, f;
  list<particle>::const_iterator ip, fp;
		
  i  = this->i_mesons();
  ip = i.begin();
  while (ip != i.end()) {
    A += ip->get3P();
    ++ip;
  }

  f  = this->f_mesons();
  fp = f.begin();
  while (fp != f.end()) {
    C += fp->get3P();
    ++fp;
  }
		
  if ((A < threeVec(1e-4, 0, 0)) || (C < threeVec(1e-4, 0, 0)))
    return threeVec(0, 0, 0);
		
  N  = A / C;
  N *= (1 / N.len());
		
  return N;
}


threeVec
event::baryonPlane() const
{
  threeVec B, D, N;
  list<particle> i, f;
  list<particle>::const_iterator ip, fp;
		
  i  = this->i_baryons();
  ip = i.begin();
  while (ip != i.end()) {
    B += ip->get3P();
    ++ip;
  }

  f  = this->f_baryons();
  fp = f.begin();
  while (fp != f.end()) {
    D += fp->get3P();
    ++fp;
  }
		
  if ((B < threeVec(1e-4, 0, 0)) || (D < threeVec(1e-4, 0, 0)))
    return threeVec(0, 0, 0);
		
  N  = B / D;
  N *= (1 / N.len());
		
  return N;
}


void
event::print() const
{
  cout << "beam: ";
  this->_beam->get4P().print();
  cout << "target: ";
  this->_target->get4P().print();
  cout << "final particles: ";
  cout << endl;
  list<particle>::const_iterator p = this->_final.begin();
  while (p != this->_final.end()) {
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


ostream&
event::write2(ostream& out) const
{
  const unsigned int nmbDigits = numeric_limits<double>::digits10 + 1;
  ostringstream      s;
  s.precision(nmbDigits);
  s.setf(ios_base::scientific, ios_base::floatfield);
  fourVec v = this->beam().get4P();
  s << "B " << name2id(this->_beam->Name(), this->_beam->Charge()) << " " 
    << this->_beam->Charge() << " " 
    << v.t() << " "
    << v.x() << " " << v.y() << " " << v.z() << " " 
    << endl;
  v = this->target().get4P();
  s << "T " << name2id(this->_target->Name(), this->_target->Charge()) << " " 
    << this->_target->Charge() << " " 
    << v.t() << " "
    << v.x() << " " << v.y() << " " << v.z() << " " 
    << endl;
  list<particle>::const_iterator part = this->_final.begin();
  while (part != this->_final.end()) {
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
event::read2(istream& is)
{
  int ptype, q;
  double px, py, pz, t;
  string name;
		
  char Tag = 0;
  this->erase();
  while (!(is >> Tag).eof()) {
    switch (Tag) {
    case 'I':
      {
	is >> ptype >> q >> t >> px >> py >> pz;
	name = id2name((Geant_ID)ptype);
	particle part(PDGtable.get(name), q);
	part.set4P(fourVec(t, threeVec(px, py, pz)));
	this->addinitial(part);
      }
      break;
    case 'F':
      {
	is >> ptype >> q >> t >> px >> py >> pz;
	name = id2name((Geant_ID)ptype);
	particle part(PDGtable.get(name), q);
	part.set4P(fourVec(t, threeVec(px, py, pz)));
	this->addfinal(part);
      }
      break;
    case 'B':
      {
	is >> ptype >> q >> t >> px >> py >> pz;
	name = id2name((Geant_ID)ptype);
	particle part(PDGtable.get(name), q);
	part.set4P(fourVec(t, threeVec(px, py, pz)));
	this->beam(part);
      }
      break;
    case 'T':
      {
	is >> ptype >> q >> t >> px >> py >> pz;
	name = id2name((Geant_ID)ptype);
	particle part(PDGtable.get(name), q);
	part.set4P(fourVec(t, threeVec(px, py, pz)));
	this->target(part);
      }
      break;
    case 'E':
      return is;
    }
  }
  return is;
}


ostream&
event::write1(ostream& out) const
{
  const unsigned int nmbDigits = numeric_limits<double>::digits10 + 1;
  ostringstream      s;
  s.precision(nmbDigits);
  s.setf(ios_base::scientific, ios_base::floatfield);
  s << this->_final.size() + 1 << endl;
  fourVec v = this->beam().get4P();
  s << name2id(this->_beam->Name(), this->_beam->Charge()) << " " 
    << this->_beam->Charge() << " " 
    << v.x() << " " << v.y() << " " << v.z() << " " 
    << v.t() << endl;
  list<particle>::const_iterator part = this->_final.begin();
  while (part != this->_final.end()) {
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


istream&
event::read1(istream& is)
{
  int nparticles = 0;
  int ptype, q;
  double px, py, pz, t;
  string name;

  this->erase();

  particle Target(PDGtable.get("p"), 1);
  Target.set4P(fourVec(Target.Mass(), threeVec(0, 0, 0)));
  this->target(Target);
		
  is >> nparticles;
  for (int i = 0; i < nparticles; ++i) {
    is >> ptype >> q >> px >> py >> pz >> t;
    name = id2name((Geant_ID)ptype);
    if (i == 0) {
      particle Beam(PDGtable.get(name), q);
      Beam.setName(name);
      Beam.set4P(fourVec(t, threeVec(px, py, pz)));
      this->beam(Beam);
    } else {
      particle part(PDGtable.get(name), q);
      part.setName(name);
      part.set4P(fourVec(t, threeVec(px, py, pz)));
      this->addfinal(part);
    }
  }
  return is;
}
		

event&
event::setIOVersion(int ver)
{
  if ((ver >= 1) && (ver <= 2))
    this->_ioversion = ver;
  else {
    cerr << "unknown io version " << ver << endl;
    throw ("UnknownIOVersion");
  }
  return *this;
}


event
operator * (const lorentzTransform& L,
	    const event&            e)
{
  event r;
  r.beam(L * e.beam());
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
