#include <cstdlib>

#include "pputil.h"
	

using namespace std;

	
complex<double>
D(const double alpha,
  const double beta,
  const double gamma,
  const int    j,
  const int    m,
  const int    n)
{
  const complex<double> i(0, 1);
  const complex<double> c = exp(-i * ((m / 2.0) * alpha + (n / 2.0) * gamma) ) * d_jmn_b(j, m, n, beta);
  return c;
}


double
d_jmn_b(int    J,
	int    M,
	int    N,
	double beta)
{
  int temp_M, k, k_low, k_hi;
  double const_term = 0.0, sum_term = 0.0, d = 1.0;
  int m_p_n, j_p_m, j_p_n, j_m_m, j_m_n;
  int kmn1, kmn2, jmnk, jmk, jnk;
  double kk;

  if (J < 0 || abs(M) > J || abs(N) > J) {
    cerr << endl;
    cerr << "d_jmn_b: you have entered an illegal number for J, M, N." << endl;
    cerr << "Must follow these rules: J >= 0, abs(M) <= J, and abs(N) <= J." << endl;
    cerr << "J = " << J <<  " M = " << M <<  " N = " << N << endl;
    exit (1);
  }

  if (beta < 0) {
    beta = fabs (beta);
    temp_M = M;
    M = N;
    N = temp_M;
  }

  m_p_n = (M + N) / 2;
  j_p_m = (J + M) / 2;
  j_m_m = (J - M) / 2;
  j_p_n = (J + N) / 2;
  j_m_n = (J - N) / 2;
	
  kk = (double) fact (j_p_m) * (double) fact (j_m_m) * (double) fact (j_p_n) * (double) fact (j_m_n) ;
  const_term = pow ((-1.0), (j_p_m)) * sqrt (kk);	
 
  k_low = MAX (0, m_p_n);
  k_hi = MIN (j_p_m, j_p_n);

  for (k = k_low; k <= k_hi; k++) {

    kmn1 = 2 * k - (M + N) / 2;
    jmnk = J + (M + N) / 2 - 2 * k;
    jmk = (J + M) / 2 - k;
    jnk = (J + N) / 2 - k;
    kmn2 = k - (M + N) / 2;

    sum_term += pow ((-1.0), (k)) *
      ((pow (cos (beta / 2.0), kmn1)) * (pow (sin (beta / 2.0), jmnk))) /
      (fact (k) * fact (jmk) * fact (jnk) * fact (kmn2));
  }

  d = const_term * sum_term;
  return d;
}


double
clebsch(const int j1,
	const int j2,
	const int j3,
	const int m1,
	const int m2,
	const int m3)
{
  int nu = 0;
  double exp;
  double n0, n1, n2, n3, n4, n5;
  double d0, d1, d2, d3, d4;
  double sum;
  double A;

  if ((m1 + m2) != m3) {
    return 0;
  }

  sum = 0;
  while ( (d3=(j1-j2-m3)/2+nu) < 0
	  || (n2=(j1-m1)/2+nu) < 0 ) { nu++;}
  while ( (d1=(j3-j1+j2)/2-nu) >= 0
	  && (d2=(j3+m3)/2-nu) >= 0
	  && (n1=(j2+j3+m1)/2-nu) >= 0 ) {
    d3=((j1-j2-m3)/2+nu);
    n2=((j1-m1)/2+nu);
    d0=dfact((double) nu);
    exp=nu+(j2+m2)/2;
    n0 = (double) pow(-1,exp);
    sum+=(n0*dfact(n1)*dfact(n2))/(d0*dfact(d1)*dfact(d2)*dfact(d3));
    nu++;
  }

  if ( sum == 0 ) {
    return 0;
  }

  n0 = j3+1;
  n1 = dfact((double) (j3+j1-j2)/2);
  n2 = dfact((double) (j3-j1+j2)/2);
  n3 = dfact((double) (j1+j2-j3)/2);
  n4 = dfact((double) (j3+m3)/2);
  n5 = dfact((j3-m3)/2);
	
  d0 = dfact((double) (j1+j2+j3)/2+1);
  d1 = dfact((double) (j1-m1)/2);
  d2 = dfact((double) (j1+m1)/2);
  d3 = dfact((double) (j2-m2)/2);
  d4 = dfact((double) (j2+m2)/2);

  A = ((double) (n0*n1*n2*n3*n4*n5))/((double) (d0*d1*d2*d3*d4));
	
  return pow(A,0.5)*sum;
}


double
dfact(const double i)
{
  if (i < 0.00001)
    return 1;
  if (i < 0)
    return 0;
  return i * dfact(i - 1);
}


double
F(const int    n,
  const double p)
{
#define Pr 0.1973 // Gev/c corresponds to 1 fermi
  double ret;
  double z = (p/Pr) * (p/Pr);
  int m = n/2;
  switch (m) {
  case 0:
    ret = 1.0;
    break;
  case 1:
    ret = sqrt((2.0 * z)/(z + 1));
    break;
  case 2:
    ret = sqrt((13.0 * z * z)/(pow(z - 3.0,2.0) + 9.0 * z));
    break;
  case 3:
    ret = sqrt( (277.0 * pow(z,3.0))/(z * pow(z - 15.0,2.0) + 9.0 * pow(2.0 * z - 5.0,2.0)));
    break;
  case 4:
    ret = sqrt( (12746.0 * pow(z,4.0))/( pow(z * z - 45.0 * z + 105.0,2.0) + 25.0 * z * pow(2.0 * z - 21.0,2.0)));
    break;
  case 5:
    ret =  sqrt(z*z*z*z*z/(893025.0 +99225.0*z +6300.0*z*z +315.0*z*z*z +15.0*z*z*z*z +z*z*z*z*z));
    break;
  default:
    cerr << "Blatt-Weisskopf called for undefined L = " << n/2 << endl;
    ret = 1.0;
    break;
  }
  return ret;
}


double
lambda(const double a,
       const double b,
       const double c)
{
  return a * a + b * b + c * c - 2.0 * (a * b + b * c + c * a);
}


complex<double>
q(const double M,
  const double m1,
  const double m2)
{
  double lam = lambda(M * M, m1 * m1, m2 * m2);
  complex<double> ret;
  if (lam < 0)
    return complex<double>(0.0, sqrt(fabs(lam / (4 * M * M))));
  return complex<double>(sqrt(lam / (4 * M * M)), 0.0 );
}


int
fact(int i)
{
  int f = 1;
  if (i == 0 || i == 1)
    f = 1;
  else {
    while (i > 0) {
      f = f * i;
      i--;
    }
  }
  return f;
}


int ntab    = 0;
int tabsize = 8;


void
addtab()
{
  ntab++;
}


void
subtab()
{
  ntab--;
  if (ntab < 0)
    ntab = 0;
}


void
ptab()
{
  for (int i = 0; i < ntab; i++)
    for (int s = 0; s < tabsize ; s++)
      cout << " ";
}


void
settab(const int nchar)
{
  tabsize = nchar;
}


string
id2name(const Geant_ID type) 
{
  switch (type) {
  case g_EtaPrime:
    return "eta'(958)";
    break;
  case g_PiMinus:
  case g_PiPlus:
    return "pi";
    break;
  case g_KMinus:
  case g_KPlus:
    return "K";
    break;
  case g_KShort:
    return "K0";
    break;
  case g_Pi0:
    return "pi0";
    break;
  case g_Eta:
    return "eta";
    break;
  case g_Rho0:
    return ("rho(770)");
    break;
  case g_RhoPlus:
    return ("rho(770)");
    break;
  case g_RhoMinus:
    return ("rho(770)");
    break;
  case g_omega:
    return ("omega(782)");
    break;
  case g_phiMeson:
    return ("phi(1020)");
    break;
  case g_Proton:
    return "p";
    break;
  case g_AntiProton:
    return("pbar");
    break;
  case g_Neutron:
    return "n";
    break;
  case g_Gamma:
    return ("gamma");
    break;
  case g_Electron:
  case g_Positron:
    return ("e");
    break;
  case g_Deuteron:
    return("d");
    break;
  case g_Lambda:
    return ("lambda");
    break;
  default:
    return "unknown";
    break;
  }
}


Geant_ID
name2id(const string& name,
	const int     q)
{
  if (name == "pi") {
    switch (q) {
    case -1:
      return g_PiMinus;
      break;
    case 1:
      return g_PiPlus;
      break;
    }
  }
  else if (name == "omega(782)") {
    return(g_omega);
  }
  else if (name == "e") {
    switch (q) {
    case -1:
      return g_Electron;
      break;
    case 1:
      return g_Positron;
      break;
    }
  }
  else if (name == "K") {
    switch (q) {
    case -1:
      return g_KMinus;
      break;
    case 1:
      return g_KPlus;
      break;
    }
  }
  else if (name == "K0") {
    return g_KShort;
  }
  else if (name == "KLong") {
    return g_KLong;
  }
  else if (name == "pi0") {
    return g_Pi0;
  }
  else if (name == "eta") {
    return g_Eta;
  }
  else if (name == "phi(1020)") {
    return g_phiMeson;
  }
  else if (name == "p") {
    return g_Proton;
  }
  else if (name == "pbar") {
    return g_AntiProton;
  }
  else if (name == "n") {
    return g_Neutron;
  }
  else if (name == "d") {
    return g_Deuteron;
  }
  else if (name == "gamma") {
    return g_Gamma;
  }
  else if (name == "eta'(958)") {
    return g_EtaPrime;
  }
  else if (name == "lambda") {
    return g_Lambda;
  }
  else
    return (g_Unknown);
  return (g_Unknown);
}


string
itos(const int i)
{
  int    digits = (int)log10((float) i) + 2;
  char*  c_s    = (char*)malloc(digits * sizeof(char));
  string s;
  sprintf(c_s, "%d", i);
  s = c_s;
  return s;
}


string
chargetos(int charge)
{
  string s;
  string c;
  if (charge) {
    c      = (charge < 0) ? "-" : "+";
    charge = abs(charge);
    while (charge--)
      s += c;
  }
  return s;
}
