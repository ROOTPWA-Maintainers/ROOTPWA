#ifndef PPUTIL_H
#define PPUTIL_H
	

#include <iostream>
#include <string>
#include <complex>
#include <cstdio>


typedef enum {
    g_Unknown = 0,
    g_Gamma = 1,
    g_Positron = 2,
    g_Electron = 3,
    g_Neutrino = 4,
    g_MuonPlus = 5,
    g_MuonMinus = 6,
    g_Pi0 = 7,
    g_PiPlus = 8,
    g_PiMinus = 9,
    g_KLong = 10,
    g_KPlus = 11,
    g_KMinus = 12,
    g_Neutron = 13,
    g_Proton = 14,
    g_AntiProton = 15,
    g_KShort = 16,
    g_Eta = 17,
    g_Lambda = 18,
    g_SigmaPlus = 19,
    g_Sigma0 = 20,
    g_SigmaMinus = 21,
    g_Xi0 = 22,
    g_XiMinus = 23,
    g_OmegaMinus = 24,
    g_AntiNeutron = 25,
    g_AntiLambda = 26,
    g_AntiSigmaMinus = 27,
    g_AntiSigma0 = 28,
    g_AntiSigmaPlus = 29,
    g_AntiXi0 = 30,
    g_AntiXiPlus = 31,
    g_AntiOmegaPlus = 32,

    g_Deuteron = 45,
    g_Triton = 49,

    g_Rho0 = 57,
    g_RhoPlus = 58,
    g_RhoMinus = 59,
    g_omega = 60,
    g_EtaPrime = 61,
    g_phiMeson = 62
} Geant_ID;


#define signof(x) (x<0 ? -1 : 1)
#define MAX(x,y) (x>y ? x : y)
#define MIN(x,y) (x<y ? x : y)


inline double tilde(const int l) { return pow(l+1,0.5); }

std::complex<double> D(const double alpha,
		       const double beta,
		       const double gamma,
		       const int    j,
		       const int    n,
		       const int    m);
double d_jmn_b(int    J,
	       int    M,
	       int    N,
	       double beta);
double clebsch(const int j1,
	       const int j2,
	       const int j3,
	       const int m1,
	       const int m2,
	       const int m3);
extern "C" void clebs_(int*,
		       int*,
		       int*,
		       int*,
		       int*,
		       int*,
		       int*,
		       int*);
int fact(int i);
double dfact(const double i);
double F(const int n,
	 const double p);
double lambda(const double a,
	      const double b,
	      const double c);
std::complex<double> q(const double M,
		       const double m1,
		       const double m2);

void addtab();
void subtab();
void ptab();
void settab(const int);
std::string itos(const int);
std::string chargetos(int charge);

Geant_ID    name2id(const std::string& name,
		    const int          q);
std::string id2name(const Geant_ID);


#endif  // PPUTIL_H
