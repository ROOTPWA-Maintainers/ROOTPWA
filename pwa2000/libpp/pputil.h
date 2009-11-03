#line 4 "pputil.nw"
#ifndef PPUTIL_H
#define PPUTIL_H
	
#line 14 "pputil.nw"
#include <iostream>
#include <string>
#include <complex>
#include <cstdio>


#line 26 "pputil.nw"
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



#line 7 "pputil.nw"
	
#line 135 "pputil.nw"
inline double tilde(int l)
{
	return pow(l+1,0.5);
}


#line 8 "pputil.nw"
	
#line 77 "pputil.nw"
#define signof(x) (x<0 ? -1 : 1)
#define MAX(x,y) (x>y ? x : y)
#define MIN(x,y) (x<y ? x : y)


#line 9 "pputil.nw"
	
#line 97 "pputil.nw"
	std::complex<double> D(double alpha,double beta,double gamma,int j,int n,int m);
	double clebsch(int j1,int j2,int j3,int m1,int m2,int m3);
	extern "C" void clebs_(int*,int*,int*,int*,int*,int*,int*,int*);
	double d_jmn_b(int J, int M, int N, double beta);
	int fact(int i);
	double dfact(double i);
	double F(int n,double p);
	double lambda(double a, double b, double c);
	std::complex<double> q(double M, double m1, double m2);


#line 116 "pputil.nw"
	void addtab() ;
	void subtab() ;
	void ptab() ;
	void settab(int) ;
	std::string itos(int);
	std::string chargetos(int);

#line 128 "pputil.nw"
	Geant_ID  name2id( std::string name,int q);
	std::string  id2name( Geant_ID );

#line 10 "pputil.nw"
#endif


