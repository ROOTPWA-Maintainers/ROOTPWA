#include <iostream>
#include <fstream>
#include <string>
#include "TFhh.h"
using namespace std;

extern ClebschGordanBox box;

Int_t debugFhh=1;

TFhh::TFhh(Int_t J_, Int_t S1, Int_t S2,
		Int_t lambda_, Int_t nu_,
		Int_t nLS, TLSAmpl* *LSampl,
		Int_t even_contr_) {

	J=J_;
	lambda= lambda_;
	nu= nu_;
	even_contraction = even_contr_;
	Int_t delta=lambda-nu;
	Nterms=0;

	name_str=new char[15];
	sprintf(name_str, "F_%lld_%lld", lambda, nu);

	for (Int_t iLS=0; iLS<nLS; iLS++) {
		if (LSampl[iLS]->Getdelta() ==  delta ||
				LSampl[iLS]->Getdelta() == -delta ) {
			// cout << "iLS=" << iLS << ", delta=" << delta << endl;
			TFracNum *CG3S = box.GetCG(LSampl[iLS]->GetS(), S1, S2);
			TFracNum SpinCouplFac = CG3S[CGIndex(S1, lambda, S2, -nu)];
			if (SpinCouplFac==TFracNum_Zero) {
				if (debugFhh==2)
					cout << "Clebsch-Gordan is zero" << endl;
			}
			else {
				TLSContrib* *newPTR = new TLSContrib*[Nterms+1];
				for (Int_t iC=0; iC<Nterms; iC++) newPTR[iC] = LSt[iC];
				newPTR[Nterms] = new TLSContrib(LSampl[iLS], delta, SpinCouplFac);
				// delete[] LSt;
				Nterms++;
				LSt=newPTR;
			}
		}
	}
	if (debugFhh) Print();
}

TFhh::TFhh(TFhh *sFhh, char flag) {
	if (flag != 'i' && flag !='m'){
		cerr << "TFhh::TFhh unknown flag "<<flag<<endl;
		return;
	}
	if (debugFhh)
		cout << "Initializing from single Amplitude" << endl;
	if ( ( flag =='i' && (sFhh->GetJ()) % 2      ) ||
			( flag =='m' && (sFhh->GetJ()) % 2 == 0 )    ) {
		cout << sFhh->GetName() << "[symm] = 0" << endl;
		Nterms=0;
	}
	else {
		name_str = new char[15];
		sprintf(name_str, "%s[symm]", sFhh->GetName());
		J        = sFhh->GetJ();
		lambda   = sFhh->GetLambda();
		nu       = sFhh->GetNu();
		even_contraction = sFhh->GetEvenContr();
		Nterms   = sFhh->GetNterms();
		LSt      = sFhh->GetLStPtr();
		Print();
	}
}

TFhh::TFhh(TFhh *sFhh, TFhh *xFhh) {
	name_str = new char[15];
	sprintf(name_str, "%s[symm]", sFhh->GetName());
	J        = sFhh->GetJ();
	lambda   = sFhh->GetLambda();
	nu       = sFhh->GetNu();
	even_contraction = sFhh->GetEvenContr();
	Nterms=0;
	LSt   =0;
	if (J != xFhh->GetJ() || even_contraction != xFhh->GetEvenContr() ||
			lambda != xFhh->GetNu() || nu != xFhh->GetLambda() ){
		cerr << "TFhh::TFhh(TFhh *, TFhh*): Something is wrong,"    << endl
			<< " source amplitudes for symmetrization do not match" << endl;
		return;
	}

	// Since some LS-terms cancel out, first a max-length array of
	// LSContrib pointers is filled, and then squeezed to only
	// non-zero contributions

	Int_t Ns = sFhh->GetNterms();
	Int_t Nx = xFhh->GetNterms();
	Nterms = Ns + Nx;
	TLSContrib* *pLSt = new TLSContrib*[Nterms];

	Int_t prevterms=0;

	for (Int_t i=0; i<Nterms; i++) {
		TLSContrib *to_be_added;
		if (i<Ns) to_be_added = sFhh->GetLStPtr()[i];
		else      to_be_added = xFhh->GetLStPtr()[i-Ns];
		Int_t found_in_prevterm=0;
		for (Int_t j=0; j<prevterms; j++) {
			if ( pLSt[j]->SameParameter(to_be_added) ) {
				found_in_prevterm=1;
				if (i<Ns) pLSt[j]->Add(to_be_added, false);
				else      pLSt[j]->Add(to_be_added, true);
			}
		}
		if (!found_in_prevterm) {
			if (i<Ns) pLSt[prevterms] = new TLSContrib(to_be_added, false);
			else      pLSt[prevterms] = new TLSContrib(to_be_added, true);
			prevterms++;
		}
	}

	//
	// Cancel zeros
	//
	Int_t non_zeros=0;
	for (Int_t i=0; i<prevterms; i++) {
		if ( pLSt[i]->GetNterms() != 0 ) non_zeros++;
	}

	Nterms=non_zeros;

	if (Nterms) {
		LSt = new TLSContrib*[Nterms];
		Int_t j=0;
		for (Int_t i=0; i<prevterms; i++) {
			if ( pLSt[i]->GetNterms() != 0 ) {
				LSt[j]=new TLSContrib(pLSt[i], false);
				j++;
			}
		}
	}
	Print();

}


Int_t TFhh::NonRelLimit(){
	NNRterms=0;
	for (Int_t i=0; i<Nterms; i++){
		if (! LSt[i]->IsPureRelativistic()) {
			Int_t jfound=-1;
			for (Int_t j=0; j<NNRterms; j++) {
				if (NRLSt[j]->CheckJLS(LSt[i])){
					jfound=j;
					break;
				}
			}
			if (jfound != -1) {
				NRLSt[jfound]->Add(LSt[i]);
			}
			else {
				NNRterms++;
				TLSNonRel **newNRLS=new TLSNonRel*[NNRterms];
				for (Int_t j=0; j<NNRterms-1; j++) {
					newNRLS[j]=NRLSt[j];
				}
				newNRLS[NNRterms-1] = new TLSNonRel(LSt[i]);
				NRLSt=newNRLS;
			}
		}
	}
	cout << name_str << " (NR) = "  << endl;
	for (Int_t j=0; j<NNRterms; j++) {
		cout << "LS=" << NRLSt[j]->GetL() << NRLSt[j]->GetS() << ": ";
		NRLSt[j]->Print();
	}
	return 0;
}

Int_t TFhh::PrintNRG(){
	cout << name_str << " (NR) = "  << endl;
	for (Int_t j=0; j<NNRterms; j++) {
		cout << "LS=" << NRLSt[j]->GetL() << NRLSt[j]->GetS() << " => ";
		NRLSt[j]->PrintG();
	}
	return 0;
}


Int_t
TFhh::Print() {
	cout << name_str << " =";
	if (even_contraction) cout << endl;
	else cout << " (iw)" << endl;
	for (Int_t iLSt=0; iLSt<Nterms; iLSt++) {
		LSt[iLSt]->Print();
	}
	if (Nterms==0) cout << " 0" << endl;
	return 0;
}

Int_t debugJSS=9999;

Int_t TJSS::CalcAmpl() {
	if (debugJSS>=2) {
		cout << "  Decay channel:   " << JMother;
		if (etaJ>0) cout << "+  ->  ";
		else        cout << "-  ->  ";
		cout << SDecay1;
		if (eta1>0) cout << "+ ";
		else        cout << "- ";
		cout << SDecay2;
		if (eta2>0) cout << "+";
		else        cout << "-";
		cout << endl;
	}

	// range for coupled Spin
	Smin=SDecay1-SDecay2; if (Smin<0) Smin=-Smin;
	Smax=SDecay1+SDecay2;

	if (debugJSS>=2) {
		cout << "possible S:";
		Int_t iS=Smin;
		while (iS<=Smax) {
			cout << " " << iS;
			iS++;
		}
	}
	cout << endl;

	Int_t intr_parity = etaJ*eta1*eta2;

	Lmax=JMother+Smax;
	Lmin=Lmax;
	for (Int_t iS=Smin; iS<=Smax; iS++) {
		Int_t Lm1=JMother-iS; if (Lm1<0) Lm1=-Lm1;
		if (Lm1<Lmin) Lmin=Lm1;
	}

	if (debugJSS>=2) cout << "possible L:";
	Int_t NL=0;
	Int_t *fL=0;
	Int_t              testL=0; // L even
	if (intr_parity<0) testL=1; // L odd
	Int_t tL=testL;

	while(testL<=Lmax) {
		if (testL>=Lmin) NL++;
		testL+=2;
	}

	if (NL) {
		fL=new Int_t[NL];
		NL=0;
		while(tL<=Lmax) {
			if (tL>=Lmin) {
				fL[NL]=tL;
				if (debugJSS>=2) cout << " " << fL[NL];
				NL++;
			}
			tL+=2;
		}
	}
	if (debugJSS>=2) cout << endl;

	Int_t even_contraction=1;
	if ((SDecay1+SDecay2+testL-JMother) % 2) even_contraction=0;

	if (debugJSS>=2) {
		if (even_contraction) cout << "contraction only with g~" << endl;
		else                  cout << "contraction with eps*p and g~" << endl;
	}

	Int_t                MaxPsiInternal=SDecay1;
	if (SDecay2<SDecay1) MaxPsiInternal=SDecay2;

	Int_t MaxPsiPhi = SDecay1+SDecay2;
	if (JMother<MaxPsiPhi) MaxPsiPhi=JMother;

	Int_t preloop=1;

	NLSAmpl=0;
	LSAmplitudes = 0;

	while (preloop==1 || preloop==0) {

		if (preloop==0) {
			LSAmplitudes = new TLSAmpl*[NLSAmpl];
			NLSAmpl=0;
			if (debugJSS>=2) cout << endl << "*************" << endl;
		}

		for (Int_t il=0; il<NL; il++) {

			Int_t L = fL[il];

			Int_t MaxPsiChi = SDecay1+SDecay2;
			if (L<MaxPsiChi) MaxPsiChi=L;

			Int_t MaxChiPhi = JMother;
			if (L<JMother) MaxChiPhi=L;

			// possible S in this case:
			Int_t SminL=JMother-L;
			if (SminL<0) SminL=-SminL;
			if (SminL<Smin) SminL=Smin;
			Int_t SmaxL=JMother+L;
			if (SmaxL>Smax) SmaxL=Smax;
			for (Int_t S_L=SminL; S_L<=SmaxL; S_L++) {
				if (debugJSS>=2) cout << "Amplitudes for L=" << L << " S="<< S_L
					<< "  Rank scheme [ "
						<< SDecay1+SDecay2 << " "
						<< L << " "
						<< JMother << "]" << endl;
				Int_t totalRank=SDecay1+SDecay2+L+JMother;
				Int_t MaxDelta=S_L;
				if (JMother<S_L) MaxDelta=JMother;

				Int_t IndexContractions=(totalRank-3)/2;
				if (even_contraction) IndexContractions= totalRank / 2;
				if (debugJSS>=2) cout << IndexContractions
					<< " Lorentz contractions." << endl;

				Int_t MaxContractionNumber=0;

				for (Int_t PsiInternal=0; PsiInternal<=MaxPsiInternal; PsiInternal++) {
					for (Int_t cChiPhi=0; cChiPhi<=MaxChiPhi; cChiPhi++) {
						for (Int_t cPsiChi=0; cPsiChi<=MaxPsiChi; cPsiChi++) {
							for (Int_t cPsiPhi=0; cPsiPhi<=MaxPsiPhi; cPsiPhi++) {
								if (debugJSS==3)
									cout << "Checking " << PsiInternal << " " << cPsiChi << " "
										<< cChiPhi << " " << cPsiPhi; // << endl;
								if (PsiInternal+cPsiChi+cChiPhi+cPsiPhi != IndexContractions){
									if (debugJSS==3) cout << " C-" << endl;
									continue;
								}
								else if (debugJSS==3) cout << " C+";
								if (even_contraction) {
									if (2*PsiInternal+cPsiChi+cPsiPhi != SDecay1+SDecay2) {
										if (debugJSS==3) cout << "S-" << endl;
										continue;
									}
									else if (debugJSS==3) cout << "S+";
									if (cPsiChi+cChiPhi != L) {
										if (debugJSS==3) cout << "L-" << endl;
										continue;
									}
									else if (debugJSS==3) cout << "L+";
									if (cChiPhi+cPsiPhi != JMother) {
										if (debugJSS==3) cout << "J-" << endl;
										continue;
									}
									else if (debugJSS==3) cout << "J+";
								}
								else {
									if (L      -cPsiChi-cChiPhi>1) {
										if (debugJSS==3) cout << "L-" << endl;
										continue;
									}
									else if (debugJSS==3) cout << "L+";
									if (JMother-cPsiPhi-cChiPhi>1) {
										if (debugJSS==3) cout << "J-" << endl;
										continue;
									}
									else if (debugJSS==3) cout << "J+";
								}
								Int_t r_ome=SDecay1-PsiInternal;
								Int_t r_eps=SDecay2-PsiInternal;

								if (!even_contraction) {
									Int_t PsiRest=SDecay1+SDecay2-2*PsiInternal-cPsiChi-cPsiPhi;
									if ( PsiRest<0 || PsiRest>2 ) {
										if (debugJSS==3) cout << "R-" << endl;
										continue;
									}
									else if (debugJSS==3) cout << "R+";
									if (PsiRest==2){
										if (r_ome>0) r_ome--;
										else {
											if (debugJSS==3) cout << "O-";
											continue;
										}
										if (r_eps>0) r_eps--;
										else {
											if (debugJSS==3) cout << "E-";
											continue;
										}
									}
								}
								//
								// Original ordering:
								//
								//for (Int_t cPhiOmega=0; cPhiOmega<=r_ome; cPhiOmega++) {
								//  for (Int_t cChiOmega=0; cChiOmega<=r_ome-cPhiOmega;
								//       cChiOmega++) {
								//
								// For agreement with paper:
								//
								// error found 14.10.07 "r_eps" replaced with "r_ome"
								if (debugJSS==3) cout << "{"<<r_ome<<"}";
								for (Int_t cChiOmega=0; cChiOmega<=r_ome; cChiOmega++) {
									for (Int_t cPhiOmega=0; cPhiOmega<=r_ome-cChiOmega;
											cPhiOmega++) {
										Int_t cPhiEps=cPsiPhi-cPhiOmega;
										Int_t cChiEps=cPsiChi-cChiOmega;
										if (debugJSS==3) cout << "[" << cPhiEps
											<< cChiEps<< r_eps<<"]";
										if (cPhiEps<0 || cChiEps<0 ||
												cPhiEps+cChiEps>r_eps) {
											continue;
										}
										else if (debugJSS==3) cout << "E+ OK" << endl;
										//
										//
										//
										if (debugJSS>=2)
											cout << "Checking PsiInt=" << PsiInternal
												<< " cPsiChi=" << cPsiChi
												<< " cChiPhi=" << cChiPhi
												<< " cPsiPhi=" << cPsiPhi << endl
												<< " cPhiOmega=" << cPhiOmega
												<< " cChiOmega=" << cChiOmega
												<< " cPhiEps=" << cPhiEps
												<< " cChiEps=" << cChiEps << endl;
										Int_t cc=0;
										if (debugJSS>=2) {
											cout << "Contraction pattern ";
											cc=cPsiPhi; while (cc--) cout << "#";
											if (cPsiPhi) cout << "(";
											cc=cPhiOmega; while (cc--) cout << "o";
											cc=cPhiEps;   while (cc--) cout << "e";
											if (cPsiPhi) cout << ")";
											cout<<" "<< SDecay1+SDecay2;
											cc=PsiInternal; while (cc--) cout << "'";
											cout << " ";
											cc=cPsiChi; while (cc--) cout << "#";
											if (cPsiChi) cout << "(";
											cc=cChiOmega; while (cc--) cout << "o";
											cc=cChiEps;   while (cc--) cout << "e";
											if (cPsiChi) cout << ")";
											cout << " "<< L               << " ";
											cc=cChiPhi; while (cc--) cout << "#";
											cout << " "<< JMother         << " ";
											cc=cPsiPhi; while (cc--) cout << "#";
											if (cPsiPhi) cout << "(";
											cc=cPhiOmega; while (cc--) cout << "o";
											cc=cPhiEps;   while (cc--) cout << "e";
											if (cPsiPhi) cout << ")";

											if (preloop) cout << " delta:";
											else cout << endl;
										}
										for (Int_t delta=0; delta<=MaxDelta; delta++) {
											if (preloop) {
												if (debugJSS>=2) cout << " " << delta;
												NLSAmpl++;
											}
											else {
												if (debugJSS>=2) cout << " Constructing LS-Amplitude "
													<< NLSAmpl << endl;

												Int_t SameContr=0;
												for (SameContr=0; SameContr<NLSAmpl; SameContr++) {
													if (LSAmplitudes[SameContr]->
															CheckContraction(L, S_L, PsiInternal, cPsiChi,
																cChiPhi, cPsiPhi,
																cPhiOmega, cChiOmega,
																cPhiEps, cChiEps)) break;
												}
												Int_t ContractionNumber=0;
												if (SameContr<NLSAmpl)
													ContractionNumber = LSAmplitudes[SameContr]->
														GetContraction();
												else
													ContractionNumber = MaxContractionNumber+1;

												LSAmplitudes[NLSAmpl]=
													new TLSAmpl(SDecay1, SDecay2, L, JMother,
															delta, S_L,
															PsiInternal, cPsiChi, cChiPhi, cPsiPhi,
															cPhiOmega, cChiOmega, cPhiEps, cChiEps,
															ContractionNumber);

												if (LSAmplitudes[NLSAmpl]->GetNterms() ) {
													NLSAmpl++;
													if (ContractionNumber > MaxContractionNumber)
														MaxContractionNumber++;
												}
												else {
													delete LSAmplitudes[NLSAmpl];
												}
											}
										}
										if (debugJSS>=2 && preloop) cout << endl;
									}
								}
							}
							}
							}
						}
					}
				}
				if (preloop) {
					if (debugJSS)
						cout << NLSAmpl << " LS-Amplitudes to be evaluated." << endl;
				}
				preloop--;
			}
			if (debugJSS) {
				cout << NLSAmpl << " LS-Amplitudes found to be non-zero." << endl;
				cout << "++++++++++++++++++++++++++++++++++++" << endl;
				cout << "+++ Helicity-coupling amplitudes +++" << endl;
				cout << "++++++++++++++++++++++++++++++++++++" << endl;
			}

			NFhhAmpl=0;
			FhhAmpl = new TFhh*[(SDecay1+1)*(SDecay2+1)];

			for (Int_t lambda=0; lambda<=SDecay1; lambda++)
				for (Int_t nu=-SDecay2; nu<=SDecay2; nu++) {
					//    for (Int_t nu=-SDecay1; nu<=SDecay2; nu++) { bug!!!! found 4.3.08
					if (lambda==0 && nu<0) continue;
					FhhAmpl[NFhhAmpl] = new TFhh(JMother, SDecay1, SDecay2,
							lambda, nu, NLSAmpl, LSAmplitudes,
							even_contraction);
					if (FhhAmpl[NFhhAmpl]->GetNterms()) {
						NFhhAmpl++;
					}
					else {
						delete FhhAmpl[NFhhAmpl];
					}
				}

				if (debugJSS)
					cout << NFhhAmpl << " non-zero helicity-coupling amplitudes" << endl;

				NFhhIdAmpl=0;
				FhhIdAmpl = new TFhh*[(SDecay1+1)*(SDecay2+1)];

				if ( SDecay1==SDecay2 && eta1==eta2 ) {
					if (debugJSS)
						cout << endl << " for identical-particle decay:" << endl;
					for (Int_t ifhh=0; ifhh<NFhhAmpl; ifhh++) {
						FhhIdAmpl[ifhh]=0;
						if (FhhAmpl[ifhh]) {
							if (FhhAmpl[ifhh]->IsNuNu()) {
								FhhIdAmpl[ifhh] = new TFhh(FhhAmpl[ifhh], 'i');
							}
							else if (FhhAmpl[ifhh]->IsNuMinusNu() ) {
								FhhIdAmpl[ifhh] = new TFhh(FhhAmpl[ifhh], 'm');
							}
							else {
								Int_t found_partner=0;
								for (Int_t jfhh=0; jfhh<NFhhAmpl; jfhh++) {
									if (FhhAmpl[ifhh]->GetLambda()==FhhAmpl[jfhh]->GetNu() &&
											FhhAmpl[ifhh]->GetNu()==FhhAmpl[jfhh]->GetLambda() ) {
										found_partner=1;
										FhhIdAmpl[ifhh] = new TFhh(FhhAmpl[ifhh], FhhAmpl[jfhh]);
										// ** continue here **
									}
								}
								if (!found_partner)
									cerr << "?!?! No partner for amplitude "<< FhhAmpl[ifhh]->GetName()
										<< endl;
							}
						}
					}


				}

				cout << NFhhAmpl<< " amplitudes: non-relativistic limit" << endl;
				for (Int_t i=0; i<NFhhAmpl; i++){
					FhhAmpl[i]->NonRelLimit();
				}

				cout << "Check non-relativistic G's" << endl;
				for (Int_t i=0; i<NFhhAmpl; i++){
					FhhAmpl[i]->PrintNRG();
				}

				return 0;
				}

			Int_t TJSS::PrintHFILE() {
				char DecayName[10];
				sprintf(DecayName,"%lld%lld%lld%c%c", JMother, SDecay1, SDecay2,
						etaJ*eta1*eta2==-1?'n':'p',' ');
				char ofname[20];
				sprintf(ofname, "CalcAmpl-%s.h", DecayName);
				ofstream ofs(ofname);
				ofs << "// CalcAmpl output for " << DecayName << endl;
				ofs << "const int FhhAmpl_" << DecayName << "[] = { " << endl;
				ofs << "  " << NFhhAmpl
					<< ",               // number of Fhh amplitudes" << endl;
				for (int i=0;i<NFhhAmpl; i++){
					ofs << "  "
						<< FhhAmpl[i]->GetJ()       <<", "
						<< FhhAmpl[i]->GetLambda()  <<", "
						<< FhhAmpl[i]->GetNu()      <<", "
						<< FhhAmpl[i]->GetEvenContr() <<",      // "
						<< "F"
						<< FhhAmpl[i]->GetLambda()
						<< FhhAmpl[i]->GetNu()
						<< ": J, lambda, nu, even_contr"
						<< endl;
					ofs << "    " << FhhAmpl[i]->GetNterms()
						<< ",             // number of contributions" << endl;
					for (int j=0; j<FhhAmpl[i]->GetNterms(); j++) {
						TLSContrib* lsa = FhhAmpl[i]->GetLStPtr()[j];
						ofs << "    "
							<< lsa->GetJ() <<", "
							<< lsa->GetL() <<", "
							<< lsa->GetS() <<", "
							<< lsa->GetDelta()<<", "
							<< lsa->GetRunningNumber()
							<< ",    // contr. F"
							<< FhhAmpl[i]->GetLambda()
							<< FhhAmpl[i]->GetNu() << "-" << j
							<< ": J, L, S, delta, #[g,f,h,...]" << endl;
						ofs << "      "
							<< lsa->GetNterms() << ", "
							<< lsa->GetNormFactor()->GetSign() << ", "
							<< lsa->GetNormFactor()->GetNumerator() << ", "
							<< lsa->GetNormFactor()->GetDenominator()
							<< ",     // number of terms, squared norm factor sign/N/D" << endl;
						for (int k=0; k<lsa->GetNterms(); k++) {
							ofs << "      "
								<< lsa->GetTermFracNum()[k].GetSign() << ", "
								<< lsa->GetTermFracNum()[k].GetNumerator() << ", "
								<< lsa->GetTermFracNum()[k].GetDenominator() << ", "
								<< lsa->GetTermg1pot()[k] << ", "
								<< lsa->GetTermg2pot()[k]
								<< ",  // squared sign/N/D, exponents of g_s and g_sigma" << endl;
						}
					}
				}
				ofs << "};" << endl;

				return 0;
			}
