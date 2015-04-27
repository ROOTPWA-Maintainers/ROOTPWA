#include <TJSS.h>

#include <iostream>
#include <fstream>
#include <cstdio>

using namespace std;

unsigned int TJSS::_debugLevel = 9999;

long TJSS::CalcAmpl() {
	if (_debugLevel >= 2) {
		cout << "  Decay channel:   " << JMother;
		if (etaJ > 0)
			cout << "+  ->  ";
		else
			cout << "-  ->  ";
		cout << SDecay1;
		if (eta1 > 0)
			cout << "+ ";
		else
			cout << "- ";
		cout << SDecay2;
		if (eta2 > 0)
			cout << "+";
		else
			cout << "-";
		cout << endl;
	}

	// range for coupled Spin
	long Smin = SDecay1 - SDecay2;
	if (Smin < 0)
		Smin = -Smin;
	long Smax = SDecay1 + SDecay2;

	if (_debugLevel >= 2) {
		cout << "possible S:";
		long iS = Smin;
		while (iS <= Smax) {
			cout << " " << iS;
			iS++;
		}
	}
	cout << endl;

	const long intr_parity = etaJ * eta1 * eta2;

	long Lmax = JMother + Smax;
	long Lmin = Lmax;
	for (long iS = Smin; iS <= Smax; iS++) {
		long Lm1 = JMother - iS;
		if (Lm1 < 0)
			Lm1 = -Lm1;
		if (Lm1 < Lmin)
			Lmin = Lm1;
	}

	if (_debugLevel >= 2)
		cout << "possible L:";
	long NL = 0;
	long *fL = 0;
	long testL = 0; // L even
	if (intr_parity < 0)
		testL = 1; // L odd
	long tL = testL;

	while (testL <= Lmax) {
		if (testL >= Lmin)
			NL++;
		testL += 2;
	}

	if (NL) {
		fL = new long[NL];
		NL = 0;
		while (tL <= Lmax) {
			if (tL >= Lmin) {
				fL[NL] = tL;
				if (_debugLevel >= 2)
					cout << " " << fL[NL];
				NL++;
			}
			tL += 2;
		}
	}
	if (_debugLevel >= 2)
		cout << endl;

	long even_contraction = 1;
	if ((SDecay1 + SDecay2 + testL - JMother) % 2)
		even_contraction = 0;

	if (_debugLevel >= 2) {
		if (even_contraction)
			cout << "contraction only with g~" << endl;
		else
			cout << "contraction with eps*p and g~" << endl;
	}

	long MaxPsiInternal = SDecay1;
	if (SDecay2 < SDecay1)
		MaxPsiInternal = SDecay2;

	long MaxPsiPhi = SDecay1 + SDecay2;
	if (JMother < MaxPsiPhi)
		MaxPsiPhi = JMother;

	long preloop = 1;

	NLSAmpl = 0;
	LSAmplitudes = 0;

	while (preloop == 1 || preloop == 0) {

		if (preloop == 0) {
			LSAmplitudes = new TLSAmpl*[NLSAmpl];
			NLSAmpl = 0;
			if (_debugLevel >= 2)
				cout << endl << "*************" << endl;
		}

		for (long il = 0; il < NL; il++) {

			long L = fL[il];

			long MaxPsiChi = SDecay1 + SDecay2;
			if (L < MaxPsiChi)
				MaxPsiChi = L;

			long MaxChiPhi = JMother;
			if (L < JMother)
				MaxChiPhi = L;

			// possible S in this case:
			long SminL = JMother - L;
			if (SminL < 0)
				SminL = -SminL;
			if (SminL < Smin)
				SminL = Smin;
			long SmaxL = JMother + L;
			if (SmaxL > Smax)
				SmaxL = Smax;
			for (long S_L = SminL; S_L <= SmaxL; S_L++) {
				if (_debugLevel >= 2)
					cout << "Amplitudes for L=" << L << " S=" << S_L
							<< "  Rank scheme [ " << SDecay1 + SDecay2 << " "
							<< L << " " << JMother << "]" << endl;
				long totalRank = SDecay1 + SDecay2 + L + JMother;
				long MaxDelta = S_L;
				if (JMother < S_L)
					MaxDelta = JMother;

				long IndexContractions = (totalRank - 3) / 2;
				if (even_contraction)
					IndexContractions = totalRank / 2;
				if (_debugLevel >= 2)
					cout << IndexContractions << " Lorentz contractions."
							<< endl;

				long MaxContractionNumber = 0;

				for (long PsiInternal = 0; PsiInternal <= MaxPsiInternal;
						PsiInternal++) {
					for (long cChiPhi = 0; cChiPhi <= MaxChiPhi; cChiPhi++) {
						for (long cPsiChi = 0; cPsiChi <= MaxPsiChi;
								cPsiChi++) {
							for (long cPsiPhi = 0; cPsiPhi <= MaxPsiPhi;
									cPsiPhi++) {
								if (_debugLevel == 3)
									cout << "Checking " << PsiInternal << " "
											<< cPsiChi << " " << cChiPhi << " "
											<< cPsiPhi; // << endl;
								if (PsiInternal + cPsiChi + cChiPhi + cPsiPhi
										!= IndexContractions) {
									if (_debugLevel == 3)
										cout << " C-" << endl;
									continue;
								} else if (_debugLevel == 3)
									cout << " C+";
								if (even_contraction) {
									if (2 * PsiInternal + cPsiChi + cPsiPhi
											!= SDecay1 + SDecay2) {
										if (_debugLevel == 3)
											cout << "S-" << endl;
										continue;
									} else if (_debugLevel == 3)
										cout << "S+";
									if (cPsiChi + cChiPhi != L) {
										if (_debugLevel == 3)
											cout << "L-" << endl;
										continue;
									} else if (_debugLevel == 3)
										cout << "L+";
									if (cChiPhi + cPsiPhi != JMother) {
										if (_debugLevel == 3)
											cout << "J-" << endl;
										continue;
									} else if (_debugLevel == 3)
										cout << "J+";
								} else {
									if (L - cPsiChi - cChiPhi > 1) {
										if (_debugLevel == 3)
											cout << "L-" << endl;
										continue;
									} else if (_debugLevel == 3)
										cout << "L+";
									if (JMother - cPsiPhi - cChiPhi > 1) {
										if (_debugLevel == 3)
											cout << "J-" << endl;
										continue;
									} else if (_debugLevel == 3)
										cout << "J+";
								}
								long r_ome = SDecay1 - PsiInternal;
								long r_eps = SDecay2 - PsiInternal;

								if (!even_contraction) {
									long PsiRest = SDecay1 + SDecay2
											- 2 * PsiInternal - cPsiChi
											- cPsiPhi;
									if (PsiRest < 0 || PsiRest > 2) {
										if (_debugLevel == 3)
											cout << "R-" << endl;
										continue;
									} else if (_debugLevel == 3)
										cout << "R+";
									if (PsiRest == 2) {
										if (r_ome > 0)
											r_ome--;
										else {
											if (_debugLevel == 3)
												cout << "O-";
											continue;
										}
										if (r_eps > 0)
											r_eps--;
										else {
											if (_debugLevel == 3)
												cout << "E-";
											continue;
										}
									}
								}
								//
								// Original ordering:
								//
								//for (long cPhiOmega=0; cPhiOmega<=r_ome; cPhiOmega++) {
								//  for (long cChiOmega=0; cChiOmega<=r_ome-cPhiOmega;
								//       cChiOmega++) {
								//
								// For agreement with paper:
								//
								// error found 14.10.07 "r_eps" replaced with "r_ome"
								if (_debugLevel == 3)
									cout << "{" << r_ome << "}";
								for (long cChiOmega = 0; cChiOmega <= r_ome;
										cChiOmega++) {
									for (long cPhiOmega = 0;
											cPhiOmega <= r_ome - cChiOmega;
											cPhiOmega++) {
										long cPhiEps = cPsiPhi - cPhiOmega;
										long cChiEps = cPsiChi - cChiOmega;
										if (_debugLevel == 3)
											cout << "[" << cPhiEps << cChiEps
													<< r_eps << "]";
										if (cPhiEps < 0 || cChiEps < 0
												|| cPhiEps + cChiEps > r_eps) {
											continue;
										} else if (_debugLevel == 3)
											cout << "E+ OK" << endl;
										//
										//
										//
										if (_debugLevel >= 2)
											cout << "Checking PsiInt="
													<< PsiInternal
													<< " cPsiChi=" << cPsiChi
													<< " cChiPhi=" << cChiPhi
													<< " cPsiPhi=" << cPsiPhi
													<< endl << " cPhiOmega="
													<< cPhiOmega
													<< " cChiOmega="
													<< cChiOmega << " cPhiEps="
													<< cPhiEps << " cChiEps="
													<< cChiEps << endl;
										long cc = 0;
										if (_debugLevel >= 2) {
											cout << "Contraction pattern ";
											cc = cPsiPhi;
											while (cc--)
												cout << "#";
											if (cPsiPhi)
												cout << "(";
											cc = cPhiOmega;
											while (cc--)
												cout << "o";
											cc = cPhiEps;
											while (cc--)
												cout << "e";
											if (cPsiPhi)
												cout << ")";
											cout << " " << SDecay1 + SDecay2;
											cc = PsiInternal;
											while (cc--)
												cout << "'";
											cout << " ";
											cc = cPsiChi;
											while (cc--)
												cout << "#";
											if (cPsiChi)
												cout << "(";
											cc = cChiOmega;
											while (cc--)
												cout << "o";
											cc = cChiEps;
											while (cc--)
												cout << "e";
											if (cPsiChi)
												cout << ")";
											cout << " " << L << " ";
											cc = cChiPhi;
											while (cc--)
												cout << "#";
											cout << " " << JMother << " ";
											cc = cPsiPhi;
											while (cc--)
												cout << "#";
											if (cPsiPhi)
												cout << "(";
											cc = cPhiOmega;
											while (cc--)
												cout << "o";
											cc = cPhiEps;
											while (cc--)
												cout << "e";
											if (cPsiPhi)
												cout << ")";

											if (preloop)
												cout << " delta:";
											else
												cout << endl;
										}
										for (long delta = 0; delta <= MaxDelta;
												delta++) {
											if (preloop) {
												if (_debugLevel >= 2)
													cout << " " << delta;
												NLSAmpl++;
											} else {
												if (_debugLevel >= 2)
													cout
															<< " Constructing LS-Amplitude "
															<< NLSAmpl << endl;

												long SameContr = 0;
												for (SameContr = 0;
														SameContr < NLSAmpl;
														SameContr++) {
													if (LSAmplitudes[SameContr]->CheckContraction(
															L, S_L, PsiInternal,
															cPsiChi, cChiPhi,
															cPsiPhi, cPhiOmega,
															cChiOmega, cPhiEps,
															cChiEps))
														break;
												}
												long ContractionNumber = 0;
												if (SameContr < NLSAmpl)
													ContractionNumber =
															LSAmplitudes[SameContr]->GetContraction();
												else
													ContractionNumber =
															MaxContractionNumber
																	+ 1;

												LSAmplitudes[NLSAmpl] =
														new TLSAmpl(SDecay1,
																SDecay2, L,
																JMother, delta,
																S_L,
																PsiInternal,
																cPsiChi,
																cChiPhi,
																cPsiPhi,
																cPhiOmega,
																cChiOmega,
																cPhiEps,
																cChiEps,
																ContractionNumber);

												if (LSAmplitudes[NLSAmpl]->GetNterms()) {
													NLSAmpl++;
													if (ContractionNumber
															> MaxContractionNumber)
														MaxContractionNumber++;
												} else {
													delete LSAmplitudes[NLSAmpl];
												}
											}
										}
										if (_debugLevel >= 2 && preloop)
											cout << endl;
									}
								}
							}
						}
					}
				}
			}
		}
		if (preloop) {
			if (_debugLevel)
				cout << NLSAmpl << " LS-Amplitudes to be evaluated." << endl;
		}
		preloop--;
	}
	if (_debugLevel) {
		cout << NLSAmpl << " LS-Amplitudes found to be non-zero." << endl;
		cout << "++++++++++++++++++++++++++++++++++++" << endl;
		cout << "+++ Helicity-coupling amplitudes +++" << endl;
		cout << "++++++++++++++++++++++++++++++++++++" << endl;
	}

	NFhhAmpl = 0;
	FhhAmpl = new TFhh*[(SDecay1 + 1) * (SDecay2 + 1)];

	for (long lambda = 0; lambda <= SDecay1; lambda++)
		for (long nu = -SDecay2; nu <= SDecay2; nu++) {
			//    for (long nu=-SDecay1; nu<=SDecay2; nu++) { bug!!!! found 4.3.08
			if (lambda == 0 && nu < 0)
				continue;
			FhhAmpl[NFhhAmpl] = new TFhh(JMother, SDecay1, SDecay2, lambda, nu,
					NLSAmpl, LSAmplitudes, even_contraction);
			if (FhhAmpl[NFhhAmpl]->GetNterms()) {
				NFhhAmpl++;
			} else {
				delete FhhAmpl[NFhhAmpl];
			}
		}

	if (_debugLevel)
		cout << NFhhAmpl << " non-zero helicity-coupling amplitudes" << endl;

	NFhhIdAmpl = 0;
	FhhIdAmpl = new TFhh*[(SDecay1 + 1) * (SDecay2 + 1)];

	if (SDecay1 == SDecay2 && eta1 == eta2) {
		if (_debugLevel)
			cout << endl << " for identical-particle decay:" << endl;
		for (long ifhh = 0; ifhh < NFhhAmpl; ifhh++) {
			FhhIdAmpl[ifhh] = 0;
			if (FhhAmpl[ifhh]) {
				if (FhhAmpl[ifhh]->IsNuNu()) {
					FhhIdAmpl[ifhh] = new TFhh(FhhAmpl[ifhh], 'i');
				} else if (FhhAmpl[ifhh]->IsNuMinusNu()) {
					FhhIdAmpl[ifhh] = new TFhh(FhhAmpl[ifhh], 'm');
				} else {
					long found_partner = 0;
					for (long jfhh = 0; jfhh < NFhhAmpl; jfhh++) {
						if (FhhAmpl[ifhh]->GetLambda() == FhhAmpl[jfhh]->GetNu()
								&& FhhAmpl[ifhh]->GetNu()
										== FhhAmpl[jfhh]->GetLambda()) {
							found_partner = 1;
							FhhIdAmpl[ifhh] = new TFhh(FhhAmpl[ifhh],
									FhhAmpl[jfhh]);
							// ** continue here **
						}
					}
					if (!found_partner)
						cerr << "?!?! No partner for amplitude "
								<< FhhAmpl[ifhh]->GetName() << endl;
				}
			}
		}

	}

	cout << NFhhAmpl << " amplitudes: non-relativistic limit" << endl;
	for (long i = 0; i < NFhhAmpl; i++) {
		FhhAmpl[i]->NonRelLimit();
	}

	cout << "Check non-relativistic G's" << endl;
	for (long i = 0; i < NFhhAmpl; i++) {
		FhhAmpl[i]->PrintNRG();
	}

	return 0;
}

#if(0)
long TJSS::PrintHFILE() {
	char DecayName[10];
	sprintf(DecayName, "%ld%ld%ld%c%c", JMother, SDecay1, SDecay2,
			etaJ * eta1 * eta2 == -1 ? 'n' : 'p', ' ');
	char ofname[20];
	sprintf(ofname, "CalcAmpl-%s.h", DecayName);
	ofstream ofs(ofname);
	ofs << "// CalcAmpl output for " << DecayName << endl;
	ofs << "const int FhhAmpl_" << DecayName << "[] = { " << endl;
	ofs << "  " << NFhhAmpl << ",               // number of Fhh amplitudes"
			<< endl;
	for (int i = 0; i < NFhhAmpl; i++) {
		ofs << "  " << FhhAmpl[i]->GetJ() << ", " << FhhAmpl[i]->GetLambda()
				<< ", " << FhhAmpl[i]->GetNu() << ", "
				<< FhhAmpl[i]->GetEvenContr() << ",      // " << "F"
				<< FhhAmpl[i]->GetLambda() << FhhAmpl[i]->GetNu()
				<< ": J, lambda, nu, even_contr" << endl;
		ofs << "    " << FhhAmpl[i]->GetNterms()
				<< ",             // number of contributions" << endl;
		for (int j = 0; j < FhhAmpl[i]->GetNterms(); j++) {
			TLSContrib* lsa = FhhAmpl[i]->GetLStPtr()[j];
			ofs << "    " << lsa->GetJ() << ", " << lsa->GetL() << ", "
					<< lsa->GetS() << ", " << lsa->GetDelta() << ", "
					<< lsa->GetRunningNumber() << ",    // contr. F"
					<< FhhAmpl[i]->GetLambda() << FhhAmpl[i]->GetNu() << "-"
					<< j << ": J, L, S, delta, #[g,f,h,...]" << endl;
			ofs << "      " << lsa->GetNterms() << ", "
					<< lsa->GetNormFactor()->GetSign() << ", "
					<< lsa->GetNormFactor()->GetNumerator() << ", "
					<< lsa->GetNormFactor()->GetDenominator()
					<< ",     // number of terms, squared norm factor sign/N/D"
					<< endl;
			for (int k = 0; k < lsa->GetNterms(); k++) {
				ofs << "      " << lsa->GetTermFracNum()[k].GetSign() << ", "
						<< lsa->GetTermFracNum()[k].GetNumerator() << ", "
						<< lsa->GetTermFracNum()[k].GetDenominator() << ", "
						<< lsa->GetTermg1pot()[k] << ", "
						<< lsa->GetTermg2pot()[k]
						<< ",  // squared sign/N/D, exponents of g_s and g_sigma"
						<< endl;
			}
		}
	}
	ofs << "};" << endl;

	return 0;
}
#endif
