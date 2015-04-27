#include <TJSS.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstdlib>

using namespace std;

unsigned int TJSS::_debugLevel = 9999;

void TJSS::CalcAmpl() {
	if (_debugLevel >= 2) {
		cout << "  Decay channel:   " << _JMother;
		if (_etaJ > 0) {
			cout << "+  ->  ";
		} else {
			cout << "-  ->  ";
		}
		cout << _SDecay1;
		if (_eta1 > 0) {
			cout << "+ ";
		} else {
			cout << "- ";
		}
		cout << _SDecay2;
		if (_eta2 > 0) {
			cout << "+";
		} else {
			cout << "-";
		}
		cout << endl;
	}

	// range for coupled Spin
	const long Smin = abs(_SDecay1 - _SDecay2);
	const long Smax = _SDecay1 + _SDecay2;

	if (_debugLevel >= 2) {
		cout << "possible S:";
		long iS = Smin;
		while (iS <= Smax) {
			cout << " " << iS;
			iS++;
		}
		cout << endl;
	}

	const long intr_parity = _etaJ * _eta1 * _eta2;

	long Lmax = _JMother + Smax;
	long Lmin = Lmax;
	for (long iS = Smin; iS <= Smax; iS++) {
		long Lm1 = abs(_JMother - iS);
		if (Lm1 < Lmin) {
			Lmin = Lm1;
		}
	}

	if (_debugLevel >= 2) {
		cout << "possible L:";
	}
	long NL = 0;
	long *fL = 0;
	long testL = 0; // L even
	if (intr_parity < 0) {
		testL = 1; // L odd
	}
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
				if (_debugLevel >= 2) {
					cout << " " << fL[NL];
				}
				NL++;
			}
			tL += 2;
		}
	}
	if (_debugLevel >= 2) {
		cout << endl;
	}

	const bool evenContraction = not ((_SDecay1 + _SDecay2 + testL - _JMother) % 2);

	if (_debugLevel >= 2) {
		if (evenContraction) {
			cout << "contraction only with g~" << endl;
		} else {
			cout << "contraction with eps*p and g~" << endl;
		}
	}

	long MaxPsiInternal = _SDecay1;
	if (_SDecay2 < _SDecay1) {
		MaxPsiInternal = _SDecay2;
	}

	long MaxPsiPhi = _SDecay1 + _SDecay2;
	if (_JMother < MaxPsiPhi) {
		MaxPsiPhi = _JMother;
	}

	for (long il = 0; il < NL; il++) {

		long L = fL[il];

		long MaxPsiChi = _SDecay1 + _SDecay2;
		if (L < MaxPsiChi) {
			MaxPsiChi = L;
		}

		long MaxChiPhi = _JMother;
		if (L < _JMother) {
			MaxChiPhi = L;
		}

		// possible S in this case:
		long SminL = _JMother - L;
		if (SminL < 0) {
			SminL = -SminL;
		}
		if (SminL < Smin) {
			SminL = Smin;
		}
		long SmaxL = _JMother + L;
		if (SmaxL > Smax) {
			SmaxL = Smax;
		}
		for (long S_L = SminL; S_L <= SmaxL; S_L++) {
			if (_debugLevel >= 2)
				cout << "Amplitudes for L=" << L << " S=" << S_L
				     << "  Rank scheme [ " << _SDecay1 + _SDecay2 << " "
				     << L << " " << _JMother << "]" << endl;
			long totalRank = _SDecay1 + _SDecay2 + L + _JMother;
			long MaxDelta = S_L;
			if (_JMother < S_L) {
				MaxDelta = _JMother;
			}

			long IndexContractions = (totalRank - 3) / 2;
			if (evenContraction) {
				IndexContractions = totalRank / 2;
			}
			if (_debugLevel >= 2) {
				cout << IndexContractions << " Lorentz contractions." << endl;
			}

			long MaxContractionNumber = 0;

			for (long PsiInternal = 0; PsiInternal <= MaxPsiInternal; PsiInternal++) {
				for (long cChiPhi = 0; cChiPhi <= MaxChiPhi; cChiPhi++) {
					for (long cPsiChi = 0; cPsiChi <= MaxPsiChi; cPsiChi++) {
						for (long cPsiPhi = 0; cPsiPhi <= MaxPsiPhi; cPsiPhi++) {

							if(initialCheck_ASK_JAN_WHAT_THAT_SHOULD_BE_CALLED(PsiInternal, cPsiChi, cChiPhi, cPsiPhi, IndexContractions, evenContraction, L)) {
								continue;
							}

							pair<long, long> rPair;
							if(not getRomeAndREps_ASK_JAN_WHAT_THAT_SHOULD_BE_CALLED(PsiInternal, cPsiChi, cPsiPhi, evenContraction, rPair)) {
								continue;
							}
							const long r_ome = rPair.first;
							const long r_eps = rPair.second;
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
							if (_debugLevel >= 3) {
								cout << "{" << r_ome << "}";
							}
							for (long cChiOmega = 0; cChiOmega <= r_ome; cChiOmega++) {
								for (long cPhiOmega = 0; cPhiOmega <= r_ome - cChiOmega; cPhiOmega++) {
									const long cPhiEps = cPsiPhi - cPhiOmega;
									const long cChiEps = cPsiChi - cChiOmega;
									if (_debugLevel >= 3) {
										cout << "[" << cPhiEps << cChiEps << r_eps << "]";
									}
									if ( (cPhiEps < 0) or (cChiEps < 0) or (cPhiEps + cChiEps > r_eps) ) {
										continue;
									} else if (_debugLevel >= 3) {
										cout << "E+ OK" << endl;
									}
									if (_debugLevel >= 2) {
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
										cout << "Contraction pattern " << getContractionPattern(cPsiPhi,
										                                                        cPhiOmega,
										                                                        cPhiEps,
										                                                        PsiInternal,
										                                                        cPsiChi,
										                                                        cChiOmega,
										                                                        cChiEps,
										                                                        cChiPhi,
										                                                        L)
										     << endl;
									}
									for (long delta = 0; delta <= MaxDelta; delta++) {
										if (_debugLevel >= 2) {
											cout << " Constructing LS-Amplitude " << _LSAmplitudes.size() << endl;
										}
										unsigned int SameContr = 0;
										for (SameContr = 0; SameContr < _LSAmplitudes.size(); SameContr++) {
											if (_LSAmplitudes[SameContr]->CheckContraction(
											    L, S_L, PsiInternal,
											    cPsiChi, cChiPhi,
											    cPsiPhi, cPhiOmega,
											    cChiOmega, cPhiEps,
											    cChiEps))
											{
												break;
											}
										}
										long ContractionNumber = 0;
										if (SameContr < _LSAmplitudes.size()) {
											ContractionNumber = _LSAmplitudes[SameContr]->GetContraction();
										} else {
											ContractionNumber = MaxContractionNumber + 1;
										}

										TLSAmpl* ampl =
												new TLSAmpl(_SDecay1,
														_SDecay2, L,
														_JMother, delta,
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

										if (ampl->GetNterms()) {
											if (ContractionNumber > MaxContractionNumber) {
												MaxContractionNumber++;
											}
											_LSAmplitudes.push_back(ampl);
										} else {
											delete ampl;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	if (_debugLevel) {
		cout << _LSAmplitudes.size() << " LS-Amplitudes found to be non-zero." << endl;
		cout << "++++++++++++++++++++++++++++++++++++" << endl;
		cout << "+++ Helicity-coupling amplitudes +++" << endl;
		cout << "++++++++++++++++++++++++++++++++++++" << endl;
	}

	if(not _FhhAmpl.empty()) {
		cerr << "_FhhAmpl not empty. Aborting..." << endl;
		throw;
	}

	for (long lambda = 0; lambda <= _SDecay1; lambda++) {
		for (long nu = -_SDecay2; nu <= _SDecay2; nu++) {
			//    for (long nu=-SDecay1; nu<=SDecay2; nu++) { bug!!!! found 4.3.08
			if (lambda == 0 && nu < 0) {
				continue;
			}
			TFhh* fhh = new TFhh(_JMother, _SDecay1, _SDecay2, lambda, nu, _LSAmplitudes, evenContraction);
			if (fhh->GetNterms()) {
				_FhhAmpl.push_back(fhh);
			} else {
				delete fhh;
			}
		}
	}

	if (_debugLevel) {
		cout << _FhhAmpl.size() << " non-zero helicity-coupling amplitudes" << endl;
	}

	_FhhIdAmpl = vector<TFhh*>((_SDecay1 + 1) * (_SDecay2 + 1), 0);

	if (_SDecay1 == _SDecay2 && _eta1 == _eta2) {
		if (_debugLevel) {
			cout << endl << " for identical-particle decay:" << endl;
		}
		for (unsigned int ifhh = 0; ifhh < _FhhAmpl.size(); ifhh++) {
			if (_FhhAmpl[ifhh]) {
				if (_FhhAmpl[ifhh]->IsNuNu()) {
					_FhhIdAmpl[ifhh] = new TFhh(_FhhAmpl[ifhh], 'i');
				} else if (_FhhAmpl[ifhh]->IsNuMinusNu()) {
					_FhhIdAmpl[ifhh] = new TFhh(_FhhAmpl[ifhh], 'm');
				} else {
					long found_partner = 0;
					for (unsigned int jfhh = 0; jfhh < _FhhAmpl.size(); jfhh++) {
						if ( (_FhhAmpl[ifhh]->GetLambda() == _FhhAmpl[jfhh]->GetNu()) and
						     (_FhhAmpl[ifhh]->GetNu() == _FhhAmpl[jfhh]->GetLambda()) )
						{
							found_partner = 1;
							_FhhIdAmpl[ifhh] = new TFhh(_FhhAmpl[ifhh], _FhhAmpl[jfhh]);
							// ** continue here **
						}
					}
					if (!found_partner) {
						cerr << "?!?! No partner for amplitude " << _FhhAmpl[ifhh]->GetName() << endl;
					}
				}
			}
		}

	}

	cout << _FhhAmpl.size() << " amplitudes: non-relativistic limit" << endl;
	for (unsigned int i = 0; i < _FhhAmpl.size(); i++) {
		_FhhAmpl[i]->NonRelLimit();
	}

	cout << "Check non-relativistic G's" << endl;
	for (unsigned int i = 0; i < _FhhAmpl.size(); i++) {
		_FhhAmpl[i]->PrintNRG();
	}

}


bool TJSS::initialCheck_ASK_JAN_WHAT_THAT_SHOULD_BE_CALLED(const long& PsiInternal,
                                                           const long& cPsiChi,
                                                           const long& cChiPhi,
                                                           const long& cPsiPhi,
                                                           const long& IndexContractions,
                                                           const long& even_contraction,
                                                           const long& L) const
{
	if (_debugLevel >= 3) {
		cout << "Checking " << PsiInternal << " "
		     << cPsiChi << " " << cChiPhi << " "
		     << cPsiPhi; // << endl;
	}
	if ( (PsiInternal + cPsiChi + cChiPhi + cPsiPhi) != IndexContractions) {
		if (_debugLevel >= 3) {
			cout << " C-" << endl;
		}
		return true;
	} else if (_debugLevel >= 3) {
		cout << " C+";
	}
	if (even_contraction) {
		if ( (2 * PsiInternal + cPsiChi + cPsiPhi) != _SDecay1 + _SDecay2) {
			if (_debugLevel >= 3) {
				cout << "S-" << endl;
			}
			return true;
		} else if (_debugLevel >= 3) {
			cout << "S+";
		}
		if (cPsiChi + cChiPhi != L) {
			if (_debugLevel >= 3) {
				cout << "L-" << endl;
			}
			return true;
		} else if (_debugLevel >= 3) {
			cout << "L+";
		}
		if (cChiPhi + cPsiPhi != _JMother) {
			if (_debugLevel >= 3) {
				cout << "J-" << endl;
			}
			return true;
		} else if (_debugLevel >= 3) {
			cout << "J+";
		}
	} else {
		if (L - cPsiChi - cChiPhi > 1) {
			if (_debugLevel >= 3) {
				cout << "L-" << endl;
			}
			return true;
		} else if (_debugLevel >= 3) {
			cout << "L+";
		}
		if (_JMother - cPsiPhi - cChiPhi > 1) {
			if (_debugLevel >= 3) {
				cout << "J-" << endl;
			}
			return true;
		} else if (_debugLevel >= 3) {
			cout << "J+";
		}
	}
	return false;
}


bool TJSS::getRomeAndREps_ASK_JAN_WHAT_THAT_SHOULD_BE_CALLED(const long& PsiInternal,
                                                             const long& cPsiChi,
                                                             const long& cPsiPhi,
                                                             const bool& evenContraction,
                                                             pair<long, long>& rPair) const
{

	long r_ome = _SDecay1 - PsiInternal;
	long r_eps = _SDecay2 - PsiInternal;

	if (!evenContraction) {
		long PsiRest = _SDecay1 + _SDecay2 - 2 * PsiInternal - cPsiChi - cPsiPhi;
		if (PsiRest < 0 or PsiRest > 2) {
			if (_debugLevel >= 3) {
				cout << "R-" << endl;
			}
			return false;
		} else if (_debugLevel >= 3) {
			cout << "R+";
		}
		if (PsiRest == 2) {
			if (r_ome > 0) {
				r_ome--;
			} else {
				if (_debugLevel >= 3) {
					cout << "O-";
				}
				return false;
			}
			if (r_eps > 0) {
				r_eps--;
			} else {
				if (_debugLevel >= 3) {
					cout << "E-";
				}
				return false;
			}
		}
	}
	rPair.first = r_ome;
	rPair.second = r_eps;
	return true;
}


string TJSS::getContractionPattern(const long& cPsiPhi,
                                   const long& cPhiOmega,
                                   const long& cPhiEps,
                                   const long& PsiInternal,
                                   const long& cPsiChi,
                                   const long& cChiOmega,
                                   const long& cChiEps,
                                   const long& cChiPhi,
                                   const long& L) const
{
	stringstream sstr;
	long cc = cPsiPhi;
	while (cc--) {
		sstr << "#";
	}
	if (cPsiPhi) {
		sstr << "(";
	}
	cc = cPhiOmega;
	while (cc--) {
		sstr << "o";
	}
	cc = cPhiEps;
	while (cc--) {
		sstr << "e";
	}
	if (cPsiPhi) {
		sstr << ")";
	}
	sstr << " " << _SDecay1 + _SDecay2;
	cc = PsiInternal;
	while (cc--) {
		sstr << "'";
	}
	sstr << " ";
	cc = cPsiChi;
	while (cc--) {
		sstr << "#";
	}
	if (cPsiChi) {
		sstr << "(";
	}
	cc = cChiOmega;
	while (cc--) {
		sstr << "o";
	}
	cc = cChiEps;
	while (cc--) {
		sstr << "e";
	}
	if (cPsiChi) {
		sstr << ")";
	}
	sstr << " " << L << " ";
	cc = cChiPhi;
	while (cc--) {
		sstr << "#";
	}
	sstr << " " << _JMother << " ";
	cc = cPsiPhi;
	while (cc--) {
		sstr << "#";
	}
	if (cPsiPhi) {
		sstr << "(";
	}
	cc = cPhiOmega;
	while (cc--) {
		sstr << "o";
	}
	cc = cPhiEps;
	while (cc--) {
		sstr << "e";
	}
	if (cPsiPhi) {
		sstr << ")";
	}
	return sstr.str();
}


#if(0)
long TJSS::PrintHFILE() {
	char DecayName[10];
	sprintf(DecayName, "%ld%ld%ld%c%c", _JMother, _SDecay1, _SDecay2,
			_etaJ * _eta1 * _eta2 == -1 ? 'n' : 'p', ' ');
	char ofname[20];
	sprintf(ofname, "CalcAmpl-%s.h", DecayName);
	ofstream ofs(ofname);
	ofs << "// CalcAmpl output for " << DecayName << endl;
	ofs << "const int FhhAmpl_" << DecayName << "[] = { " << endl;
	ofs << "  " << _NFhhAmpl << ",               // number of Fhh amplitudes"
			<< endl;
	for (int i = 0; i < _NFhhAmpl; i++) {
		ofs << "  " << _FhhAmpl[i]->GetJ() << ", " << _FhhAmpl[i]->GetLambda()
				<< ", " << _FhhAmpl[i]->GetNu() << ", "
				<< _FhhAmpl[i]->GetEvenContr() << ",      // " << "F"
				<< _FhhAmpl[i]->GetLambda() << _FhhAmpl[i]->GetNu()
				<< ": J, lambda, nu, even_contr" << endl;
		ofs << "    " << _FhhAmpl[i]->GetNterms()
				<< ",             // number of contributions" << endl;
		for (int j = 0; j < _FhhAmpl[i]->GetNterms(); j++) {
			TLSContrib* lsa = _FhhAmpl[i]->GetLStPtr()[j];
			ofs << "    " << lsa->GetJ() << ", " << lsa->GetL() << ", "
					<< lsa->GetS() << ", " << lsa->GetDelta() << ", "
					<< lsa->GetRunningNumber() << ",    // contr. F"
					<< _FhhAmpl[i]->GetLambda() << _FhhAmpl[i]->GetNu() << "-"
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
