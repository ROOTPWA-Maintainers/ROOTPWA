#include "TLSAmpl.h"

#include <iostream>
#include <string>

#include <reportingUtils.hpp>

using namespace std;

unsigned int TLSAmpl::_debugLevel = 1;

TLSAmpl::TLSAmpl(long RankS1,
                 long RankS2,
                 long RankL,
                 long RankJ,
                 long delta,
                 long S_L,
                 long cPsiInt,
                 long cPsiChi,
                 long cChiPhi,
                 long cPsiPhi,
                 long cPhiOme,
                 long cChiOme,
                 long cPhiEps,
                 long cChiEps,
                 long contractionNumber)
	: _J(RankJ),
	  _L(RankL),
	  _S(S_L),
	  _delta(delta),
	  _contractionNumber(contractionNumber),
	  _cPsiInt(cPsiInt),
	  _cPsiChi(cPsiChi),
	  _cChiPhi(cChiPhi),
	  _cPsiPhi(cPsiPhi),
	  _cPhiOme(cPhiOme),
	  _cChiOme(cChiOme),
	  _cPhiEps(cPhiEps),
	  _cChiEps(cChiEps)
{
	const long totalR = RankS1 + RankS2 + RankL + RankJ;
	const long contractions = 2 * (cPsiInt + cPsiChi + cChiPhi + cPsiPhi);
	const bool evenContraction = (not (totalR % 2));

	if ( (evenContraction == 1 and contractions != totalR) or
	     (evenContraction == 0 and contractions != totalR - 3))
	{
		printErr << "Invalid contraction occurred. Aborting..." << endl;
		throw;
	}

	if (_debugLevel) {
		cout << "LSAmpl: " << RankS1 << " " << RankS2 << " L=" << RankL << " "
		     << RankJ << " d=" << _delta << " S=" << S_L << " c: " << cPsiInt
		     << " " << cPsiChi << " " << cChiPhi << " " << cPsiPhi << " s: "
		     << cPhiOme << " " << cChiOme << " " << cPhiEps << " "
		     << cChiEps;
		if (not evenContraction) {
			cout << " (iw)";
		}
		cout << endl;
	}

	TSpinWaveFunction WFS1(RankS1, 's');
	TSpinWaveFunction WFS2(RankS2, 's');

	TTensorSum* TSS = WFS1.GetSpinCoupledTensorSum(&WFS2, _delta, S_L);
	if (_debugLevel >= 2) {
		TSS->Print();
	}

	if (cPsiInt) {
		long ires = TSS->SpinInnerContraction(cPsiInt);
		if (ires == 0) {
			if (_debugLevel) {
				cout << "Inner contraction is zero." << endl;
			}
			_Nterms = 0;
			return;
		}
		if (_debugLevel >= 2) {
			TSS->Print();
		}
	}

	TSpinWaveFunction WFL(RankL, 'l');
	TTensorSum *TSL = WFL.GetTensorSum('c', 0);
	if (_debugLevel >= 2) {
		TSL->Print();
	}

	TTensorSum *TSLS = TSS->LSContraction(TSL, cPsiChi, cChiOme, cChiEps, 'c');
	if (TSLS->GetNterms() == 0) {
		if (_debugLevel) {
			cout << "LS contraction is zero." << endl;
		}
		_Nterms = 0;
		return;
	}
	if (_debugLevel >= 2) {
		TSLS->Print();
	}

	TSpinWaveFunction WFJ(RankJ, 'c');
	TTensorSum* TSJ = WFJ.GetTensorSum('p', _delta);
	if (_debugLevel >= 2) {
		TSJ->Print();
	}

	TTensorSum* TSLSJ = TSLS->LSContraction(TSJ, cPsiPhi, cPhiOme, cPhiEps, 'p');

	if (TSLSJ->GetNterms() == 0) {
		if (_debugLevel) {
			cout << "JS contraction is zero." << endl;
		}
		_Nterms = 0;
		return;
	}

	if (_debugLevel >= 2) {
		TSLSJ->Print();
	}

	_TSScalar = TSLSJ->LJContraction(cChiPhi, evenContraction);

	if (_debugLevel >= 2) {
		TSLSJ->Print();
	}

	_Nterms = _TSScalar->GetNterms();

	if (_Nterms == 0) {
		if (_debugLevel) {
			cout << "LJ contraction is zero." << endl;
		}
		return;
	}

	if (_debugLevel) {
		//TSScalar->Print();
		_TSScalar->Print('s');
	}
}


bool TLSAmpl::CheckContraction(long L_,
                               long S_,
                               long cPsI_,
                               long cPsC_,
                               long cCP_,
                               long cPsP_,
                               long cPO_,
                               long cCO_,
                               long cPE_,
                               long cCE_) const
{
	if (_L != L_          or
	    _S != S_          or
	    _cPsiInt != cPsI_ or
	    _cPsiChi != cPsC_ or
	    _cChiPhi != cCP_  or
	    _cPsiPhi != cPsP_ or
	    _cPhiOme != cPO_  or
	    _cChiOme != cCO_  or
	    _cPhiEps != cPE_  or
	    _cChiEps != cCE_)
	{
		return false;
	}
	return true;
}
