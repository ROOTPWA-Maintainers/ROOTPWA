#include "TLSAmpl.h"

#include <iostream>
#include <string>

#include <reportingUtils.hpp>

using namespace std;

unsigned int TLSAmpl::_debugLevel = 1;

TLSAmpl::TLSAmpl(const long& RankS1,
                 const long& RankS2,
                 const long& RankL,
                 const long& RankJ,
                 const long& delta,
                 const long& S_L,
                 const long& cPsiInt,
                 const long& cPsiChi,
                 const long& cChiPhi,
                 const long& cPsiPhi,
                 const long& cPhiOme,
                 const long& cChiOme,
                 const long& cPhiEps,
                 const long& cChiEps,
                 const long& contractionNumber)
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

	if ( (    evenContraction and contractions != totalR) or
	     (not evenContraction and contractions != totalR - 3))
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

	TTensorSum TSS = WFS1.GetSpinCoupledTensorSum(WFS2, _delta, S_L);
	if (_debugLevel >= 2) {
		TSS.Print();
	}

	if (cPsiInt) {
		long ires = TSS.SpinInnerContraction(cPsiInt);
		if (ires == 0) {
			if (_debugLevel) {
				cout << "Inner contraction is zero." << endl;
			}
			_Nterms = 0;
			return;
		}
		if (_debugLevel >= 2) {
			TSS.Print();
		}
	}

	TSpinWaveFunction WFL(RankL, 'l');
	TTensorSum TSL = WFL.GetTensorSum('c', 0);
	if (_debugLevel >= 2) {
		TSL.Print();
	}

	TTensorSum TSLS = TSS.LSContraction(TSL, cPsiChi, cChiOme, cChiEps, 'c');
	if (TSLS.GetNterms() == 0) {
		if (_debugLevel) {
			cout << "LS contraction is zero." << endl;
		}
		_Nterms = 0;
		return;
	}
	if (_debugLevel >= 2) {
		TSLS.Print();
	}

	TSpinWaveFunction WFJ(RankJ, 'c');
	TTensorSum TSJ = WFJ.GetTensorSum('p', _delta);
	if (_debugLevel >= 2) {
		TSJ.Print();
	}

	TTensorSum TSLSJ = TSLS.LSContraction(TSJ, cPsiPhi, cPhiOme, cPhiEps, 'p');

	if (TSLSJ.GetNterms() == 0) {
		if (_debugLevel) {
			cout << "JS contraction is zero." << endl;
		}
		_Nterms = 0;
		return;
	}

	if (_debugLevel >= 2) {
		TSLSJ.Print();
	}

	_TSScalar = TSLSJ.LJContraction(cChiPhi, evenContraction);

	if (_debugLevel >= 2) {
		TSLSJ.Print();
	}

	_Nterms = _TSScalar.GetNterms();

	if (_Nterms == 0) {
		if (_debugLevel) {
			cout << "LJ contraction is zero." << endl;
		}
		return;
	}

	if (_debugLevel) {
		//TSScalar.Print();
		_TSScalar.Print('s');
	}
}


bool TLSAmpl::CheckContraction(const long& L,
                               const long& S,
                               const long& cPsI,
                               const long& cPsC,
                               const long& cCP,
                               const long& cPsP,
                               const long& cPO,
                               const long& cCO,
                               const long& cPE,
                               const long& cCE) const
{
	if (_L != L          or
	    _S != S          or
	    _cPsiInt != cPsI or
	    _cPsiChi != cPsC or
	    _cChiPhi != cCP  or
	    _cPsiPhi != cPsP or
	    _cPhiOme != cPO  or
	    _cChiOme != cCO  or
	    _cPhiEps != cPE  or
	    _cChiEps != cCE)
	{
		return false;
	}
	return true;
}
