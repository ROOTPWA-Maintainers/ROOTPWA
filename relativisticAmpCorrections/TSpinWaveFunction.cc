#include "TSpinWaveFunction.h"

#include <iostream>
#include <string>

#include <reportingUtils.hpp>

#include "ClebschGordanBox.h"

using namespace std;

unsigned int TSpinWaveFunction::_debugSpinWave = 0;

TSpinWaveFunction::TSpinWaveFunction(long J_, char type_) {

	_J = J_;
	_type = type_;

	_max_pzm = 1;
	_pot3 = new long[_J];
	for (long i = 0; i < _J; i++) {
		_pot3[i] = _max_pzm;
		_max_pzm *= 3;
	}
	_mi = new long[_max_pzm * _J];
	_M = new long[_max_pzm];

	_coeff = new TFracNum[_max_pzm];
	for (long pzm = 0; pzm < _max_pzm; pzm++) {
		_coeff[pzm] = TFracNum::Zero;
	}

	for (long pzm = 0; pzm < _max_pzm; pzm++) {
		long rest = pzm;
		long i = _J - 1;
		_M[pzm] = 0;
		while (i >= 0) {
			long ipot3 = rest / _pot3[i];
			_mi[_J * pzm + i] = ipot3 - 1;
			_M[pzm] += _mi[_J * pzm + i];
			//cout << "Setting mi["<<J<<"*"<<pzm<<"+"<<i<<"]="<<mi[J*pzm+i]<<endl;
			rest -= ipot3 * _pot3[i];
			i--;
		}
	}

	for (long MM = _J; MM >= -_J; MM--) {
		if (_debugSpinWave >= 2) {
			cout << "Phi(" << _J << "," << MM << ")=" << endl;
		}
		for (long pzm = _max_pzm - 1; pzm >= 0; pzm--) {
			if (((_type == 's' or  _type == 'c') and _M[pzm] == MM) or
			     (_type == 'l' and _M[pzm] == MM and MM == 0)) {
				long m0 = 0;
				for (long i = 0; i < _J; i++) {
					if (_debugSpinWave >= 2) {
						if (_mi[_J * pzm + i] == 1) {
							cout << "+";
						}
						if (_mi[_J * pzm + i] == 0) {
							cout << "0";
						}
						if (_mi[_J * pzm + i] == -1) {
							cout << "-";
						}
					}
					if (_mi[_J * pzm + i] == 0) {
						m0++;
					}
				}
				if (_type == 's' || _type == 'c') {
					_coeff[pzm] = TFracNum::am0_to_J(_J, MM, m0);
				}
				if (_type == 'l') {
					_coeff[pzm] = TFracNum::cm0_sub_ell_2(_J, m0);
				}
				if (_debugSpinWave >= 2) {
					cout << " * sqrt("
					     << _coeff[pzm].FracString()
					     << ") (eq. 4.1) [" << pzm << "]" << endl;
				}
			}
		}
	}
}

TTensorSum
TSpinWaveFunction::GetTensorSum(char name, long delta) {

	TTensorSum ts;

	for (long pzm = _max_pzm - 1; pzm >= 0; pzm--) {
		if (_M[pzm] == delta) {
			vector<long> pzm_field(_J);
			for (long i = 0; i < _J; i++) {
				pzm_field[i] = _mi[_J * pzm + i];
			}
			if (not (_coeff[pzm] == TFracNum::Zero)) {
				ts.AddTerm(TTensorTerm(name, pzm_field, _coeff[pzm]));
			}
		}
	}
	return ts;
}

TTensorSum
TSpinWaveFunction::GetSpinCoupledTensorSum(const TSpinWaveFunction& E,
                                           const long& delta,
                                           const long& S) {

	if (!(_type == 's' and E._type == 's')) {
		printErr << "GetSpinCoupledTensorSum only for spin-type wave functions!" << endl;
		throw;
	}

	long twoS1 = 2 * _J;
	long twoS2 = 2 * E._J;
	long twoS = 2 * S;

	long twoSmin = twoS1 - twoS2;
	if (twoSmin < 0) {
		twoSmin = -twoSmin;
	}

	if (twoS < twoSmin || twoS > twoS1 + twoS2) {
		printErr << "GetSpinCoupledTensorSum: no coupling "
		         << _J << " + " << E._J << " -> " << S << endl;
		throw;
	}

	//TFracNum *S1S2 = ClebschGordan(twoS, twoS1, twoS2);
	TFracNum* S1S2 = ClebschGordanBox::instance()->GetCG(S, _J, E._J);
	if (_debugSpinWave >= 1) {
		cout << "Clebsch-Gordans calculated for "
		     << twoS << "," << twoS1 << "," << twoS2 << ": " << S1S2 << endl;
	}

	TTensorSum ts;

	for (long MM1 = _J; MM1 >= -_J; MM1--) {
		for (long pzm1 = _max_pzm - 1; pzm1 >= 0; pzm1--) {
			if (_M[pzm1] == MM1 and not (_coeff[pzm1] == TFracNum::Zero)) {

				for (long MM2 = E._J; MM2 >= -_J; MM2--) {
					for (long pzm2 = E._max_pzm - 1; pzm2 >= 0; pzm2--) {
						if (E._M[pzm2] == MM2 and not (E._coeff[pzm2] == TFracNum::Zero)) {

							if (MM1 + MM2 == delta) {

								vector<long> pzm1_field(_J);
								for (long i = 0; i < _J; i++) {
									pzm1_field[i] = _mi[_J * pzm1 + i];
								}

								vector<long> pzm2_field(E._J);
								for (long i = 0; i < E._J; i++) {
									pzm2_field[i] = E._mi[E._J * pzm2 + i];
								}

								TTensorTerm newterm('o', pzm1_field, _coeff[pzm1]);

								TFracNum coeff2_CG = S1S2[ClebschGordanBox::CGIndex(_J, MM1, E._J, MM2)];
								coeff2_CG = coeff2_CG * E._coeff[pzm2];

								newterm.Multiply('e', pzm2_field, coeff2_CG);

								ts.AddTerm(newterm);
							}
						}
					}
				}
			}
		}
	}
	return ts;
}

long
TSpinWaveFunction::CheckCGFormula() {

	TFracNum **J1J = new TFracNum*[_J - 1];
	for (long i = 0; i < _J - 1; i++) {
		// long twoJ=2*(i+2);
		// J1J[i] = ClebschGordan(twoJ, 2, twoJ-2);
		J1J[i] = ClebschGordanBox::instance()->GetCG(i + 2, 1, i + 1);
	}

	TFracNum *coeffCG = new TFracNum[_max_pzm];

	for (long pzm = 0; pzm < _max_pzm; pzm++) {
		coeffCG[pzm] = TFracNum::One;
		long m = _mi[_J * pzm];
		if (_debugSpinWave >= 2) {
			cout << m << " x ";
		}
		for (long jj = 1; jj < _J; jj++) {
			long mj = _mi[_J * pzm + jj];
			// m1=i/max2-J1, m2=i%max2-J2 ==> i=(m1+J1)*max2
			if (_debugSpinWave >= 2) {
				cout << mj << ": * CG[" << jj - 1 << "]["
				     << (mj + 1) * (2 * jj + 1) + jj + m
				     << "]:" << J1J[jj - 1][(mj + 1) * (2 * jj + 1) + jj + m].Dval()
				     << endl;
			}
			coeffCG[pzm] = coeffCG[pzm] * J1J[jj - 1][(mj + 1) * (2 * jj + 1) + jj + m];
			m += mj;
			if (_debugSpinWave >= 2) {
				cout << m << " x ";
			}
		}
		if (_debugSpinWave >= 2) {
			cout << "done." << endl;
		}
	}

	for (long MM = _J; MM >= 0; MM--) {
		cout << "Phi(" << _J << "," << MM << ")=" << endl;
		for (long pzm = _max_pzm - 1; pzm >= 0; pzm--) {
			if (_M[pzm] == MM) {
				cout << coeffCG[pzm].FracString() << " * (";
				long m0 = 0;
				for (long i = 0; i < _J; i++) {
					if (_mi[_J * pzm + i] == 1) {
						cout << "+";
					}
					if (_mi[_J * pzm + i] == 0) {
						cout << "0";
						m0++;
					}
					if (_mi[_J * pzm + i] == -1) {
						cout << "-";
					}
				}
				cout << ")   [ (4.1): " << _coeff[pzm].FracString() << " ]" << endl;
			}
		}
	}
	return 0;
}
