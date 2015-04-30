#include "TSpinWaveFunction.h"

#include <iostream>
#include <string>

#include <reportingUtils.hpp>

#include "ClebschGordanBox.h"

using namespace std;

unsigned int TSpinWaveFunction::_debugSpinWave = 0;

TSpinWaveFunction::TSpinWaveFunction(const size_t& J, const char& type)
	: _J(J),
	  _type(type)
{

	vector<long> pot3(_J);
	{
		long max_pzm = 1;
		for (size_t i = 0; i < _J; i++) {
			pot3[i] = max_pzm;
			max_pzm *= 3;
		}
		_mi = vector<long>(max_pzm * _J);
		_M_and_coeff = vector<std::pair<long, TFracNum> >(max_pzm);
	}

	for(size_t pzm = 0; pzm < _M_and_coeff.size(); ++pzm) {
		long rest = pzm;
		long i = _J - 1;
		long M = 0;
		while (i >= 0) {
			long ipot3 = rest / pot3[i];
			_mi[_J * pzm + i] = ipot3 - 1;
			M += _mi[_J * pzm + i];
			//cout << "Setting mi["<<J<<"*"<<pzm<<"+"<<i<<"]="<<mi[J*pzm+i]<<endl;
			rest -= ipot3 * pot3[i];
			i--;
		}
		_M_and_coeff[pzm].first  = M;
		_M_and_coeff[pzm].second = TFracNum::Zero;
	}

	for (long MM = _J; MM >= -(long)_J; MM--) {
		if (_debugSpinWave >= 2) {
			cout << "Phi(" << _J << "," << MM << ")=" << endl;
		}
		for (long pzm = (long)_M_and_coeff.size() - 1; pzm >= 0; pzm--) {
			if (((_type == 's' or  _type == 'c') and _M_and_coeff[pzm].first == MM) or
			     (_type == 'l' and _M_and_coeff[pzm].first == MM and MM == 0)) {
				long m0 = 0;
				for (size_t i = 0; i < _J; i++) {
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
					_M_and_coeff[pzm].second = TFracNum::am0_to_J(_J, MM, m0);
				}
				if (_type == 'l') {
					_M_and_coeff[pzm].second = TFracNum::cm0_sub_ell_2(_J, m0);
				}
				if (_debugSpinWave >= 2) {
					cout << " * sqrt("
					     << _M_and_coeff[pzm].second.FracString()
					     << ") (eq. 4.1) [" << pzm << "]" << endl;
				}
			}
		}
	}
}

TTensorSum
TSpinWaveFunction::GetTensorSum(const char& name, const long& delta) {
	TTensorSum ts;
	for (long pzm = (long)_M_and_coeff.size() - 1; pzm >= 0; pzm--) {
		if (_M_and_coeff[pzm].first == delta) {
			vector<long> pzm_field(_J);
			for (size_t i = 0; i < _J; i++) {
				pzm_field[i] = _mi[_J * pzm + i];
			}
			if (not (_M_and_coeff[pzm].second == TFracNum::Zero)) {
				ts.AddTerm(TTensorTerm(name, pzm_field, _M_and_coeff[pzm].second));
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

	for (long MM1 = _J; MM1 >= -(long)_J; MM1--) {
		for (long pzm1 = (long)_M_and_coeff.size() - 1; pzm1 >= 0; pzm1--) {
			if (_M_and_coeff[pzm1].first == MM1 and not (_M_and_coeff[pzm1].second == TFracNum::Zero)) {
				for (long MM2 = E._J; MM2 >= -(long)_J; MM2--) {
					for (long pzm2 = (long)E._M_and_coeff.size() - 1; pzm2 >= 0; pzm2--) {
						if (E._M_and_coeff[pzm2].first == MM2 and not (E._M_and_coeff[pzm2].second == TFracNum::Zero)) {

							if (MM1 + MM2 == delta) {

								vector<long> pzm1_field(_J);
								for (size_t i = 0; i < _J; i++) {
									pzm1_field[i] = _mi[_J * pzm1 + i];
								}

								vector<long> pzm2_field(E._J);
								for (size_t i = 0; i < E._J; i++) {
									pzm2_field[i] = E._mi[E._J * pzm2 + i];
								}

								TTensorTerm newterm('o', pzm1_field, _M_and_coeff[pzm1].second);

								TFracNum coeff2_CG = S1S2[ClebschGordanBox::CGIndex(_J, MM1, E._J, MM2)];
								coeff2_CG = coeff2_CG * E._M_and_coeff[pzm2].second;

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

	if(_J == 0) {
		printErr << "J = 0, that is bad..." << endl;
		throw;
	}

	TFracNum **J1J = new TFracNum*[_J - 1];
	for (size_t i = 0; i < _J - 1; i++) {
		// long twoJ=2*(i+2);
		// J1J[i] = ClebschGordan(twoJ, 2, twoJ-2);
		J1J[i] = ClebschGordanBox::instance()->GetCG(i + 2, 1, i + 1);
	}

	TFracNum *coeffCG = new TFracNum[_M_and_coeff.size()];

	for (size_t pzm = 0; pzm < _M_and_coeff.size(); pzm++) {
		coeffCG[pzm] = TFracNum::One;
		long m = _mi[_J * pzm];
		if (_debugSpinWave >= 2) {
			cout << m << " x ";
		}
		for (size_t jj = 1; jj < _J; jj++) {
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
		for (size_t pzm = _M_and_coeff.size() - 1; pzm >= 0; pzm--) {
			if (_M_and_coeff[pzm].first == MM) {
				cout << coeffCG[pzm].FracString() << " * (";
				long m0 = 0;
				for (size_t i = 0; i < _J; i++) {
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
				cout << ")   [ (4.1): " << _M_and_coeff[pzm].second.FracString() << " ]" << endl;
			}
		}
	}
	return 0;
}
