#include <iostream>
#include <string>
#include "ClebschGordanBox.h"
using namespace std;

unsigned int ClebschGordanBox::debugCG = 0;
unsigned int ClebschGordanBox::debugCGBox = 0;

ClebschGordanBox* ClebschGordanBox::_instance = 0;

ClebschGordanBox* ClebschGordanBox::instance()
{
	if(not _instance) {
		_instance = new ClebschGordanBox();
	}
	return _instance;
}


TFracNum*
ClebschGordanBox::GetCG(const long& J, const long& J1, const long& J2) {
	for (long i = 0; i < NCG; i++) {
		if (J == CGJ[i] and J1 == CGJ1[i] and J2 == CGJ2[i]) {
			return CG[i];
		}
	}
	if (debugCGBox) {
		cout << "Registering Clebsch-Gordans for J=" << J << " J1=" << J1
		     << " J2=" << J2 << endl;
	}
	NCG++;
	long *newJ = new long[NCG];
	long *newJ1 = new long[NCG];
	long *newJ2 = new long[NCG];
	TFracNum* *newCG = new TFracNum*[NCG];
	for (long i = 0; i < NCG - 1; i++) {
		newJ[i] = CGJ[i];
		newJ1[i] = CGJ1[i];
		newJ2[i] = CGJ2[i];
		newCG[i] = CG[i];
	}
	newJ[NCG - 1] = J;
	newJ1[NCG - 1] = J1;
	newJ2[NCG - 1] = J2;
	newCG[NCG - 1] = ClebschGordan(2 * J, 2 * J1, 2 * J2);

	if (NCG - 1) {
		delete[] CGJ;    // wrong delete corrected 27.10.08
		delete[] CGJ1;
		delete[] CGJ2;
		delete[] CG;
	}

	CGJ = newJ;
	CGJ1 = newJ1;
	CGJ2 = newJ2;
	CG = newCG;

	return newCG[NCG - 1];
}

TFracNum*
ClebschGordanBox::ClebschGordan(long twoJ, long twoJ1, long twoJ2) {

	if (twoJ < 0 or twoJ1 < 0 or twoJ2 < 0 or twoJ < twoJ1 - twoJ2 or
	    twoJ < twoJ2 - twoJ1  or twoJ > twoJ1 + twoJ2)
	{
		cerr << "Clebsch-Gordan-Coefficient only for positive J's"
		     << " with |J1-J2|<=J<=J1+J2" << endl;
	}

	const long max1 = twoJ1 + 1;
	const long max2 = twoJ2 + 1;

	TFracNum* CG = new TFracNum[max1 * max2];

	// first non-zero element
	long mintwom1 = -twoJ1;
	if (twoJ + twoJ2 < twoJ1) {
		mintwom1 = -twoJ - twoJ2;
	}
	const long minIndex = (mintwom1 + twoJ1) / 2 * max2 + max2 - 1;
	CG[minIndex] = TFracNum::One;
	if (debugCG == 2) {
		cout << "#############################################" << endl
		     << " first row CG[" << minIndex << "]:"            << endl
		     << CG[minIndex];
	}

	for (long twom1 = -twoJ1; twom1 <= twoJ1; twom1 += 2) {
		const long twoi1 = twom1 + twoJ1;
		long twom2 = 0;
		long twom = 0;
		long twoi2 = 0;
		long twoi = 0;
		long i = 0;
		for (twom2 = twoJ2 - 2; twom2 >= -twoJ2; twom2 -= 2) {
			if (debugCG == 2) {
				cout << " Setting twom1=" << twom1 << ", twom2=" << twom2 << endl;
			}
			twom = twom1 + twom2;
			twoi2 = twom2 + twoJ2;
			twoi = twoi1 * max2 + twoi2;
			i = twoi / 2;
			if (twom1 + twom2 < -twoJ or twom1 + twom2 > twoJ) {
				CG[i] = TFracNum::Zero;
				if (debugCG == 2) {
					cout << "#############################################" << endl
					     << "zero CG[" << i << "]:"                         << endl
					     << CG[i];
				}
			} else if (twom1 == mintwom1) {
				CG[i] = CG[i + 1] * TFracNum((twoJ - twom) * (twoJ + twom + 2),
				                             (twoJ2 - twom2) * (twoJ2 + twom2 + 2));
				if (debugCG == 2) {
					cout << "#############################################" << endl
					     << " first row CG[" << i << "]:"                   << endl
					     << CG[i];
				}
			} else if (twom1 + twom2 != -twoJ) {
				const long Denom = (twoJ1 + twom1) * (twoJ1 - twom1 + 2);
				TFracNum cgp = CG[i - max2]     * TFracNum((twoJ + twom)   * (twoJ - twom + 2),   Denom);
				TFracNum cgm = CG[i - max2 + 1] * TFracNum((twoJ2 - twom2) * (twoJ2 + twom2 + 2), Denom);
				if (debugCG == 2) {
					cout << "cgp:" << endl;
					cout << cgp;
					cout << "cgm:" << endl;
					cout << cgm;
				}

				TFracNum cgmixed = cgp * cgm;
				if (cgmixed.Sqrt()) {
					const bool flipsign = cgm > cgp;
					cgp.Abs();
					cgm.Abs();
					CG[i] = cgp + cgm + TFracNum::mTwo * cgmixed;
					if (flipsign) {
						CG[i].FlipSign();
					}
					if (debugCG == 2) {
						cout << "#############################################" << endl
						     << " middle CG Works fine!!! CG[" << i << "]:"     << endl
						     << CG[i];
					}
				} else {
					if (debugCG == 2) {
						cout << "2*J=" << twoJ << ",2*J1=" << twoJ1 << ",2*J2="
						     << twoJ2 << ",2*m=" << twom << ",2*m1=" << twom1
						     << ",2*m2=" << twom2 << endl
						     << cgp << cgm;
					}
					return 0;
				}
			} else {
				CG[i] = CG[i - max2 + 1] * TFracNum(-(twoJ2 - twom2) * (twoJ2 + twom2 + 2),
				                                     (twoJ1 + twom1) * (twoJ1 - twom1 + 2));
				if (debugCG == 2) {
					cout << "#############################################" << endl
					     << " first row CG[" << i << "]:"                   << endl
					     << CG[i];
				}
			}
		}
		//
		// special case maximum m2=J2
		//
		twom2 = twoJ2;
		if (debugCG == 2) {
			cout << " special setting twom1=" << twom1 << ", twom2=" << twom2 << endl;
		}
		twom = twom1 + twom2;
		twoi2 = twom2 + twoJ2;
		twoi = twoi1 * max2 + twoi2;
		i = twoi / 2;
		if (twom1 + twom2 < -twoJ or twom1 + twom2 > twoJ) {
			CG[i] = TFracNum::Zero;
			if (debugCG == 2) {
				cout << "#############################################" << endl
				     << "first line Simple setting: CG[" << i << "]:"   << endl
				     << CG[i];
			}
		} else if (twom1 > mintwom1 and twom1 <= -mintwom1) {
			const long Denom = (twoJ - twom + 2) * (twoJ + twom);
			TFracNum cg1 = CG[i - max2] * TFracNum((twoJ1 + twom1) * (twoJ1 - twom1 + 2), Denom);
			TFracNum cg2 = CG[i - 1]    * TFracNum((twoJ2 + twom2) * (twoJ2 - twom2 + 2), Denom);
			TFracNum cgmixed = cg1 * cg2;
			if (debugCG == 2) {
				cout << "cg1:" << endl;
				cout << cg1;
				cout << "cg2:" << endl;
				cout << cg2;
				cout << "cgmixed:" << endl;
				cout << cgmixed;
			}
			if (cgmixed.Sqrt()) {
				cg1.Abs();
				cg2.Abs();
				CG[i] = cg1 + cg2 + TFracNum::Two * cgmixed;
				if (debugCG == 2) {
					cout << "#############################################" << endl;
					cout << " first line Works fine!!! CG[" << i << "]:"    << endl;
					cout << CG[i];
				}
			} else {
				if (debugCG == 2) {
					cout << "2*J=" << twoJ << ",2*J1=" << twoJ1 << ",2*J2="
					     << twoJ2 << ",2*m=" << twom << ",2*m1=" << twom1
					     << ",2*m2=" << twom2 << endl
					     << cg1 << cg2;
				}
				return 0;
			}
		}
		//
		// special case m2=J-m1 (when CG(m1-1,m2)=0)
		//
		if (twoJ < twoJ1 + twoJ2   and
		    twom1 > mintwom1       and
		    twom1 <= -mintwom1     and
		    twoJ + twom1 <= twoJ2)
		{
			twom2 = -twoJ - twom1;
			twom = -twoJ;
			twoi2 = twom2 + twoJ2;
			twoi = twoi1 * max2 + twoi2;
			i = twoi / 2;
			const long Denom = (twoJ2 - twom2) * (twoJ2 + twom2 + 2);
			TFracNum cg1 = CG[i + 1]        * TFracNum((twoJ - twom) * (twoJ + twom + 2),     Denom);
			TFracNum cg2 = CG[i - max2 + 1] * TFracNum((twoJ1 + twom1) * (twoJ1 - twom1 + 2), Denom);
			if (debugCG == 2) {
				cout << "cg1:" << endl;
				cout << cg1;
				cout << "cg2:" << endl;
				cout << cg2;
			}
			TFracNum cgmixed = cg1 * cg2;
			if (cgmixed.Sqrt()) {
				const bool flipsign = cg2 > cg1;
				cg1.Abs();
				cg2.Abs();
				CG[i] = cg1 + cg2 + TFracNum::mTwo * cgmixed;
				if (flipsign) {
					CG[i].FlipSign();
				}
				if (debugCG == 2) {
					cout << "#############################################" << endl;
					cout << "lower triangle Works fine!!! CG[" << i << "]:" << endl;
					cout << CG[i];
				}
			} else {
				if (debugCG == 2) {
					cout << "2*J=" << twoJ << ",2*J1=" << twoJ1 << ",2*J2="
					     << twoJ2 << ",2*m=" << twom << ",2*m1=" << twom1
					     << ",2*m2=" << twom2 << endl
					     << cg1 << cg2;
				}
				return 0;
			}
		}
	}

	long overall_sign = 0;

	for (long i = 0; i < max1 * max2; i++) {
		const long twom1 = 2 * (i / max2) - twoJ1;
		const long twom2 = 2 * (i % max2) - twoJ2;
		if (twom1 + twom2 == twoJ and !(CG[i] == TFracNum::Zero)) {
			overall_sign = CG[i].GetSign();
		}
	}

	for (long twom = -twoJ; twom <= twoJ; twom += 2) {
		TFracNum Sum = TFracNum::Zero;
		for (long i = 0; i < max1 * max2; i++) {
			// i=i1*max2+i2 => i1=i/max2, i2=i%max2
			// i1=m1+J1, i2=m2+J2 => m1=i/max2-J1, m2=i%max2-J2
			const long twom1 = 2 * (i / max2) - twoJ1;
			const long twom2 = 2 * (i % max2) - twoJ2;
			if (twom1 + twom2 == twom) {
				TFracNum CGAbs = CG[i];
				CGAbs.Abs();
				Sum += CGAbs;
				if (debugCG == 2) {
					cout << "2m=" << twom << ": adding " << CGAbs.FracString()
					     << " -> " << Sum.FracString() << endl;
				}
			}
		}
		TFracNum SumInv = Sum;
		SumInv.Invert();
		for (long i = 0; i < max1 * max2; i++) {
			const long twom1 = 2 * (i / max2) - twoJ1;
			const long twom2 = 2 * (i % max2) - twoJ2;
			if (twom1 + twom2 == twom) {
				CG[i] = CG[i] * SumInv;
				if (overall_sign == -1) {
					CG[i].FlipSign();
				}
			}
		}
	}
	return CG;
}
