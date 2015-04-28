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

#if(0)
TFracNum* ClebschGordanInteger(long J, long J1, long J2) {

	long nom_ptr[1] = { 1 };
	TFracNum two(1, 0, nom_ptr, 0, 1);
	TFracNum mtwo(1, 0, nom_ptr, 0, -1);

	long max1 = 2 * J1 + 1;
	long max2 = 2 * J2 + 1;
	TFracNum *CG = new TFracNum[max1 * max2];

	CG[max2 - 1] = TFracNum(0, 0, 0, 0, 1);

	for (long m1 = -J1; m1 <= J1; m1++) {
		long i1 = m1 + J1;
		long m2 = 0;
		long m = 0;
		long i2 = 0;
		long i = 0;
		for (m2 = J2 - 1; m2 >= -J2; m2--) {
			if (debugCG == 2)
				cout << " Setting m1=" << m1 << ", m2=" << m2 << endl;
			m = m1 + m2;
			i2 = m2 + J2;
			i = i1 * max2 + i2;
			if (m1 + m2 < -J || m1 + m2 > J) {
				CG[i] = TFracNum(0, 0, 0, 0, 0);
				if (debugCG == 2) {
					cout << "#############################################"
							<< endl;
					cout << "zero CG[" << i << "]:" << endl;
					CG[i].Print();
				}
			} else if (m1 == -J1) {
				CG[i] = CG[i + 1]
						* TFracNum((J - m) * (J + m + 1),
								(J2 - m2) * (J2 + m2 + 1));
				if (debugCG == 2) {
					cout << "#############################################"
							<< endl;
					cout << " first row CG[" << i << "]:" << endl;
					CG[i].Print();
				}
			} else if (m1 + m2 != -J) {
				long Denom = (J1 + m1) * (J1 - m1 + 1);
				TFracNum cgp = CG[i - max2]
						* TFracNum((J + m) * (J - m + 1), Denom);
				TFracNum cgm = CG[i - max2 + 1]
						* TFracNum((J2 - m2) * (J2 + m2 + 1), Denom);
				if (debugCG == 2) {
					cout << "cgp:" << endl;
					cgp.Print();
					cout << "cgm:" << endl;
					cgm.Print();
				}

				TFracNum cgmixed = cgp * cgm;
				if (cgmixed.Sqrt()) {
					bool flipsign = cgm > cgp;
					cgp.Abs();
					cgm.Abs();
					CG[i] = cgp + cgm + mtwo * cgmixed;
					if (flipsign)
						CG[i].FlipSign();
					if (debugCG == 2) {
						cout << "#############################################"
								<< endl;
						cout << " middle CG Works fine!!! CG[" << i << "]:"
								<< endl;
						CG[i].Print();
					}
				} else {
					if (debugCG == 2) {
						cout << "J=" << J << ",J1=" << J1 << ",J2=" << J2
								<< ",m=" << m << ",m1=" << m1 << ",m2=" << m2
								<< endl;
						cgp.Print();
						cgm.Print();
					}
					return 0;
				}
			}
		}
		//
		// special case maximum m2=J2
		//
		m2 = J2;
		m = m1 + m2;
		i2 = m2 + J2;
		i = i1 * max2 + i2;
		if (m1 + m2 > J) {
			CG[i] = TFracNum(0, 0, 0, 0, 0);
			if (debugCG == 2) {
				cout << "#############################################" << endl;
				cout << "first line Simple setting: CG[" << i << "]:" << endl;
				CG[i].Print();
			}
		} else if (m1 != -J1) {
			long Denom = (J - m + 1) * (J + m);
			TFracNum cg1 = CG[i - max2]
					* TFracNum((J1 + m1) * (J1 - m1 + 1), Denom);
			TFracNum cg2 = CG[i - 1]
					* TFracNum((J2 + m2) * (J2 - m2 + 1), Denom);
			TFracNum cgmixed = cg1 * cg2;
			if (debugCG == 2) {
				cout << "cg1:" << endl;
				cg1.Print();
				cout << "cg2:" << endl;
				cg2.Print();
				cout << "cgmixed:" << endl;
				cgmixed.Print();
			}
			if (cgmixed.Sqrt()) {
				cg1.Abs();
				cg2.Abs();
				CG[i] = cg1 + cg2 + two * cgmixed;
				if (debugCG == 2) {
					cout << "#############################################"
							<< endl;
					cout << " first line Works fine!!! CG[" << i << "]:"
							<< endl;
					CG[i].Print();
				}
			} else {
				if (debugCG == 2) {
					cout << "J=" << J << ",J1=" << J1 << ",J2=" << J2 << ",m="
							<< m << ",m1=" << m1 << ",m2=" << m2 << endl;
					cg1.Print();
					cg2.Print();
				}
				return 0;
			}
		}
		//
		// special case m2=J-m1 (when CG(m1-1,m2)=0)
		//
		if (J < J1 + J2 && m1 != -J1 && J + m1 <= J2) {
			m2 = -J - m1;
			m = -J;
			i2 = m2 + J2;
			i = i1 * max2 + i2;
			long Denom = (J2 - m2) * (J2 + m2 + 1);
			TFracNum cg1 = CG[i + 1] * TFracNum((J - m) * (J + m + 1), Denom);
			TFracNum cg2 = CG[i - max2 + 1]
					* TFracNum((J1 + m1) * (J1 - m1 + 1), Denom);
			if (debugCG == 2) {
				cout << "cg1:" << endl;
				cg1.Print();
				cout << "cg2:" << endl;
				cg2.Print();
			}
			TFracNum cgmixed = cg1 * cg2;
			if (cgmixed.Sqrt()) {
				bool flipsign = cg2 > cg1;
				cg1.Abs();
				cg2.Abs();
				CG[i] = cg1 + cg2 + mtwo * cgmixed;
				if (flipsign)
					CG[i].FlipSign();
				if (debugCG == 2) {
					cout << "#############################################"
							<< endl;
					cout << "lower triangle Works fine!!! CG[" << i << "]:"
							<< endl;
					CG[i].Print();
				}
			} else {
				if (debugCG == 2) {
					cout << "J=" << J << ",J1=" << J1 << ",J2=" << J2 << ",m="
							<< m << ",m1=" << m1 << ",m2=" << m2 << endl;
					cg1.Print();
					cg2.Print();
				}
				return 0;
			}
		}
	}

	long overall_sign = 0;
	TFracNum ZERO(0, 0, 0, 0, 0);

	for (long i = 0; i < max1 * max2; i++) {
		long m1 = i / max2 - J1;
		long m2 = i % max2 - J2;
		if (m1 + m2 == J && !(CG[i] == ZERO)) {
			overall_sign = CG[i].GetSign();
		}
	}

	for (long m = -J; m <= J; m++) {
		TFracNum Sum(0, 0, 0, 0, 0);
		for (long i = 0; i < max1 * max2; i++) {
			// i=i1*max2+i2 => i1=i/max2, i2=i%max2
			// i1=m1+J1, i2=m2+J2 => m1=i/max2-J1, m2=i%max2-J2
			long m1 = i / max2 - J1;
			long m2 = i % max2 - J2;
			if (m1 + m2 == m) {
				TFracNum CGAbs = CG[i];
				CGAbs.Abs();
				Sum = Sum + CGAbs;
			}
		}
		TFracNum SumInv = Sum;
		SumInv.Invert();
		for (long i = 0; i < max1 * max2; i++) {
			long m1 = i / max2 - J1;
			long m2 = i % max2 - J2;
			if (m1 + m2 == m) {
				CG[i] = CG[i] * SumInv;
				if (overall_sign == -1)
					CG[i].FlipSign();
			}
		}
	}

	return CG;
}
#endif


TFracNum*
ClebschGordanBox::GetCG(long J, long J1, long J2) {
	for (long i = 0; i < NCG; i++) {
		if (J == CGJ[i] && J1 == CGJ1[i] && J2 == CGJ2[i])
			return CG[i];
	}
	if (debugCGBox)
		cout << "Registering Clebsch-Gordans for J=" << J << " J1=" << J1
				<< " J2=" << J2 << endl;
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

	if (twoJ < 0 || twoJ1 < 0 || twoJ2 < 0 || twoJ < twoJ1 - twoJ2
			|| twoJ < twoJ2 - twoJ1 || twoJ > twoJ1 + twoJ2) {
		cerr << "Clebsch-Gordan-Coefficient only for positive J's"
				<< " with |J1-J2|<=J<=J1+J2" << endl;
	}

	long nom_ptr[1] = { 1 };
	TFracNum two(1, 0, nom_ptr, 0, 1);
	TFracNum mtwo(1, 0, nom_ptr, 0, -1);

	long max1 = twoJ1 + 1;
	long max2 = twoJ2 + 1;

	TFracNum *CG = new TFracNum[max1 * max2];

	// first non-zero element
	long mintwom1 = -twoJ1;
	if (twoJ + twoJ2 < twoJ1)
		mintwom1 = -twoJ - twoJ2;
	long minIndex = (mintwom1 + twoJ1) / 2 * max2 + max2 - 1;
	CG[minIndex] = TFracNum(0, 0, 0, 0, 1);
	if (debugCG == 2) {
		cout << "#############################################" << endl;
		cout << " first row CG[" << minIndex << "]:" << endl;
		CG[minIndex].Print();
	}

	for (long twom1 = -twoJ1; twom1 <= twoJ1; twom1 += 2) {
		long twoi1 = twom1 + twoJ1;
		long twom2 = 0;
		long twom = 0;
		long twoi2 = 0;
		long twoi = 0;
		long i = 0;
		for (twom2 = twoJ2 - 2; twom2 >= -twoJ2; twom2 -= 2) {
			if (debugCG == 2)
				cout << " Setting twom1=" << twom1 << ", twom2=" << twom2
						<< endl;
			twom = twom1 + twom2;
			twoi2 = twom2 + twoJ2;
			twoi = twoi1 * max2 + twoi2;
			i = twoi / 2;
			if (twom1 + twom2 < -twoJ || twom1 + twom2 > twoJ) {
				CG[i] = TFracNum(0, 0, 0, 0, 0);
				if (debugCG == 2) {
					cout << "#############################################"
							<< endl;
					cout << "zero CG[" << i << "]:" << endl;
					CG[i].Print();
				}
			} else if (twom1 == mintwom1) {
				CG[i] = CG[i + 1]
						* TFracNum((twoJ - twom) * (twoJ + twom + 2),
								(twoJ2 - twom2) * (twoJ2 + twom2 + 2));
				if (debugCG == 2) {
					cout << "#############################################"
							<< endl;
					cout << " first row CG[" << i << "]:" << endl;
					CG[i].Print();
				}
			} else if (twom1 + twom2 != -twoJ) {
				long Denom = (twoJ1 + twom1) * (twoJ1 - twom1 + 2);
				TFracNum cgp = CG[i - max2]
						* TFracNum((twoJ + twom) * (twoJ - twom + 2), Denom);
				TFracNum cgm = CG[i - max2 + 1]
						* TFracNum((twoJ2 - twom2) * (twoJ2 + twom2 + 2),
								Denom);
				if (debugCG == 2) {
					cout << "cgp:" << endl;
					cgp.Print();
					cout << "cgm:" << endl;
					cgm.Print();
				}

				TFracNum cgmixed = cgp * cgm;
				if (cgmixed.Sqrt()) {
					bool flipsign = cgm > cgp;
					cgp.Abs();
					cgm.Abs();
					CG[i] = cgp + cgm + mtwo * cgmixed;
					if (flipsign)
						CG[i].FlipSign();
					if (debugCG == 2) {
						cout << "#############################################"
								<< endl;
						cout << " middle CG Works fine!!! CG[" << i << "]:"
								<< endl;
						CG[i].Print();
					}
				} else {
					if (debugCG == 2) {
						cout << "2*J=" << twoJ << ",2*J1=" << twoJ1 << ",2*J2="
								<< twoJ2 << ",2*m=" << twom << ",2*m1=" << twom1
								<< ",2*m2=" << twom2 << endl;
						cgp.Print();
						cgm.Print();
					}
					return 0;
				}
			} else {
				CG[i] = CG[i - max2 + 1]
						* TFracNum(-(twoJ2 - twom2) * (twoJ2 + twom2 + 2),
								(twoJ1 + twom1) * (twoJ1 - twom1 + 2));
				if (debugCG == 2) {
					cout << "#############################################"
							<< endl;
					cout << " first row CG[" << i << "]:" << endl;
					CG[i].Print();
				}
			}
		}
		//
		// special case maximum m2=J2
		//
		twom2 = twoJ2;
		if (debugCG == 2)
			cout << " special setting twom1=" << twom1 << ", twom2=" << twom2
					<< endl;
		twom = twom1 + twom2;
		twoi2 = twom2 + twoJ2;
		twoi = twoi1 * max2 + twoi2;
		i = twoi / 2;
		if (twom1 + twom2 < -twoJ || twom1 + twom2 > twoJ) {
			CG[i] = TFracNum(0, 0, 0, 0, 0);
			if (debugCG == 2) {
				cout << "#############################################" << endl;
				cout << "first line Simple setting: CG[" << i << "]:" << endl;
				CG[i].Print();
			}
		} else if (twom1 > mintwom1 && twom1 <= -mintwom1) {
			long Denom = (twoJ - twom + 2) * (twoJ + twom);
			TFracNum cg1 = CG[i - max2]
					* TFracNum((twoJ1 + twom1) * (twoJ1 - twom1 + 2), Denom);
			TFracNum cg2 = CG[i - 1]
					* TFracNum((twoJ2 + twom2) * (twoJ2 - twom2 + 2), Denom);
			TFracNum cgmixed = cg1 * cg2;
			if (debugCG == 2) {
				cout << "cg1:" << endl;
				cg1.Print();
				cout << "cg2:" << endl;
				cg2.Print();
				cout << "cgmixed:" << endl;
				cgmixed.Print();
			}
			if (cgmixed.Sqrt()) {
				cg1.Abs();
				cg2.Abs();
				CG[i] = cg1 + cg2 + two * cgmixed;
				if (debugCG == 2) {
					cout << "#############################################"
							<< endl;
					cout << " first line Works fine!!! CG[" << i << "]:"
							<< endl;
					CG[i].Print();
				}
			} else {
				if (debugCG == 2) {
					cout << "2*J=" << twoJ << ",2*J1=" << twoJ1 << ",2*J2="
							<< twoJ2 << ",2*m=" << twom << ",2*m1=" << twom1
							<< ",2*m2=" << twom2 << endl;
					cg1.Print();
					cg2.Print();
				}
				return 0;
			}
		}
		//
		// special case m2=J-m1 (when CG(m1-1,m2)=0)
		//
		if (twoJ < twoJ1 + twoJ2 && twom1 > mintwom1 && twom1 <= -mintwom1
				&& twoJ + twom1 <= twoJ2) {
			twom2 = -twoJ - twom1;
			twom = -twoJ;
			twoi2 = twom2 + twoJ2;
			twoi = twoi1 * max2 + twoi2;
			i = twoi / 2;
			long Denom = (twoJ2 - twom2) * (twoJ2 + twom2 + 2);
			TFracNum cg1 = CG[i + 1]
					* TFracNum((twoJ - twom) * (twoJ + twom + 2), Denom);
			TFracNum cg2 = CG[i - max2 + 1]
					* TFracNum((twoJ1 + twom1) * (twoJ1 - twom1 + 2), Denom);
			if (debugCG == 2) {
				cout << "cg1:" << endl;
				cg1.Print();
				cout << "cg2:" << endl;
				cg2.Print();
			}
			TFracNum cgmixed = cg1 * cg2;
			if (cgmixed.Sqrt()) {
				bool flipsign = cg2 > cg1;
				cg1.Abs();
				cg2.Abs();
				CG[i] = cg1 + cg2 + mtwo * cgmixed;
				if (flipsign)
					CG[i].FlipSign();
				if (debugCG == 2) {
					cout << "#############################################"
							<< endl;
					cout << "lower triangle Works fine!!! CG[" << i << "]:"
							<< endl;
					CG[i].Print();
				}
			} else {
				if (debugCG == 2) {
					cout << "2*J=" << twoJ << ",2*J1=" << twoJ1 << ",2*J2="
							<< twoJ2 << ",2*m=" << twom << ",2*m1=" << twom1
							<< ",2*m2=" << twom2 << endl;
					cg1.Print();
					cg2.Print();
				}
				return 0;
			}
		}
	}

	long overall_sign = 0;
	TFracNum ZERO(0, 0, 0, 0, 0);

	for (long i = 0; i < max1 * max2; i++) {
		long twom1 = 2 * (i / max2) - twoJ1;
		long twom2 = 2 * (i % max2) - twoJ2;
		if (twom1 + twom2 == twoJ && !(CG[i] == ZERO)) {
			overall_sign = CG[i].GetSign();
		}
	}

	for (long twom = -twoJ; twom <= twoJ; twom += 2) {
		TFracNum Sum(0, 0, 0, 0, 0);
		for (long i = 0; i < max1 * max2; i++) {
			// i=i1*max2+i2 => i1=i/max2, i2=i%max2
			// i1=m1+J1, i2=m2+J2 => m1=i/max2-J1, m2=i%max2-J2
			long twom1 = 2 * (i / max2) - twoJ1;
			long twom2 = 2 * (i % max2) - twoJ2;
			if (twom1 + twom2 == twom) {
				TFracNum CGAbs = CG[i];
				CGAbs.Abs();
				Sum = Sum + CGAbs;
				if (debugCG == 2)
					cout << "2m=" << twom << ": adding " << CGAbs.FracString()
							<< " -> " << Sum.FracString() << endl;
			}
		}
		TFracNum SumInv = Sum;
		SumInv.Invert();
		for (long i = 0; i < max1 * max2; i++) {
			long twom1 = 2 * (i / max2) - twoJ1;
			long twom2 = 2 * (i % max2) - twoJ2;
			if (twom1 + twom2 == twom) {
				CG[i] = CG[i] * SumInv;
				if (overall_sign == -1)
					CG[i].FlipSign();
				CG[i].SetINTs();
			}
		}
	}
	return CG;
}
