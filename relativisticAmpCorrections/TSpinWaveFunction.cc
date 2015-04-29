#include <iostream>
#include <string>
#include "TSpinWaveFunction.h"
using namespace std;


Int_t debugSpinWave=0;


TSpinWaveFunction::TSpinWaveFunction(Int_t J_, char type_) {

	J=J_;
	type=type_;

	max_pzm=1;
	pot3 = new Int_t[J];
	for (Int_t i=0; i<J; i++) {
		pot3[i]=max_pzm;
		max_pzm*=3;
	}
	mi = new Int_t[max_pzm*J];
	M  = new Int_t[max_pzm];

	coeff = new TFracNum[max_pzm];
	for (Int_t pzm=0; pzm<max_pzm; pzm++) {
		coeff[pzm] = TFracNum::Zero;
	}

	for (Int_t pzm=0; pzm<max_pzm; pzm++) {
		Int_t rest=pzm;
		Int_t i=J-1;
		M[pzm]=0;
		while (i>=0) {
			Int_t ipot3=rest/pot3[i];
			mi[J*pzm+i]=ipot3-1;
			M[pzm]+=mi[J*pzm+i];
			//cout << "Setting mi["<<J<<"*"<<pzm<<"+"<<i<<"]="<<mi[J*pzm+i]<<endl;
			rest-=ipot3*pot3[i];
			i--;
		}
	}

	for (Int_t MM=J; MM>=-J; MM--) {
		if (debugSpinWave==2) cout << "Phi("<<J<<","<<MM<<")="<<endl;
		for (Int_t pzm=max_pzm-1; pzm>=0; pzm--) {
			if ( ((type=='s' || type=='c') && M[pzm]==MM) ||
					(type=='l' && M[pzm]==MM && MM==0) ) {
				Int_t m0=0;
				for (Int_t i=0; i<J; i++) {
					if (debugSpinWave==2) {
						if (mi[J*pzm+i]== 1) cout <<"+";
						if (mi[J*pzm+i]== 0) cout <<"0";
						if (mi[J*pzm+i]==-1) cout <<"-";
					}
					if (mi[J*pzm+i]== 0) m0++;
				}
				if (type=='s' || type=='c')
					coeff[pzm] = TFracNum::am0_to_J(J,MM,m0);
				if (type=='l')
					coeff[pzm] = TFracNum::cm0_sub_ell_2(J,m0);
				if (debugSpinWave==2)
					cout << " * sqrt("
						<< coeff[pzm].FracString()
						<<") (eq. 4.1) ["<<pzm<<"]" <<endl;
			}
		}
	}
}

TTensorSum*
TSpinWaveFunction::GetTensorSum(char name, Int_t delta) {

	TTensorSum *ts = new TTensorSum();

	for (Int_t pzm=max_pzm-1; pzm>=0; pzm--) {
		if ( M[pzm]==delta ) {
			Int_t *pzm_field = new Int_t[J];
			for (Int_t i=0; i<J; i++) pzm_field[i] = mi[J*pzm+i];

			if ( !(coeff[pzm]==TFracNum::Zero) )
				ts->AddTerm(new TTensorTerm(name, J, pzm_field, &(coeff[pzm])));
		}
	}
	return ts;
}

TTensorSum*
TSpinWaveFunction::GetSpinCoupledTensorSum(TSpinWaveFunction* E,
		Int_t delta, Int_t S) {

	if ( !(type=='s' && E->type=='s') ) {
		cerr<<"GetSpinCoupledTensorSum only for spin-type wave functions!"<<endl;
		return 0;
	}

	Int_t twoS1=2*J;
	Int_t twoS2=2*E->J;
	Int_t twoS=2*S;

	Int_t twoSmin=twoS1-twoS2;
	if (twoSmin<0) twoSmin=-twoSmin;

	if (twoS<twoSmin || twoS>twoS1+twoS2){
		cerr<<"GetSpinCoupledTensorSum: no coupling "
			<< J << " + " << E->J << " -> "<< S << endl;
		return 0;
	}

	//TFracNum *S1S2 = ClebschGordan(twoS, twoS1, twoS2);
	TFracNum* S1S2 = ClebschGordanBox::instance()->GetCG(S, J, E->J);
	if (debugSpinWave==1){
		cout << "Clebsch-Gordans calculated for "
			<<twoS<<","<<twoS1<<","<<twoS2<<": "<< S1S2 << endl;
	}

	TTensorSum *ts = new TTensorSum();

	for (Int_t MM1=J; MM1>=-J; MM1--) {
		for (Int_t pzm1=max_pzm-1; pzm1>=0; pzm1--) {
			if ( M[pzm1]==MM1 && !(coeff[pzm1]==TFracNum::Zero) ) {

				for (Int_t MM2=E->J; MM2>=-J; MM2--) {
					for (Int_t pzm2=E->max_pzm-1; pzm2>=0; pzm2--) {
						if ( E->M[pzm2]==MM2 && !(E->coeff[pzm2]==TFracNum::Zero) ) {

							if (MM1+MM2==delta) {

								Int_t *pzm1_field = new Int_t[J];
								for (Int_t i=0; i<J; i++)
									pzm1_field[i] = mi[J*pzm1+i];

								Int_t *pzm2_field = new Int_t[E->J];
								for (Int_t i=0; i<E->J; i++)
									pzm2_field[i] = E->mi[E->J*pzm2+i];

								TTensorTerm* newterm =
									new TTensorTerm('o', J, pzm1_field, &(coeff[pzm1]));

								TFracNum coeff2_CG = S1S2[ClebschGordanBox::CGIndex(J, MM1, E->J, MM2)];
								coeff2_CG = coeff2_CG * E->coeff[pzm2];

								newterm->Multiply('e', E->J, pzm2_field, &coeff2_CG);

								ts->AddTerm(newterm);
							}
						}
					}
				}
			}
		}
	}
	return ts;
}

Int_t
TSpinWaveFunction::CheckCGFormula() {

	TFracNum **J1J=new TFracNum*[J-1];
	for (Int_t i=0; i<J-1; i++){
		// Int_t twoJ=2*(i+2);
		// J1J[i] = ClebschGordan(twoJ, 2, twoJ-2);
		J1J[i] = ClebschGordanBox::instance()->GetCG(i+2, 1, i+1);
	}

	TFracNum *coeffCG = new TFracNum[max_pzm];

	for (Int_t pzm=0; pzm<max_pzm; pzm++) {
		coeffCG[pzm]=TFracNum::One;
		Int_t m=mi[J*pzm];
		if (debugSpinWave==2) cout << m << " x ";
		for (Int_t jj=1; jj<J; jj++) {
			Int_t mj=mi[J*pzm+jj];
			// m1=i/max2-J1, m2=i%max2-J2 ==> i=(m1+J1)*max2
			if (debugSpinWave==2) {
				cout << mj << ": * CG["<<jj-1<<"]["
					<<                            (mj+1)*(2*jj+1) + jj+m
					<<"]:"<<           J1J[jj-1][ (mj+1)*(2*jj+1) + jj+m ].Dval()
					<<endl;
			}
			coeffCG[pzm]=coeffCG[pzm] * J1J[jj-1][ (mj+1)*(2*jj+1) + jj+m ];
			m+=mj;
			if (debugSpinWave==2) cout <<  m << " x ";
		}
		if (debugSpinWave==2) {
			cout <<"done."<<endl;
		}
	}

	for (Int_t MM=J; MM>=0; MM--) {
		cout << "Phi("<<J<<","<<MM<<")="<<endl;
		for (Int_t pzm=max_pzm-1; pzm>=0; pzm--) {
			if (M[pzm]==MM) {
				cout << coeffCG[pzm].FracString()<< " * (";
				Int_t m0=0;
				for (Int_t i=0; i<J; i++) {
					if (mi[J*pzm+i]== 1)  cout <<"+";
					if (mi[J*pzm+i]== 0) {cout <<"0"; m0++;}
					if (mi[J*pzm+i]==-1)  cout <<"-";
				}
				cout << ")   [ (4.1): "<< coeff[pzm].FracString()<<" ]" <<endl;
			}
		}
	}
	return 0;
}
