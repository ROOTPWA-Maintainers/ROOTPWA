#include <iostream>
#include <string>
#include "TJwfTensor.h"
using namespace std;

TTensorTerm::TTensorTerm(char name, Int_t RJ, Int_t* pzm_field,
		TFracNum* prefac_){
	Rome=0; ome_pzm=0;
	Reps=0; eps_pzm=0;
	Rchi=0; chi_pzm=0;
	Rphi=0; phi_pzm=0;
	gam_s_pot=0;
	gam_sig_pot=0;
	prefac=*prefac_;
	if (name == 'o') {Rome=RJ; ome_pzm=pzm_field;}
	if (name == 'e') {Reps=RJ; eps_pzm=pzm_field;}
	if (name == 'c') {Rchi=RJ; chi_pzm=pzm_field;}
	if (name == 'p') {Rphi=RJ; phi_pzm=pzm_field;}
}

TTensorTerm::TTensorTerm(TTensorTerm *S, TTensorTerm *L,
		Int_t contractions,
		Int_t o_share, Int_t e_share,
		char con_type) {

	if ( L->Rome || L->Reps ) return;
	Rome = S->Rome;
	Reps = S->Reps;
	Rchi = S->Rchi + L->Rchi;
	Rphi = S->Rphi + L->Rphi;
	gam_s_pot   = S->gam_s_pot  +L->gam_s_pot;
	gam_sig_pot = S->gam_sig_pot+L->gam_sig_pot;
	prefac=S->prefac*L->prefac;

	Int_t ocount=o_share;
	Int_t ecount=e_share;
	for (Int_t con=0; con<contractions; con++) {
		Int_t cp = 0;
		if (con_type=='c') { cp=L->chi_pzm[Rchi-1]; Rchi--; }
		else               { cp=L->phi_pzm[Rphi-1]; Rphi--; }
		Int_t oe = 0;
		if (ocount) {              // o is to be contracted
			oe = S->ome_pzm[Rome-1];
			Rome--;
		}
		else {                     // e is to be contracted
			oe = S->eps_pzm[Reps-1];
			ecount--;
			Reps--;
		}
		if ( con_type == 'c' && ( (oe==1&&cp==-1) || (oe==-1&&cp==1) ) )
			prefac.FlipSign();
		else if
			(  con_type == 'p' && ( (oe==1&&cp==1) || (oe==-1&&cp==-1) ) ) ;
		else if (oe==0&&cp==0 ) {
			if (ocount) gam_s_pot++;
			else        gam_sig_pot++;
		}
		else {
			prefac=TFracNum(0,0,0,0,0);
			return;
		}
		if (ocount) ocount--;
	}

	ome_pzm = new Int_t[Rome];
	for (Int_t i=0; i<Rome; i++) ome_pzm[i] = S->ome_pzm[i];
	eps_pzm = new Int_t[Reps];
	for (Int_t i=0; i<Reps; i++) eps_pzm[i] = S->eps_pzm[i];
	chi_pzm = new Int_t[Rchi];
	for (Int_t i=0; i<Rchi; i++) {
		if (i < L->Rchi) chi_pzm[i] = L->chi_pzm[i];
		else             chi_pzm[i] = S->chi_pzm[i-(L->Rchi)];
	}
	phi_pzm = new Int_t[Rphi];
	for (Int_t i=0; i<Rphi; i++){
		if (i < L->Rphi) phi_pzm[i] = L->phi_pzm[i];
		else             phi_pzm[i] = S->phi_pzm[i-(L->Rphi)];
	}
}

Int_t
TTensorTerm::LJContraction(Int_t ncon, Int_t even) {

	for (Int_t con=0; con<ncon; con++) {
		Int_t c=chi_pzm[Rchi-1];
		Int_t p=phi_pzm[Rphi-1];
		if ( c != p ) {
			prefac=TFracNum(0,0,0,0,0);
			return 0;
		}
		Rchi--;
		Rphi--;
	}

	Int_t trouble=0;

	if (even==0) {
		if (Rome+Reps+Rchi+Rphi!=3) {
			cerr << "TTensorTerm::LJContraction:"
				<<" Contraction ended with wrong number of indices!!" <<endl;
			return 0;
		}
		else { // eq. (5.21) - (5.23)
			Int_t os=-100; if (Rome==1) os=ome_pzm[0];
			Int_t es=-100; if (Reps==1) es=eps_pzm[0];
			Int_t cs=-100; if (Rchi==1) cs=chi_pzm[0];
			Int_t ps=-100; if (Rphi==1) ps=phi_pzm[0];
			if ( Rome==0 ) {
				Reps--; Rchi--; Rphi--;
				if      (es== 1 && cs==-1 && ps== 0 )                  ;
				else if (es==-1 && cs== 1 && ps== 0 ) prefac.FlipSign();
				else if (es== 1 && cs== 0 && ps== 1 )                  ;
				else if (es==-1 && cs== 0 && ps==-1 ) prefac.FlipSign();
				else if (es== 0 && cs== 1 && ps== 1 ){prefac.FlipSign();gam_sig_pot++;}
				else if (es== 0 && cs==-1 && ps==-1 ){                  gam_sig_pot++;}
				else { prefac=TFracNum(0,0,0,0,0); return 1;}
			}
			else if ( Reps==0 ) {
				Rome--; Rchi--; Rphi--;
				if      (os== 1 && cs==-1 && ps== 0 )                  ;
				else if (os==-1 && cs== 1 && ps== 0 ) prefac.FlipSign();
				else if (os== 1 && cs== 0 && ps== 1 )                  ;
				else if (os==-1 && cs== 0 && ps==-1 ) prefac.FlipSign();
				else if (os== 0 && cs== 1 && ps== 1 ){prefac.FlipSign();gam_s_pot++;}
				else if (os== 0 && cs==-1 && ps==-1 ){                  gam_s_pot++;}
				else { prefac=TFracNum(0,0,0,0,0); return 1;}
			}
			else if ( Rchi==0 ) {
				Rome--; Reps--; Rphi--;
				if      (os== 1 && es==-1 && ps== 0 )                  ;
				else if (os==-1 && es== 1 && ps== 0 ) prefac.FlipSign();
				else if (os== 1 && es== 0 && ps== 1 ){                  gam_sig_pot++;}
				else if (os==-1 && es== 0 && ps==-1 ){prefac.FlipSign();gam_sig_pot++;}
				else if (os== 0 && es== 1 && ps== 1 ){prefac.FlipSign();gam_s_pot++;}
				else if (os== 0 && es==-1 && ps==-1 ){                  gam_s_pot++;}
				else { prefac=TFracNum(0,0,0,0,0); return 1;}
			}
			else if ( Rphi==0 ) {
				Rome--; Reps--; Rchi--;
				if      (os== 1 && es==-1 && cs== 0 )                  ;
				else if (os==-1 && es== 1 && cs== 0 ) prefac.FlipSign();
				else if (os== 1 && es== 0 && cs==-1 ){prefac.FlipSign();gam_sig_pot++;}
				else if (os==-1 && es== 0 && cs== 1 ){                  gam_sig_pot++;}
				else if (os== 0 && es== 1 && cs==-1 ){                  gam_s_pot++;}
				else if (os== 0 && es==-1 && cs== 1 ){prefac.FlipSign();gam_s_pot++;}
				else { prefac=TFracNum(0,0,0,0,0); return 1;}
			}
			else { trouble=5;}
		}
	}

	else {
		if (Rome!=0 || Reps!=0 || Rchi!=0 || Rphi!=0 ) trouble=1;
	}
	if (trouble){
		cerr << "TTensorTerm::LJContraction: "
			<< "Invalid espilon-contraction occurred "<< trouble << endl;
		return 0;
	}
	return 1;
}

Int_t
TTensorTerm::Multiply(char name, Int_t RJ, Int_t* pzm_field,
		TFracNum* prefac_){

	prefac = prefac * *prefac_;

	Int_t Merr=0;
	if (name == 'o') {if (Rome) Merr=1; else {Rome=RJ; ome_pzm=pzm_field;}}
	if (name == 'e') {if (Reps) Merr=1; else {Reps=RJ; eps_pzm=pzm_field;}}
	if (name == 'c') {if (Rchi) Merr=1; else {Rchi=RJ; chi_pzm=pzm_field;}}
	if (name == 'p') {if (Rphi) Merr=1; else {Rphi=RJ; phi_pzm=pzm_field;}}
	if (Merr) {
		cerr << "TTensorTerm::Multiply: Each type can be multiplied only once!"
			<<endl;
		return 0;
	}
	return 1;
}

Int_t
TTensorTerm::SpinInnerContraction(Int_t cPsiInt) {
	Int_t res=0;
	for (Int_t ic=0; ic<cPsiInt; ic++) {
		res=0;
		Int_t o = ome_pzm[Rome-1];
		Int_t e = eps_pzm[Reps-1];
		if ( (o==1&&e==-1) || (o==-1&&e==1) ) {
			prefac.FlipSign();
			res=1;
		}
		if   (o==0&&e==0) {
			gam_s_pot   += 1;
			gam_sig_pot += 1;
			res=1;
		}
		Rome--;
		Reps--;
	}
	return res;
}

Int_t
TTensorTerm::SameStructure(TTensorTerm *other) {
	if (Rome==other->Rome &&
			Reps==other->Reps &&
			Rchi==other->Rchi &&
			Rphi==other->Rphi &&
			gam_s_pot==other->gam_s_pot &&
			gam_sig_pot==other->gam_sig_pot) return 1;
	return 0;
}

Int_t
TTensorTerm::AddTwoTerms(TTensorTerm *other) {
	if (!SameStructure(other)) {
		cerr << "NO NO NO these terms cannot be added!" << endl;
		return 0;
	}
	else {
		TFracNum *sum = prefac.SumSignedRoots(&(other->prefac));
		if (sum) {
			prefac= *sum;
			return 1;
		}
	}
	return 0;
}

Int_t
TTensorTerm::Print(char flag) {
	if (flag=='s') cout << prefac.FracStringSqrt()<<" ";
	else           cout << prefac.FracString()<<" ";
	if (Rome) {
		cout << "o(";
		for (Int_t i=0; i<Rome; i++) {
			if (ome_pzm[i]== 1) cout <<"+";
			if (ome_pzm[i]== 0) cout <<"0";
			if (ome_pzm[i]==-1) cout <<"-";
		}
		cout << ")";
	}
	if (Reps) {
		cout << "e(";
		for (Int_t i=0; i<Reps; i++) {
			if (eps_pzm[i]== 1) cout <<"+";
			if (eps_pzm[i]== 0) cout <<"0";
			if (eps_pzm[i]==-1) cout <<"-";
		}
		cout << ")";
	}
	if (Rchi) {
		cout << "c(";
		for (Int_t i=0; i<Rchi; i++) {
			if (chi_pzm[i]== 1) cout <<"+";
			if (chi_pzm[i]== 0) cout <<"0";
			if (chi_pzm[i]==-1) cout <<"-";
		}
		cout << ")";
	}
	if (Rphi) {
		cout << "p(";
		for (Int_t i=0; i<Rphi; i++) {
			if (phi_pzm[i]== 1) cout <<"+";
			if (phi_pzm[i]== 0) cout <<"0";
			if (phi_pzm[i]==-1) cout <<"-";
		}
		cout << ")";
	}
	if (gam_s_pot) {
		if (gam_s_pot==1) cout << " gs";
		else              cout << " gs^" << gam_s_pot;
	}
	if (gam_sig_pot) {
		if (gam_sig_pot==1) cout << " gsig";
		else                cout << " gsig^" << gam_sig_pot;
	}
	return 0;
}

//Int_t
//TTensorSum::Print(char flag='n'){ // CINT limitation for overloading
Int_t
TTensorSum::Print(char flag){
	for (Int_t i=0; i<Nterms; i++) {
		terms[i].Print(flag);
		if (i<Nterms-1) cout << " ";
	}
	cout << endl;
	return 0;
};

Int_t
TTensorSum::AddTerm (TTensorTerm* addt) {
	TTensorTerm* nt=new TTensorTerm[Nterms+1];
	for (Int_t i=0; i<Nterms; i++)
		nt[i]=terms[i];
	nt[Nterms]= *addt;
	delete[] terms;
	terms=nt;
	Nterms++;
	return 0;
}

Int_t
TTensorSum::SpinInnerContraction(Int_t cPsiInt) {
	Int_t non_zero_terms=0;
	Int_t val[Nterms];
	for (Int_t i=0; i<Nterms; i++) {
		val[i] = terms[i].SpinInnerContraction(cPsiInt);
		if (val[i]!=0) non_zero_terms++;
	}
	if (non_zero_terms<Nterms) {
		TTensorTerm* nt=new TTensorTerm[non_zero_terms];
		Int_t j=0;
		for (Int_t i=0; i<Nterms; i++)
			if (val[i]) {
				nt[j]=terms[i];
				j++;
			}
		delete[] terms;
		terms=nt;
		Nterms=non_zero_terms;
	}
	return Nterms;
}

TTensorSum*
TTensorSum::LSContraction(TTensorSum *L, Int_t contr,
		Int_t co, Int_t ce, char con_type) {

	TTensorSum *tls = new TTensorSum();

	for (Int_t i=0; i<Nterms; i++) {
		for (Int_t j=0; j<L->Nterms; j++) {
			TTensorTerm *nt = new TTensorTerm(&(terms[i]), &(L->terms[j]),
					contr, co, ce, con_type);
			if ( nt->IsNonZero() ) tls->AddTerm(nt);
		}
	}

	return tls;
}

TTensorSum*
TTensorSum::LJContraction(Int_t cChiPhi, Int_t even) {

	for (Int_t i=0; i<Nterms; i++) terms[i].LJContraction(cChiPhi, even);

	TTensorSum *tls = new TTensorSum();

	for (Int_t i=0; i<Nterms; i++) {
		Int_t found_same=0;
		for (Int_t j=0; j<tls->Nterms; j++) {
			if (terms[i].SameStructure(&(tls->terms[j]))) {
				if ( (tls->terms[j]).AddTwoTerms(&(terms[i])) ) {
					found_same=1;
					break;
				}
			}
		}
		if (!found_same) {
			TTensorTerm *nt = new TTensorTerm(terms[i]);
			tls->AddTerm(nt);
		}
	}

	TTensorSum *tls_nonzero = new TTensorSum();

	for (Int_t i=0; i<tls->Nterms; i++) {
		if ( tls->terms[i].IsNonZero() ) {
			TTensorTerm *nt = new TTensorTerm(tls->terms[i]);
			tls_nonzero->AddTerm(nt);
		}
	}

	return tls_nonzero;
}
