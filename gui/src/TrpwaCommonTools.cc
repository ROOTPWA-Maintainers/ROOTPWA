/*
 * TrpwaCommonTools.cc
 *
 *  Created on: Dec 8, 2010
 *      Author: Promme
 */

#include "TrpwaCommonTools.h"

// check if character is a sign and return the result
bool TrpwaCommonTools::IsSign(char character, int& result){
	result = 0;
	if (character == '-' || character == '+'){
		if (character == '-') result = -1; else result = +1;
		return true;
	}
	else return false;
}

// check if a number is coded to the text and return it
// the number must start at the given position
// pos is set to the next character after the number
bool TrpwaCommonTools::IsNumber(string text, int& result, unsigned int& pos){
	result = 0;
	if (pos < 0) return false;
	unsigned int posend;
	for (posend = pos; posend < text.size(); posend++){
		if (!isdigit(text[posend])) break;
	}
	if (posend == pos) return false;
	result = atoi(text.substr(pos,posend-pos).c_str());
	pos = posend;
	return true;
}

string TrpwaCommonTools::RemovePathExtension(string filename, string extension){
	// remove the extension if available
	if (extension != "" && filename.find(extension.c_str())!=string::npos){
		int pos_ext = filename.find(extension.c_str());
		filename.erase(pos_ext, filename.size()-pos_ext);
	}
	// remove the path if given
	int slashpos = filename.rfind('/');
	if (slashpos != (int) string::npos){
		filename.erase(0, slashpos+1);
	}
	return filename;
}

void TrpwaCommonTools::GetJPCMreflISO1lsISO2(string wavename, int& J, int& P, int& C, int& M, int& refl, string& iso1, string& iso2, int& l, int& s){
	// remove path and/or extension that might be .amp or .key if existent
	string key = RemovePathExtension(wavename, ".amp");
	key = RemovePathExtension(key, ".key");
	J = 0;
	P = 0;
	C = 0;
	M = 0;
	refl = 0;
	iso1 = "";
	iso2 = "";
	l = 0;
	s = 0;

	// IGJPCMeIso1_ls_iso2

	//cout << " TrpwaSessionManager::GetJPCMreflISO1lsISO2() not implemented yet" << endl;
	unsigned int charposlow(0);
	unsigned int charposhigh(0);

	unsigned int stringsize = key.size();
	// expecting the isospin
	int I;
	if (!IsNumber(key,I,charposlow)){
		cout << " error decoding wavename " << key << ", parity expected at position " << charposlow << endl;
		return;
	}
//cout << "G" << endl;
	// a sign for G parity
	int G;
	if (!IsSign(key[charposlow],G)){
		cout << " error decoding wavename " << key << ", G parity not found " << endl;
		return;
	}
	// do nothing with G parity information
//cout << "J" << endl;
	charposlow++;
	// expecting an integer for J
	if (!IsNumber(key,J,charposlow)){
		cout << " error decoding wavename " << key << ", J not found " << endl;
		return;
	}
//cout << "P" << endl;
	if (!IsSign(key[charposlow],P)){
		cout << " error decoding wavename " << key << ", parity not found " << endl;
		return;
	}
//cout << "C" << endl;
    charposlow++;
	// a sign for charge coniugation
	if (!IsSign(key[charposlow],C)){
		cout << " error decoding wavename " << key << ", C parity found " << endl;
		return;
	}
//cout << "M" << endl;
	charposlow++;
	// a digit for the M quantum number expected
	if (!IsNumber(key,M,charposlow)){
		cout << " error decoding wavename " << key << ", M not found " << endl;
		return;
	}
//cout << "e" << endl;
	// a sign expected for reflectivity
	if (!IsSign(key[charposlow],refl)){
		cout << " error decoding wavename " << key << ", reflectivity not found " << endl;
		return;
	}
//cout << "iso1" << endl;
	charposlow++;
	charposhigh = charposlow;
	// reading isobar 1
	while(charposhigh < stringsize && key[charposhigh] != '_'){
		charposhigh++;
	}
	iso1 = key.substr(charposlow, charposhigh-charposlow);
	// skip the underscore
//cout << iso1 << "l" << endl;
	charposlow=charposhigh;
	charposlow++;
	// decode the l
	if (!isdigit(key[charposlow])){
		cout << " error decoding wavename " << key << ", l expected at position " << charposlow << endl;
		return;
	}
	char _l[] = {key[charposlow],'\0'};
	l = atoi(_l);
//cout << "s" << endl;
	charposlow++;
	if (!isdigit(key[charposlow])){
		cout << " error decoding wavename " << key << ", s expected at position " << charposlow << endl;
		return;
	}
	char _s[] = {key[charposlow],'\0'};
	s = atoi(_s);
//cout << "iso2" << endl;
	charposlow++;
	if (key[charposlow] != '_'){
		cout << " error decoding wavename " << key << ", _ expected at position " << charposlow << endl;
		return;
	}
	charposlow++;
	// read the second isobar
	iso2 = key.substr(charposlow, stringsize-charposlow);
}
