/*
 * TrpwaCommonTools.cc
 *
 *  Created on: Dec 8, 2010
 *      Author: Promme
 */

#ifndef TRPWACOMMONTOOLS_CC_
#define TRPWACOMMONTOOLS_CC_

#include "TrpwaCommonTools.h"
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <fstream>
#include <sys/types.h>
#include <sys/stat.h>

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

// check whether a file exists
bool TrpwaCommonTools::FileExists(string filename){
	  ifstream ifile(filename.c_str());
	  return ifile;
}

// check whether a directory exists
bool TrpwaCommonTools::DirExists(string dirname){
	struct stat st;
	if(stat(dirname.c_str(),&st) == 0 && S_ISDIR(st.st_mode))
		return true;
	else
		return false;
}

int TrpwaCommonTools::GetDir (string path,
		vector<string> &files, string filterext, bool rmext){
    DIR *dp;
    struct dirent *dirp;
    files.clear();
    if((dp  = opendir(path.c_str())) == NULL) {
        cout << "Error(" << errno << ") opening " << path << endl;
        return errno;
    }

    while ((dirp = readdir(dp)) != NULL) {
    	string _filename = dirp->d_name;
    	if ((_filename.size() <= filterext.size())||
    		(filterext != "" && _filename.compare(_filename.size()-filterext.size(), filterext.size(), filterext) != 0)) continue;
    	if (rmext) _filename.erase(_filename.size()-filterext.size(), filterext.size());
    	files.push_back(_filename);
    }
    closedir(dp);
    return (signed) files.size();
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

#endif
