// small macro converting a given wavelist to a latex file
// note: this script does not work for sequential decays
// author: Promme@web.de
// date: 16.12.10
//
#include <string>
#include <iostream>

using namespace std;

// taken from TrpwaCommonTools.cc in the gui/src directory 16.12.10

// check if character is a sign and return the result
bool IsSign(char character, int& result){
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
bool IsNumber(string text, int& result, unsigned int& pos){
	result = 0;
	//if (pos < 0) return false;
	unsigned int posend;
	for (posend = pos; posend < text.size(); posend++){
		if (!isdigit(text[posend])) break;
	}
	if (posend == pos) return false;
	result = atoi(text.substr(pos,posend-pos).c_str());
	pos = posend;
	return true;
}

string RemovePathExtension(string filename, string extension){
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

void GetJPCMreflISO1lsISO2(string wavename, int& J, int& P, int& C, int& M, int& refl, string& iso1, string& iso2, int& l, int& s){
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

#include <fstream>
#include <map>
#include <sstream>

void wavelist_to_latex(string wavelistname, string latexfilename = "wavelist.tex"){
	ifstream wavelistfile(wavelistname.c_str());
	if (!wavelistfile){
		cout << " could not open wavelist " << wavelistname << endl;
		return;
	}
	ofstream latexfile(latexfilename.c_str());
	if (!latexfile){
		cout << " could not open latexfile " << latexfilename << endl;
		return;
	}
	map<string,string> isobars; // map of existing isobar names in the naming sheme and an index for it

	int nlines(0);
	int nisobars(0);
	latexfile << "\\begin{table}" << endl;
	latexfile << "  \\myfloatalign " << endl;
	latexfile << " \\begin{tabularx}{\\textwidth}{XXX} \\toprule " << endl;
	latexfile << " \\tableheadline{$J^{PC}$} &";
	latexfile << " \\tableheadline{$M\\epsilon$} &";
	latexfile << " \\tableheadline{$iso1 \\left[ ^l _s \\right] iso2$} \\\\ \\midrule" << endl;
	while (1){
		char oneline[1024];
		wavelistfile.getline(oneline, 1024);
		if (!wavelistfile.good()) break;
		// skip excluded entries
		if (oneline[0] == '#') continue;
		int J,P,C,M,refl,l,s;
		string iso1, iso2, signP, signC, signrefl;
		string wavename = oneline;
		GetJPCMreflISO1lsISO2(wavename, J, P, C, M, refl, iso1, iso2, l, s);
		// search if this isobars exist already and create/retrieve the index of it		
		map<string,string>::iterator it = isobars.find(iso1);
		if (it != isobars.end()){
		} else {
			// construct a counting sheme out of 26 letters starting from #65 (A-Z)
			stringstream iso;
			int rest = nisobars;
			int counter[] = {65,65,65};
			while(rest){
				counter[0]++;
				if (counter[0] == 91){ // reached Z
					counter[1]++;
					counter[0] = 65; // restart from A
				}
				if (counter[1] == 91){
					counter[2]++;
					counter[1] = 65;
				}
					rest--;
			}
			iso << "\\isobar" << (char) counter[2] << (char) counter[1] << (char) counter[0];
			isobars[iso1]=iso.str();
			nisobars++;
		}
		it = isobars.find(iso2);
		if (it != isobars.end()){
		} else {
			// construct a counting sheme out of 26 letters starting from #65 (A-Z)
			stringstream iso;
			int rest = nisobars;
			int counter[] = {65,65,65};
			while(rest){
				counter[0]++;
				if (counter[0] == 91){ // reached Z
					counter[1]++;
					counter[0] = 65; // restart from A
				}
				if (counter[1] == 91){
					counter[2]++;
					counter[1] = 65;
				}
					rest--;
			}
			iso << "\\isobar" << (char) counter[2] << (char) counter[1] << (char) counter[0];
			isobars[iso2]=iso.str();
			nisobars++;
		}
		// give the isobars a latex command name
		iso1 = isobars[iso1];
		iso2 = isobars[iso2]; 
		signP = P < 0 ? "-" : "+";
		signC = C < 0 ? "-" : "+";
		signrefl = refl < 0 ? "-" : "+";
		latexfile << " $" << J << "^{" << signP << signC << "}$ & $";
		latexfile << M << signrefl << "$ & $";
		latexfile << " " << iso1;
		latexfile << " \\left[ ^" << l << " _" << s << " \\right]";
		latexfile << " " << iso2 << " $ \\\\ " << endl; 
		nlines++;
	}
	latexfile << " \\bottomrule " << endl;
	latexfile << " \\end{tabularx}" << endl;
	latexfile << "   \\caption[final partial wave set]{The final partial wave set used for fitting.}" << endl;
    latexfile << " \\label{tab:pwset}" << endl;
	latexfile << "\\end{table}" << endl;
	latexfile << endl;

	latexfile << " % specification of latex commands to be filled by the user and put to proper place " << endl;
	// create a definition of latex commands to be specified by the user
	for (map<string,string>::iterator it = isobars.begin(); it != isobars.end(); it++){
		latexfile << "\\newcommand{";
		latexfile << it->second <<"}{"<< it->first <<"}" << endl;
	}

	cout << " converted " << nlines << " waves " << endl;
	latexfile.close();
	wavelistfile.close();
}
