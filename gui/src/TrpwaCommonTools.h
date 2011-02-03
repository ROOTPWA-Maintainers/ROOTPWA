/*
 * TrpwaCommonTools.h
 *
 *  Created on: Dec 8, 2010
 *      Author: Promme
 */

#include <string>
#include <iostream>
#include <cstdlib>
#include <vector>

using namespace std;

#ifndef TRPWACOMMONTOOLS_H_
#define TRPWACOMMONTOOLS_H_

namespace TrpwaCommonTools{
	// check if character is a sign and return the result
	bool IsSign(char character, int& result);
	// check if a number is coded to the text and return it
	// the number must start at the given position
	// pos is set to the next character after the number
	bool IsNumber(string text, int& result, unsigned int& pos);
	// get the corresponding variables to the coded wave name
	void GetJPCMreflISO1lsISO2(string wavename, int& J, int& P, int& C, int& M, int& refl, string& iso1, string& iso2, int& l, int& s);
	// remove a given extension and the path (if it exists)
	string RemovePathExtension(string filename, string extension = "");
	// check whether a file exists
	bool FileExists(string filename);
	// check whether a directory exists
	bool DirExists(string dirname);
	// get all files in a given path
	// returns the number of entries
	// files is filled with the filenames in this directory
	// filterext is the extension to be specified for filtering
	// if rmext the files will be delivered without the extension
	int GetDir (string path, vector<string> &files, string filterext = "", bool rmext = false);
	// simple statsus bar
	// percent is 0. to 1.
	void DrawProgressBar(int len, double percent);
}

#endif
