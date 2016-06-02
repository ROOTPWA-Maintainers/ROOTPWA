///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
//
//    This file is part of rootpwa
//
//    rootpwa is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    rootpwa is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with rootpwa. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
//
// Description:
//      class that reads/writes wave description from/to keyfiles,
//      constructs decay topologies and amplitudes
//      can be saved into .root files
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <fstream>
#include <map>
#include <cstdio>

#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <boost/assign.hpp>

#include "TClass.h"

#include "amplitudeMetadata.h"
#include "libConfigUtils.hpp"
#include "conversionUtils.hpp"
#include "diffractiveDissVertex.h"
#include "leptoProductionVertex.h"
#include "isobarHelicityAmplitude.h"
#include "isobarCanonicalAmplitude.h"
#include "waveDescription.h"


using namespace std;
using namespace boost;
using namespace boost::assign;
using namespace libconfig;
using namespace rpwa;


namespace {


	typedef boost::shared_ptr<Config> configPtr;


	const string binBorders     = "bins";
	const string degreeName     = "degree";
	const string mMinName       = "mMin";
	const string mMaxName       = "mMax";
	const string realCoeffName  = "realCoefficients";
	const string imagCoeffName  = "imagCoefficients";
	const string expandName     = "expand";
	const string expandTypeName = "type";
	const string expandNameName = "name";


	string
	getExpandType(const Setting* setting)
	{
		string expandType;
		if (not setting->lookupValue(expandTypeName, expandType)) {
			printErr << "no expand type given." << endl;
			return "";
		}
		return expandType;
	}


	string
	getExpandName(const Setting* setting)
	{
		string binningName;
		if (not setting->lookupValue(expandNameName, binningName)) {
			printErr << "no binning name given." << endl;
			return "";
		}
		return binningName;
	}


	vector<double>
	getChebyshevCoefficients(const size_t n)
	{
		vector<double> coefficientsNMinus2(1,1.);
		if (n == 0) {
			return coefficientsNMinus2;
		}
		vector<double> coefficientsNMinus1(1,0.);
		coefficientsNMinus1.push_back(1.);
		if (n == 1) {
			return coefficientsNMinus1;
		}
		vector<double> coefficientsN;
		for (size_t degree = 2; degree < n+1; ++degree) {
			coefficientsN = vector<double>(degree+1,0.); // A polynomial of degree n has n+1 coefficients
			coefficientsN[0] = -coefficientsNMinus2[0];
			for (size_t i = 1; i < degree-1; ++i) {
				coefficientsN[i] = 2*coefficientsNMinus1[i-1] - coefficientsNMinus2[i];
			}
			coefficientsN[degree-1] = 2*coefficientsNMinus1[degree-2];
			coefficientsN[degree  ]  = 2*coefficientsNMinus1[degree-1];
			coefficientsNMinus2 = coefficientsNMinus1;
			coefficientsNMinus1 = coefficientsN;
		}
		return coefficientsN;
	}


	vector<double>
	getBinning(const Setting* setting)
	{
		const Setting* binningSetting = findLibConfigList(*setting, binBorders, true);
		if (not binningSetting) {
			printErr << "no binning given to expand '" << expandName << "' "
			         << "in path '" << setting->getPath() << "'." << endl;
			return vector<double>();
		}
		const int length = binningSetting->getLength();
		vector<double> binningVector(length);
		for (int bin=0; bin<length; ++bin) {
			binningVector[bin] = (*binningSetting)[bin];
		}
		return binningVector;
	}


	Setting*
	findExpand(const Setting &setting)
	{
		Setting* expand  = findLibConfigGroup(setting, expandName, false);
		if (expand) {
			return expand;
		}
		const int length = setting.getLength();
		for (int i=0; i<length; ++i) {
			Setting &actSetting = setting[i];
			const int actLength = actSetting.getLength();
			if (actLength > 0) {
				Setting* recurse = findExpand(actSetting);
				if (recurse) {
					return recurse;
				}
			}
		}
		return NULL;
	}


	vector<configPtr>
	expand(const configPtr& config)
	{
		vector<configPtr> expanded;

		Setting& root     = config->getRoot();
		Setting* toExpand = findExpand(root);
		if (not toExpand) {
			expanded.push_back(config);
			return expanded;
		}
		Setting& parent = toExpand->getParent();

		const string expandType = getExpandType(toExpand);
		if (expandType == "binning") {
			const string         binningName = getExpandName(toExpand);
			const vector<double> binning     = getBinning(toExpand);

			parent.remove(expandName);

			for (size_t bin=0; bin<binning.size()-1; ++bin) {
				parent.add(binningName,Setting::TypeList);
				parent[binningName.c_str()].add(Setting::TypeFloat) = binning[bin];
				parent[binningName.c_str()].add(Setting::TypeFloat) = binning[bin+1];

				configPtr newConfig(new Config);
				copyConfig(*config, *newConfig, waveDescription::debug());

				parent.remove(binningName);

				const vector<configPtr> newExpanded = expand(newConfig);
				expanded.insert(expanded.end(), newExpanded.begin(), newExpanded.end());
			}
		} else if (expandType == "taylor") {
			int degree = 0;
			if (not toExpand->lookupValue(degreeName, degree)) {
				printErr << "no degree given. Aborting..." << endl;
				throw;
			}
			parent.remove(expandName);
			for (int i = 0; i < degree+1; ++i) {
				parent.add(realCoeffName, Setting::TypeList);
				parent.add(imagCoeffName, Setting::TypeList);
				for (int j = 0; j < i; ++j) {
					parent[realCoeffName.c_str()].add(Setting::TypeFloat) = 0.;
					parent[imagCoeffName.c_str()].add(Setting::TypeFloat) = 0.;
				}
				parent[realCoeffName.c_str()].add(Setting::TypeFloat) = 1.;
				parent[imagCoeffName.c_str()].add(Setting::TypeFloat) = 0.;
				configPtr newConfig(new Config);
				copyConfig(*config, *newConfig, waveDescription::debug());

				parent.remove(realCoeffName);
				parent.remove(imagCoeffName);

				const vector<configPtr> newExpanded = expand(newConfig);
				expanded.insert(expanded.end(), newExpanded.begin(), newExpanded.end());
			}
		} else if (expandType == "fourier") {
			int degree = 0;
			if (not toExpand->lookupValue(degreeName, degree)) {
				printErr << "no degree given. Aborting..." << endl;
				throw;
			}
			parent.remove(expandName);
			for (int i = -degree; i < degree+1; ++i) {
				parent.add(degreeName.c_str(), Setting::TypeInt) = i;

				configPtr newConfig(new Config);
				copyConfig(*config, *newConfig, waveDescription::debug());

				parent.remove(degreeName);

				const vector<configPtr> newExpanded = expand(newConfig);
				expanded.insert(expanded.end(), newExpanded.begin(), newExpanded.end());
			}
		} else if (expandType == "chebyshev") {
			int degree = 0;
			if (not toExpand->lookupValue(degreeName, degree)) {
				printErr << "no degree given. Aborting..." << endl;
				throw;
			}
			parent.remove(expandName);
			for (int i = 0; i < degree; ++i) {
				parent.add(realCoeffName, Setting::TypeList);
				parent.add(imagCoeffName, Setting::TypeList);
				vector<double> coefficients = getChebyshevCoefficients(i);
				for (int j = 0; j < i+1; ++j) {
					parent[realCoeffName.c_str()].add(Setting::TypeFloat) = coefficients[j];
					parent[imagCoeffName.c_str()].add(Setting::TypeFloat) = 0.;
				}
				configPtr newConfig(new Config);
				copyConfig(*config, *newConfig, waveDescription::debug());

				parent.remove(realCoeffName);
				parent.remove(imagCoeffName);

				const vector<configPtr> newExpanded = expand(newConfig);
				expanded.insert(expanded.end(), newExpanded.begin(), newExpanded.end());
			}
		} else {
			printErr << "expand type '" << expandType << "' unknown." << endl;
		}

		return expanded;
	}


}


ClassImp(waveDescription);


bool waveDescription::_debug = false;


map<string,string> waveDescription::isobars = map_list_of
	("pi0",         "\\pi^0")
	("pi+",         "\\pi^+")
	("pi-",         "\\pi^-")
	("pi+-",        "\\pi^\\pm")
	("pi-+",        "\\pi^\\mp")
	("eta0",        "\\eta")
	("sigma0",      "\\sigma")
	("rho(770)0",   "\\rho^0(770)")
	("a1(1260)-",   "a_1^-(1260)")
	("a2(1320)-",   "a_2^-(1320)")
	("rho(1450)0",  "\\rho^0(1450)")
	("rho(1700)0",  "\\rho^0(1700)")
	("pi(1300)-",   "\\pi^-(1300)")
	("pi(1800)-",   "\\pi^-(1800)")
	("pi2(1670)-",  "\\pi^-_2(1670)")
	("f0(1370)0",   "f_0^0(1370)")
	("f0(1500)0",   "f_0^0(1500)")
	("f0(1700)0",   "f_0^0(1700)")
	("f1(1285)0",   "f_1^0(1285)")
	("f1(1420)0",   "f_1^0(1420)")
	("b1(1235)0",   "b_1^0(1235)")
	("b1(1800)0",   "b_1^0(1800)")
	("b0(1800)0",   "b_0^0(1800)")
	("b2(1800)0",   "b_2^0(1800)")
	("b1(1500)0",   "b_1^0(1500)")
	("f2(1270)0",   "f_2^0(1270)")
	("f2(1950)0",   "f_2^0(1950)")
	("f2(1565)0",   "f_2^0(1565)")
	("f2(2010)0",   "f_2^0(2010)")
	("eta(1440)0" , "\\eta^0(1420)")
	("eta2(1645)0", "\\eta_2^0(1645)")
	("eta1(1600)0", "\\eta_1^0(1600)")
	("pi1(1600)-",  "\\pi_1^-(1600)")
	("rho3(1690)0", "\\rho_3^0(1690)")
	("rho(1600)0",  "\\rho^0(1600)");


waveDescription::waveDescription()
	: TObject          (),
	  _key             (0),
	  _keyFileParsed   (false),
	  _keyFileLocalCopy("")
{
	//waveDescription::Class()->IgnoreTObjectStreamer();  // don't store TObject's fBits and fUniqueID
}


waveDescription::waveDescription(const waveDescription& waveDesc)
	: TObject          (waveDesc),
	  _key             (0),
	  _keyFileParsed   (false),
	  _keyFileLocalCopy(waveDesc._keyFileLocalCopy)
{
	//waveDescription::Class()->IgnoreTObjectStreamer();  // don't store TObject's fBits and fUniqueID
	if (waveDesc._keyFileParsed)
		parseKeyFileLocalCopy();
}

waveDescription::waveDescription(const amplitudeMetadata* amplitudeMeta)
	: TObject          (),
	  _key             (0),
	  _keyFileParsed   (false),
	  _keyFileLocalCopy(amplitudeMeta->keyfileContent())
{
	//waveDescription::Class()->IgnoreTObjectStreamer();  // don't store TObject's fBits and fUniqueID
	parseKeyFileLocalCopy();
}


waveDescription::~waveDescription()
{
	if (_key)
		delete _key;
}


void
waveDescription::clear()
{
	if (_key)
		delete _key;
	_key              = 0;
	_keyFileParsed    = false;
	_keyFileLocalCopy = "";
}


waveDescription&
waveDescription::operator =(const waveDescription& waveDesc)
{
	if (this != &waveDesc) {
		TObject::operator =(waveDesc);
		_key              = 0;
		_keyFileParsed    = false;
		_keyFileLocalCopy = waveDesc._keyFileLocalCopy;
		if (waveDesc._keyFileParsed)
			parseKeyFileLocalCopy();
	}
	return *this;
}


vector<waveDescriptionPtr>
waveDescription::parseKeyFile(const string& keyFileName)
{
	configPtr config(new Config);
	if (not parseLibConfigFile(keyFileName, *config, _debug)) {
		printWarn << "problems reading key file '" << keyFileName << "'. "
		          << "cannot construct wave description." << endl;
		return vector<waveDescriptionPtr>();
	}
	const vector<configPtr> configs = expand(config);

	vector<waveDescriptionPtr> waveDescriptions;
	for (size_t c = 0; c < configs.size(); ++c) {
		const string confString = getConfigString(*configs[c]);

		waveDescriptionPtr waveDesc(new waveDescription);
		waveDesc->_keyFileLocalCopy = confString;
		const bool result = waveDesc->parseKeyFileLocalCopy();
		if (!result) // error message already printed in parseKeyFileLocalCopy
			return vector<waveDescriptionPtr>();

		waveDescriptions.push_back(waveDesc);
	}
	return waveDescriptions;
}


vector<waveDescriptionPtr>
waveDescription::parseKeyFileContent(const string& keyFileContent)
{
	if (keyFileContent == "") {
		printWarn << "empty key file content string. cannot construct wave description." << endl;
		return vector<waveDescriptionPtr>();
	}
	configPtr config(new Config);
	if (not parseLibConfigString(keyFileContent, *config, _debug)) {
		printWarn << "problems reading key file content string:" << endl;
		printKeyFileContent(cout, keyFileContent);
		cout << "    cannot construct wave description." << endl;
		return vector<waveDescriptionPtr>();
	}
	const vector<configPtr> configs = expand(config);

	vector<waveDescriptionPtr> waveDescriptions;
	for (size_t c = 0; c < configs.size(); ++c) {
		const string confString = getConfigString(*configs[c]);

		waveDescriptionPtr waveDesc(new waveDescription);
		waveDesc->_keyFileLocalCopy = confString;
		const bool result = waveDesc->parseKeyFileLocalCopy();
		if (!result) // error message already printed in parseKeyFileLocalCopy
			return vector<waveDescriptionPtr>();

		waveDescriptions.push_back(waveDesc);
	}
	return waveDescriptions;
}


bool
waveDescription::parseKeyFileLocalCopy()
{
	_keyFileParsed = false;
	if (_keyFileLocalCopy == "") {
		printWarn << "empty key file content string. cannot construct wave description." << endl;
		return false;
	}
	if (not _key)
		_key = new Config();
	if (not parseLibConfigString(_keyFileLocalCopy, *_key, _debug)) {
		printWarn << "problems reading key file content string:" << endl;
		printKeyFileContent(cout);
		cout << "    cannot construct wave description." << endl;
		delete _key;
		_key = 0;
		return false;
	}
	_keyFileParsed = true;
	return true;
}


ostream&
waveDescription::printKeyFileContent(ostream& out, const string& keyFileContent)
{
	if (keyFileContent != "") {
		typedef tokenizer<char_separator<char> > tokenizer;
		char_separator<char> separator("\n");
		tokenizer            keyFileLines(keyFileContent, separator);
		unsigned int         lineNumber = 0;
		for (tokenizer::iterator i = keyFileLines.begin(); i != keyFileLines.end(); ++i)
			out << setw(5) << ++lineNumber << "  " << *i << endl;
	} else
		out << "key file content string is empty" << endl;
	return out;
}


ostream&
waveDescription::printKeyFileContent(ostream& out) const
{
	return printKeyFileContent(out, _keyFileLocalCopy);
}


bool
waveDescription::constructDecayTopology(isobarDecayTopologyPtr& topo,
                                        const bool              fromTemplate) const
{
	if (not _key or not _keyFileParsed) {
		printWarn << "parsing was not successful. cannot construct decay topology." << endl;
		return false;
	}

	if (topo)
		topo.reset();
	topo = isobarDecayTopologyPtr();  // null pointer

	const Setting& rootKey = _key->getRoot();
	printInfo << "constructing decay topology from key file" << endl;

	// find wave group
	const Setting* decayVertKey = findLibConfigGroup(rootKey, "decayVertex");
	if (not decayVertKey) {
		printWarn << "cannot find 'wave' group. cannot construct decay topology." << endl;
		return false;
	}

	// find  X quantum numbers group
	const Setting* XQnKey = findLibConfigGroup(*decayVertKey, "XQuantumNumbers", not fromTemplate);
	particlePtr    X;
	if (not XQnKey)
		if (not fromTemplate) {
			printWarn << "cannot find 'XQuantumNumbers' group. cannot construct decay topology." << endl;
			return false;
		} else {
			// create default X particle
			if (_debug)
				printDebug << "cannot find 'XQuantumNumbers' group. "
				           << "creating X particle with default properties." << endl;
			X = createParticle("X", not fromTemplate);
		}
	else
		// create X particle from key
		if (not constructXParticle(*XQnKey, X)) {
			printWarn << "problems constructing X particle. cannot construct decay topology." << endl;
			return false;
		}

	// create production vertex
	productionVertexPtr prodVert = productionVertexPtr();
	if (not constructProductionVertex(rootKey, X, prodVert)) {
		printWarn << "problems constructing production vertex. cannot construct decay topology." << endl;
		return false;
	}

	// find X decay group
	const Setting* XDecayKey = findLibConfigGroup(*decayVertKey, "XDecay");
	if (not XDecayKey) {
		printWarn << "cannot find 'XDecay' group. cannot construct decay topology." << endl;
		return false;
	}
	// traverse decay chain and create final state particles and isobar decay vertices
	vector<isobarDecayVertexPtr> decayVertices;
	vector<particlePtr>          fsParticles;
	if (not constructDecayVertex(*XDecayKey, X, decayVertices, fsParticles, fromTemplate)) {
		printWarn << "problems constructing decay chain. cannot construct decay topology." << endl;
		return false;
	}

	// create decay isobar decay topology
	topo = createIsobarDecayTopology(prodVert, decayVertices, fsParticles);
	// backward compatibility: allow sloppy key files, where charges of
	// isobars are not explicitely defined
	topo->calcIsobarCharges();

	printSucc << "constructed decay topology from key file" << endl;
	return true;
}


bool
waveDescription::constructAmplitude(isobarAmplitudePtr& amplitude) const
{
	isobarDecayTopologyPtr topo;
	if (not constructDecayTopology(topo)) {
		printWarn << "problems constructing decay topology. cannot construct decay amplitude." << endl;
		return false;
	}
	return constructAmplitude(amplitude, topo);
}


bool
waveDescription::constructAmplitude(isobarAmplitudePtr&           amplitude,
                                    const isobarDecayTopologyPtr& topo) const
{
	if (not topo) {
		printWarn << "null pointer to decay topology. cannot construct decay amplitude." << endl;
		return false;
	}
	if (not topo->checkTopology() or not topo->checkConsistency()) {
		printWarn << "decay topology has issues. cannot construct decay amplitude." << endl;
		return false;
	}
	if (amplitude)
		amplitude.reset();
	// get amplitude parameters from key file
	// default values
	string formalism            = "helicity";
	bool   boseSymmetrize       = true;
	bool   isospinSymmetrize    = true;
	bool   useReflectivityBasis = true;
	// find amplitude group
	const Setting&            rootKey      = _key->getRoot();
	const Setting* amplitudeKey = findLibConfigGroup(rootKey, "amplitude", false);
	if (amplitudeKey) {
		if (amplitudeKey->lookupValue("formalism", formalism) and _debug)
			printDebug << "setting amplitude formalism to '" << formalism << "'" << endl;
		if (amplitudeKey->lookupValue("boseSymmetrize", boseSymmetrize) and _debug)
			printDebug << "setting amplitude option 'boseSymmetrize' to "
			           << trueFalse(boseSymmetrize) << endl;
		if (amplitudeKey->lookupValue("isospinSymmetrize", isospinSymmetrize) and _debug)
			printDebug << "setting amplitude option 'isospinSymmetrize' to "
			           << trueFalse(isospinSymmetrize) << endl;
		if (amplitudeKey->lookupValue("useReflectivityBasis", useReflectivityBasis) and _debug)
			printDebug << "setting amplitude option 'useReflectivityBasis' to "
			           << trueFalse(useReflectivityBasis) << endl;
	}
	// construct amplitude
	amplitude = mapAmplitudeType(formalism, topo);
	if (_debug)
		printDebug << "constructed amplitude '"<< amplitude->name() << "': "
		           << enDisabled(boseSymmetrize      ) << " Bose symmetrization, "
		           << enDisabled(isospinSymmetrize   ) << " isospin symmetrization, "
		           << enDisabled(useReflectivityBasis) << " reflectivity basis" << endl;
	if (not amplitude) {
		printWarn << "problems constructing decay amplitude." << endl;
		return false;
	}
	amplitude->enableBoseSymmetrization   (boseSymmetrize      );
	amplitude->enableIsospinSymmetrization(isospinSymmetrize   );
	amplitude->enableReflectivityBasis    (useReflectivityBasis);
	return true;
}


namespace {


	string
	nameForMassDependence(const particlePtr          P,
	                      const isobarDecayVertexPtr vertex)
	{
		ostringstream name;
		if (vertex->massDependence() && vertex->massDependence()->name() == "binned") {
			binnedMassDependencePtr massDep = static_pointer_cast<binnedMassDependence>(vertex->massDependence());
			name << "binned[" << spinQn(P->isospin()) << parityQn(P->G()) << ","
			     << spinQn(P->J()) << parityQn(P->P()) << parityQn(P->C()) << ","
			     << massDep->getMassMin() << "," << massDep->getMassMax() << "]";
		} else if(vertex->massDependence() && vertex->massDependence()->name() == "sawtooth") {
			sawtoothMassDependencePtr massDep = static_pointer_cast<sawtoothMassDependence>(vertex->massDependence());
			name << "sawtooth[" << spinQn(P->isospin()) << parityQn(P->G()) << ","
			     << spinQn(P->J()) << parityQn(P->P()) << parityQn(P->C()) << ","
			     << massDep->getMassMin() << "," << massDep->getMassMax() << "]";
		} else if (vertex->massDependence() && vertex->massDependence()->name() == "polynomial") {
			polynomialMassDependencePtr massDep = static_pointer_cast<polynomialMassDependence>(vertex->massDependence());
			name << "polynomial[" << spinQn(P->isospin()) << parityQn(P->G()) << ","
			     << spinQn(P->J()) << parityQn(P->P()) << parityQn(P->C()) << ","
			     << massDep->getMassMin() << "," << massDep->getMassMax();
			const vector<complex<double> >& coefficients = massDep->getCoefficients();
			for (size_t i = 0; i < coefficients.size(); ++i) {
				name << ",[" << coefficients[i].real() << "," << coefficients[i].imag() << "]";
			}
			name << "]";
		} else if (vertex->massDependence() && vertex->massDependence()->name() == "complexExponential") {
			complexExponentialMassDependencePtr massDep = static_pointer_cast<complexExponentialMassDependence>(vertex->massDependence());
			name << "complexExponential[" << spinQn(P->isospin()) << parityQn(P->G()) << ","
			     << spinQn(P->J()) << parityQn(P->P()) << parityQn(P->C()) << ","
			     << massDep->getDegree() << "," << massDep->getMassMin() << "," << massDep->getMassMax() << "]";
		} else if (vertex->massDependence() && vertex->massDependence()->name() == "arbitraryFunction") {
			arbitraryFunctionMassDependencePtr massDep = static_pointer_cast<arbitraryFunctionMassDependence>(vertex->massDependence());
			name << "arbitrary[" << spinQn(P->isospin()) << parityQn(P->G()) << ","
			     << spinQn(P->J()) << parityQn(P->P()) << parityQn(P->C()) << ","
			     << massDep->getRealFunctionString() << "," << massDep->getImagFunctionString() << "]";
		} else
			name << P->name();

		return name.str();
	}


}


string
waveDescription::waveNameFromTopology(isobarDecayTopology         topo,
                                      const isobarDecayVertexPtr& currentVertex)
{
	ostringstream waveName;
	if (currentVertex == interactionVertexPtr()) {
		if (not topo.checkTopology() or not topo.checkConsistency()) {
			printWarn << "decay topology has issues. cannot construct wave name." << endl;
			return "";
		}
		// X quantum numbers
		const particle& X = *(topo.XParticle());
		waveName << "[" << spinQn(X.isospin()) << parityQn(X.G()) << ","
		         << spinQn(X.J()) << parityQn(X.P()) << parityQn(X.C()) << ","
		         << spinQn(X.spinProj()) << parityQn(X.reflectivity()) << "]"
		         << waveNameFromTopology(topo, topo.XIsobarDecayVertex());
	} else {
		// recurse down decay chain
		// first daughter
		waveName << "=[";
		if (topo.isFsParticle(currentVertex->daughter1()))
			waveName << currentVertex->daughter1()->name();
		else {
			particlePtr P = currentVertex->daughter1();
			isobarDecayVertexPtr vertex = static_pointer_cast<isobarDecayVertex>(topo.toVertex(P));
			waveName << nameForMassDependence(P, vertex);
			waveName << waveNameFromTopology(topo, vertex);
		}

		// L, S
		waveName << "[" << spinQn(currentVertex->L()) << "," << spinQn(currentVertex->S()) << "]";

		// second daughter
		if (topo.isFsParticle(currentVertex->daughter2()))
			waveName << currentVertex->daughter2()->name();
		else {
			particlePtr P = currentVertex->daughter2();
			isobarDecayVertexPtr vertex = static_pointer_cast<isobarDecayVertex>(topo.toVertex(P));
			waveName << nameForMassDependence(P, vertex);
			waveName << waveNameFromTopology(topo, vertex);
		}
		waveName << "]";
	}
	{
		string name = waveName.str();
		replace_all(name, "(", "_");
		replace_all(name, ")", "_");
		return name;
	}
}


string
waveDescription::waveLaTeXFromTopology(isobarDecayTopology         topo,
                                       const isobarDecayVertexPtr& currentVertex)
{
	ostringstream waveLaTeX;
	if (currentVertex == interactionVertexPtr()) {
		if (not topo.checkTopology() or not topo.checkConsistency()) {
			printWarn << "decay topology has issues. cannot construct wave LaTeX." << endl;
			return "";
		}
		// X quantum numbers
		const particle& X = *(topo.XParticle());
		waveLaTeX << spinQn(X.isospin()) << "^{"<< parityQn(X.G()) << "}"
		          << spinQn(X.J()) << "^{" << parityQn(X.P()) << parityQn(X.C()) << "}"
		          << spinQn(X.spinProj()) << "^{" << parityQn(X.reflectivity()) << "}\\quad & "
		          << waveLaTeXFromTopology(topo, topo.XIsobarDecayVertex());
	}
	else if (not (topo.isFsParticle(currentVertex->daughter1())
	             and topo.isFsParticle(currentVertex->daughter2()))) {
		// recurse down decay chain
		// do this only if not both daughters are fs partiles

		const bool isXdecay = ( currentVertex == topo.XIsobarDecayVertex() );

		// first daughter
		string dau1 = isobars[currentVertex->daughter1()->name()];
		if (dau1.length()<2) {
			dau1  = "{\\bf ";
			dau1 += currentVertex->daughter1()->name();
			dau1 += "}";
		}
		if (!isXdecay)
			waveLaTeX << "\\rightarrow\\left\\{ ";
		waveLaTeX << dau1;
		if (not topo.isFsParticle(currentVertex->daughter1()))
			waveLaTeX << waveLaTeXFromTopology(topo, static_pointer_cast<isobarDecayVertex>(topo.toVertex(currentVertex->daughter1())));

		// L, S
		waveLaTeX << "\\ells{" << spinQn(currentVertex->L()) << "}{" << spinQn(currentVertex->S()) << "}";

		// second daughter
		string dau2=isobars[currentVertex->daughter2()->name()];
		if (dau2.length()<2) {
			dau2  = "{\\bf ";
			dau2 += currentVertex->daughter2()->name();
			dau2 += "}";
		}
		waveLaTeX << dau2;
		if (not topo.isFsParticle(currentVertex->daughter2()))
			waveLaTeX << waveLaTeXFromTopology(topo, static_pointer_cast<isobarDecayVertex>(topo.toVertex(currentVertex->daughter2())));
		if (!isXdecay)
			waveLaTeX << "\\right\\} ";
	}

	return waveLaTeX.str();
}


bool
waveDescription::readKeyFileIntoLocalCopy(const string& keyFileName)
{
	// read key file content into string
	ifstream keyFile(keyFileName.c_str());
	if (not keyFile or not keyFile.good()) {
		printWarn << "cannot read from file '"  << keyFileName << "'" << endl;
		return false;
	}
	_keyFileLocalCopy = "";
	string line;
	while(getline(keyFile, line))
		_keyFileLocalCopy += line + "\n";
	return true;
}


bool
waveDescription::constructXParticle(const Setting& XQnKey,
                                    particlePtr&   X)
{
	if (_debug)
		printDebug << "reading X quantum numbers from '" << XQnKey.getPath() << "':" << endl;
	// get X quantum numbers
	map<string, int> mandatoryXQn, optionalXQn;
	mandatoryXQn["isospin"    ];
	optionalXQn ["G"          ];
	mandatoryXQn["J"          ];
	mandatoryXQn["P"          ];
	optionalXQn ["C"          ];
	mandatoryXQn["M"          ];
	optionalXQn ["refl"       ];
	optionalXQn ["baryonNmb"  ];
	optionalXQn ["strangeness"];
	optionalXQn ["charm"      ];
	optionalXQn ["beauty"     ];
	bool success = true;
	for (map<string, int>::iterator i = mandatoryXQn.begin(); i != mandatoryXQn.end(); ++i)
		if (not XQnKey.lookupValue(i->first, i->second)) {
			printWarn << "cannot find integer field '" << i->first << "in "
			          << "'" << XQnKey.getPath() << "'" << endl;
			success = false;
		}
	for (map<string, int>::iterator i = optionalXQn.begin(); i != optionalXQn.end(); ++i)
		if (not XQnKey.lookupValue(i->first, i->second))
			i->second = 0;
	if (not success) {
		printWarn << "cannot find X quantum numbers. cannot construct decay topology." << endl;
		return false;
	}
	// create X particle
	X = createParticle("X",
	                   mandatoryXQn["isospin"], optionalXQn["G"],
	                   mandatoryXQn["J"], mandatoryXQn["P"], optionalXQn["C"],
	                   mandatoryXQn["M"], optionalXQn["refl"]);
	X->setBaryonNmb(optionalXQn["baryonNmb"]);
	X->setSCB(optionalXQn["strangeness"], optionalXQn["charm"], optionalXQn["beauty"]);
	if (_debug)
		printDebug << "constructed X particle: " << X->qnSummary() << endl;
	return true;
}


productionVertexPtr
waveDescription::mapProductionVertexType(const Setting&       prodVertKey,
                                         const string&        vertType,
                                         const particlePtr&   X)
{
	productionVertexPtr prodVert;
	bool                success = true;
	if ((vertType == "diffractiveDissVertex") or (vertType == "leptoProductionVertex")) {
		const string             prodKinParticleNames[] = {"beam", "target", "recoil"};
		map<string, particlePtr> prodKinParticles;
		for (unsigned int i = 0; i < sizeof(prodKinParticleNames) / sizeof(string); ++i)
			prodKinParticles[prodKinParticleNames[i]] = particlePtr();
		// construct particles
		const Setting* beamParticleKey = 0;  // needed for leptoproduction vertex
		for (map<string, particlePtr>::iterator i = prodKinParticles.begin();
		     i != prodKinParticles.end(); ++i) {
			const Setting* particleKey = findLibConfigGroup(prodVertKey, i->first,
			                                                (i->first != "recoil") ? true : false);
			if (particleKey) {
				if (i->first == "beam")
					beamParticleKey = particleKey;
				if (not constructParticle(*particleKey, i->second))
					success = false;
			} else if (i->first != "recoil")
				success = false;
		}
		if (success) {
			if (vertType == "diffractiveDissVertex")
				prodVert = createDiffractiveDissVertex(prodKinParticles["beam"], prodKinParticles["target"],
				                                       X, prodKinParticles["recoil"]);
			else if (vertType == "leptoProductionVertex" ) {
				prodVert = createLeptoProductionVertex(prodKinParticles["beam"], prodKinParticles["target"],
				                                       X, prodKinParticles["recoil"]);
				double beamLongPol = 0;
				if ((beamParticleKey->lookupValue("longPol", beamLongPol))) {
					if (_debug)
						printDebug << "setting polarization of beam " << prodKinParticles["beam"]->qnSummary()
						           << " to " << beamLongPol << endl;
				} else
					printWarn << "no polarization is given for beam " << prodKinParticles["beam"]->qnSummary()
					          << ". assuming unpolarized beam." << endl;
				static_pointer_cast<leptoProductionVertex>(prodVert)->setBeamPol(beamLongPol);
			}
		}
	} else
		printWarn << "unknown production vertex type '" << vertType << "'" << endl;
	if (_debug) {
		if (success and prodVert)
			printDebug << "constructed " << *prodVert << endl;
		else
			printWarn << "problems constructing production vertex of type '" << vertType << "'" << endl;
	}
	return prodVert;
}


bool
waveDescription::constructProductionVertex(const Setting&       rootKey,
                                           const particlePtr&   X,
                                           productionVertexPtr& prodVert)
{
	// find production vertex group
	const Setting* prodVertKey = findLibConfigGroup(rootKey, "productionVertex");
	if (not prodVertKey) {
		printWarn << "cannot find 'productionVertex' group" << endl;
		return false;
	}
	if (_debug)
		printDebug << "reading production vertex from '" << prodVertKey->getPath() << "':" << endl;
	//bool success = true;
	// get vertex type
	string vertType;
	if (not prodVertKey->lookupValue("type", vertType)) {
		printWarn << "cannot find 'type' entry in '" << prodVertKey->getPath() << "'" << endl;
		//success = false;
	}
	// create production vertex
	prodVert = mapProductionVertexType(*prodVertKey, vertType, X);
	return (prodVert) ? true : false;
}


bool
waveDescription::constructParticle(const Setting& particleKey,
                                   particlePtr&   particle,
                                   const bool     requirePartInTable)
{
	if (_debug)
		printDebug << "reading particle information from '" << particleKey.getPath() << "':" << endl;
	string name;
	if (particleKey.lookupValue("name", name)) {
		particle = createParticle(name, requirePartInTable);
		int spinProj;
		if (particleKey.lookupValue("spinProj", spinProj))
			particle->setSpinProj(spinProj);
		if (_debug)
			printDebug << "created particle " << particle->qnSummary() << flush;
		int index;
		if (particleKey.lookupValue("index", index)) {
			particle->setIndex(index);
			if (_debug)
				cout << " with index [" << index << "]" << flush;
		}
		if (_debug)
			cout << endl;
		return true;
	} else {
		printWarn << "cannot construct particle from entry '" << particleKey.getPath() << "'. "
		          << "name field is missing." << endl;
		particle = particlePtr();
		return false;
	}
}


massDependencePtr
waveDescription::mapMassDependenceType(const Setting* massDepKey)
{
	string massDepType = "";
	if (massDepKey) {
		// NULL is a valid value for massDepKey, in this case the
		// 'massDep' group is not present in the keyfile, and the
		// default 'BreitWigner' is to be used.
		massDepKey->lookupValue("name", massDepType);
	}

	massDependencePtr massDep;
	if (   (massDepType == "relativisticBreitWigner")
	    or (massDepType == ""))  // default mass dependence
		massDep = createRelativisticBreitWigner();
	else if (massDepType == "constWidthBreitWigner")
		massDep = createConstWidthBreitWigner();
	else if (massDepType == "flat")
		massDep = createFlatMassDependence();
	else if (massDepType == "f_0(980)")
		massDep = createF0980BreitWigner();
	else if (massDepType == "f_0(980)Flatte")
		massDep = createF0980Flatte();
	else if (massDepType == "piPiSWaveAuMorganPenningtonM")
		massDep = createPiPiSWaveAuMorganPenningtonM();
	else if (massDepType == "piPiSWaveAuMorganPenningtonVes")
		massDep = createPiPiSWaveAuMorganPenningtonVes();
	else if (massDepType == "piPiSWaveAuMorganPenningtonKachaev")
		massDep = createPiPiSWaveAuMorganPenningtonKachaev();
	else if (massDepType == "rhoPrime")
		massDep = createRhoPrimeMassDep();
	else if (massDepType == "binned") {
		const Setting* bounds = findLibConfigList(*massDepKey, "bounds" , false);
		if (not bounds) {
			printErr << "no bounds given for binned mass dependence." << endl;
			throw;
		}
		const int length = bounds->getLength();
		if (length != 2) {
			printErr << "bounds do not have the required length (expected: 2, found: " << length << ")." << endl;
			throw;
		}
		const double mMin = (*bounds)[0];
		const double mMax = (*bounds)[1];
		if (mMin > mMax) {
			printErr << "bounds are not ordered: mMin(" << mMin << ") > mMax(" << mMax << ")." << endl;
			throw;
		}
		massDep = createBinnedMassDependence(mMin, mMax);
	}
	else if (massDepType == "sawtooth") {
		const Setting* bounds = findLibConfigList(*massDepKey, "bounds" , false);
		if (not bounds) {
			printErr << "no bounds given for sawtooth mass dependence." << endl;
			throw;
		}
		const int length = bounds->getLength();
		if (length != 2) {
			printErr << "bounds do not have the required length (expected: 2, found: " << length << ")." << endl;
			throw;
		}
		const double mMin = (*bounds)[0];
		const double mMax = (*bounds)[1];
		if (mMin > mMax) {
			printErr << "bounds are not ordered: mMin(" << mMin << ") > mMax(" << mMax << ")." << endl;
			throw;
		}
		massDep = createSawtoothMassDependence(mMin, mMax);
	}
	else if (massDepType == "polynomial") {
		double mMin = 0.;
		bool mMinFound = massDepKey->lookupValue("mMin", mMin);
		double mMax = 0.;
		bool mMaxFound = massDepKey->lookupValue("mMax", mMax);
		if (mMinFound != mMaxFound) {
			printErr << "something wrong with the setting of the mass boundaries. Aborting..." << endl;
			throw;
		}
		if (not mMinFound) {
			mMin = -1.; // No mass boundaries were given. So no scaling is set.
			mMax =  1.; // This is achieved by setting mMin = -1. and mMax = 1.
		}
		const Setting* realCoefficients = findLibConfigList(*massDepKey, "realCoefficients" , false);
		if (not realCoefficients) {
			printErr << "no real coefficients given for polynomial mass dependence" << endl;
			throw;
		}
		const int length = realCoefficients->getLength();
		bool hasImagCoefficients = true;
		const Setting* imagCoefficients = findLibConfigList(*massDepKey, "imagCoefficients" , false);
		if (not imagCoefficients) {
			hasImagCoefficients = false;
		} else {
			if (length != imagCoefficients->getLength()) {
				printErr << "different number of real and imag coefficients given." << endl;
				throw;
			}
		}
		vector<complex<double> > coefficients;
		for (int i = 0; i < length; ++i) {
			if (hasImagCoefficients) {
				coefficients.push_back(complex<double>((*realCoefficients)[i], (*imagCoefficients)[i]));
			} else {
				coefficients.push_back(complex<double>((*realCoefficients)[i], 0.));
			}
		}
		massDep = createPolynomialMassDependencePtr(coefficients, mMin, mMax);
	}
	else if (massDepType == "complexExponential") {
		int degree = 0;
		if (not massDepKey->lookupValue("degree", degree)) {
			printErr << "no degree specified" << endl;
			throw;
		}
		double mMin = 0.;
		if (not massDepKey->lookupValue("mMin", mMin)) {
			printErr << "no mMin specified" << endl;
			throw;

		}
		double mMax = 0.;
		if (not massDepKey->lookupValue("mMax", mMax)) {
			printErr << "no mMax specified" << endl;
			throw;
		}
		massDep = createComplexExponentialMassDependencePtr(degree, mMin, mMax);
	}
	else if (massDepType == "arbitraryFunction") {
		string functionName = "";
		if (not massDepKey->lookupValue("functionName", functionName)) {
			printErr << "no functionName specified" << endl;
			throw;
		}
		string realFunctionString = "";
		if (not massDepKey->lookupValue("realFunction", realFunctionString)) {
			printErr << "no realFuntion specified" << endl;
			throw;
		}
		string imagFunctionString = "";
		if (not massDepKey->lookupValue("imagFunction", imagFunctionString)) {
			imagFunctionString = "0.";
		}
		massDep = createArbitraryFunctionMassDependencePtr(functionName, realFunctionString, imagFunctionString);
	}
	else {
		printWarn << "unknown mass dependence '" << massDepType << "'. using Breit-Wigner." << endl;
		massDep = createRelativisticBreitWigner();
	}
	return massDep;
}


bool
waveDescription::constructDecayVertex(const Setting&                parentKey,
                                      const particlePtr&            parentParticle,
                                      vector<isobarDecayVertexPtr>& decayVertices,
                                      vector<particlePtr>&          fsParticles,
                                      const bool                    fromTemplate)
{
	if (_debug)
		printDebug << "reading decay vertex information from '" << parentKey.getPath() << "':" << endl;
	bool success = true;

	const Setting*      fsPartKeys = findLibConfigList(parentKey, "fsParticles", false);
	vector<particlePtr> fsDaughters;
	if (fsPartKeys)
		for (int i = 0; i < fsPartKeys->getLength(); ++i) {
			particlePtr p;
			if (constructParticle((*fsPartKeys)[i], p)) {
				fsDaughters.push_back(p);
				fsParticles.push_back(p);
			} else
				success = false;
		}

	const Setting*      isobarKeys = findLibConfigList(parentKey, "isobars", false);
	vector<particlePtr> isobarDaughters;
	if (isobarKeys)
		for (int i = 0; i < isobarKeys->getLength(); ++i) {
			particlePtr p;
			if (constructParticle((*isobarKeys)[i], p, not fromTemplate))
				isobarDaughters.push_back(p);
			else
				success = false;
			success &= constructDecayVertex((*isobarKeys)[i], isobarDaughters.back(),
			                                decayVertices, fsParticles, fromTemplate);
		}

	const unsigned int nmbDaughters = fsDaughters.size() + isobarDaughters.size();
	if (nmbDaughters != 2) {
		printWarn << "cannot construct isobar vertex, because number of daughters "
		          << "(" << nmbDaughters << ") in '" << parentKey.getPath() << "' "
		          << "is not 2" << endl;
		success = false;
	}
	if (not success)
		return false;

	// get isobar vertex parameters
	int L = 0, S = 0;
	if ((   not parentKey.lookupValue("L", L)
	     or not parentKey.lookupValue("S", S)) and not fromTemplate)
		printWarn << "Either L or S are not specified in '" << parentKey.getPath() << "'. "
		          << "using zero." << endl;

	// get mass dependence
	massDependencePtr massDep;
	if (parentParticle->bareName() != "X") {
		const Setting* massDepKey  = findLibConfigGroup(parentKey, "massDep", false);
		massDep = mapMassDependenceType(massDepKey);
	}

	// if there is 1 final state particle and 1 isobar put them in the
	// same order as in the key file
	vector<particlePtr> daughters(2, particlePtr());
	if ((fsDaughters.size() == 1) and (isobarDaughters.size() == 1)) {
		if (fsPartKeys->getIndex() < isobarKeys->getIndex()) {
			daughters[0] = fsDaughters    [0];
			daughters[1] = isobarDaughters[0];
		} else {
			daughters[0] = isobarDaughters[0];
			daughters[1] = fsDaughters    [0];
		}
	} else if (fsDaughters.size() == 2)
		daughters = fsDaughters;
	else if (isobarDaughters.size() == 2)
		daughters = isobarDaughters;

	// construct isobar decay vertex
	decayVertices.push_back(createIsobarDecayVertex(parentParticle, daughters[0],
	                                                daughters[1], L, S, massDep));
	if (_debug)
		printDebug << "constructed " << *decayVertices.back() << endl;
	return true;
}


isobarAmplitudePtr
waveDescription::mapAmplitudeType(const string&                 formalismType,
                                  const isobarDecayTopologyPtr& topo)
{
	isobarAmplitudePtr amplitude;
	if (formalismType == "helicity")
		amplitude = createIsobarHelicityAmplitude(topo);
	else if (formalismType == "canonical")
		amplitude = createIsobarCanonicalAmplitude(topo);
	else {
		printWarn << "unknown amplitude formalism '" << formalismType << "'. "
		          << "using helicity formalism." << endl;
		amplitude = createIsobarHelicityAmplitude(topo);
	}
	return amplitude;
}


bool
waveDescription::setProductionVertexKeys(Setting&                   prodVertKey,
                                         const productionVertexPtr& prodVert)
{
	const string prodVertName = prodVert->name();
	if ((prodVertName == "diffractiveDissVertex") or (prodVertName == "leptoProductionVertex")) {
		if (_debug)
			printDebug << "setting keys for " << *prodVert << endl;
		prodVertKey.add("type", Setting::TypeString) = prodVertName;
		map<string, particlePtr> prodKinParticles;
		if (prodVertName == "diffractiveDissVertex") {
			const diffractiveDissVertexPtr& vert = static_pointer_cast<diffractiveDissVertex>(prodVert);
			prodKinParticles["beam"  ] = vert->beam();
			prodKinParticles["target"] = vert->target();
			if (vert->recoil()->name() != vert->target()->name())
				prodKinParticles["recoil"] = vert->recoil();
		} else if (prodVertName == "leptoProductionVertex") {
			const leptoProductionVertexPtr& vert = static_pointer_cast<leptoProductionVertex>(prodVert);
			prodKinParticles["beam"  ] = vert->beamLepton();
			prodKinParticles["target"] = vert->target();
			if (vert->recoil()->name() != vert->target()->name())
				prodKinParticles["recoil"] = vert->recoil();
		}
		// write particles
		Setting* beamParticleKey = 0;  // needed for leptoproduction vertex
		for (map<string, particlePtr>::iterator i = prodKinParticles.begin();
		     i != prodKinParticles.end(); ++i) {
			Setting& particleKey = prodVertKey.add(i->first, Setting::TypeGroup);
			particleKey.add("name", Setting::TypeString) = i->second->name();
			if (i->first == "beam")
				beamParticleKey = &particleKey;
		}
		if ((prodVertName == "leptoProductionVertex") and beamParticleKey)
			beamParticleKey->add("longPol", Setting::TypeFloat)
				= static_pointer_cast<leptoProductionVertex>(prodVert)->beamPol();
		return true;
	} else {
		printWarn << "setting of keys for production vertex of this type is not yet implemented:"
		          << *prodVert << endl;
		return false;
	}
}


bool
waveDescription::setXQuantumNumbersKeys(Setting&        XQnKey,
                                        const particle& X)
{
	if (_debug)
		printDebug << "setting quantum number keys for " << X.qnSummary() << endl;
	XQnKey.add("isospin",       Setting::TypeInt) = X.isospin();
	if (X.G() != 0)
		XQnKey.add("G",           Setting::TypeInt) = X.G();
	XQnKey.add("J",             Setting::TypeInt) = X.J();
	XQnKey.add("P",             Setting::TypeInt) = X.P();
	if (X.C() != 0)
		XQnKey.add("C",           Setting::TypeInt) = X.C();
	XQnKey.add("M",             Setting::TypeInt) = X.spinProj();
	if (X.reflectivity() != 0)
		XQnKey.add("refl",        Setting::TypeInt) = X.reflectivity();
	if (X.baryonNmb() != 0)
		XQnKey.add("baryonNmb",   Setting::TypeInt) = X.baryonNmb();
	if (X.strangeness() != 0)
		XQnKey.add("strangeness", Setting::TypeInt) = X.strangeness();
	if (X.charm() != 0)
		XQnKey.add("charm",       Setting::TypeInt) = X.charm();
	if (X.beauty() != 0)
		XQnKey.add("beauty",      Setting::TypeInt) = X.beauty();
	return true;
}


bool
waveDescription::setMassDependence(Setting&              isobarDecayKey,
                                   const massDependence& massDep,
                                   const bool            XDecay)
{
	const string massDepName = massDep.name();
	if (XDecay && massDepName == "flat")
		// default for X
		return true;
	else if ((not XDecay) && massDepName == "relativisticBreitWigner")
		// default mass dependence for isobars
		return true;
	else {
		if (_debug)
			printDebug << "setting key for '" << massDep << "' mass dependence to "
			           << "'" << massDepName << "'" << endl;
		Setting& massDepKey = isobarDecayKey.add("massDep", Setting::TypeGroup);
		massDepKey.add("name", Setting::TypeString) = massDepName;

		if (massDepName == "binned") {
			// for this mass dependence additionally the mass bounds
			// have to be stored in the keyfile.
			const binnedMassDependence& binned = dynamic_cast<const binnedMassDependence&>(massDep);

			Setting& bounds = massDepKey.add("bounds", Setting::TypeList);
			bounds.add(Setting::TypeFloat) = binned.getMassMin();
			bounds.add(Setting::TypeFloat) = binned.getMassMax();
		}
		if (massDepName == "sawtooth") {
			// for this mass dependence additionally the mass bounds
			// have to be stored in the keyfile.
			const sawtoothMassDependence& binned = dynamic_cast<const sawtoothMassDependence&>(massDep);

			Setting& bounds = massDepKey.add("bounds", Setting::TypeList);
			bounds.add(Setting::TypeFloat) = binned.getMassMin();
			bounds.add(Setting::TypeFloat) = binned.getMassMax();
		}
		if (massDepName == "polynomial") {
			// for this mass dependence additionally the polynomial coefficients have to be stored
			const polynomialMassDependence& polynomial = dynamic_cast<const polynomialMassDependence&>(massDep);
			const vector<complex<double> >& coefficients = polynomial.getCoefficients();
			Setting& realCoefficients = massDepKey.add("realCoefficients", Setting::TypeList);
			Setting& imagCoefficients = massDepKey.add("imagCoefficients", Setting::TypeList);
			for (size_t i = 0; i < coefficients.size(); ++i) {
				realCoefficients.add(Setting::TypeFloat) = coefficients[i].real();
				imagCoefficients.add(Setting::TypeFloat) = coefficients[i].imag();
			}
			massDepKey.add("mMin", Setting::TypeFloat) = polynomial.getMassMin();
			massDepKey.add("mMax", Setting::TypeFloat) = polynomial.getMassMax();
		}
		if (massDepName == "complexExponential") {
			// for this mass dependence additionally the degree and range have to be stored
			const complexExponentialMassDependence& exponential = dynamic_cast<const complexExponentialMassDependence&>(massDep);
			massDepKey.add("degree", Setting::TypeInt) = exponential.getDegree();
			massDepKey.add("mMin", Setting::TypeFloat) = exponential.getMassMin();
			massDepKey.add("mMax", Setting::TypeFloat) = exponential.getMassMax();
		}
		if (massDepName == "arbitraryFunction") {
			// For this mass dependence, the function definitions have to be stored
			const arbitraryFunctionMassDependence& arbitrary = dynamic_cast<const arbitraryFunctionMassDependence&>(massDep);
			massDepKey.add("functionName", Setting::TypeString) = arbitrary.getName();
			massDepKey.add("realFunction", Setting::TypeString) = arbitrary.getRealFunctionString();
			massDepKey.add("imagFunction", Setting::TypeString) = arbitrary.getImagFunctionString();

		}
	}
	return true;
}


bool
waveDescription::setXDecayKeys(Setting&                   parentDecayKey,
                               const isobarDecayTopology& topo,
                               const isobarDecayVertex&   vert)
{
	if (_debug)
		printDebug << "setting keys for decay " << vert << endl;
	setMassDependence(parentDecayKey, *vert.massDependence(), *topo.XIsobarDecayVertex() == vert);
	vector<particlePtr> fsParticles;
	vector<particlePtr> isobars;
	for (unsigned int i = 0; i < vert.nmbOutParticles(); ++i) {
		const particlePtr& part = vert.outParticles()[i];
		if (topo.isFsParticle(part))
			fsParticles.push_back(part);
		else
			isobars.push_back(part);
	}
	bool success = true;
	if (isobars.size() > 0) {
		Setting& isobarsKey = parentDecayKey.add("isobars", Setting::TypeList);
		for (unsigned int i = 0; i < isobars.size(); ++i) {
			Setting& isobarKey = isobarsKey.add(Setting::TypeGroup);
			isobarKey.add("name", Setting::TypeString) = isobars[i]->name();
			if (not setXDecayKeys(isobarKey, topo,
			                      *static_pointer_cast<isobarDecayVertex>(topo.toVertex(isobars[i]))))
				success = false;
		}
	}
	parentDecayKey.add("L", Setting::TypeInt) = (int)vert.L();
	parentDecayKey.add("S", Setting::TypeInt) = (int)vert.S();
	if (fsParticles.size() > 0) {
		Setting& fsParticlesKey = parentDecayKey.add("fsParticles", Setting::TypeList);
		for (unsigned int i = 0; i < fsParticles.size(); ++i) {
			Setting& fsParticleKey = fsParticlesKey.add(Setting::TypeGroup);
			fsParticleKey.add("name", Setting::TypeString) = fsParticles[i]->name();
		}
	}
	return success;
}


bool
waveDescription::setKeysFromTopology(Setting&                   rootKey,
                                     const isobarDecayTopology& topo,
                                     const bool                 setProdVert)
{
	if (not topo.checkTopology() or not topo.checkConsistency()) {
		printWarn << "decay topology has issues. cannot write key file." << endl;
		return false;
	}
	if (_debug)
		printDebug << "setting keys for " << topo;
	if (setProdVert) {
		Setting& prodVertKey = rootKey.add("productionVertex", Setting::TypeGroup);
		setProductionVertexKeys(prodVertKey, topo.productionVertex());
	}
	Setting& decayVertKey = rootKey.add     ("decayVertex",     Setting::TypeGroup);
	Setting& XQnKey       = decayVertKey.add("XQuantumNumbers", Setting::TypeGroup);
	Setting& XDecayKey    = decayVertKey.add("XDecay",          Setting::TypeGroup);
	if (not setXQuantumNumbersKeys(XQnKey, *(topo.XParticle()))) {
		printWarn << "problems setting X quantum numbers" << endl;
		return false;
	}
	if (not setXDecayKeys(XDecayKey, topo, *(topo.XIsobarDecayVertex()))) {
		printWarn << "problems setting X decay" << endl;
		return false;
	}
	return true;
}


bool
waveDescription::setAmplitude(Setting&               amplitudeKey,
                              const isobarAmplitude& amplitude)
{
	if (_debug)
		printDebug << "setting amplitude keys." << endl;
	string       formalism     = "";
	const string amplitudeName = amplitude.name();
	if (amplitudeName == "isobarCanonicalAmplitude") {
		formalism = "canonical";
	} else if (amplitudeName != "isobarHelicityAmplitude") {
		printWarn << "unknown amplitude type '" << amplitudeName << "'" << endl;
		return false;
	}
	if (formalism != "") {
		if (_debug)
			printDebug << "setting 'formalism' key for '" << amplitudeName << "' to "
			           << "'" << formalism << "'" << endl;
		amplitudeKey.add("formalism", Setting::TypeString) = formalism;
	}
	if (not amplitude.boseSymmetrization())
		amplitudeKey.add("boseSymmetrize", Setting::TypeBoolean) = false;
	if (not amplitude.isospinSymmetrization())
		amplitudeKey.add("isospinSymmetrize", Setting::TypeBoolean) = false;
	if (not amplitude.reflectivityBasis())
		amplitudeKey.add("useReflectivityBasis", Setting::TypeBoolean) = false;
	return true;
}


bool
waveDescription::setKeysFromAmplitude(Setting&               rootKey,
                                      const isobarAmplitude& amplitude,
                                      const bool             setProdVert)
{
	if (_debug)
		printInfo << "setting keys for "<< amplitude;
	if (not setKeysFromTopology(rootKey, *amplitude.decayTopology(), setProdVert)) {
		printWarn << "problems setting keys for decay topology" << endl;
		return false;
	}
	Setting& amplitudeKey = rootKey.add("amplitude", Setting::TypeGroup);
	if (not setAmplitude(amplitudeKey, amplitude)) {
		printWarn << "problems setting amplitude paramaters" << endl;
		return false;
	}
	return true;
}


bool
waveDescription::writeKeyFile(FILE&                      outStream,
                              const isobarDecayTopology& topo,
                              const bool                 writeProdVert)
{
	Config   key;
	Setting& rootKey = key.getRoot();
	if (not setKeysFromTopology(rootKey, topo, writeProdVert)) {
		printWarn << "problems writing keys for decay topology. cannot write key file." << endl;
		return false;
	}
	try {
		key.write(&outStream);
	} catch (const FileIOException& ioEx) {
		printWarn << "I/O error while writing key file" << endl;
		return false;
	}
	return true;
}


bool
waveDescription::writeKeyFile(FILE&                  outStream,
                              const isobarAmplitude& amplitude,
                              const bool             writeProdVert)
{
	Config   key;
	Setting& rootKey = key.getRoot();
	if (not setKeysFromAmplitude(rootKey, amplitude, writeProdVert)) {
		printWarn << "problems writing keys for amplitude. cannot write key file." << endl;
		return false;
	}
	try {
		key.write(&outStream);
	} catch (const FileIOException& ioEx) {
		printWarn << "I/O error while writing key file" << endl;
		return false;
	}
	return true;
}
