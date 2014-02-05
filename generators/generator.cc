
#include<limits>

#include<TRandom3.h>

#include "generator.h"
#include "randomNumberGenerator.h"
#include <reportingUtils.hpp>

using namespace rpwa;
using namespace std;


ostream& generator::convertEventToAscii(ostream&  out,
                                     const particle& beam,
                                     const vector<particle>& finalState)
{
	if(!out) {
		printErr << "output stream is not writable." << endl;
		throw;
	}
	unsigned int nmbDaughters = finalState.size();
	const TLorentzVector& beamLorentzVector = beam.lzVec();
	out << nmbDaughters + 1 << endl;
	// beam particle: geant ID, charge, p_x, p_y, p_z, E
	out << beam.geantId() << " " << beam.charge() << " "
	    << maxPrecision(beamLorentzVector.Px()) << " "
	    << maxPrecision(beamLorentzVector.Py()) << " "
	    << maxPrecision(beamLorentzVector.Pz()) << " "
	    << maxPrecision(beamLorentzVector.E()) << endl;
	for(unsigned int i = 0; i < nmbDaughters; ++i) {
		const TLorentzVector& hadron = finalState[i].lzVec();
		// hadron: geant ID, charge, p_x, p_y, p_z, E
		out << finalState[i].geantId() << " " << finalState[i].charge() << " "
		    << maxPrecision(hadron.Px()) << " "
		    << maxPrecision(hadron.Py()) << " "
		    << maxPrecision(hadron.Pz()) << " "
		    << maxPrecision(hadron.E()) << endl;
	}
	return out;
}

ostream& generator::convertEventToComgeant(ostream& out,
                                           const particle& beam,
                                           const particle& recoil,
                                           const TVector3& vertex,
                                           const vector<rpwa::particle>& finalState,
                                           bool writeBinary)
{
	if(!out) {
		cerr << "Output stream is not writable." << endl;
		throw;
	}

	unsigned int nmbDaughters = finalState.size();
	const TLorentzVector& beamLorentzVector = beam.lzVec();
	const TLorentzVector& recoilLorentzVector = recoil.lzVec();

	if(not writeBinary) { // Write text file.
		// total number of particles including recoil proton and beam particle
		out << nmbDaughters + 1 + 1 << endl;
		// vertex position in cm
		// note that Comgeant's coordinate system is different
		out << maxPrecision(vertex.Z()) << " "
		    << maxPrecision(vertex.X()) << " "
		    << maxPrecision(vertex.Y()) << endl;
		// beam particle: geant ID , -p_z, -p_x, -p_y must go the opposite direction upstream and should be defined as mulike with PID 44 in Comgeant
		out << "44 " << maxPrecision(-beamLorentzVector.Pz()) << " "
		             << maxPrecision(-beamLorentzVector.Px()) << " "
		             << maxPrecision(-beamLorentzVector.Py()) << endl;
		// the recoil proton
		out << "14 " << maxPrecision(recoilLorentzVector.Pz()) << " "
		             << maxPrecision(recoilLorentzVector.Px()) << " "
		             << maxPrecision(recoilLorentzVector.Py()) << endl;
		for (unsigned int i = 0; i < nmbDaughters; ++i) {
			const TLorentzVector& hadron = finalState[i].lzVec();
			// hadron: geant ID, p_z, p_x, p_y
			out << finalState[i].geantId() << " "
			    << maxPrecision(hadron.Pz()) << " "
			    << maxPrecision(hadron.Px()) << " "
			    << maxPrecision(hadron.Py()) << endl;// << " " << hadron->E() << endl;
		}
	} else {
		int intval;
		float floatval;
		intval = (int)nmbDaughters+1+1; out.write((char*)&intval,4);
		// vertex position in cm
		// note that Comgeant's coordinate system is different
		floatval = (float)vertex.Z(); out.write((char*)&floatval,4);
		floatval = (float)vertex.X(); out.write((char*)&floatval,4);
		floatval = (float)vertex.Y(); out.write((char*)&floatval,4);
		// beam particle: geant ID , -p_z, -p_x, -p_y must go the opposite direction upstream and should be defined as mulike with PID 44 in Comgeant
		intval = 44; out.write((char*)&intval,4);
		floatval = (float)-beamLorentzVector.Pz(); out.write((char*)&floatval,4);
		floatval = (float)-beamLorentzVector.Px(); out.write((char*)&floatval,4);
		floatval = (float)-beamLorentzVector.Py(); out.write((char*)&floatval,4);
		// the recoil proton
		intval = 14; out.write((char*)&intval,4);
		floatval = (float)recoilLorentzVector.Pz(); out.write((char*)&floatval,4);
		floatval = (float)recoilLorentzVector.Px(); out.write((char*)&floatval,4);
		floatval = (float)recoilLorentzVector.Py(); out.write((char*)&floatval,4);
		for (unsigned int i = 0; i < nmbDaughters; ++i) {
			const TLorentzVector& hadron = finalState[i].lzVec();
			// hadron: geant ID, p_z, p_x, p_y
			intval = (int)finalState[i].geantId(); out.write((char*)&intval,4);
			floatval = (float)hadron.Pz(); out.write((char*)&floatval,4);
			floatval = (float)hadron.Px(); out.write((char*)&floatval,4);
			floatval = (float)hadron.Py(); out.write((char*)&floatval,4);
		}
	}
	return out;
}
