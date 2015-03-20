#include <complex>
#include <vector>

using namespace rpwa;

namespace cupwa {
	
	struct vertexdata {
		int JX, Lambda, J1, L1, S1, J2, L2, S2;
		double theta1, phi1, theta2, phi2, wX, wf2, massf2, massX, wpi, qX, qf2;
	};
	
	void DMAccess(vertexdata* data, int n);
	
	void gpuAmpcall(vertexdata* data, int n, const std::vector<std::complex<double> >& cpudata);
	
}
