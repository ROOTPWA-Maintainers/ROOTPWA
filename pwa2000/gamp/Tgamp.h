#include <complex>
#include <vector>
#include <string>
#include <pputil.h>
#include <keyfile.h>



class Tgamp {
public:
  
  std::complex<double> Amp(unsigned int i, event& e);
  void addWave(const string& keyfilename){m_waves.push_back(keyfilename);}

private:
  std::vector<std::string> m_waves;

};
