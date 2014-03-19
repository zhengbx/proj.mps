#include <vector>
#include <string>
#include <cstdlib>

#define BOOST_NO_CXX11_SCOPED_ENUMS
#include "newmat10/newmatap.h"
#include "newmat10/newmatio.h"

using std::vector;
using std::string;

void permute(Matrix& orbs, vector<int>& order);
void read_config(string file, double& thr1p, double& thrnp, int& M, bool& calc_spectra, bool& savemps);
Matrix read_orbitals(string file);
string mktmpdir(const string& prefix);
