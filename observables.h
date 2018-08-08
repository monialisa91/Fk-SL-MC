#include <armadillo>
#include <math.h>
#include <cstdlib>
#include <ctime>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace arma;
using namespace std;

double ipr(vec eigenvectors) {
	int n = eigenvectors.n_elem;
	double ipr1 = 0;
	double ipr2 = 0;

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			ipr1 += eigenvectors(i, j) * eigenvectors(i, j);
			ipr2 += pow(eigenvectors(i, j), 4);
		}
	}

	return ipr2 / ipr1;
}


double correlation() {

}



