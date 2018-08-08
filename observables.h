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

double ipr(mat eigenvectors) {
	int n = eigenvectors.n_rows;
	double ipr1 = 0;
	double ipr2 = 0;

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			ipr2 += eigenvectors(i, j) * eigenvectors(i, j);
			ipr1 += pow(eigenvectors(i, j), 4);
		}
	}

	return ipr1 / ipr2;
}

double correlation() {

}

// Fermi -- Dirac statistics
double nferm(vec eigenvalues, double beta, double cp) {
	double nf = 0;
	int size = eigenvalues.n_elem;

	for (int i = 0; i < size; i++) {
		nf += 0.5 * (1 - tanh(eigenvalues(i) - cp) * beta / 2);
	}

	return nf;
}

// Density of states (DOS)

vec dos(vec eigenvalues, double U) {
	int size = eigenvalues.n_elem;
	int dos_len = 500;
	int k;
	double dos_min = -4.2;
	double dos_max = 4.2 + U;
	double dos_step = (dos_max - dos_min) / (dos_len);
	double En;
	vec dos_vec(dos_len);

	for (int i = 0; i < size; i++) {
		En = eigenvalues(i);
		En = En - dos_min;
		En = En / dos_step;
		k = (int) En;

		if (k < dos_len) {
			dos_vec(k) += 1;
		}

	}

	return dos_vec;

}

vec ipr_omega(mat eigenvectors, vec eigenvalues, double U, int ms) {
	int n = eigenvalues.n_elem;
	int k;
	double dos_min = -4.2 - U / 2;
	double dos_max = 4.2 + U / 2;
	double dos_step = (dos_max - dos_min) / (ms);
	double E, EA, ipr1, ipr2, ipr;
	vec count_t(ms);
	vec ipr_tt(ms);
	vec ipr_t(ms);

	for (int i = 0; i < n; i++) {
		ipr1 = ipr2 = 0;
		E = eigenvalues(i) - U / 2;
		EA = abs(E) / dos_step;
		k = (int) E;
		if (E < 0)
			k = -k;
		k += int(ms / 2.0) + 1;

		if (k >= 0 && k < ms) {
			for (int j = 0; j < n; j++) {
				ipr1 += pow(eigenvectors(i, j), 4);
				ipr2 += eigenvectors(i, j) * eigenvectors(i, j);
			}
			ipr = ipr1 / ipr2;
			count_t(k) += 1;
			ipr_tt(k) += ipr;

		}

	}

	for(int i = 0; i<ms; i++)  {
		if(count_t(i) > 0) ipr_t(i) += ipr_tt(i)/count_t(i);
	}

	return ipr_tt;
}


