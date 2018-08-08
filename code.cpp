#include <armadillo>
#include <math.h>
#include <cstdlib>
#include <ctime>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "observables.h" // calculation of observables
#include "some_functions.h" // helping functions

using namespace arma;
using namespace std;

int main() {

	arma_rng::set_seed_random();

	// PARAMETERS
	int L = 4; // linear size of the lattice
	int MC_steps = 4000;
	int t = -1;
	int n = 1; // the multiple of space between the independent configurations
	int space = L * L * n; // space between the independent configurations
	int mes = 100; // number of measurements to average over
	double U = 4.0; // potential
	double cp = U / 2; //chemical potential
	double T = 0.1;
	double beta = 1.0 / T;
	double r, delta, E_new, E0, corr, corr2, ipr, n_electr;
	double delta_sub, delta_sub2, energies, energies2, cv; // heat capacity
	double c_sub, cc, E_norm;
	mat new_hamiltonian;
	mat initial_hamiltonian = random_hamiltonian(U, L, t);
	vec eigval;
	mat eigvec;
	eig_sym(eigval, eigvec, initial_hamiltonian);
	E0 = sum_eigvals(eigval, beta, cp, t);

// thermalisation

	for (int i = 0; i < MC_steps; i++) {
		new_hamiltonian = Swap_sites(initial_hamiltonian);
		eig_sym(eigval, eigvec, new_hamiltonian);
		E_new = sum_eigvals(eigval, beta, cp, t);
		delta = E_new - E0;
		r = ((double) rand() / (RAND_MAX));
		if (delta < 0 || exp(-delta * beta) > r) {
			initial_hamiltonian = new_hamiltonian;
			E0 = E_new;
		}

	}

// measurement
	int n_conf = 0; // iterator over independent configurations
	corr = 0;
	corr2 = 0;
	ipr = 0;
	n_electr = 0;
	delta_sub = 0;
	delta_sub2 = 0;
	energies = 0;
	energies2 = 0;

	for (int k = 0; k < mes; k++) {
		while (n_conf < space) {
			new_hamiltonian = Swap_sites(initial_hamiltonian);
			eig_sym(eigval, eigvec, new_hamiltonian);
			E_new = sum_eigvals(eigval, beta, cp, t);
			delta = E_new - E0;
			r = ((double) rand() / (RAND_MAX));
			if (delta < 0 || exp(-delta * beta) > r) {
				initial_hamiltonian = new_hamiltonian;
				E0 = E_new;
				n_conf++;
			}
		}
		// correlation
		corr += correlation(new_hamiltonian, U, L);
		corr2 += correlation(new_hamiltonian, U, L) * correlation(new_hamiltonian, U, L);
		// IPR
		ipr += ipr(eigvec);
		// NFERM
		n_electr += nferm(eigval, beta, cp);
		//sub_diff
		delta_sub += sub_lattice(new_hamiltonian);
		delta_sub2 += sub_lattice(new_hamiltonian) * sub_lattice(new_hamiltonian);
		// energies
		E_norm = E_new/(L * L);
		energies += E_norm;
		energies2 += E_norm * E_norm;

	}
	// AVERAGES

	energies = energies/mes ;
	energies2 = energies2/mes;
	cv = energies2 - energies * energies;// heat capacity
	cv = cv * beta * beta;

	n_electr = n_electr/mes;

	delta_sub = delta_sub/mes;
	delta_sub2 = delta_sub2/mes;
	c_sub = delta_sub2 - delta_sub * delta_sub;
	c_sub = beta * c_sub;

	corr = corr/mes;
	corr2 = corr2/mes;
	cc = corr2 - corr * corr;
	cc = beta * cc;


}
