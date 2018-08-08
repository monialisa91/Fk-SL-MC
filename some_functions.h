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


/* Pick random non--zero diagonal element of the matrix lattice and swap with
  the random zero diagonal element */

mat Swap_sites(mat lattice) {
	double r;
	int lattice_n1, lattice_n2;
	int size = lattice.n_rows;
	double swap_var;

	do {
		r  = ((double) rand() / (RAND_MAX));
		lattice_n1 = int(r * size);
	} while(lattice(lattice_n1, lattice_n1) < 0.00001);

	do {
		r = ((double) rand() / (RAND_MAX));
		lattice_n2 = int(r * size);
	} while(lattice(lattice_n2, lattice_n2) > 0.00001);


	swap_var = lattice(lattice_n1, lattice_n1);
	lattice(lattice_n1, lattice_n1) = lattice(lattice_n2, lattice_n2);
	lattice(lattice_n2, lattice_n2) = swap_var;

	return lattice;

}

double sum_eigvals(vec eigenvalues, double beta, double cp, int t) {
	int N = eigenvalues.n_elem;
	double E = 0;
	for(int j=0; j<N; j++){
			E += log(1+exp(-beta*eigenvalues(j)-cp));
		}
	return -t*E;
}


/* initial random matrix: the interaction terms are on the diagonal cp = 0.5 */
/* hopping integral -1 included : kinetic term */
/* L--original lattice size */
mat random_hamiltonian(double U, int L, int t) {
	int size_matrix = L*L;
	vec lista(size_matrix);
	mat lattice(size_matrix, size_matrix);
	lattice.zeros();
	// ordered list of terms with U and without U
	for(int i=0; i<size_matrix/2; i++){
		lista(i) = U;
	}
	for(int i=size_matrix/2; i<size_matrix; i++){
			lista(i) = 0;
	}
	vec matrix_shuffled = shuffle(lista); // shuffled elements

	for(int i=0; i<size_matrix; i++) lattice(i, i) = matrix_shuffled(i);
	// hopping integral

	for(int i = 0; i<size_matrix-L; i++) {
		lattice(i, i+L) = -t;
		lattice(i+L, i) = -t;
	}

	for(int i=0; i<L; i++) {
			lattice(i, i+size_matrix-L) = -t;
			lattice(i+size_matrix-L, i) = -t;
	}

	for(int i=0; i<=L-1; i++){
		lattice(i*L, i*L+L-1) = -t;
		lattice(i*L+L-1,  i*L) = -t;
	}

	for(int i=1; i<=L-1; i++) {
		for(int j=0; j<=L-1; j++) {
			lattice(i+j*L-1, i+j*L) = -t;
			lattice(i+j*L, i+j*L-1) = -t;
		}
	}
	return lattice;
}



void print_matrix(mat lattice) {
	int N = lattice.n_rows;
	int M = lattice.n_cols;
	for(int i = 0; i<N; i++) {
		for(int j = 0; j<M; j++){
			cout<<lattice(i,j)<<" ";
		}
		cout<<endl;
	}
}

void print_vector(vec eigvals) {
	int N = eigvals.n_elem;
	for(int  j=0; j<N; j++){
		cout<<eigvals(j)<<" ";
		}
	}

void save_conf(mat hamiltonian, ofstream file) {
	int size = hamiltonian.n_rows;
	for(int i = 0; i<size; i++) {
		if(hamiltonian(i, i) >0.0001) {
			file << 1 << " ";
		}
		else {
			file << 0 << " ";
		}
	}
	file << "\n";
}








