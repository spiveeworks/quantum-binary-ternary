#include "matrix_quadratic.h"

// generates a two-level matrix to eliminate an entry from a given matrix
void eliminator(int n, int i, int j, num *out, num *x) {
	num a = x[j*n+j];
	num b = x[i*n+j];
	num modsq = num_add(num_mod_sq(a), num_mod_sq(b));
	num mod = num_inv_sqrt(modsq);
	// printf("|"); num_print(a); printf("|^2 + |"); num_print(b); printf("|^2 = "); num_print(modsq); printf("\n");
	// printf("1/sqrt("); num_print(modsq); printf(") = "); num_print(mod); printf("\n");
	a = num_mul(a, mod);
	b = num_mul(b, mod);
	range (k, n) {
		range (l, n) {
			out[k*n+l] = num_from_int(0);
		}
	}
	range (k, n) {
		out[k*n+k] = num_from_int(1);
	}
	out[j*n+j] = num_conj_i(a);
	// @Performance make num_mul_int and use it here
	out[i*n+i] = num_mul(a, num_from_int(-1));
	out[i*n+j] = b;
	out[j*n+i] = num_conj_i(b);
}

int main() {
	num phase[4] = {};
	phase[0].c[0][0][0] = 1;
	phase[1].c[1][0][0] = 1;
	phase[2].c[0][0][0] = -1;
	phase[3].c[1][0][0] = -1;
	phase[0].de = 1;
	phase[1].de = 1;
	phase[2].de = 1;
	phase[3].de = 1;
	num fourier[4][4];
	mat_fourier(4, &fourier[0][0], phase, 4);
	printf("fourier: \n");
	mat_print(4, &fourier[0][0], ", ", "\n");

	const int REMOVE = (4*(4-1))/2-1;
	num u[REMOVE][4][4];
	num fred[REMOVE][4][4];
	int j = 0;
	int i = j+1;
	range(ui, REMOVE) {
		num *fprev = ui > 0 ? &fred[ui-1][0][0] : &fourier[0][0];
		num *ucurr = &u[ui][0][0];
		num *fcurr = &fred[ui][0][0];

		eliminator(4, i, j, ucurr, fprev);
		mat_reduce(4, ucurr);
		printf("\nu[%lu]: \n", ui);
		mat_print(4, ucurr, ", ", "\n");

		mat_mul(4, fcurr, ucurr, fprev);
		mat_reduce(4, fcurr);
		printf("\nfred[%lu]: \n", ui);
		mat_print(4, fcurr, ", ", "\n");

		if (i == 3) {
			j++;
			i = j+1;
		} else {
			i++;
		}
	}
}

