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
			out[k*n+l] = (num){{0, 1}, {0, 1}, 1};
		}
	}
	range (k, n) {
		out[k*n+k] = (num){{1, 1}, {0, 1}, 1};
	}
	a.im.whole = -a.im.whole;
	out[j*n+j] = a;
	a.re.whole = -a.re.whole;
	out[i*n+i] = a;
	out[i*n+j] = b;
	b.im.whole = -b.im.whole;
	out[j*n+i] = b;
}

int main() {
	num fourier[4][4];
	range(i, 4) {
		range (j, 4) {
			switch ((i*j)%4) {
				case 0:
					fourier[i][j] = (num) {{1, 1}, {0, 1}, 2};
					break;
				case 1:
					fourier[i][j] = (num) {{0, 1}, {1, 1}, 2};
					break;
				case 2:
					fourier[i][j] = (num) {{-1, 1}, {0, 1}, 2};
					break;
				case 3:
					fourier[i][j] = (num) {{0, 1}, {-1, 1}, 2};
			}
		}
	}
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

