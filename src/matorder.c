#include "matrix_quadratic.h"
#include <stdint.h>

typedef uint8_t u8;

#define MAT_DIM 6

void test_perm_fourier_order() {
	num h2[MAT_DIM][MAT_DIM];
	initialize(2, 3, &h2[0][0], GATE_FOURIER, GATE_IDENT);

	num h3[MAT_DIM][MAT_DIM];
	initialize(2, 3, &h3[0][0], GATE_IDENT, GATE_FOURIER);

	u8 m[MAT_DIM] = {};
	int total_ph2 = 0;
	int total_ph3 = 0;
	int total_ph23 = 0;
	bool done = false;
	while (!done) {
		bool invertible = true;
		range(i, MAT_DIM) {
			range(j, i) {
				if (m[i] == m[j]) {
					invertible = false;
				}
			}
		}
		if (invertible) {
			num m_mat[MAT_DIM][MAT_DIM];
			mat_perm(MAT_DIM, &m_mat[0][0], m);
			num result[MAT_DIM][MAT_DIM];

			mat_mul(MAT_DIM, &result[0][0], &m_mat[0][0], &h2[0][0]);
			int order_ph2 = mat_calc_order(MAT_DIM, &result[0][0], false);
			mat_mul(MAT_DIM, &result[0][0], &h2[0][0], &m_mat[0][0]);
			int order_h2p = mat_calc_order(MAT_DIM, &result[0][0], false);
			if (order_ph2 != -1 && order_h2p != -1) {
				/*
				printf("P H2 has order %d, and H2 P has order %d, where P:\n", order_ph2, order_h2p);
				mat_print(MAT_DIM, &m_mat[0][0], ", ", "\n");
				printf("\n");
				*/
				total_ph2 += 1;
			}

			mat_mul(MAT_DIM, &result[0][0], &m_mat[0][0], &h3[0][0]);
			int order_ph3 = mat_calc_order(MAT_DIM, &result[0][0], false);
			mat_mul(MAT_DIM, &result[0][0], &h3[0][0], &m_mat[0][0]);
			int order_h3p = mat_calc_order(MAT_DIM, &result[0][0], false);
			if (order_ph3 != -1 && order_h3p != -1) {
				/*
				printf("P H3 has order %d, and H3 P has order %d, where P:\n", order_ph3, order_h3p);
				mat_print(MAT_DIM, &m_mat[0][0], ", ", "\n");
				printf("\n");
				*/
				total_ph3 += 1;
				if (order_ph2 != -1 && order_h2p != -1) {
					printf("ord(P H2) = %d, ord(H2 P) = %d, ord(P H3) = %d, prd(H3 P) = %d\n",
						order_ph2, order_h2p, order_ph3, order_h3p);
					mat_print(MAT_DIM, &m_mat[0][0], ", ", "\n");
					printf("\n");
					total_ph23 += 1;
				}
			}
		}

		for (u8 x = 0; x < MAT_DIM; x++) {
			m[x] += 1;
			if (m[x] < MAT_DIM) { break; }
			m[x] = 0;
			if (x+1 == MAT_DIM) { done = true; }
		}
	}
	printf("P H2 and H2 P were both finite order for %d distinct permutations P\n\n", total_ph2);
	printf("P H3 and H3 P were both finite order for %d distinct permutations P\n\n", total_ph3);
	printf("All four were finite order for %d distinct permutations P\n", total_ph23);
}

int main() {
	test_perm_fourier_order();
}

