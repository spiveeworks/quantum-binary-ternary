#include "generator.h"
#include "matrix_quadratic.h"

#define MAT_DIM_A 2
#define MAT_DIM_B 3
#define MAT_DIM (MAT_DIM_A * MAT_DIM_B)
const size_t MAT_SIZE = (sizeof(num)*MAT_DIM*MAT_DIM);

void mat_mul_reduce(void *out, void *x, void *y) {
	mat_mul(MAT_DIM, out, x, y);
	mat_reduce(MAT_DIM, out);
}

// not sure what we should do for 0/0
bool num_is_zero(num x) {
	range(i, 2) {
		range(j, 2) {
			range(k, 2) {
				if (x.c[i][j][k] != 0) {
					return false;
				}
			}
		}
	}
	return true;
}

void mat_mul_phase(void *out_v, void *x, void *y) {
	mat_mul(MAT_DIM, out_v, x, y);
	// think of this as a pointer to num[MAT_DIM][MAT_DIM]
	num (*out)[MAT_DIM] = out_v;
	size_t pi = 0, pj = 0;
	mat_reduce(MAT_DIM, out_v);
	while (num_is_zero(out[pi][pj])) {
		if (pi < MAT_DIM) {
			pi += 1;
		} else if (pj < MAT_DIM) {
			pj += 1;
			pi = 0;
		} else {
			break;
		}
	}
	num phase = out[pi][pj];
	if (!num_is_zero(phase)) {
		// normalise
		phase = num_mul(phase, num_inv_sqrt(num_mod_sq(phase)));
		// invert
		phase = num_reduce(num_inv(phase));
		range(i, MAT_DIM) {
			range(j, MAT_DIM) {
				out[i][j] = num_mul(out[i][j], phase);
			}
		}
		mat_reduce(MAT_DIM, out_v);
	}
}

num omega12() {
	num out = {};
	out.c[0][0][1] = 1;
	out.c[1][0][0] = 1;
	out.de = 2;
	return out;
}

enum {
	GATE_X,
	GATE_Z,
	GATE_H,
	GATE_D,
	GATE_COUNT,
};

void clifford_gen(int n, num *out) {
	if (12 % n != 0) {
		printf("Tried to generate Clifford matrices for n = %d\n", n);
		exit(1);
	}
	num roots_12[12];
	num roots_2n[2*n];
	num roots[n];
	{
		num curr = num_from_int(1);
		num w = omega12();
		range(i, 12) {
			roots_12[i] = curr;
			curr = num_reduce(num_mul(curr, w));
		}
		range(i, 2*n) {
			roots_2n[i] = roots_12[i*6/n];
		}
		range(i, n) {
			roots[i] = roots_12[i*12/n];
		}
	}
	mat_shift(n, out);
	mat_clock(n, out + n*n, roots, n);
	mat_fourier(n, out + 2*n*n, roots, n);
	mat_cliff_diag(n, out + 3*n*n, roots_2n, 2*n);
}

void find_group() {
	num ident_a[MAT_DIM_A][MAT_DIM_A];
	mat_ident(MAT_DIM_A, &ident_a[0][0]);
	num gates_a[GATE_COUNT][MAT_DIM_A][MAT_DIM_A];
	clifford_gen(MAT_DIM_A, &gates_a[0][0][0]);

	num ident_b[MAT_DIM_B][MAT_DIM_B];
	mat_ident(MAT_DIM_B, &ident_b[0][0]);
	num gates_b[GATE_COUNT][MAT_DIM_B][MAT_DIM_B];
	clifford_gen(MAT_DIM_B, &gates_b[0][0][0]);

	char gate_base_names[GATE_COUNT][2] = {"X", "Z", "H", "D"};
	char gate_names[2*GATE_COUNT][8];
	num gates[2*GATE_COUNT][MAT_DIM][MAT_DIM];
	range(i, GATE_COUNT) {
		mat_kronecker(MAT_DIM_A, MAT_DIM_B,
				&gates[i][0][0], &gates_a[i][0][0], &ident_b[0][0]);
		mat_reduce(MAT_DIM, &gates[i][0][0]);
		sprintf(gate_names[i], "%s%d", gate_base_names[i], MAT_DIM_A);

		size_t j = i + GATE_COUNT;
		mat_kronecker(MAT_DIM_A, MAT_DIM_B,
				&gates[j][0][0], &ident_a[0][0], &gates_b[i][0][0]);
		mat_reduce(MAT_DIM, &gates[j][0][0]);
		sprintf(gate_names[j], "%s%d", gate_base_names[i], MAT_DIM_B);
	}

	void *gen[GATE_COUNT];
	range(i, 2*GATE_COUNT) {
		gen[i] = gates[i];
	/* uncomment if debugging gen_paths
		printf("%s:\n", gate_names[i]);
		mat_print(MAT_DIM, (num*)gates[i], ", ", "\n");
		printf("\n");
	*/
	}

	int gate_count = MAT_DIM_B == 1 ? GATE_COUNT : 2*GATE_COUNT;
	PathList paths = gen_paths(gate_count, gen, MAT_SIZE, mat_mul_phase);
	range(i, paths.path_count) {
		PathNode *curr = &paths.paths[i];
		while (curr) {
			printf(" %s", gate_names[curr->last.i]);
			curr = curr->pred;
		}
		printf(":\n");
		mat_print(MAT_DIM, paths.paths[i].result, ", ", "\n");
		printf("\n");
	}
	printf("Total of %lu gates generated\n", paths.path_count);
}

int main() {
	find_group();
}
