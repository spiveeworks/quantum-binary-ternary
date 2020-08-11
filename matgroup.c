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

void mat_mul_phase(void *out_v, void *x, void *y) {
	mat_mul(MAT_DIM, out_v, x, y);
	// think of this as a pointer to num[MAT_DIM][MAT_DIM]
	num (*out)[MAT_DIM] = out_v;
	size_t pi = 0, pj = 0;
	mat_reduce(MAT_DIM, out_v);
	while (num_eq(out[pi][pj], num_from_int(0))) {
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
	if (!num_eq(phase, num_from_int(0))) {
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

void mat_print_default(void *x) {
	mat_print(MAT_DIM, x, ", ", "\n");
}

bool mat_is_diag(void *x_v) {
	num (*x)[MAT_DIM] = x_v;
	range(i, MAT_DIM) {
		range(j, MAT_DIM) {
			if (i != j && !num_eq(x[i][j], num_from_int(0))) {
				return false;
			}
		}
	}
	return true;
}

bool mat_is_ident(num (*x)[MAT_DIM]) {
	range(i, MAT_DIM) {
		range(j, MAT_DIM) {
			if (i != j && !num_eq(x[i][j], num_from_int(0))) {
				return false;
			}
			if (i == j && !num_eq(x[i][j], num_from_int(1))) {
				return false;
			}
		}
	}
	return true;
}

const int max_order = 28;

int mat_calc_order(num (*x)[MAT_DIM]) {
	if (mat_is_ident(x)) { return 0; }
	num y[MAT_DIM][MAT_DIM];
	range(i, MAT_DIM) {
		range(j, MAT_DIM) {
			y[i][j] = x[i][j];
		}
	}
	int order = 1;
	while (!mat_is_ident(y)) {
		num z[MAT_DIM][MAT_DIM];
		mat_mul(MAT_DIM, &z[0][0], &x[0][0], &y[0][0]);
		range(i, MAT_DIM) {
			range(j, MAT_DIM) {
				y[i][j] = num_reduce(z[i][j]);
			}
		}
		order += 1;
		if (order >= max_order) {
			printf("Matrix has order >= %d:\n", max_order);
			mat_print(MAT_DIM, &x[0][0], ", ", "\n");
			printf("\nraised to power %d:\n", order);
			mat_print(MAT_DIM, &y[0][0], ", ", "\n");
			return order;
		}
	}
	return order;
}

bool mat_print_order_ret_false(void *x) {
	int ord = mat_calc_order(x);
	if (ord >= max_order) {
		printf("Order = %d\n", ord);
		return true;
	}
	return false;
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
	mat_fourier(n, out + 0*n*n, roots, n);
	mat_cliff_diag(n, out + 1*n*n, roots_2n, 2*n);
	mat_clock(n, out + 2*n*n, roots, n);
	mat_shift(n, out + 3*n*n);
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

	char gate_base_names[GATE_COUNT][2] = {"H", "D", "Z", "X"};
	char gate_name_buffs[2*GATE_COUNT+2][8];
	char *gate_names[2*GATE_COUNT+2];
	range(i, 2*GATE_COUNT+2) {
		gate_names[i] = gate_name_buffs[i];
	}
	num gates[2*GATE_COUNT+2][MAT_DIM][MAT_DIM];
	range(i, GATE_COUNT) {
		mat_kronecker(MAT_DIM_A, MAT_DIM_B,
				&gates[2*i][0][0], &gates_a[i][0][0], &ident_b[0][0]);
		mat_reduce(MAT_DIM, &gates[2*i][0][0]);
		sprintf(gate_names[2*i], "%s%d", gate_base_names[i], MAT_DIM_A);

		mat_kronecker(MAT_DIM_A, MAT_DIM_B,
				&gates[2*i+1][0][0], &ident_a[0][0], &gates_b[i][0][0]);
		mat_reduce(MAT_DIM, &gates[2*i+1][0][0]);
		sprintf(gate_names[2*i+1], "%s%d", gate_base_names[i], MAT_DIM_B);
	}
	int gate_count;
	if (MAT_DIM_B == 1) {
		gate_count = GATE_COUNT;
	} else {
		gate_count = 2*GATE_COUNT;
		mat_sum_right(MAT_DIM_A, MAT_DIM_B, &gates[gate_count][0][0]);
		sprintf(gate_names[gate_count], "C%d(X%d)", MAT_DIM_A, MAT_DIM_B);
		gate_count += 1;
		mat_sum_left(MAT_DIM_A, MAT_DIM_B, &gates[gate_count][0][0]);
		sprintf(gate_names[gate_count], "C%d(X%d)", MAT_DIM_B, MAT_DIM_A);
		gate_count += 1;
	}

	void *gen[2*GATE_COUNT+2];
	range(i, gate_count) {
		gen[i] = gates[i];
		printf("%s:\n", gate_names[i]);
		mat_print(MAT_DIM, (num*)gates[i], ", ", "\n");
		printf("\n");
	}

	PathList paths = gen_paths(
		gate_count, gen, gate_names, MAT_SIZE,
		mat_mul_reduce, mat_print_default, mat_print_order_ret_false, SEARCH_BREADTH_FIRST);
	printf("Total of %lu gates generated\n", paths.path_count);
}

int main() {
	find_group();
}
