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

bool mat_print_order_ret_true(void *x) {
	int ord = mat_calc_order(MAT_DIM, x, true);
	if (ord == -1) {
		return true;
		//exit(0);
	} else {
		//printf("Order = %d\n", ord);
		return false;
	}
	return true;
}

void find_group() {
	char *gate_names[10] = {};
	num gates[10][MAT_DIM][MAT_DIM];
	int gate_count = 0;

	gate_names[gate_count] = "X2";
	initialize(MAT_DIM_A, MAT_DIM_B, &gates[gate_count][0][0], GATE_SHIFT, GATE_IDENT);
	gate_count += 1;

	gate_names[gate_count] = "X3";
	initialize(MAT_DIM_A, MAT_DIM_B, &gates[gate_count][0][0], GATE_IDENT, GATE_SHIFT);
	gate_count += 1;

	gate_names[gate_count] = "Z2";
	initialize(MAT_DIM_A, MAT_DIM_B, &gates[gate_count][0][0], GATE_CLOCK, GATE_IDENT);
	gate_count += 1;

	gate_names[gate_count] = "Z3";
	initialize(MAT_DIM_A, MAT_DIM_B, &gates[gate_count][0][0], GATE_IDENT, GATE_CLOCK);
	gate_count += 1;

	gate_names[gate_count] = "H2";
	initialize(MAT_DIM_A, MAT_DIM_B, &gates[gate_count][0][0], GATE_FOURIER, GATE_IDENT);
	gate_count += 1;

	gate_names[gate_count] = "H3";
	initialize(MAT_DIM_A, MAT_DIM_B, &gates[gate_count][0][0], GATE_IDENT, GATE_FOURIER);
	gate_count += 1;

	gate_names[gate_count] = "D2";
	initialize(MAT_DIM_A, MAT_DIM_B, &gates[gate_count][0][0], GATE_CLIFF_DIAG, GATE_IDENT);
	gate_count += 1;

	gate_names[gate_count] = "D3";
	initialize(MAT_DIM_A, MAT_DIM_B, &gates[gate_count][0][0], GATE_IDENT, GATE_CLIFF_DIAG);
	gate_count += 1;

	gate_names[gate_count] = "C2X3";
	initialize(MAT_DIM_A, MAT_DIM_B, &gates[gate_count][0][0], GATE_CONTROL_1, GATE_SHIFT);
	//gate_count += 1;

	gate_names[gate_count] = "C3X2";
	initialize(MAT_DIM_A, MAT_DIM_B, &gates[gate_count][0][0], GATE_SHIFT, GATE_CONTROL_1);
	//gate_count += 1;

	gate_names[gate_count] = "C2S3";
	initialize(MAT_DIM_A, MAT_DIM_B, &gates[gate_count][0][0], GATE_CONTROL_1, GATE_TRANSPOSE);
	gate_count += 1;

	void *gen[gate_count];
	range(i, gate_count) {
		gen[i] = gates[i];
		printf("%s:\n", gate_names[i]);
		mat_print(MAT_DIM, (num*)gates[i], ", ", "\n");
		printf("\n");
	}

	PathList paths = gen_paths(
		gate_count, gen, gate_names, MAT_SIZE,
		mat_mul_phase, mat_print_default, mat_print_order_ret_true, SEARCH_BREADTH_FIRST
	);
	printf("Total of %lu gates generated\n", paths.path_count);
}

int main() {
	find_group();
}


