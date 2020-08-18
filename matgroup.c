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

const int max_precision = 100;

int mat_calc_order(num (*x)[MAT_DIM], bool print) {
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
		bool overprecise = false;
		range(i, MAT_DIM) {
			range(j, MAT_DIM) {
				y[i][j] = num_reduce(z[i][j]);
				if (y[i][j].de > max_precision) {
					overprecise = true;
				}
				range(a, 2) {
					range(b, 2) {
						range(c, 2) {
							int val = y[i][j].c[a][b][c];
							if (val > max_precision || -val > max_precision) {
								overprecise = true;
							}
						}
					}
				}
			}
		}
		order += 1;
		if (overprecise) {
			if (print) {
				mat_print(MAT_DIM, &x[0][0], ", ", "\n");
				printf("\nraised to power %d:\n", order);
				mat_print(MAT_DIM, &y[0][0], ", ", "\n");
			}
			return -1;
		}
	}
	return order;
}

bool mat_print_order_ret_true(void *x) {
	int ord = mat_calc_order(x, true);
	if (ord == -1) {
		return true;
		//exit(0);
	} else {
		//printf("Order = %d\n", ord);
		return false;
	}
	return true;
}

#define SMALLEST_ROOT 24
bool roots_initialized = false;
num roots[SMALLEST_ROOT];

num root_of_unity(int r, int d) {
	if (SMALLEST_ROOT * r % d != 0) {
		printf("Cannot calculate root of unity, %d does not divide %d\n",
			d, SMALLEST_ROOT * r);
		exit(1);
	}

	if (!roots_initialized) {
		num x = num_from_int(1);
		num y = num_from_int(0);
		y.c[0][1][0] = 1;
		y.c[0][1][1] = 1;
		y.c[1][1][0] = -1;
		y.c[1][1][1] = 1;
		y.de = 4;

		range(i, SMALLEST_ROOT) {
			if (i > 0 && num_eq(x, num_from_int(1))) {
				printf("Error while initializing roots:\n");
				num_print(y);
				printf("raised to %lu is already 1\n", i);
				exit(1);
			}
			roots[i] = x;
			x = num_reduce(num_mul(x, y));
		}
		if (!num_eq(x, num_from_int(1))) {
			printf("Error while initializing roots:\n");
			num_print(y);
			printf("raised to %d is only ", SMALLEST_ROOT);
			num_print(x);
			printf("\n");
			exit(1);
		}

		roots_initialized = true;
	}
	return roots[SMALLEST_ROOT * r / d % SMALLEST_ROOT];
}

enum Gate {
	GATE_IDENT,
	GATE_SHIFT,
	GATE_CLOCK,
	GATE_FOURIER,
	GATE_CLIFF_DIAG,
	GATE_TRANSPOSE,
	GATE_CONTROL_0,
	GATE_CONTROL_1,
	GATE_CONTROL_N,
};

void initialize_single(int n, num *out, enum Gate gate) {
	switch(gate) {
		case GATE_IDENT:
		{
			mat_ident(n, out);
			break;
		}
		case GATE_SHIFT:
		{
			mat_shift(n, out);
			break;
		}
		case GATE_CLOCK:
		{
			num phase[n];
			range(i, n) { phase[i] = root_of_unity(i, n); }
			mat_clock(n, out, phase, n);
			break;
		}
		case GATE_FOURIER:
		{
			num phase[n];
			range(i, n) { phase[i] = root_of_unity(i, n); }
			mat_fourier(n, out, phase, n);
			break;
		}
		case GATE_CLIFF_DIAG:
		{
			num phase[2*n];
			range(i, 2*n) { phase[i] = root_of_unity(i, 2*n); }
			mat_cliff_diag(n, out, phase, 2*n);
			break;
		}
		case GATE_TRANSPOSE:
		{
			unsigned char p[n];
			range(i, n) {
				p[i] = i;
			}
			p[0] = 1;
			p[1] = 0;
			mat_perm(n, out, p);
			break;
		}
		case GATE_CONTROL_0:
		case GATE_CONTROL_1:
		case GATE_CONTROL_N:
		{
			printf("Tried to initialize a single gate with a control specifier\n");
			exit(1);
		}
		default:
		{
			printf("Unknown gate specifier %d\n", gate);
			exit(1);
		}
	}
}

void initialize(num (*out)[MAT_DIM], enum Gate gate_a, enum Gate gate_b) {
	num a[MAT_DIM_A][MAT_DIM_A];
	if (gate_a < GATE_CONTROL_0) {
		initialize_single(MAT_DIM_A, &a[0][0], gate_a);
	}
	num b[MAT_DIM_B][MAT_DIM_B];
	if (gate_b < GATE_CONTROL_0) {
		initialize_single(MAT_DIM_B, &b[0][0], gate_b);
	}
	if (gate_a < GATE_CONTROL_0 && gate_b < GATE_CONTROL_0) {
		mat_kronecker(MAT_DIM_A, MAT_DIM_B, &out[0][0], &a[0][0], &b[0][0]);
	} else if (gate_b < GATE_CONTROL_0) {
		int control;
		switch(gate_a) {
			case GATE_CONTROL_0:
				control = 0;
				break;
			case GATE_CONTROL_1:
				control = 1;
				break;
			case GATE_CONTROL_N:
				control = MAT_DIM_A - 1;
				break;
			default:
				printf("Invalid control specifier %d\n", gate_a);
				exit(1);
		}
		mat_control_target(MAT_DIM_A, MAT_DIM_B, &out[0][0], control, &b[0][0]);
	} else if (gate_a < GATE_CONTROL_0) {
		int control;
		switch(gate_b) {
			case GATE_CONTROL_0:
				control = 0;
				break;
			case GATE_CONTROL_1:
				control = 1;
				break;
			case GATE_CONTROL_N:
				control = MAT_DIM_B - 1;
				break;
			default:
				printf("Invalid control specifier %d\n", gate_b);
				exit(1);
		}
		mat_target_control(MAT_DIM_A, MAT_DIM_B, &out[0][0], &a[0][0], control);
	} else {
		printf("Tried to initialize a gate with only control specifiers\n");
		exit(1);
	}
	mat_reduce(MAT_DIM, &out[0][0]);
}

void find_group() {
	char *gate_names[10] = {};
	num gates[10][MAT_DIM][MAT_DIM];
	int gate_count = 0;

	gate_names[gate_count] = "X2";
	initialize(gates[gate_count], GATE_SHIFT, GATE_IDENT);
	gate_count += 1;

	gate_names[gate_count] = "X3";
	initialize(gates[gate_count], GATE_IDENT, GATE_SHIFT);
	gate_count += 1;

	gate_names[gate_count] = "Z2";
	initialize(gates[gate_count], GATE_CLOCK, GATE_IDENT);
	gate_count += 1;

	gate_names[gate_count] = "Z3";
	initialize(gates[gate_count], GATE_IDENT, GATE_CLOCK);
	gate_count += 1;

	gate_names[gate_count] = "H2";
	initialize(gates[gate_count], GATE_FOURIER, GATE_IDENT);
	gate_count += 1;

	gate_names[gate_count] = "H3";
	initialize(gates[gate_count], GATE_IDENT, GATE_FOURIER);
	gate_count += 1;

	gate_names[gate_count] = "D2";
	initialize(gates[gate_count], GATE_CLIFF_DIAG, GATE_IDENT);
	gate_count += 1;

	gate_names[gate_count] = "D3";
	initialize(gates[gate_count], GATE_IDENT, GATE_CLIFF_DIAG);
	gate_count += 1;

	gate_names[gate_count] = "C2X3";
	initialize(gates[gate_count], GATE_CONTROL_1, GATE_SHIFT);
	//gate_count += 1;

	gate_names[gate_count] = "C3X2";
	initialize(gates[gate_count], GATE_SHIFT, GATE_CONTROL_1);
	//gate_count += 1;

	gate_names[gate_count] = "C2S3";
	initialize(gates[gate_count], GATE_CONTROL_1, GATE_TRANSPOSE);
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

void test_perm_fourier_order() {
	num h2[MAT_DIM][MAT_DIM];
	initialize(h2, GATE_FOURIER, GATE_IDENT);

	num h3[MAT_DIM][MAT_DIM];
	initialize(h3, GATE_IDENT, GATE_FOURIER);

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
			int order_ph2 = mat_calc_order(result, false);
			mat_mul(MAT_DIM, &result[0][0], &h2[0][0], &m_mat[0][0]);
			int order_h2p = mat_calc_order(result, false);
			if (order_ph2 != -1 && order_h2p != -1) {
				/*
				printf("P H2 has order %d, and H2 P has order %d, where P:\n", order_ph2, order_h2p);
				mat_print(MAT_DIM, &m_mat[0][0], ", ", "\n");
				printf("\n");
				*/
				total_ph2 += 1;
			}

			mat_mul(MAT_DIM, &result[0][0], &m_mat[0][0], &h3[0][0]);
			int order_ph3 = mat_calc_order(result, false);
			mat_mul(MAT_DIM, &result[0][0], &h3[0][0], &m_mat[0][0]);
			int order_h3p = mat_calc_order(result, false);
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
	//find_group();
	test_perm_fourier_order();
}
