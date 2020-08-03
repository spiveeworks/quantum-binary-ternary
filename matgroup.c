#include "generator.h"
#include "matrix_quadratic.h"

#define MAT_DIM 2
const size_t MAT_SIZE = (sizeof(num)*MAT_DIM*MAT_DIM);

typedef num Mat[MAT_DIM][MAT_DIM];
// variables of type Mat can be used as MatPtr, but MatPtr makes for a more
// convenient parameter
typedef num (*MatPtr)[MAT_DIM];

void mat_mul_reduce(void *out, void *x, void *y) {
	mat_mul(MAT_DIM, out, x, y);
	mat_reduce(MAT_DIM, out);
}

// not robust to weird values of x, e.g. sqrt(0)/0
bool num_is_zero(num x) {
	return x.re.whole == 0 && x.im.whole == 0;
}

void mat_mul_phase(void *out_v, void *x, void *y) {
	mat_mul(MAT_DIM, out_v, x, y);
	MatPtr out = out_v;
	size_t pi = 0, pj = 0;
	while (num_is_zero(out[pi][pj])) {
		if (pj < MAT_DIM) {
			pj += 1;
		} else if (pi < MAT_DIM) {
			pi += 1;
			pj = 0;
		} else {
			break;
		}
	}
	num phase = out[pi][pj];
	if (!num_is_zero(phase)) {
		// normalise
		phase = num_mul(phase, num_inv_sqrt(num_mod_sq(phase)));
		// invert
		phase = num_inv(phase);
		range(i, MAT_DIM) {
			range(j, MAT_DIM) {
				out[i][j] = num_mul(out[i][j], phase);
			}
		}
	}
	mat_reduce(MAT_DIM, out_v);
}

void shift(MatPtr out) {
	mat_zero(MAT_DIM, &out[0][0]);
	range(i, MAT_DIM) {
		out[(i+1)%MAT_DIM][i].re.whole = 1;
	}
}

void clock(MatPtr out, num *phase, size_t phase_len) {
	range(i, MAT_DIM) {
		out[i][i] = phase[i%phase_len];
	}
}

void fourier(MatPtr out, num *phase, size_t phase_len) {
	range(i, MAT_DIM) {
		range(j, MAT_DIM) {
			num coord = phase[i*j%phase_len];
			coord.re.qu *= MAT_DIM;
			coord.im.qu *= MAT_DIM;
			coord.de *= MAT_DIM;
			out[i][j] = num_reduce(coord);
		}
	}
}

void cliff_diag(MatPtr out, num *phase, size_t phase_len) {
	range(i, MAT_DIM) {
		out[i][i] = phase[i*(MAT_DIM-i)%phase_len];
	}
}

void find_group() {
#define MAX_ROOTS 7
	num principle_roots[MAX_ROOTS] = {
		{{1, 1}, {0, 1}, 1}, // arbitrary
		{{1, 1}, {0, 1}, 1}, // 1
		{{-1, 1}, {0, 1}, 1}, // -1
		{{-1, 1}, {1, 3}, 2}, // (-1 + sqrt(3)i)/2
		{{0, 1}, {1, 1}, 1}, // i
		{{1, 1}, {0, 1}, 1}, // not representable
		{{1, 1}, {1, 3}, 2}, // (1 + sqrt(3)i)/2
	};
	num roots[MAX_ROOTS][MAX_ROOTS];
	range(i, MAX_ROOTS) {
		num curr = {{1,1},{0,1},1};
		range(j, MAX_ROOTS) {
			roots[i][j] = curr;
			curr = num_reduce(num_mul(curr, principle_roots[i]));
		}
	}
#define GATES_MAX 4
	Mat gates[GATES_MAX];
	char *gate_names[GATES_MAX];
	int gates_count = 0;
	range(i, GATES_MAX) {
		mat_zero(MAT_DIM, &gates[i][0][0]);
	}
	gate_names[gates_count] = "X2";
	shift(gates[gates_count]);
	gates_count+=1;

	gate_names[gates_count] = "Z2";
	clock(gates[gates_count], roots[2], 2);
	gates_count+=1;

	gate_names[gates_count] = "H2";
	fourier(gates[gates_count], roots[2], 2);
	gates_count+=1;

	gate_names[gates_count] = "D2";
	cliff_diag(gates[gates_count], roots[4], 4);
	gates_count+=1;

	void *gen[gates_count];
	range(i, gates_count) {
		gen[i] = gates[i];
	}

	PathList paths = gen_paths(gates_count, gen, MAT_SIZE, mat_mul_phase);
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
}

int main() {
	find_group();
}
