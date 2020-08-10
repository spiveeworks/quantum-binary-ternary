#include "generator.h"

typedef u8 *Map;
typedef u8 const *ConstMap;

void print_cycle_decomposition(ConstMap m, u8 num) {
	bool printed[num];
	for (u8 i = 0; i < num; i++) {
		printed[i] = false;
	}
	bool err = false;
	bool any_printed = false;
	for (u8 i = 0; i < num; i++) {
		if (m[i] == i || printed[i]) { continue; }
		any_printed = true;
		printf("(%u", i);
		printed[i] = true;
		u8 d = i;
		while (!printed[m[d]]) {
			d = m[d];
			if (num > 10) {
				printf(",");
			}
			printf("%u", d);
			printed[d] = true;
		}
		if (m[d] != i) {
			err = true;
			if (num > 10) {
				printf(",");
			}
			printf("%u", m[d]);
		}
		printf(")");
	}
	if (!any_printed) {
		printf("()");
	}
	if (err) {
		printf("WARNING non-cycle was generated\n");
	}
}

#define NUM 6

void compose_permutations(void *out_v, void *x_v, void *y_v) {
	u8 *out = out_v;
	u8 *x = x_v;
	u8 *y = y_v;
	for (size_t i = 0; i < NUM; i++) {
		out[i] = x[y[i]];
	}
}

void print_permutation(void *x) {
	print_cycle_decomposition(x, NUM);
}

void findgen() {
#define MAPS_LEN 8
	struct {
		char *name;
		u8 map[6];
		bool include;
	} maps[MAPS_LEN] = {
		{"cx3", {0,1,2,4,5,3}, false},
		{"cx2", {0,4,2,3,1,5}, false},
		{"x3", {1,2,0,4,5,3}, true},
		{"s3", {1,0,2,4,3,5}, true},
		{"x2", {3,4,5,0,1,2}, true},
		{"x2x3", {4,5,3,1,2,0}, false},
		{"g6", {1,2,3,4,5,0}, false},
		{"cs3", {0,1,2,4,3,5}, false},
	};
	size_t gen_len = 0;
	void* gen[MAPS_LEN];
	char* gen_names[MAPS_LEN];
	for (size_t i = 0; i < MAPS_LEN; i++) {
		if (maps[i].include) {
			gen[gen_len] = maps[i].map;
			gen_names[gen_len] = maps[i].name;
			gen_len++;
		}
	}
	PathList paths = gen_paths(gen_len, gen, gen_names, NUM, compose_permutations, print_permutation, NULL);
	size_t sn = 1;
	for (size_t i = 1; i <= NUM; i++) {
		sn *= i;
	}
	printf("Hit %u out of %u permutations\n", paths.path_count, sn);
	printf("Max word length was %d\n", paths.paths[paths.path_count - 1].len.n);
}

void findcliffordish() {
	u8 shift[NUM][NUM];
	u8 x2x3[NUM] = {4,5,3,1,2,0};
	for (u8 x = 0; x < NUM; x++) {
		shift[0][x] = x;
	}
	for (u8 i = 1; i < NUM; i++) {
		for (u8 x = 0; x < NUM; x++) {
			shift[i][x] = x2x3[shift[i-1][x]];
		}
		printf("x2x3^%u = ", i);
		print_cycle_decomposition(shift[i], NUM);
		printf("\n");
	}
	u8 m[NUM] = {};
	bool done = false;
	while (!done) {
		u8 minv[NUM] = {255, 255, 255, 255, 255, 255};
		for (u8 x = 0; x < NUM; x++) {
			minv[m[x]] = x;
		}
		bool invertible = true;
		for (u8 y = 0; y < NUM; y++) {
			if (minv[y] == 255) {
				invertible = false;
			}
		}
		if (invertible) {
			u8 result[NUM];
			for (u8 x = 0; x < NUM; x++) {
				result[x] = minv[x2x3[m[x]]];
			}
			for (u8 i = 0; i < NUM; i++) {
				bool match = true;
				for (u8 x = 0; x < NUM; x++) {
					if (result[x] != shift[i][x]) {
						match = false;
					}
				}
				if (match) {
					print_cycle_decomposition(minv, NUM);
					printf(" x2x3 ");
					print_cycle_decomposition(m, NUM);
					printf(" = x2x3^%u\n", i);
				}
			}
		}

		for (u8 x = 0; x < NUM; x++) {
			m[x] += 1;
			if (m[x] < NUM) { break; }
			m[x] = 0;
			if (x+1 == NUM) { done = true; }
		}
	}
}

int main() {
	findgen();
	findcliffordish();
}

