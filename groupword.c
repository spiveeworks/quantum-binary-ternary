#include<stdbool.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

typedef unsigned char u8;
typedef u8 *Map;
typedef u8 const *ConstMap;

typedef struct { size_t i; } GenIndex;
typedef struct { size_t n; } PathLength;

typedef struct PathNode {
	PathLength len;
	struct PathNode *pred;
	// we're using the convention x . y = x after y
	// so by last we mean leftmost
	GenIndex last;
	void *result;
} PathNode;

typedef struct {
	size_t path_count;
	PathNode *paths;
} PathList;


PathList gen_paths(
	size_t gen_len, void** gen, size_t elem_size,
	void (*compose)(void*,void*,void*)
) {
#define MAX_PATH_COUNT 1000
	static PathNode paths[MAX_PATH_COUNT];
	PathNode *next_node = paths;
	size_t path_count = 0;

#define DEQ_CAP 1024
	PathNode *deq[DEQ_CAP];
	size_t deq_start = 0;

	for (size_t i = 0; i < gen_len; i++) {
		next_node->len.n = 1;
		next_node->pred = NULL;
		next_node->last.i = i;
		next_node->result = gen[i];
		deq[i] = next_node;
		next_node++;
		path_count++;
	}

	size_t deq_end = gen_len;

	while (deq_end > deq_start) {
		PathNode* curr = deq[deq_start % DEQ_CAP];
		deq_start++;
		for (size_t i = 0; i < gen_len; i++) {
			u8 result[elem_size];
			compose(result, gen[i], curr->result);
			bool new = true;
			for (PathNode *it = paths; it < next_node; it++) {
				if (memcmp(result, it->result, elem_size) == 0) {
					new = false;
				}
			}
			if (!new) { continue; }
			next_node->len.n = curr->len.n + 1;
			next_node->pred = curr;
			next_node->last.i = i;
			next_node->result = malloc(elem_size);
			memcpy(next_node->result, result, elem_size);
			if (deq_end - deq_start >= DEQ_CAP) {
				printf("Queue is full!\n");
				exit(1);
			}
			deq[deq_end % DEQ_CAP] = next_node;
			deq_end++;
			next_node++;
			path_count++;
		}
	}

	return (PathList){path_count, paths};
}

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

void compose_permutations(u8 *out, const u8 *x, const u8 *y) {
	for (size_t i = 0; i < NUM; i++) {
		out[i] = x[y[i]];
	}
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
	PathList paths = gen_paths(gen_len, gen, NUM, compose_permutations);
	for (size_t i = 0; i < paths.path_count; i++) {
		print_cycle_decomposition(paths.paths[i].result, NUM);
		printf(" =");
		PathNode *curr = &paths.paths[i];
		while (curr) {
			printf(" %s", gen_names[curr->last.i]);
			curr = curr->pred;
		}
		printf("\n");
	}
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

