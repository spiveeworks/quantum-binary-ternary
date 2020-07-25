#include<stdbool.h>
#include<stdio.h>
#include<stdlib.h>

typedef unsigned char u8;
typedef u8 *Map;
typedef u8 const *ConstMap;

typedef struct { size_t i; } MapIndex;

MapIndex map_index(ConstMap m, u8 num) {
	size_t mul = 1;
	MapIndex result = {};
	for (size_t i = 0; i < num; i++) {
		result.i += m[i]*mul;
		mul *= num;
	}
	return result;
}

void from_map_index(Map m, u8 num, MapIndex mi) {
	for (size_t i = 0; i < num; i++) {
		m[i] = mi.i%num;
		mi.i /= num;
	}
}

typedef struct { size_t i; } PathIndex;
typedef struct { size_t n; } PathLength;

#define MAX_PATH_STEPS 100
typedef struct {
	PathLength len;
	MapIndex result;
	PathIndex pred;
	// we're using the convention x . y = x after y
	// so by last we mean leftmost
	MapIndex last;
} PathNode;

typedef struct {
	size_t path_count;
	PathNode *paths;
} PathList;


PathList gen_paths(size_t gen_len, MapIndex *gen_mi, u8 num) {
	u8 gen[gen_len][num];
	for (size_t i = 0; i < gen_len; i++) {
		from_map_index(gen[i], num, gen_mi[i]);
	}
#define MAX_PATH_COUNT 1000
	static PathNode paths[MAX_PATH_COUNT];

	u8 ident[num];
	for (u8 i = 0; i < num; i++) { ident[i] = i; }
	paths[0].len.n = 0;
	paths[0].result = map_index(ident, num);
	paths[0].pred.i = ~0U;
	paths[0].last.i = ~0U;
	PathList out = {1,paths};

#define DEQ_CAP 1024
	PathIndex deq[DEQ_CAP];
	size_t deq_start = 0;
	deq[0].i = 0;
	size_t deq_end = 1;

	while (deq_end > deq_start) {
		PathIndex curr_pi = deq[deq_start % DEQ_CAP];
		deq_start++;
		PathLength curr_len = out.paths[curr_pi.i].len;
		MapIndex curr_mi = out.paths[curr_pi.i].result;
		u8 curr[num];
		from_map_index(curr, num, curr_mi);
		for (size_t i = 0; i < gen_len; i++) {
			// result = gen[i] . curr
			// i.e. gen[i] after curr
			u8 result[num];
			for (u8 x = 0; x < num; x++) {
				result[x] = gen[i][curr[x]];
			}
			MapIndex result_mi = map_index(result, num);
			bool new = true;
			for (PathIndex pi = {0}; new && pi.i < out.path_count; pi.i++) {
				if (out.paths[pi.i].result.i == result_mi.i) {
					new = false;
				}
			}
			if (!new) { continue; }
			PathIndex next_pi = {out.path_count++};
			out.paths[next_pi.i].len.n = curr_len.n + 1;
			out.paths[next_pi.i].result = result_mi;
			out.paths[next_pi.i].pred = curr_pi;
			out.paths[next_pi.i].last = gen_mi[i];
			if (deq_end - deq_start >= DEQ_CAP) {
				printf("Queue is full!\n");
				exit(1);
			}
			deq[deq_end % DEQ_CAP] = next_pi;
			deq_end++;
		}
	}

	return out;
}

void print_cycle_decomposition(ConstMap m, u8 num) {
	bool printed[num];
	for (u8 i = 0; i < num; i++) {
		printed[i] = false;
	}
	bool err = false;
	for (u8 i = 0; i < num; i++) {
		if (m[i] == i || printed[i]) { continue; }
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
	if (err) {
		printf("WARNING non-cycle was generated\n");
	}
}

int main() {
#define NUM 6
	u8 cx3[NUM] = {0,1,2,4,5,3};
	u8 cx2[NUM] = {0,4,2,3,1,5};
	u8 x3[NUM] = {1,2,0,4,5,3};
	u8 x2[NUM] = {3,4,5,0,1,2};
#define GEN_LEN 4
	ConstMap gen[GEN_LEN] = {cx3, cx2, x3, x2};
	MapIndex gen_mi[GEN_LEN];
	for (size_t i = 0; i < GEN_LEN; i++) {
		gen_mi[i] = map_index(gen[i], NUM);
	}
	PathList paths = gen_paths(GEN_LEN, gen_mi, NUM);
	for (size_t i = 1; i < paths.path_count; i++) {
		u8 buff[NUM];
		from_map_index(buff, NUM, paths.paths[i].result);
		print_cycle_decomposition(buff, NUM);
		printf(" =");
		PathIndex pi = {i};
		while (pi.i != 0) {
			printf(" ");
			from_map_index(buff, NUM, paths.paths[pi.i].last);
			print_cycle_decomposition(buff, NUM);
			pi = paths.paths[pi.i].pred;
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

