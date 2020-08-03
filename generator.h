#pragma once

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
			bool isnew = true;
			for (PathNode *it = paths; it < next_node; it++) {
				if (memcmp(result, it->result, elem_size) == 0) {
					isnew = false;
				}
			}
			if (!isnew) { continue; }
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

