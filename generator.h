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
	size_t gen_len, void** gen, char **names, size_t elem_size,
	void (*compose)(void*,void*,void*), void (*print_elem)(void*)
) {
#define MAX_PATH_COUNT 10000
	static PathNode paths[MAX_PATH_COUNT];
	PathNode *next_node = paths;
	size_t path_count = 0;

	size_t searched_count = 0;

	for (size_t i = 0; i < gen_len; i++) {
		next_node->len.n = 1;
		next_node->pred = NULL;
		next_node->last.i = i;
		next_node->result = gen[i];
		next_node++;
		path_count++;
	}

	while (searched_count < path_count) {
		PathNode* curr = &paths[searched_count];
		searched_count++;
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
			if (path_count >= MAX_PATH_COUNT) {
				printf("Path List is full!\n");
				exit(1);
			}
			next_node->len.n = curr->len.n + 1;
			next_node->pred = curr;
			next_node->last.i = i;
			next_node->result = malloc(elem_size);
			memcpy(next_node->result, result, elem_size);
			{
				PathNode *curr = next_node;
				while (curr) {
					printf(" %s", names[curr->last.i]);
					curr = curr->pred;
				}
				printf(":\n");
				print_elem(next_node->result);
				printf("\n");
			}
			next_node++;
			path_count++;
		}
	}

	return (PathList){path_count, paths};
}

