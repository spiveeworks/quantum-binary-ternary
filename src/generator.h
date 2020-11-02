#pragma once

#include<stdbool.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#define MAX_PATH_COUNT 1000000
#include "hashmap.h"

typedef unsigned char u8;

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

enum SearchStrategy {
	SEARCH_BREADTH_FIRST,
	SEARCH_DEPTH_FIRST,
};

PathList gen_paths(
	size_t gen_len, void** gen, char **names, size_t elem_size,
	void (*compose)(void*,void*,void*), void (*print_elem)(void*),
	bool (*print_cond)(void*), enum SearchStrategy strategy
) {
	static PathNode paths[MAX_PATH_COUNT];
	PathNode *next_node = paths;
	size_t path_count = 0;

	struct PathNode *curr = paths;
	GenIndex next_generator = {0};

	for (size_t i = 0; i < gen_len; i++) {
		next_node->len.n = 1;
		next_node->pred = NULL;
		next_node->last.i = i;
		next_node->result = gen[i];
		hash_insert(next_node->result, elem_size, next_node);
		next_node++;
		path_count++;
	}

	while (true) {
		bool done = false;
		while (!done && next_generator.i == gen_len) {
			switch (strategy) {
				case SEARCH_BREADTH_FIRST:
					next_generator.i = 0;
					curr += 1;
					done = curr == next_node;
					break;
				case SEARCH_DEPTH_FIRST:
					if (curr->pred == NULL) {
						done = true;
					} else {
						next_generator.i = curr->last.i + 1;
						curr = curr->pred;
					}
					break;
			}
		}
		if (done) { break; }

		u8 result[elem_size];
		compose(result, gen[next_generator.i], curr->result);
		if (hash_lookup(result, elem_size) != NULL) {
			next_generator.i += 1;
			continue;
		}
		if (path_count >= MAX_PATH_COUNT) {
			printf("Path List is full!\n");
			exit(1);
		}
		next_node->len.n = curr->len.n + 1;
		next_node->pred = curr;
		next_node->last = next_generator;
		next_node->result = malloc(elem_size);
		memcpy(next_node->result, result, elem_size);
		hash_insert(next_node->result, elem_size, next_node);
		if (print_cond == NULL || print_cond(next_node->result)){
			PathNode *printing = next_node;
			while (printing) {
				printf(" %s", names[printing->last.i]);
				printing = printing->pred;
			}
			printf(":\n");
			print_elem(next_node->result);
			printf("\n");
		}
		switch (strategy) {
			case SEARCH_BREADTH_FIRST:
				next_generator.i += 1;
				break;
			case SEARCH_DEPTH_FIRST:
				curr = next_node;
				next_generator.i = 0;
				break;
		}
		next_node++;
		path_count++;
	}

	return (PathList){path_count, paths};
}

