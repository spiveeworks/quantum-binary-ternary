#pragma once

#include<stdbool.h>
#include<stdint.h>
#include<string.h>

#ifdef MAX_PATH_COUNT
#define TABLE_LEN MAX_PATH_COUNT
#else
#define TABLE_LEN 100000
#endif

struct HashNode {
    uint64_t hash;
    void *key;
    void *val;
} table[TABLE_LEN] = {};

uint64_t hash_calc(void *val_v, size_t size) {
    char *val = (char*)val_v;
    uint64_t result = 0xFEDCBA987654321;
    for (size_t i = 0; i < size; i++) {
        char x = val[i];
        // just put the bits in different places, nothing special happening
        result = (result << 5U) + (result >> 2U) + (x << 17U) + x;
    }
    // final touch to give higher bits more effect, and to get correct range
    return (result + (result >> 32U)) % TABLE_LEN;
}

// Robin Hood Hashing: swap new_hash with curr_hash if new_hash has been
// probing from a (cyclically) smaller index than curr_hash
#define swap_cond(ind, curr_hash, new_hash) (((new_hash)+TABLE_LEN-(ind)-1)%TABLE_LEN < ((curr_hash)+TABLE_LEN-(ind)-1)%TABLE_LEN)

void* hash_lookup(void *key, size_t key_size) {
    uint64_t hash = hash_calc(key, key_size);
    size_t index = hash % TABLE_LEN;
    
    // probe until either we find the correct key, or we find a node that
    // would have been replaced by Robin Hood Hashing 
    while (true) {
        if (table[index].key == NULL || swap_cond(index, table[index].hash, hash)) {
            return NULL;
        }
        if (memcmp(table[index].key, key, key_size) == 0) {
            return table[index].val;
        }
        index = (index + 1) % TABLE_LEN;
    }
}

void hash_insert(void *key, size_t key_size, void *val) {
    struct HashNode new_node;
    new_node.hash = hash_calc(key, key_size);
    new_node.key = key;
    new_node.val = val;
    size_t index = new_node.hash;
    bool done = false;
    // rearrange nodes according to Robin Hood Hashing
    while (!done) {
        if (table[index].key == NULL) {
            table[index] = new_node;
            done = true;
        } else if (swap_cond(index, table[index].hash, new_node.hash)) {
            struct HashNode swap = table[index];
            table[index] = new_node;
            new_node = swap;
        } else if (memcmp(table[index].key, new_node.key, key_size) == 0) {
            table[index] = new_node;
            done = true;
        }
        index = (index + 1) % TABLE_LEN;
    }
}

