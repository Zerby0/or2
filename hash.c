#include "hash.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>

#define RETURN_FAIL \
	if (CUCKOO_EXIT_ON_FAIL) { \
		fprintf(stderr, "Failed: %s at %s:%d\n", __func__, __FILE__, __LINE__); \
		exit(EXIT_FAILURE); \
		return false; \
	} else { \
		return false; \
	}

// Fast logarithm base 2 calculation
int fast_log2(unsigned long long x) {
    return ((unsigned) (8*sizeof (unsigned long long) - __builtin_clzll((x)) - 1));
}

// Forward declarations
typedef struct CuckooHash CuckooHash;

// Function pointer types for hash functions and comparison
typedef size_t (*HashFunction)(const void* key, size_t table_size);
typedef bool (*CompareFunction)(const void* a, const void* b);

// Helper functions to access elements in the tables
static inline void* get_element(void* table, size_t index, size_t obj_size) {
    return (char*)table + (index * obj_size);
}
static inline void set_element(void* table, size_t index, const void* value, size_t obj_size) {
    memcpy((char*)table + (index * obj_size), value, obj_size);
}

// Calculate stash capacity based on table size
static inline size_t calculate_stash_capacity(size_t table_size) {
    int capacity  = fast_log2(table_size > 0 ? table_size : 1);
    return capacity > 4 ? capacity : 4;
}

// Calculate MAX_LOOPS based on table size
static inline int calculate_max_loops(size_t table_size) {
    int max_loops = 2 * fast_log2(table_size > 0 ? table_size : 1);
    return max_loops > 8 ? max_loops : 8;
}

// Create a new cuckoo hash table
bool cuckoo_create(CuckooHash* hash,
                   size_t initial_size, double max_load_factor, size_t obj_size,
				   HashFunction hash1, HashFunction hash2, CompareFunction compare) {
    // Validate parameters
    if (initial_size < 1 || max_load_factor <= 0 || max_load_factor >= 1 || obj_size == 0
			|| hash1 == NULL || hash2 == NULL || compare == NULL) {
		RETURN_FAIL;
    }
    
    // Initialize tables with actual object size
    hash->table1 = malloc(initial_size * obj_size);
    hash->table2 = malloc(initial_size * obj_size);
    
    // Initialize stash size based on log2 of table size
    hash->stash_capacity = calculate_stash_capacity(initial_size);
    hash->stash = malloc(hash->stash_capacity * obj_size);
    
    // Initialize occupation status arrays
    hash->occupied1 = (bool*)calloc(initial_size, sizeof(bool));
    hash->occupied2 = (bool*)calloc(initial_size, sizeof(bool));
    
    if (hash->table1 == NULL || hash->table2 == NULL || hash->stash == NULL || hash->occupied1 == NULL || hash->occupied2 == NULL) {
        free(hash->table1);
        free(hash->table2);
        free(hash->stash);
        free(hash->occupied1);
        free(hash->occupied2);
        free(hash);
		RETURN_FAIL;
    }
    
    hash->table_size = initial_size;
    hash->stash_count = 0;
    hash->count = 0;
    hash->obj_size = obj_size;
    hash->max_load_factor = max_load_factor;
    hash->hash1 = hash1;
    hash->hash2 = hash2;
    hash->compare = compare;
    
    return hash;
}

// Free the hash table
void cuckoo_free(CuckooHash* hash) {
    if (hash == NULL) return;
    
    // Free tables and stash
    free(hash->table1);
    free(hash->table2);
    free(hash->stash);
    free(hash->occupied1);
    free(hash->occupied2);
}

// Check if an object exists in the hash table
int cuckoo_contains(const CuckooHash* hash, const void* data) {
    if (hash == NULL || data == NULL) return false;
    
    // Check table1
    size_t index1 = hash->hash1(data, hash->table_size);
    if (hash->occupied1[index1]) {
        void* element = get_element(hash->table1, index1, hash->obj_size);
        if (hash->compare(element, data)) {
            return 1;
        }
    }
    
    // Check table2
    size_t index2 = hash->hash2(data, hash->table_size);
    if (hash->occupied2[index2]) {
        void* element = get_element(hash->table2, index2, hash->obj_size);
        if (hash->compare(element, data)) {
            return 2;
        }
    }
    
    // Check stash (only up to stash_count elements)
    for (size_t i = 0; i < hash->stash_count; i++) {
        void* element = get_element(hash->stash, i, hash->obj_size);
        if (hash->compare(element, data)) {
            return 3;
        }
    }
    
    return false;
}

// Forward declaration for insert function (needed for resize)
bool cuckoo_insert(CuckooHash* hash, const void* data);

// Resize the hash table
bool cuckoo_resize(CuckooHash* hash) {
    if (hash == NULL) return false;
    
    size_t new_size = hash->table_size * 2;
    size_t new_stash_capacity = calculate_stash_capacity(new_size);
    
    // Allocate new tables
    void* new_table1 = malloc(new_size * hash->obj_size);
    void* new_table2 = malloc(new_size * hash->obj_size);
    void* new_stash = malloc(new_stash_capacity * hash->obj_size);
    bool* new_occupied1 = (bool*)calloc(new_size, sizeof(bool));
    bool* new_occupied2 = (bool*)calloc(new_size, sizeof(bool));
    
    if (new_table1 == NULL || new_table2 == NULL || new_stash == NULL || new_occupied1 == NULL || new_occupied2 == NULL) {
        free(new_table1);
        free(new_table2);
        free(new_stash);
        free(new_occupied1);
        free(new_occupied2);
		RETURN_FAIL;
    }
    
    // Save old tables
    void* old_table1 = hash->table1;
    void* old_table2 = hash->table2;
    void* old_stash = hash->stash;
    bool* old_occupied1 = hash->occupied1;
    bool* old_occupied2 = hash->occupied2;
    size_t old_size = hash->table_size;
    size_t old_stash_count = hash->stash_count;
    
    // Update hash table with new tables
    hash->table1 = new_table1;
    hash->table2 = new_table2;
    hash->stash = new_stash;
    hash->occupied1 = new_occupied1;
    hash->occupied2 = new_occupied2;
    hash->table_size = new_size;
    hash->stash_capacity = new_stash_capacity;
    hash->count = 0;
    hash->stash_count = 0;
    
    // Reinsert all objects
    bool failed = false;
    
    // From table1
    for (size_t i = 0; i < old_size && !failed; i++) {
        if (old_occupied1[i]) {
            void* element = get_element(old_table1, i, hash->obj_size);
            if (!cuckoo_insert(hash, element)) {
                failed = true;
            }
        }
    }
    
    // From table2
    for (size_t i = 0; i < old_size && !failed; i++) {
        if (old_occupied2[i]) {
            void* element = get_element(old_table2, i, hash->obj_size);
            if (!cuckoo_insert(hash, element)) {
                failed = true;
            }
        }
    }
    
    // From stash (only up to stash_count elements)
    for (size_t i = 0; i < old_stash_count && !failed; i++) {
        void* element = get_element(old_stash, i, hash->obj_size);
        if (!cuckoo_insert(hash, element)) {
            failed = true;
        }
    }
    
    if (failed) {
        // Clean up and revert to old state
        free(new_table1);
        free(new_table2);
        free(new_stash);
        free(new_occupied1);
        free(new_occupied2);
        hash->table1 = old_table1;
        hash->table2 = old_table2;
        hash->stash = old_stash;
        hash->occupied1 = old_occupied1;
        hash->occupied2 = old_occupied2;
        hash->table_size = old_size;
        hash->stash_count = old_stash_count;
		RETURN_FAIL;
    }
    
    // Free old tables
    free(old_table1);
    free(old_table2);
    free(old_stash);
    free(old_occupied1);
    free(old_occupied2);
    
    return true;
}

// Insert an object into the hash table
bool cuckoo_insert(CuckooHash* hash, const void* data) {
    if (hash == NULL || data == NULL) return false;
    
    // Check if object already exists
    if (cuckoo_contains(hash, data)) {
        return true; // Object already in the hash table
    }
    
    // Check load factor and resize if needed
    double load_factor = (double)(hash->count + 1) / (hash->table_size * 2);
    if (load_factor > hash->max_load_factor) {
        if (!cuckoo_resize(hash)) {
            return false; // Resize failed
        }
    }
    
    // Calculate MAX_LOOPS based on current table size
    int MAX_LOOP = calculate_max_loops(hash->table_size);
    int loop_count = 0;
    
    // Try table1 first
    size_t index1 = hash->hash1(data, hash->table_size);
    if (!hash->occupied1[index1]) {
        set_element(hash->table1, index1, data, hash->obj_size);
        hash->occupied1[index1] = true;
        hash->count++;
        return true;
    }
    
    // Try table2 next
    size_t index2 = hash->hash2(data, hash->table_size);
    if (!hash->occupied2[index2]) {
        set_element(hash->table2, index2, data, hash->obj_size);
        hash->occupied2[index2] = true;
        hash->count++;
        return true;
    }
    
    // Make a copy of the data for cuckoo displacement
    char curr_data[hash->obj_size];
    memcpy(curr_data, data, hash->obj_size);

    // Start cuckoo hashing
    bool kick_slot1 = true; // Whether we're moving out table1 or table2
    
    while (loop_count < MAX_LOOP) {
        // Get index for current object
		index1 = hash->hash1(curr_data, hash->table_size);
		index2 = hash->hash2(curr_data, hash->table_size);
        
        // Swap objects
        char tmp[hash->obj_size];
        
		if (!hash->occupied1[index1]) {
			// Empty slot found
			set_element(hash->table1, index1, curr_data, hash->obj_size);
			hash->occupied1[index1] = true;
			hash->count++;
			return true;
		} else if (!hash->occupied2[index2]) {
			// Empty slot found
			set_element(hash->table2, index2, curr_data, hash->obj_size);
			hash->occupied2[index2] = true;
			hash->count++;
			return true;
        } else {
			void* table = kick_slot1 ? hash->table1 : hash->table2;
			int idx = kick_slot1 ? index1 : index2;
            // Save existing element
            memcpy(tmp, get_element(table, idx, hash->obj_size), hash->obj_size);
            // Replace with new element
            set_element(table, idx, curr_data, hash->obj_size);
			// Continue with the displaced object
			memcpy(curr_data, tmp, hash->obj_size);
			kick_slot1 = !kick_slot1; // Toggle between tables
			loop_count++;
        }
    }

    // If we get here, we couldn't find a slot after MAX_LOOP iterations
    // Try to add to the stash if there's room
    if (hash->stash_count < hash->stash_capacity) {
        // Add directly to the end of the stash
        set_element(hash->stash, hash->stash_count, curr_data, hash->obj_size);
        hash->stash_count++;
        hash->count++;
        return true;
    }
    
    // If stash is full, we fail
	RETURN_FAIL;
}

// Remove an object from the hash table
bool cuckoo_remove(CuckooHash* hash, const void* data) {
    if (hash == NULL || data == NULL) return false;
    
    // Check table1
    size_t index1 = hash->hash1(data, hash->table_size);
    if (hash->occupied1[index1]) {
        void* element = get_element(hash->table1, index1, hash->obj_size);
        if (hash->compare(element, data)) {
            hash->occupied1[index1] = false;
            hash->count--;
            return true;
        }
    }
    
    // Check table2
    size_t index2 = hash->hash2(data, hash->table_size);
    if (hash->occupied2[index2]) {
        void* element = get_element(hash->table2, index2, hash->obj_size);
        if (hash->compare(element, data)) {
            hash->occupied2[index2] = false;
            hash->count--;
            return true;
        }
    }
    
    // Check stash
    for (size_t i = 0; i < hash->stash_count; i++) {
        void* element = get_element(hash->stash, i, hash->obj_size);
        if (hash->compare(element, data)) {
            // Item found - remove it and shift remaining elements to keep stash compact
            if (i < hash->stash_count - 1) {
                // Shift all elements after the removed one to the left
                memmove(
                    get_element(hash->stash, i, hash->obj_size),
                    get_element(hash->stash, i + 1, hash->obj_size),
                    (hash->stash_count - i - 1) * hash->obj_size
                );
            }
            hash->stash_count--;
            hash->count--;
            return true;
        }
    }
    
    return false; // Not found
}

// Get the current number of elements in the hash table
size_t cuckoo_size(const CuckooHash* hash) {
    return hash ? hash->count : 0;
}

// Get the current load factor
double cuckoo_load_factor(const CuckooHash* hash) {
    if (hash == NULL) return 0.0;
    return (double)hash->count / (hash->table_size * 2);
}

// Get stash count
size_t cuckoo_stash_count(const CuckooHash* hash) {
    return hash ? hash->stash_count : 0;
}

// Get stash capacity
size_t cuckoo_stash_capacity(const CuckooHash* hash) {
    return hash ? hash->stash_capacity : 0;
}


// ---------------------------------- \\
// test code for cuckoo hashing below \\
// ---------------------------------- \\


/* 
bool int_compare(const void* a, const void* b) {
    return *(int*)a == *(int*)b;
}
size_t int_hash1(const void* key, size_t table_size) {
    int x = *(int*)key;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    return x % table_size;
}
size_t int_hash2(const void* key, size_t table_size) {
	int a = *(int*)key;
   a = (a+0x7ed55d16) + (a<<12);
   a = (a^0xc761c23c) ^ (a>>19);
   a = (a+0x165667b1) + (a<<5);
   a = (a+0xd3a2646c) ^ (a<<9);
   a = (a+0xfd7046c5) + (a<<3);
   a = (a^0xb55a4f09) ^ (a>>16);
   return a % table_size;
}

// these hash functions are actually not very good for cuckoo hashing
//size_t int_hash1(const void* key, size_t table_size) {
    //int value = *(int*)key;
    //return ((uint64_t)value * 2654435761ULL) % table_size;
//}
//size_t int_hash2(const void* key, size_t table_size) {
    //int value = *(int*)key;
    //return ((uint64_t)value * 1103515245ULL) % table_size;
//}

int main() {
    // Create a cuckoo hash table for integers
    CuckooHash* hash = cuckoo_create(16, 0.3, sizeof(int), int_hash1, int_hash2, int_compare);
    if (hash == NULL) {
        printf("Failed to create hash table\n");
        return 1;
    }
    
    // Print initial stash capacity and MAX_LOOPS
    printf("Initial table size: %zu\n", hash->table_size);
    printf("Initial stash capacity: %zu (log2(%zu))\n", hash->stash_capacity, hash->table_size);
    printf("Initial MAX_LOOPS: %d (2*log2(%zu))\n", calculate_max_loops(hash->table_size), hash->table_size);
    
    // Insert some data
    int data1 = 42;
    int data2 = 100;
    int data3 = 999;
    
    cuckoo_insert(hash, &data1);
    cuckoo_insert(hash, &data2);
    cuckoo_insert(hash, &data3);
    
    // Check if objects exist
    printf("Contains 42: %d\n", cuckoo_contains(hash, &data1));
    printf("Contains 100: %d\n", cuckoo_contains(hash, &data2));
    printf("Contains 999: %d\n", cuckoo_contains(hash, &data3));
    
    // Test a value that should not be in the table
    int not_found = 123;
    printf("Contains 123: %d\n", cuckoo_contains(hash, &not_found));
    
    // Remove an object
    cuckoo_remove(hash, &data2);
    printf("After removal, contains 100: %d\n", cuckoo_contains(hash, &data2));
    
    // Print some stats
    printf("Size: %zu\n", cuckoo_size(hash));
    printf("Load factor: %.2f\n", cuckoo_load_factor(hash));
    printf("Stash count: %zu/%zu\n", cuckoo_stash_count(hash), cuckoo_stash_capacity(hash));
    
    // Insert many elements to trigger resize
    printf("Inserting many elements...\n");
	const int MIN = 0, MAX = 10000000;
    for (int i = MIN; i < MAX; i++) {
        if (!cuckoo_insert(hash, &i)) {
			printf("Failed to insert %d\n", i);
			printf("\nChecking...\n");
			int tc[3];
			for (int i = MIN; i < MAX; i++) {
				int c = cuckoo_contains(hash, &i);
				if (c) tc[c-1]++;
			}
			printf("Table1: %d, Table2: %d, Stash: %d\n", tc[0], tc[1], tc[2]);
			exit(1);
		}
    }
    
    printf("New table size: %zu\n", hash->table_size);
    printf("New stash capacity: %zu (log2(%zu))\n", hash->stash_capacity, hash->table_size);
    printf("New MAX_LOOPS: %d (2*log2(%zu))\n", calculate_max_loops(hash->table_size), hash->table_size);
    printf("New size: %zu\n", cuckoo_size(hash));
    printf("New load factor: %.2f\n", cuckoo_load_factor(hash));
    printf("Stash count: %zu/%zu\n", cuckoo_stash_count(hash), cuckoo_stash_capacity(hash));

    printf("\nChecking...\n");
	int tc[3];
	for (int i = MIN; i < MAX; i++) {
		int c = cuckoo_contains(hash, &i);
		if (!c) {
			printf("Failed to find %d\n", i);
			exit(1);
		}
		tc[c-1]++;
	}
	printf("Table1: %d, Table2: %d, Stash: %d\n", tc[0], tc[1], tc[2]);
    
    printf("\nRemoving...\n");
    int remove_val = 1005;
    cuckoo_remove(hash, &remove_val);
    remove_val = 1010;
    cuckoo_remove(hash, &remove_val);
    
    // Verify items are still accessible
    printf("Checking some items are still accessible:\n");
    for (int i = 1000; i < 1020; i++) {
		printf("Contains %d: %d\n", i, cuckoo_contains(hash, &i));
    }
    
    printf("Stash count after removals: %zu/%zu\n", cuckoo_stash_count(hash), cuckoo_stash_capacity(hash));
    
    cuckoo_free(hash);
    
    printf("Test completed successfully\n");
    return 0;
}

*/
