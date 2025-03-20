#ifndef CUCKOO_HASH_H
#define CUCKOO_HASH_H

#define CUCKOO_EXIT_ON_FAIL 1

#include <stddef.h>
#include <stdbool.h>

/**
 * @brief Cuckoo hash table structure
 */
typedef struct CuckooHash CuckooHash;

/**
 * @brief Function pointer type for hash functions
 * @param key The key to hash
 * @param table_size The size of the hash table
 * @return The hash value
 */
typedef size_t (*HashFunction)(const void* key, size_t table_size);

/**
 * @brief Function pointer type for key comparison
 * @param a First key
 * @param b Second key
 * @return true if keys are equal, false otherwise
 */
typedef bool (*CompareFunction)(const void* a, const void* b);

// Cuckoo hash table structure
struct CuckooHash {
    void* table1;            // First hash table - raw memory block
    void* table2;            // Second hash table - raw memory block
    void* stash;             // Stash for failed insertions - raw memory block
    bool* occupied1;         // Tracks which slots are occupied in table1
    bool* occupied2;         // Tracks which slots are occupied in table2
    size_t table_size;       // Size of each table
    size_t stash_capacity;   // Capacity of the stash
    size_t stash_count;      // Current number of elements in stash
    size_t count;            // Total number of elements in the hash table
    size_t obj_size;         // Size of each object
    double max_load_factor;  // Maximum load factor before resize
    HashFunction hash1;      // First hash function
    HashFunction hash2;      // Second hash function
    CompareFunction compare; // Function to compare two objects
};

/**
 * @brief Create a new cuckoo hash table
 * @param initial_size Initial size of each hash table
 * @param max_load_factor Maximum load factor before resizing (0 < max_load_factor < 1)
 * @param obj_size Size of each object in bytes
 * @param hash1 First hash function
 * @param hash2 Second hash function
 * @param compare Function to compare two keys
 * @return true if successful, false otherwise
 */
bool cuckoo_create(CuckooHash* hash,
                   size_t initial_size, double max_load_factor, 
                   size_t obj_size, HashFunction hash1, HashFunction hash2, 
                   CompareFunction compare);

/**
 * @brief Free the hash table and all associated memory
 * @param hash The hash table to free
 */
void cuckoo_free(CuckooHash* hash);

/**
 * @brief Check if an object exists in the hash table
 * @param hash The hash table
 * @param data The object to look for
 * @return 0 if not found, 1 if found in table1, 2 if found in table2, 3 if found in stash
 */
int cuckoo_contains(const CuckooHash* hash, const void* data);

/**
 * @brief Insert an object into the hash table
 * @param hash The hash table
 * @param data The object to insert
 * @return true if insertion was successful, false otherwise
 */
bool cuckoo_insert(CuckooHash* hash, const void* data);

/**
 * @brief Remove an object from the hash table
 * @param hash The hash table
 * @param data The object to remove
 * @return true if object was found and removed, false otherwise
 */
bool cuckoo_remove(CuckooHash* hash, const void* data);

/**
 * @brief Get the current number of elements in the hash table
 * @param hash The hash table
 * @return Number of elements
 */
size_t cuckoo_size(const CuckooHash* hash);

/**
 * @brief Get the current load factor of the hash table
 * @param hash The hash table
 * @return Current load factor (elements / (2 * table_size))
 */
double cuckoo_load_factor(const CuckooHash* hash);

/**
 * @brief Get the number of elements stored in the stash
 * @param hash The hash table
 * @return Number of elements in stash
 */
size_t cuckoo_stash_count(const CuckooHash* hash);

/**
 * @brief Get the maximum capacity of the stash
 * @param hash The hash table
 * @return Stash capacity
 */
size_t cuckoo_stash_capacity(const CuckooHash* hash);

#endif /* CUCKOO_HASH_H */
