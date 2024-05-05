#ifndef INTERFAISS_H
#define INTERFAISS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

#define INTERFAISS_READ_BUFFER_SIZE 1048576
#define INTERFAISS_INIT_BUFFER_SIZE 8

typedef enum {
    INDEX_FLAT=0,
    INDEX_PQ=1,
    INDEX_HNSW=2,
    INDEX_NUM_IDX=3,
} index_type_t;

typedef struct {
    index_type_t type;
    void* index;
    int d;
} index_t;

typedef float dim_t;
typedef int64_t idx_t;

typedef struct {
    // PQ index parameters
    int pq_m;       // Number of subquantizers
    int pq_nbits;   // Bits per subquantizer
    // HNSW index parameters
    int hnsw_m;               // Hyperparameter M
    //int hnsw_efConstruction;  // Hyperparameter efConstruction for building, not used yet
} index_parameters_t;

index_t* interfaiss_create_index_from_file(const char* filename, index_type_t type, index_parameters_t params);
index_t* interfaiss_create_index(index_type_t type, int d, index_parameters_t params);
void interfaiss_query_index(index_t* idx, idx_t n, const dim_t* queries, int d, int k, float* distances, idx_t* indices);
void interfaiss_add_data_to_index(index_t* idx, idx_t n, const dim_t* data, int d);
void interfaiss_train_index(index_t *idx, idx_t n, const dim_t *data);
void interfaiss_free_index(index_t* idx);

// Utils
void parse_embeddings(const char *filename, dim_t **dims, idx_t *no_embeddings, int *d, char ***data_names);

#ifdef __cplusplus
}
#endif

#endif // INTERFAISS_H
