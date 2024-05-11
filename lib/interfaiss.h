/*
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

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

index_t* interfaiss_create_flat_index(int d);
// PQ index parameters are: number of subquantizers, bits per subquantizer
index_t* interfaiss_create_PQ_index(int d, int m, int n_bits);
// HNSW index parameters are: hyperparameter M, hyperparameter efConstruction (efConstruction not used yet)
index_t* interfaiss_create_HNSW_index(int d, int m, int ef_construction);

void interfaiss_query_index(index_t* idx, int d, idx_t n, const dim_t* queries, idx_t* k, float** distances, idx_t** indices);
void interfaiss_add_data_to_index(index_t* idx, int d, idx_t n, const dim_t* data);
void interfaiss_train_index(index_t *idx, int d, idx_t n, const dim_t *data);
void interfaiss_free_index(index_t* idx);

#ifdef __cplusplus
}
#endif

#endif // INTERFAISS_H

