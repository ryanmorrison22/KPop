#include "interfaiss.h"
#include <faiss/IndexFlat.h>
#include <faiss/IndexPQ.h>
#include <faiss/IndexHNSW.h>
#include <stdexcept>

extern "C" {

index_t* interfaiss_create_index(index_type_t type, int d, index_parameters_t params) {
    index_t* idx = new index_t;
    idx->type = type;
    idx->d = d;

    switch (type) {
        case INDEX_FLAT:
            idx->index = new faiss::IndexFlatL2(d);
            break;
        case INDEX_PQ:
            // Check that m * nbits <= d to ensure each dimension is quantized
            if (params.pq_m * params.pq_nbits > d) {
                fprintf(stderr, "Error: Invalid parameters for PQ index. Ensure that pq_m * pq_nbits <= dimension d.\n");
                delete idx;
                return NULL;
            }
            idx->index = new faiss::IndexPQ(d, params.pq_m, params.pq_nbits);
            break;
        case INDEX_HNSW:
            idx->index = new faiss::IndexHNSWFlat(d, params.hnsw_m);
            break;
        default:
            fprintf(stderr, "Error: Invalid index type specified.\n");
            delete idx;
            return NULL;
    }
    return idx;
}

index_t* interfaiss_create_index_from_file(const char* filename, index_type_t type, index_parameters_t params) {
    dim_t* values = nullptr;
    idx_t no_embeddings = 0;
    int d = 0;
    char **data_names = nullptr;
    parse_embeddings(const_cast<char*>(filename), &values, &no_embeddings, &d, &data_names);

    index_t* idx = interfaiss_create_index(type, d, params);
    interfaiss_train_index(idx, no_embeddings, values);
    interfaiss_add_data_to_index(idx, no_embeddings, values, d);
    // Return later on if desired
    for(idx_t i = 0; i < no_embeddings; i++) {
        free(data_names[i]);
    }
    free(data_names);
    return idx;
}

void interfaiss_query_index(index_t* idx, idx_t n, const dim_t* queries, int d, int k, float* distances, idx_t* indices) {
    if (idx->d != d) {
        throw std::invalid_argument("Dimension mismatch between index and queries");
    }
    reinterpret_cast<faiss::Index*>(idx->index)->search(n, queries, k, distances, reinterpret_cast<faiss::idx_t*>(indices));
}

void interfaiss_add_data_to_index(index_t* idx, idx_t n, const dim_t* data, int d) {
    if (idx->d != d) {
        throw std::invalid_argument("Dimension mismatch between index and data");
    }
    reinterpret_cast<faiss::Index*>(idx->index)->add(n, data);
}

void interfaiss_train_index(index_t* idx, idx_t n, const dim_t* data) {
    if (!idx || !data) {
        fprintf(stderr,"Error: Invalid index or data pointer provided for training.\n");
        return;
    }

    switch (idx->type) {
        case INDEX_PQ: {
            // Training a PQ index
            auto pq_index = dynamic_cast<faiss::IndexPQ*>(reinterpret_cast<faiss::Index*>(idx->index));
            if (!pq_index->is_trained) {
                pq_index->train(n, data);
            }
            break;
        }
        case INDEX_FLAT:
        case INDEX_HNSW:
        default:
            fprintf(stderr,"Training not required or applicable for this index type.\n");
            break;
    }
}

void interfaiss_free_index(index_t* idx) {
    delete reinterpret_cast<faiss::Index*>(idx->index);
    delete idx;
}

}

void parse_embeddings(const char* filename,
                      dim_t** dims,
                      idx_t* no_embeddings,
                      int* d,
                      char*** data_names) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        perror("Failed to open file");
        return;
    }

    char line[INTERFAISS_READ_BUFFER_SIZE];
    size_t linecap = INTERFAISS_READ_BUFFER_SIZE;

    // Read the header line to determine the dimension
    if (!fgets(line, linecap, file)) {
        fclose(file);
        return;
    }

    // Count number of tabs to find dimension count
    *d = 0;
    char* p = line;
    while (*p) {
        if (*p == '\t') (*d)++;
        p++;
    }
    // Adjust dimension count (excluding name column)
    *d -= 1;

    // Initialize the embedding storage and data name storage
    *dims = NULL;
    *data_names = NULL;
    *no_embeddings = 0;
    size_t cap = 0;

    while (fgets(line, linecap, file)) {
        if (*no_embeddings >= cap) {
            cap = cap > 0 ? cap * 2 : 1;
            *dims = (dim_t*)realloc(*dims, cap * (*d) * sizeof(dim_t));
            *data_names = (char**)realloc(*data_names, cap * sizeof(char*));
        }

        p = line;
        char* end = p;
        while (*end && *end != '\t') end++; // Skip to first tab (end of name token)

        // Allocate and store the name
        (*data_names)[*no_embeddings] = (char*)malloc(end - p + 1);
        memcpy((*data_names)[*no_embeddings], p, end - p);
        (*data_names)[*no_embeddings][end - p] = '\0';

        p = end + 1; // Move past the tab
        idx_t idx = 0;

        // Parse dimensions
        while (*p && idx < *d) {
            while (*p == '\t') p++; // Skip any extra tabs
            dim_t value = strtof(p, &end); // Convert and get next start
            if (p == end) break; // No conversion happened, break
            (*dims)[(*no_embeddings) * (*d) + idx++] = value;
            p = end;
        }

        if (idx == *d) {
            (*no_embeddings)++;
        }
    }

    fclose(file);
}
