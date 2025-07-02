#ifndef CLUSTER_H
#define CLUSTER_H

#include "../../bundled/klib/include/kvec.h"
#include "ioc_types.h"
#include "minimizer.h"

typedef kvec_t(read_id_t) idv_t;

typedef struct {
    idv_t ids;
    mmv_t minimizers;
} cluster_t; /* clusters are (ordered) vectors of minimizers */
typedef kvec_t(cluster_t) clusters_t;

int cluster_reads_weighted(
    char const *const index_filename, 
    const double similarity_threshold,
    const double weight_threshold,
    clusters_t *const clusters
    /* void **mm2clusters // for post-processing */
);

/* TODO
int cluster_postprocessing(
    const double merge_threhsold,
    clusters_t *const clusters,
    void *const mm2clusters
);
*/

int cluster_save(
    clusters_t const *const clusters, 
    char const *const output_folder
);

int cluster_print(clusters_t const *const clusters);

/* 
    Use this function to deallocate the hash table (minimizers -> cluster ids).
    It must use void* since table types are bounded to the .c source file
    (klib's macros don't work well in header files). 
    The actual parameter type is (mm2cluster_t*) and the function must be called after postprocessing.
*/
int minimizers_to_clusters_map_destroy(void *const ptr);

#endif /* CLUSTER_H */
