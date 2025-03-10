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

int cluster_reads(
    char const *const index_filename, 
    const double similarity_threshold, 
    const double merge_threhsold,
    clusters_t* clusters
);

int cluster_save(
    clusters_t const *const clusters, 
    char const *const output_folder
);

#endif /* CLUSTER_H */
