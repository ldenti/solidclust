
#include <assert.h>
#include "../include/exception.h"
#include "../include/cluster.h"

int cluster_reads(char const *const index_filename) {
    int err;
    assert(index_filename);
    err = OK;
    /* memory map index file */
    /* keep 2 pointers: cumulative sum of lengths and current sketch */
    /* implement isonclust3 algorithm */
    return err;
}

int cluster_save(char const *const output_filename) {
    int err;
    assert(output_filename);
    err = OK;
    return err;
}