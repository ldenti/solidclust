#ifndef SKETCH_READS_H
#define SKETCH_READS_H

#include <stddef.h>
#include <stdint.h>
#include "ioc_types.h"

typedef struct {
    read_id_t id;
    len_t size;
} sketch_metadata_t;

int sketch_reads_from_fastq(
    char const *const input_fastq, 
    const uint8_t k, 
    const uint8_t w, 
    const unsigned char canonical,
    const uint64_t seed, 
    const double quality_threshold,
    const char *tmp_index_filename
);

#endif /* SKETCH_READS_H */
