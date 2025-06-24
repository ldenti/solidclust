#include <stddef.h>
#include <stdlib.h>
#include "../include/option.h"
#include "../../lib/include/exception.h"
#include <assert.h>

void init_options(option_t *const opts) {
    assert(opts);
    opts->input_fastq = NULL;
    opts->output_mapping = NULL;
    opts->tmp_filename = NULL;
    opts->k = 14;
    opts->w = 30;
    opts->canonical = 1;
    opts->weighted = FALSE;
    opts->quality_threshold = 0.97;
    opts->similarity_threshold = 0.5;
    opts->weight_threshold = 0.2;
    opts->post_cluster = 0.0;
    opts->seed = 42;
    opts->tmp_filename_allocated = 0;
}

void destroy_options(option_t *const opts) {
    assert(opts);
    if (opts->tmp_filename_allocated) free(opts->tmp_filename);
}
