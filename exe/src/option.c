#include <stddef.h>
#include <stdlib.h>
#include <assert.h>
#include "../include/option.h"

void init_options(option_t *const opts) {
    assert(opts);
    opts->input_fastq = NULL;
    opts->output_dir = NULL;
    opts->tmp_filename = NULL;
    opts->k = 0;
    opts->w = 0;
    opts->canonical = 0;
    opts->quality_threshold = 0.95;
    opts->similarity_threshold = 0.5;
    opts->post_cluster = 0.0;
    opts->seed = 42;
    opts->tmp_filename_allocated = 0;
}

void destroy_options(option_t *const opts) {
    assert(opts);
    if (opts->tmp_filename_allocated) free(opts->tmp_filename);
}