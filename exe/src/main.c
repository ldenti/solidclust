#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include "../../bundled/xoshiro/include/xoshiro256++.h"
#include "../../bundled/klib/include/ketopt.h"
#include "../../lib/include/exception.h"
#include "../include/option.h"
#include "../../lib/include/sketch_reads.h"
#include "../../lib/include/cluster.h"

int print_usage(FILE *ostrm);
int parse_options(int argc, char **argv, option_t *const opts);
int generate_tmp_filename(char *const buffer, size_t buffer_len);

int main(int argc, char **argv) {
    int err;
    option_t opts;
    err = OK;
    init_options(&opts);
    if (!err) err = parse_options(argc, argv, &opts);
    if (!err) err = sketch_reads_from_fastq(opts.input_fastq, opts.k, opts.w, opts.canonical, opts.seed, opts.threshold, opts.tmp_filename);
    if (!err) err = cluster_reads(opts.tmp_filename);
    if (!err) err = cluster_save(opts.output_dir);
    destroy_options(&opts);
    return err;
}

int parse_options(int argc, char **argv, option_t *const opts) {
    static ko_longopt_t longopts[] = { /*
		{ "foo", ko_no_argument,       301 },
		{ "bar", ko_required_argument, 302 },
		{ "opt", ko_optional_argument, 303 }, */
		{ NULL, 0, 0 }
	};
    static char const* const shortopts = "i:o:m:k:w:s:d:ch";
    int c;
    unsigned long buffer;
    unsigned char k_explicit;
    unsigned char w_explicit;
    ketopt_t opt;
    
    assert(argv);
    assert(opts);

    opt = KETOPT_INIT;
    k_explicit = FALSE;
    w_explicit = FALSE;
    while((c = ketopt(&opt, argc, argv, 1, shortopts, longopts)) >= 0) {
        switch (c)
        {
            case 'i': opts->input_fastq = opt.arg; break;
            case 'o': opts->output_dir = opt.arg; break;
            case 'm': 
                if (strcmp(opt.arg, "ont") == 0) {
                    if (!k_explicit) opts->k = 13;
                    if (!w_explicit) opts->w = 21;
                } else if (strcmp(opt.arg, "pacbio") == 0) {
                    if (!k_explicit) opts->k = 15;
                    if (!w_explicit) opts->w = 51;
                } else {
                    fprintf(stderr, "Invalid mode, available options are [ont, pacbio]\n");
                    return ERR_PARSE;
                }
                break;
            case 'k': 
                buffer = strtoul(opt.arg, NULL, 10);
                if (buffer > UINT8_MAX) {
                    fprintf(stderr, "k must be < 256\n");
                    return ERR_PARSE;
                }
                opts->k = (uint8_t)buffer;
                k_explicit = TRUE;
                break;
            case 'w':
                buffer = strtoul(opt.arg, NULL, 10);
                if (buffer > UINT8_MAX) {
                    fprintf(stderr, "w must be < 256\n");
                    return ERR_PARSE;
                }
                opts->w = (uint8_t)buffer;
                w_explicit = TRUE;
                break;
            case 'c':
                opts->canonical = TRUE;
                break;
            case 's':
                opts->seed = strtoul(opt.arg, NULL, 10);
                break;
            case 'd':
                opts->tmp_filename = opt.arg;
                break;
            case 'h':
                print_usage(stdout);
                exit(EXIT_SUCCESS);
            default:
                fprintf(stderr, "Unrecognized option\n");
                print_usage(stderr);
                return ERR_PARSE;
        }
    }
    if (!opts->input_fastq) {
        fprintf(stderr, "Input fastq is required\n");
        return ERR_PARSE;
    }
    if (!opts->output_dir) {
        fprintf(stderr, "Output folder is required\n");
        return ERR_PARSE;
    }
    if (!opts->tmp_filename) {
        if (!(opts->tmp_filename = malloc(TMP_FILE_SIZE + 5))) return ERR_MALLOC;
        if (generate_tmp_filename(opts->tmp_filename, TMP_FILE_SIZE + 5)) return ERR_LOGIC;
        opts->tmp_filename_allocated = TRUE;
    }
    return OK;
}

int generate_tmp_filename(char *const buffer, size_t buffer_len) {
    size_t i;
    xoshiro256pp_state_t state;
    if (buffer_len < 6) return ERR_LOGIC;
    xoshiro256pp_seed_state(&state, time(NULL));
    for (i = 0; i < buffer_len - 5; ++i) {
        buffer[i] = 'A' + (state.s[0] % 26);
        xoshiro256pp_next(&state);
    }
    if (memcpy(&buffer[buffer_len - 5], ".bin", 5) != &buffer[buffer_len - 5]) return ERR_LOGIC; /* 5 for \0 too */
    return OK;
}

int print_usage(FILE *ostrm) {
    fprintf(ostrm, "TODO\n");
    return OK;
}

