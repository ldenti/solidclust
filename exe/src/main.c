#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#define __USE_POSIX199309
#include <time.h>
#include "../../bundled/xoshiro/include/xoshiro256++.h"
#include "../../bundled/klib/include/ketopt.h"
#include "../../lib/include/exception.h"
#include "../include/option.h"
#include "../../lib/include/sketch_reads.h"
#include "../../lib/include/cluster.h"

#include <assert.h>

#define NS_IN_SEC (1000 * 1000 * 1000)
#define NS_IN_USEC 1000
#define RESOLUTION NS_IN_USEC

int print_usage(FILE *ostrm);
int parse_options(int argc, char **argv, option_t *const opts);
int generate_tmp_filename(char *const buffer, size_t buffer_len);

int main(int argc, char **argv) {
    int err;
    size_t i;
    option_t opts;
    clusters_t clusters;
    struct timespec tstart, tstop;
    unsigned long long wallclock_elapsed;
    void *mm2clusters; /*.auxiliary information used by  */

    err = OK;
    mm2clusters = NULL;
    init_options(&opts);
    kv_init(clusters);
    if (!err) {
        /* clock_gettime(CLOCK_MONOTONIC, &tstart); */
        err = parse_options(argc, argv, &opts);
        /* 
        clock_gettime(CLOCK_MONOTONIC, &tstop);
        wallclock_elapsed = (tstop.tv_sec - tstart.tv_sec) * NS_IN_SEC;
        wallclock_elapsed += tstop.tv_nsec - tstart.tv_nsec;
        fprintf(stderr, "Options parsed in %llu us\n", wallclock_elapsed / RESOLUTION); 
        */
    }

    fprintf(stderr, "\n=== Parameters:\n");
    /* fprintf(stderr, "---     k: %ld\n", opts.k); */
    /* fprintf(stderr, "---     w: %ld\n", opts.w); */
    fprintf(stderr, "---     q: %lf\n", opts.quality_threshold);
    fprintf(stderr, "---     t: %lf\n", opts.similarity_threshold);
    fprintf(stderr, "---     w: %lf\n", opts.weight_threshold);
    fprintf(stderr, "===\n\n");
    if (!err) {
        clock_gettime(CLOCK_MONOTONIC, &tstart);
        err = sketch_reads_from_fastq(opts.input_fastq, opts.k, opts.w, opts.canonical, opts.seed, opts.quality_threshold, opts.tmp_filename);
        clock_gettime(CLOCK_MONOTONIC, &tstop);
        wallclock_elapsed = (tstop.tv_sec - tstart.tv_sec) * NS_IN_SEC;
        wallclock_elapsed += tstop.tv_nsec - tstart.tv_nsec;
        fprintf(stderr, "Reads sketched in %llu us\n", wallclock_elapsed / RESOLUTION);
    }
    if (!err) {
        clock_gettime(CLOCK_MONOTONIC, &tstart);
        if (opts.weighted) err = cluster_reads_weighted(opts.tmp_filename, opts.similarity_threshold, opts.weight_threshold, &clusters);
        else err = cluster_reads(opts.tmp_filename, opts.similarity_threshold, &clusters, &mm2clusters);
        clock_gettime(CLOCK_MONOTONIC, &tstop);
        wallclock_elapsed = (tstop.tv_sec - tstart.tv_sec) * NS_IN_SEC;
        wallclock_elapsed += tstop.tv_nsec - tstart.tv_nsec;
        fprintf(stderr, "Reads clustered in %llu us\n", wallclock_elapsed / RESOLUTION);
    }
    if (!err) {
        clock_gettime(CLOCK_MONOTONIC, &tstart);
        if (opts.weighted) err = cluster_slow_postprocessing(opts.post_cluster, &clusters); /* TODO remove in final version */
        else err = cluster_postprocessing(opts.post_cluster, &clusters, mm2clusters);
        if (!err) err = minimizers_to_clusters_map_destroy(mm2clusters);
        mm2clusters = NULL;
        clock_gettime(CLOCK_MONOTONIC, &tstop);
        wallclock_elapsed = (tstop.tv_sec - tstart.tv_sec) * NS_IN_SEC;
        wallclock_elapsed += tstop.tv_nsec - tstart.tv_nsec;
        fprintf(stderr, "Postprocessing took %llu us\n", wallclock_elapsed / RESOLUTION);
    }
    if (!err) {
        clock_gettime(CLOCK_MONOTONIC, &tstart);
        err = cluster_save(&clusters, opts.output_mapping);
        clock_gettime(CLOCK_MONOTONIC, &tstop);
        wallclock_elapsed = (tstop.tv_sec - tstart.tv_sec) * NS_IN_SEC;
        wallclock_elapsed += tstop.tv_nsec - tstart.tv_nsec;
        fprintf(stderr, "Cluster written in %llu us\n", wallclock_elapsed / RESOLUTION);
    }
    /* if (!err) err = cluster_print(&clusters); */
    for (i = 0; i < clusters.n; ++i) {
        kv_destroy(clusters.a[i].ids);
        kv_destroy(clusters.a[i].minimizers);
    }
    kv_destroy(clusters);
    destroy_options(&opts);
    return err;
}

int parse_options(int argc, char **argv, option_t *const opts) {
    static ko_longopt_t longopts[] = { /*
		{ "foo", ko_no_argument,       301 },
		{ "bar", ko_required_argument, 302 },
		{ "opt", ko_optional_argument, 303 }, */
        {"post-cluster", ko_required_argument, 301},
        {"weighted", ko_required_argument, 302},
		{NULL, 0, 0}
	};
    static char const* const shortopts = "i:o:m:k:w:s:d:q:t:ch";
    int c;
    unsigned long buffer;
    unsigned char k_explicit;
    unsigned char w_explicit;
    unsigned char t1_explicit;
    unsigned char t2_explicit;
    unsigned char post_explicit;
    ketopt_t opt;
    
    assert(argv);
    assert(opts);

    opt = KETOPT_INIT;
    k_explicit = FALSE;
    w_explicit = FALSE;
    t1_explicit = FALSE;
    t2_explicit = FALSE;
    post_explicit = FALSE;
    while((c = ketopt(&opt, argc, argv, 1, shortopts, longopts)) >= 0) {
        switch (c) {
            case 'i': opts->input_fastq = opt.arg; break;
            case 'o': opts->output_mapping = opt.arg; break;
            case 'm': 
                if (strcmp(opt.arg, "ont") == 0) {
                    if (!k_explicit) opts->k = 13;
                    if (!w_explicit) opts->w = 21;
                    if (!t1_explicit) opts->quality_threshold = pow(0.95, opts->k);
                    if (!t2_explicit) opts->similarity_threshold = 0.5;
                    if (!post_explicit) opts->post_cluster = 0.5;
                } else if (strcmp(opt.arg, "pacbio") == 0) {
                    if (!k_explicit) opts->k = 15;
                    if (!w_explicit) opts->w = 51;
                    if (!t1_explicit) opts->quality_threshold = pow(0.98, opts->k);
                    if (!t2_explicit) opts->similarity_threshold = 0.5;
                    if (!post_explicit) opts->post_cluster = 0.8;
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
            case 'q':
                opts->quality_threshold = strtod(opt.arg, NULL);
                if (opts->quality_threshold < 0 || opts->quality_threshold > 1) {
                    fprintf(stderr, "quality threshold must be in [0, 1]");
                    return ERR_PARSE;
                }
                break;
            case 't':
                opts->similarity_threshold = strtod(opt.arg, NULL);
                if (opts->similarity_threshold < 0 || opts->similarity_threshold > 1) {
                    fprintf(stderr, "similarity threshold for clustering must be in [0, 1]");
                    return ERR_PARSE;
                }
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
            case 301:
                opts->post_cluster = strtod(opt.arg, NULL);
                post_explicit = TRUE;
                break;
            case 302:
                opts->weighted = TRUE;
                opts->weight_threshold = strtod(opt.arg, NULL);
                break;
            case 'h':
                print_usage(stdout);
                exit(EXIT_SUCCESS);
            default:
                fprintf(stderr, "Unrecognized option: %d -> %s\n", c, opt.arg);
                print_usage(stderr);
                return ERR_PARSE;
        }
    }
    if (!opts->input_fastq) {
        fprintf(stderr, "Input fastq is required\n");
        return ERR_PARSE;
    }
    if (!opts->tmp_filename) {
        if (!(opts->tmp_filename = malloc(TMP_FILE_SIZE + 5))) return ERR_MALLOC;
        if (generate_tmp_filename(opts->tmp_filename, TMP_FILE_SIZE + 5)) return ERR_LOGIC;
        opts->tmp_filename_allocated = TRUE;
    }
    if (opts->w >= opts->k) {
        opts->w = opts->w - opts->k + 1; /* convert w from number of bases to number of k-mers */
    } else {
        fprintf(stderr, "window size cannot be smaller than k\n");
        return ERR_PARSE;
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
    fprintf(ostrm, "isONClust3, rewritten in C\n");
    fprintf(ostrm, "Options:\n");
    fprintf(ostrm, "\t-i : Fastq file in input\n");
    fprintf(ostrm, "\t-o : output CSV file mapping reads to clusters (by index)\n");
    fprintf(ostrm, "\t-d : temporary filename [random name]\n");
    fprintf(ostrm, "\t-m : preset mode [ont, pacbio]. This option sets k, w, q, t and post-cluster with default values, if not explicitly defined\n");
    fprintf(ostrm, "\t-k : k-mer length [14]\n");
    fprintf(ostrm, "\t-w : window length for computing minimizers (in number of bases > k) [30]\n");
    fprintf(ostrm, "\t-q : quality threshold [0.97]\n");
    fprintf(ostrm, "\t-t : similarity threshold for clustering [0.5]\n");
    fprintf(ostrm, "\t-c : disable canonical minimizers [enabled]\n");
    fprintf(ostrm, "\t--post-cluster : threshold used for merging clusters together [no merging]\n");
    fprintf(ostrm, "\t-s : random seed [42]\n");
    fprintf(ostrm, "\t-h : print this help\n");
    return OK;
}

