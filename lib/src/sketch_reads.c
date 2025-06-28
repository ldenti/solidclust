#include <stdio.h>
#include <assert.h>
#include <math.h> 
#include <zlib.h>
#include "../../bundled/klib/include/kseq.h"
#include "../../bundled/klib/include/ksort.h"
#include "../include/exception.h"
#include "../include/sketch_reads.h"
#include "../include/minimizer.h"
#include "../include/ioc_types.h"

KSEQ_INIT(gzFile, gzread)

typedef struct {
    sketch_metadata_t metadata;
    size_t sketch_offset;
} size_id_handle;
typedef kvec_t(size_id_handle) lv_t;

#define SSIZE(handle) ((handle).metadata.size)
KRADIX_SORT_INIT(sketch_size, size_id_handle, SSIZE, sizeof(len_t))

#define cmp(x, y) (x < y)
KSORT_INIT(minimizer, mm_t, cmp)

typedef kvec_t(unsigned char) qfilter_t;

typedef kvec_t(double) buffer_t;

double error_prob(char Q) {
    double tmp;
    tmp = ((int)Q - 33);
    return pow(10.0, - ((double)tmp) / 10.0);
}

int make_qual_filter(char const *const qual, const size_t seq_len, const uint8_t k, const double threshold, buffer_t *const buffer, qfilter_t *const filter) {
    int err;
    size_t i;
    double kmer_quality;
    assert(qual);
    assert(filter);
    assert(error_prob('!') == 1.0);
    err = OK;
    kv_resize(double, *buffer, seq_len);
    for (i = 0; i < seq_len; ++i) {
        kv_A(*buffer, i) = 1.0 - error_prob(qual[i]);
    }

    kv_reserve(unsigned char, *filter, seq_len - k + 1);
    for (filter->n = 0; filter->n < seq_len - k + 1; ++filter->n) {
        kmer_quality = 1.0;
        for (i = 0; i < k; ++i) {
            kmer_quality *= kv_A(*buffer, filter->n + i);
        }
        filter->a[filter->n] = kmer_quality > threshold;
    }
    assert(seq_len < k || filter->n == seq_len - k + 1);
    return err;
}

int sketch_reads_from_fastq(
    char const *const input_fastq, 
    const uint8_t k, 
    const uint8_t w, 
    const unsigned char canonical,
    const uint64_t seed, 
    const double quality_threshold,
    const char *tmp_index_filename
) {
    FILE *oh;
    gzFile fp;
    kseq_t *seq;
    mmv_t mmzers; /* concatenation of minimizer representation of each read */
    lv_t lengths;
    qfilter_t quality_filter;
    buffer_t buffer;
    size_t i, j, read_id;
    uint64_t cumulative_count;
    size_id_handle handle;
    int err;
    /* */
    /* 
    char **qnames = malloc(16384 * sizeof(char *));
    int qnames_m = 16384; 
    */

    assert(input_fastq);

    err = OK;
    if (!err && !(fp = gzopen(input_fastq, "r"))) err = ERR_FILE;
    if (!err && !(seq = kseq_init(fp))) err = ERR_MALLOC;
    kv_init(mmzers);
    kv_init(quality_filter);
    kv_init(lengths);
    kv_init(buffer);
    read_id = 0;
    cumulative_count = 0;
    while(!err && kseq_read(seq) >= 0) {
        /* 
        if (read_id == qnames_m) { 	
            qnames = realloc(qnames, 2*qnames_m*sizeof(char *));
            qnames_m *= 2; 
        } 
        qnames[read_id] = malloc((seq->name.l + 1) * sizeof(char));
        strncpy(qnames[read_id], seq->name.s, seq->name.l);
        qnames[read_id][seq->name.l] = '\0'; 
        */

        if (seq->seq.l < k) { /* TODO: add min read length to CLI */
            ++read_id;
            continue;
        }
        if (!err && seq->seq.l != seq->qual.l) err = ERR_RUNTIME;
        if (!err) err = make_qual_filter(seq->qual.s, seq->qual.l, k, quality_threshold, &buffer, &quality_filter);
        if (!err) err = minimizer_from_string(seq->seq.s, quality_filter.a, seq->seq.l, k, w, canonical, seed, &mmzers); /* this appends new minimizers to end of mmzers */
        handle.metadata.id = read_id;
        handle.metadata.size = mmzers.n - cumulative_count;
        handle.sketch_offset = cumulative_count;
        if (!err && handle.metadata.size > 0) {
            ks_introsort(minimizer, handle.metadata.size, mmzers.a + cumulative_count);
            for (i = 0, j = 0; i < handle.metadata.size; ++i) {
                if (mmzers.a[cumulative_count + i] != mmzers.a[cumulative_count + j]) {
                    ++j;
                    assert(j <= i);
                    mmzers.a[cumulative_count + j] = mmzers.a[cumulative_count + i];
                }
            }
            handle.metadata.size = ++j; /* j points to the most recently inserted item, so size is j + 1 */
        }
        kv_push(size_id_handle, lengths, handle);
        mmzers.n = cumulative_count + handle.metadata.size; /* resize mmzers to its new length without duplicates */
        cumulative_count = mmzers.n;
        ++read_id;
    }
    kv_destroy(quality_filter);
    kv_destroy(buffer);
    if (seq) kseq_destroy(seq);
    if (fp) gzclose(fp);

    /* assert(lengths.n == read_id); */

    radix_sort_sketch_size(lengths.a, lengths.a + lengths.n); /* Sort by increasing length sizes. Decreasing order is done when saving to disk */
    /*
    for(i = 0; i < lengths.n -1; ++i) {
       fprintf(stderr, "SKETCH - %s %d\n", qnames[kv_A(lengths, i).metadata.id], kv_A(lengths, i).metadata.size);
    }
    for (i = 0; i < read_id; ++i) free(qnames[i]);
    free(qnames);
    */

#ifndef NDEBUG
    if (lengths.n > 0) for (i = 0; !err && i < lengths.n - 1; ++i) {
        if (kv_A(lengths, i).metadata.size > kv_A(lengths, i + 1).metadata.size) {
            fprintf(stderr, "[sketch_reads] lengths are not stored in increasing order\n");
            err = ERR_LOGIC;
        }
    }
#endif
    /* fprintf(stderr, "%llu total sketches, %llu total minimizers\n", lengths.n, mmzers.n); */
    oh = NULL;
    if (!err && !(oh = fopen(tmp_index_filename, "wb"))) err = ERR_FILE;
    cumulative_count = lengths.n;
    if (!err && fwrite(&cumulative_count, sizeof(cumulative_count), 1, oh) != 1) err = ERR_IO;
    /* fprintf(stderr, "[write] nsketches: %llu\n", cumulative_count); */
    cumulative_count = 0;

    /* ATTENTION: the following loop is in reverse order to save sketches from LONGEST to SHORTEST */
    for (i = lengths.n - 1; !err && i != SIZE_MAX; --i) { /* write cumulative lengths */
        if (fwrite(&lengths.a[i].metadata, sizeof(sketch_metadata_t), 1, oh) != 1) err = ERR_IO;
    }
    assert(ftell(oh) == (lengths.n * sizeof(sketch_metadata_t) + 8));
    for (i = lengths.n - 1; !err && i != SIZE_MAX; --i) { /* write sketches in the order given by lengths (from longest to shortest) */
        size_id_handle const *const record = &lengths.a[i]; /* DO NOT merge the two for loops since actual sketches come after their lengths */
        if (fwrite(&mmzers.a[record->sketch_offset], sizeof(mm_t), record->metadata.size, oh) != record->metadata.size) err = ERR_IO;
    }
    if (oh) {
        fflush(oh);
        fclose(oh);
        oh = NULL;
    }
    kv_destroy(lengths);
    kv_destroy(mmzers);
    return err;
}
