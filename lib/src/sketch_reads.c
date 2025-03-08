#include <stdio.h>
#include <assert.h>
#include <math.h> 
#include <zlib.h>
#include "../../bundled/klib/include/kseq.h"
#include "../../bundled/klib/include/ksort.h"
#include "../include/exception.h"
#include "../include/sketch_reads.h"
#include "../include/minimizer.h"

KSEQ_INIT(gzFile, gzread)

typedef struct {
    size_t sketch_size;
    size_t read_id;
    mm_t *sketch;
} size_id_handle;
typedef kvec_t(size_id_handle) lv_t;

#define SSIZE(handle) ((handle).sketch_size)
KRADIX_SORT_INIT(sketch_size, size_id_handle, SSIZE, sizeof(size_t))

#define cmp(x, y) (x < y)
KSORT_INIT(minimizer, mm_t, cmp)

typedef kvec_t(unsigned char) qfilter_t;

int make_qual_filter(char const *const qual, const size_t seq_len, const uint8_t k, const double threshold, qfilter_t *const filter) {
    int err;
    size_t i, buf_idx;
    double *buffer, q, kmer_quality;
    assert(qual);
    assert(filter);
    err = OK;
    filter->n = 0; /* clean */
    kmer_quality = 1;
    buf_idx = 0;
    if (!(buffer = malloc(k * sizeof(double)))) err = ERR_MALLOC;
    if (!err) for (i = 0; i < k; ++i) buffer[i] = 1;
    if (!err) for (i = 0; i < seq_len; ++i) {
        q = buffer[buf_idx];
        buffer[buf_idx] = 1.0 - pow(10.0, -((int)qual[i] - 33) / 10);
        kmer_quality /= q; /* no need to check if i >= k since buffer is all 1s at the beginning */
        kmer_quality *= buffer[buf_idx];
        kv_push(unsigned char, *filter, kmer_quality > threshold);
    }
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
    size_t i, read_id, minimizer_count;
    uint64_t cumulative_count;
    size_id_handle handle;
    int err;
    assert(input_fastq);
    err = OK;
    kv_init(mmzers);
    if (!err && !(fp = gzopen(input_fastq, "r"))) err = ERR_FILE;
    if (!err && !(seq = kseq_init(fp))) err = ERR_MALLOC;
    kv_init(quality_filter);
    kv_init(lengths);
    read_id = 0;
    cumulative_count = 0;
    while(!err && kseq_read(seq) >= 0) {
        if (!err && seq->seq.l != seq->qual.l) err = ERR_RUNTIME;
        if (!err) err = make_qual_filter(seq->qual.s, seq->qual.l, k, quality_threshold, &quality_filter);
        if (!err) err = minimizer_from_string(seq->seq.s, quality_filter.a, seq->seq.l, k, w, canonical, seed, &minimizer_count, &mmzers); /* this appends new minimizers to end of mmzers */
        ks_introsort(minimizer, minimizer_count, mmzers.a + cumulative_count);
        handle.read_id = read_id;
        handle.sketch_size = minimizer_count;
        handle.sketch = mmzers.a + cumulative_count;
        kv_push(size_id_handle, lengths, handle);
        cumulative_count += minimizer_count;
        ++read_id;
    }
    kv_destroy(quality_filter);
    if (seq) kseq_destroy(seq);
    if (fp) gzclose(fp);
    radix_sort_sketch_size(lengths.a, lengths.a + lengths.n); /* sort by decreasing length sizes */
#ifndef NDEBUG
    if (lengths.n > 0) for (i = 0; !err && i < lengths.n - 1; ++i) { /* check radix sort output */
        if (kv_A(lengths, i).sketch_size < kv_A(lengths, i + 1).sketch_size) {
            fprintf(stderr, "[sketch_reads] lengths are not stored in decreasing order\n");
            err = ERR_LOGIC;
        }
    }
#endif
    oh = NULL;
    if (!err && !(oh = fopen(tmp_index_filename, "wb"))) err = ERR_FILE;
    cumulative_count = lengths.n;
    if (!err && fwrite(&cumulative_count, sizeof(cumulative_count), 1, oh) != 1) err = ERR_IO;
    cumulative_count = 0;
    for (i = 0; i < lengths.n; ++i) { /* write cumulative lengths */
        cumulative_count += kv_A(lengths, i).sketch_size;
        if (fwrite(&cumulative_count, sizeof(cumulative_count), 1, oh) != 1) err = ERR_IO;
    }
    for (i = 0; i < lengths.n; ++i) { /* write sketches in the order given by lengths (from longest to shortest) */
        size_id_handle *record = &kv_A(lengths, i);
        if (fwrite(record->sketch, sizeof(mm_t), record->sketch_size, oh) != record->sketch_size) err = ERR_IO;
    }
    fflush(oh);
    fclose(oh);
    kv_destroy(lengths);
    kv_destroy(mmzers);
    return err;
}