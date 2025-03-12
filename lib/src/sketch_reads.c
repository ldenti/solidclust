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
    if (!err) for (i = 0; i < seq_len; ++i) { /* FIXME find a numerically stable algorithm, rn if kmer_quality goes to 0, then everything is 0 */
        q = buffer[buf_idx];
        buffer[buf_idx] = 1.0 - pow(10.0, -((int)qual[i] - 33) / 10);
        kmer_quality /= q; /* no need to check if i >= k since buffer is all 1s at the beginning */
        kmer_quality *= buffer[buf_idx];
        /* fprintf(stderr, "%d\n", kmer_quality); */
        if ((i + 1) >= k) kv_push(unsigned char, *filter, kmer_quality > threshold);
        buf_idx = (buf_idx + 1) % k;
    }
    fprintf(stderr, "seq len = %llu, k = %llu, filter length = %llu\n", seq_len, k, filter->n);
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
    size_t i, read_id;
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
        if (!err) err = minimizer_from_string(seq->seq.s, quality_filter.a, seq->seq.l, k, w, canonical, seed, &mmzers); /* this appends new minimizers to end of mmzers */
        handle.metadata.id = read_id;
        handle.metadata.size = mmzers.n - cumulative_count;
        handle.sketch_offset = cumulative_count;
        fprintf(stderr, "found %llu new minimzers from read %llu\n", handle.metadata.size, read_id);
        if (!err && handle.metadata.size > 0) ks_introsort(minimizer, handle.metadata.size, mmzers.a + cumulative_count);
        kv_push(size_id_handle, lengths, handle);
        cumulative_count = mmzers.n;
        ++read_id;
    }
    kv_destroy(quality_filter);
    if (seq) kseq_destroy(seq);
    if (fp) gzclose(fp);
    assert(lengths.n == read_id);
    radix_sort_sketch_size(lengths.a, lengths.a + lengths.n); /* sort by increasing length sizes */
#ifndef NDEBUG
    if (lengths.n > 0) for (i = 0; !err && i < lengths.n - 1; ++i) { /* check radix sort output */
        if (kv_A(lengths, i).metadata.size > kv_A(lengths, i + 1).metadata.size) {
            fprintf(stderr, "[sketch_reads] lengths are not stored in increasing order\n");
            err = ERR_LOGIC;
        }
    }
#endif
    oh = NULL;
    if (!err && !(oh = fopen(tmp_index_filename, "wb"))) err = ERR_FILE;
    cumulative_count = lengths.n;
    if (!err && fwrite(&cumulative_count, sizeof(cumulative_count), 1, oh) != 1) err = ERR_IO;
    cumulative_count = 0;
    for (i = lengths.n - 1; !err && i != SIZE_MAX; --i) { /* write cumulative lengths */
        if (fwrite(&lengths.a[i].metadata, sizeof(sketch_metadata_t), 1, oh) != 1) err = ERR_IO;
    }
    for (i = lengths.n - 1; !err && i != SIZE_MAX; --i) { /* write sketches in the order given by lengths (from longest to shortest) */
        size_id_handle *record = &kv_A(lengths, i); /* DO NOT merge the two for loops since actual sketches come after their lengths */
        if (fwrite(mmzers.a + record->sketch_offset, sizeof(mm_t), record->metadata.size, oh) != record->metadata.size) err = ERR_IO;
    }
    if (oh) {
        fflush(oh);
        fclose(oh);
    }
    kv_destroy(lengths);
    kv_destroy(mmzers);
    return err;
}

/*
0, 4, 0x105101ea0 = 4379909792
1, 3, 0x1050013c0 = 4378858432
4, 0, 0x1050013d8
3, 0, 0x1050013d8
2, 0, 0x1050013d8
*/