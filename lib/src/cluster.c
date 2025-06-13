
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include "../../bundled/klib/include/ksort.h"
#include "../../bundled/khashl/include/khashl.h"
#include "../include/exception.h"
#include "../include/minimizer.h"
#include "../include/sketch_reads.h"
#include "../include/cluster.h"

#define cmp(x, y) ((x).minimizers.n > (y).minimizers.n)
KSORT_INIT(cluster, cluster_t, cmp)

typedef kvec_t(uint32_t) clustersidv_t;
KHASHL_MAP_INIT(KH_LOCAL, mm2cluster_t, mm2cls, mm_t, clustersidv_t, kh_hash_uint64, kh_eq_generic)

size_t seek_done(unsigned char const *const flags, const size_t size, const size_t start) {
    size_t i;
    assert(flags);
    for (i = start; i < size && !flags[i]; ++i) {}
    return i;
}

size_t seek_not_done(unsigned char const *const flags, const size_t size, const size_t rstart) {
    size_t i;
    assert(flags);
    for (i = rstart; i != 0 && flags[i]; --i) {}
    return i;
}

size_t my_lower_bound(uint32_t a[], size_t n, uint32_t threshold_value) {
    size_t l, h;
    assert(a);
    l = 0;
    h = n; /* Not n - 1 */
    while (l < h) {
        size_t mid =  l + (h - l) / 2;
        if (threshold_value <= a[mid]) {
            h = mid;
        } else {
            l = mid + 1;
        }
    }
    return l;
}

typedef kvec_t(size_t) hitv_t;

int cluster_reads_fast(
    char const *const index_filename, 
    const double similarity_threshold, 
    const double merge_threhsold,
    clusters_t* clusters
) {
    int err, fd, absent;
    struct stat filestat;
    void *index;
    sketch_metadata_t *len_id;
    mm_t *mm, *mm_itr, minimizer;
    uint64_t nsketches;
    cluster_t empty_cluster;
    size_t i, j, k, best_cluster_idx, best_big_cluster_idx, back_idx, insertion_idx;
    double best_similarity, similarity;
    mmv_t buffer;
    hitv_t hits;
    mm2cluster_t *mm2clusters;
    khint_t map_itr;

    assert(index_filename);
    assert(clusters);

    err = OK;
    kv_init(empty_cluster.minimizers);
    kv_init(empty_cluster.ids);
    kv_init(buffer);
    kv_init(hits);
    mm2clusters = mm2cls_init();

    /* memory map index file */
    if (!err && (fd = open(index_filename, O_RDWR)) < 0) {
        fprintf(stderr, "opening temporary index file failed\n");
        err = ERR_FILE;
    }
    if (!err && fstat(fd, &filestat) != 0) {
        fprintf(stderr, "stat failed\n");
        err = ERR_FILE;
    }
    if (!err && (index = mmap(NULL, filestat.st_size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0)) == MAP_FAILED) err = ERR_FILE;
    /* keep 2 pointers: cumulative sum of lengths and current sketch */
    
    if (!err) {
        nsketches = *((uint64_t*)index);
        len_id = (sketch_metadata_t*)(index + sizeof(nsketches)); /* init iterator to start of metadata */
        mm = (mm_t*)(index + sizeof(nsketches) + nsketches * sizeof(sketch_metadata_t)); /* jump to start of sketches */
        /* load first read as first cluster */
        kv_init(*clusters);
        kv_push(cluster_t, *clusters, empty_cluster);
        kv_reserve(mm_t, kv_A(*clusters, 0).minimizers, len_id->size);
        /* if (!err && memcpy(kv_A(*clusters, 0).minimizers.a, mm, len_id->size * sizeof(mm_t)) != kv_A(*clusters, 0).minimizers.a) err = ERR_RUNTIME; */
        /* if (!err) kv_A(*clusters, 0).minimizers.n = len_id->size; */
        mm_itr = mm;
        for (i = 0; i < len_id->size; ++i) {
            kv_push(mm_t, kv_A(*clusters, 0).minimizers, *mm_itr);
            map_itr = mm2cls_put(mm2clusters, *mm_itr, &absent);
            kv_init(kh_val(mm2clusters, map_itr));
            kv_push(uint32_t, kh_val(mm2clusters, map_itr), 0); /* minimizers in first cluster at index 0 */
            ++mm_itr;
        }
        if (!err) kv_push(read_id_t, kv_A(*clusters, 0).ids, len_id->id);
        mm += len_id->size;
        ++len_id;
    }
    for (i = 1; !err && i < nsketches; ++i) {
        kv_resize(size_t, hits, clusters->n);
        kv_clear(hits, 0);
        mm_itr = mm;
        kv_clear(buffer, 0); /* buffer is to avoid reading multiple times from mapped memory */
        for (j = 0; j < len_id->size; ++j) { /* for each minimizer in read */
            map_itr = mm2cls_get(mm2clusters, *mm_itr);
            if (map_itr < kh_end(mm2clusters)) { /* new minimizers not seen before */
                /* do nothing, add to cluster later */
            } else { /* minimizer already seen before */
                clustersidv_t cluster_ids = kh_val(mm2clusters, map_itr);
                for (k = 0; k < cluster_ids.n; ++k) ++kv_A(hits, kv_A(cluster_ids, k));
            }
            kv_push(mm_t, buffer, *mm_itr);
            ++mm_itr;
        }
        best_cluster_idx = 0;
        if (kv_size(hits) != 0) best_big_cluster_idx = kv_A(buffer, best_cluster_idx); /* reusing variable for tmp computation */
        for (j = 0; j < kv_size(hits); ++j) { /* equivalent to best_cluster_idx = std::max_element(hits.begin(), hits.end()); */
            if (kv_A(hits, j) > best_big_cluster_idx) {
                best_cluster_idx = j;
                best_big_cluster_idx = kv_A(hits, j);
            }
        }
        
        if ((double)kv_A(hits, best_cluster_idx) / len_id->size > similarity_threshold) { /* merge read into best cluster */
            kv_push(read_id_t, kv_A(*clusters, best_cluster_idx).ids, len_id->id); /* add read to cluster */
            for (j = 0; j < kv_size(buffer); ++j) { /* for the minimizers in read */
                minimizer = kv_A(buffer, j);
                map_itr = mm2cls_put(mm2clusters, minimizer, &absent);
                if (absent) { /* new minimizers not seen before */
                    kv_init(kh_val(mm2clusters, map_itr));
                    kv_push(uint32_t, kh_val(mm2clusters, map_itr), best_cluster_idx); /* minimizers in first cluster at index 0 */
                } else { /* minimizer already seen before */
                    clustersidv_t *clsids;
                    clsids = &kh_val(mm2clusters, map_itr);
                    insertion_idx = my_lower_bound(clsids->a, clsids->n, best_cluster_idx); /* O(log2(N)) complexity, if lists are short maybe O(N) is better */
                    if (clsids->a[insertion_idx] != best_cluster_idx) kv_insert(uint32_t, *clsids, insertion_idx, best_cluster_idx); /* insert cluster idx to preserve order and set property */
                }
                /* add minimizers to clusters by preserving set property (no duplicates) and order */
                mmv_t *cluster_mms = &kv_A(*clusters, best_cluster_idx).minimizers;
                insertion_idx = my_lower_bound(cluster_mms->a, cluster_mms->n, minimizer);
                if (cluster_mms->a[insertion_idx] != minimizer) kv_insert(mm_t, *cluster_mms, insertion_idx, minimizer);
            }
        } else { /* new cluster */
            kv_push(cluster_t, *clusters, empty_cluster);
            back_idx = kv_size(*clusters) - 1;
            assert(back_idx != SIZE_MAX);
            kv_reserve(mm_t, kv_A(*clusters, back_idx).minimizers, len_id->size);
            assert(clusters->a[back_idx].minimizers.a);
            if (!err && memcpy(kv_A(*clusters, back_idx).minimizers.a, mm, len_id->size * sizeof(mm_t)) != kv_A(*clusters, back_idx).minimizers.a) err = ERR_RUNTIME;
            if (!err) kv_A(*clusters, back_idx).minimizers.n = len_id->size;
            if (!err) kv_push(read_id_t, kv_A(*clusters, back_idx).ids, len_id->id);

            for (j = 0; j < kv_size(buffer); ++j) { /* for the minimizers in read */
                minimizer = kv_A(buffer, j);
                map_itr = mm2cls_put(mm2clusters, minimizer, &absent);
                if (absent) { /* new minimizers not seen before */
                    kv_init(kh_val(mm2clusters, map_itr));
                    kv_push(uint32_t, kh_val(mm2clusters, map_itr), kv_size(*clusters) - 1); /* minimizers in first cluster at index 0 */
                } else { // minimizer already seen before
                    kv_push(uint32_t, kh_val(mm2clusters, map_itr), kv_size(*clusters) - 1);
                }
            }
        }
        mm += len_id->size;
        ++len_id;
    }
    
    kv_destroy(empty_cluster.minimizers);
    kv_destroy(empty_cluster.ids);
    kv_destroy(hits);
    mm2cls_destroy(mm2clusters);
    if (!err && munmap(index, filestat.st_size) != 0) err = ERR_RUNTIME;
    if (!err && close(fd)) return ERR_FILE;
    mm2clusters = NULL;
    len_id = NULL;
    mm = NULL;

    if (!err && merge_threhsold != 0) {
        unsigned char* done;
        ks_mergesort_cluster(clusters->n, clusters->a, NULL);
        done = NULL;
        if (!err && (done = calloc(clusters->n, 1)) == NULL) err = ERR_MALLOC;
        do {
            best_similarity = merge_threhsold;
            best_big_cluster_idx = best_cluster_idx = SIZE_MAX;
            for (i = 0; i < kv_size(*clusters); ++i) {
                cluster_t big_cluster = kv_A(*clusters, i);
                for (j = i + 1; j < kv_size(*clusters) && !done[j]; ++j) {
                    cluster_t* small_cluster = &kv_A(*clusters, j);
                    if (get_shared_count(big_cluster.minimizers.a, big_cluster.minimizers.n, small_cluster) / small_cluster->minimizers.n >= best_similarity) {
                        best_cluster_idx = j;
                        best_big_cluster_idx = i;
                    }
                }
            }
            if (best_cluster_idx != SIZE_MAX) { /* merge cluster at index best_cluster_idx into cluster at index best_big_cluster_idx */
                assert(best_big_cluster_idx != SIZE_MAX);
                if (!err) err = two_way_merge(kv_A(*clusters, best_cluster_idx).minimizers.a, kv_A(*clusters, best_cluster_idx).minimizers.n, &buffer, &kv_A(*clusters, best_big_cluster_idx));
                if (!err) { /* merge ids */
                    len_t big_size, small_size;
                    big_size = kv_A(*clusters, best_big_cluster_idx).ids.n;
                    small_size = kv_A(*clusters, best_cluster_idx).ids.n;
                    kv_reserve(read_id_t, kv_A(*clusters, best_big_cluster_idx).ids, big_size + small_size);
                    memcpy(kv_A(*clusters, best_big_cluster_idx).ids.a + big_size, kv_A(*clusters, best_cluster_idx).ids.a, small_size * sizeof(read_id_t));
                    kv_A(*clusters, best_big_cluster_idx).ids.n = big_size + small_size;
                    done[best_cluster_idx] = TRUE;
                }
            }
        } while (best_cluster_idx != SIZE_MAX);
        assert(!done[0]); /* due to sorting the first cluster is the bigger one and is always active */
        i = seek_done(done, clusters->n, 0);
        j = seek_not_done(done, clusters->n, clusters->n - 1);
        while(i < j) { /* swap to fill hole */
            cluster_t tmp = clusters->a[i];
            clusters->a[i] = clusters->a[j];
            clusters->a[j] = tmp;
            done[i] = FALSE;
            done[j] = TRUE;
            i = seek_done(done, j, i + 1);
            j = seek_not_done(done, clusters->n, j - 1);
        }
        for (j = i; j < clusters->n; ++j) { /* deallocate merged clusters which are now the tail */
            kv_destroy(kv_A(*clusters, j).minimizers);
            kv_destroy(kv_A(*clusters, j).ids); 
        }
        clusters->n = i; /* resize vector to true clusters */
        if (done) free(done);
    }
    kv_destroy(buffer);
    for (i = 0; i < kv_size(*clusters); ++i) kv_destroy(kv_A(*clusters, i).minimizers); /* deallocate all remaining minimizers */
    return err;
}

/* ----------------------------- slow version ------------------------------ */

double get_shared_count(mm_t const *const sketch, const size_t sketch_size, cluster_t const *const cluster) {
    double shared;
    size_t i, j;
    assert(sketch);
    assert(cluster);
    shared = 0;
    i = j = 0;
    while (i < sketch_size && j < cluster->minimizers.n) {
        assert(cluster->minimizers.a);
        mm_t smm, cmm;
        smm = sketch[i];
        cmm = kv_A(cluster->minimizers, j);
        if (smm == cmm) { /* advance both */
            shared += 1;
            ++i;
            ++j;
        } else if (smm > cmm) {
            ++j;
        } else {
            ++i;
        }
    }
    return shared;
}

int two_way_merge(mm_t const *const sketch, const size_t sketch_size, mmv_t *const buffer, cluster_t* cluster) { 
    int err;
    size_t i, j;
    assert(sketch);
    assert(buffer);
    assert(cluster);
    err = OK;
    buffer->n = 0; /* clear buffer */
    i = j = 0;
    while (i < sketch_size && j < kv_size(cluster->minimizers)) {
        mm_t smm, cmm;
        smm = sketch[i];
        cmm = kv_A(cluster->minimizers, j);
        if (smm == cmm) {
            kv_push(mm_t, *buffer, smm);
            ++i;
            ++j;
        } else if (smm < cmm) {
            kv_push(mm_t, *buffer, smm);
            ++i;
        } else {
            kv_push(mm_t, *buffer, cmm);
            ++j;
        }
    }
    while (i < sketch_size) {
        kv_push(mm_t, *buffer, sketch[i]);
        ++i;
    }
    while (j < kv_size(cluster->minimizers)) {
        kv_push(mm_t, *buffer, kv_A(cluster->minimizers, j));
        ++j;
    }
    kv_reserve(mm_t, cluster->minimizers, kv_size(*buffer));
    if(memcpy(cluster->minimizers.a, buffer->a, buffer->n * sizeof(mm_t)) != cluster->minimizers.a) err = ERR_RUNTIME;
    if (!err) cluster->minimizers.n = buffer->n;
    return err;
}

int cluster_reads_slow(
    char const *const index_filename, 
    const double similarity_threshold, 
    const double merge_threhsold,
    clusters_t* clusters
) {
    int err, fd;
    struct stat filestat;
    void* index;
    sketch_metadata_t* len_id;
    mm_t* mm;
    uint64_t nsketches;
    cluster_t empty_cluster;
    size_t i, j, best_cluster_idx, best_big_cluster_idx, back_idx;
    double best_similarity, similarity;
    mmv_t buffer;
    assert(index_filename);
    assert(clusters);
    err = OK;
    kv_init(empty_cluster.minimizers);
    kv_init(empty_cluster.ids);
    kv_init(buffer);
    /* memory map index file */
    if (!err && (fd = open(index_filename, O_RDWR)) < 0) {
        fprintf(stderr, "opening temporary index file failed\n");
        err = ERR_FILE;
    }
    if (!err && fstat(fd, &filestat) != 0) {
        fprintf(stderr, "stat failed\n");
        err = ERR_FILE;
    }
    if (!err && (index = mmap(NULL, filestat.st_size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0)) == MAP_FAILED) err = ERR_FILE;
    /* keep 2 pointers: cumulative sum of lengths and current sketch */
    if (!err) {
        nsketches = *((uint64_t*)index);
        len_id = (sketch_metadata_t*)(index + sizeof(nsketches)); /* init iterator to start of metadata */
        mm = (mm_t*)(index + sizeof(nsketches) + nsketches * sizeof(sketch_metadata_t)); /* jump to start of sketches */
#ifndef NDEBUG 
        {
            mm_t* mm_itr = mm;
            for (i = 0; i < nsketches; ++i) {
                sketch_metadata_t* ld_itr = len_id + i;
                fprintf(stderr, "[read] metadata %llu, %llu\n", ld_itr->id, ld_itr->size);
                for (j = 0; j < ld_itr->size; ++j) {
                    fprintf(stderr, "%llu, ", *(mm_itr + j));
                }
                fprintf(stderr, "\n");
                mm_itr += ld_itr->size;
            }
        }
#endif
        /* implement isonclust3 algorithm */
        kv_init(*clusters);
        kv_push(cluster_t, *clusters, empty_cluster);
        kv_reserve(mm_t, kv_A(*clusters, 0).minimizers, len_id->size);
        if (!err && memcpy(kv_A(*clusters, 0).minimizers.a, mm, len_id->size * sizeof(mm_t)) != kv_A(*clusters, 0).minimizers.a) err = ERR_RUNTIME;
        if (!err) kv_A(*clusters, 0).minimizers.n = len_id->size;
        if (!err) kv_push(read_id_t, kv_A(*clusters, 0).ids, len_id->id);
        mm += len_id->size;
        ++len_id;
    }
    for (i = 1; !err && i < nsketches; ++i) {
        best_cluster_idx = SIZE_MAX;
        best_similarity = similarity_threshold;
        for (j = 0; j < kv_size(*clusters); ++j) {
            /* fprintf(stderr, "%llu, %llu, %llu, %llu\n", i, j, len_id->id, len_id->size); */
            if ((similarity = get_shared_count(mm, len_id->size, &kv_A(*clusters, j)) / len_id->size) > best_similarity) {
                best_similarity = similarity;
                best_cluster_idx = j;
            }
        }
        assert(len_id->id < nsketches);
        if (best_cluster_idx != SIZE_MAX) { /* merge read and cluster */
            if (!err) err = two_way_merge(mm, len_id->size, &buffer, &kv_A(*clusters, best_cluster_idx));
            if (!err) kv_push(read_id_t, kv_A(*clusters, best_cluster_idx).ids, len_id->id);
        } else { /* append read as new cluster */
            kv_push(cluster_t, *clusters, empty_cluster);
            back_idx = kv_size(*clusters) - 1;
            assert(back_idx != SIZE_MAX);
            kv_reserve(mm_t, kv_A(*clusters, back_idx).minimizers, len_id->size);
            assert(clusters->a[back_idx].minimizers.a);
            if (!err && memcpy(kv_A(*clusters, back_idx).minimizers.a, mm, len_id->size * sizeof(mm_t)) != kv_A(*clusters, back_idx).minimizers.a) err = ERR_RUNTIME;
            if (!err) kv_A(*clusters, back_idx).minimizers.n = len_id->size;
            if (!err) kv_push(read_id_t, kv_A(*clusters, back_idx).ids, len_id->id);
        }
        mm += len_id->size;
        ++len_id;
    }
    len_id = NULL;
    mm = NULL;
    kv_destroy(empty_cluster.minimizers);
    kv_destroy(empty_cluster.ids);
    if (!err && munmap(index, filestat.st_size) != 0) err = ERR_RUNTIME;
    if (!err && close(fd)) return ERR_FILE;
    if (!err && merge_threhsold != 0) {
        unsigned char* done;
        ks_mergesort_cluster(clusters->n, clusters->a, NULL);
        done = NULL;
        if (!err && (done = calloc(clusters->n, 1)) == NULL) err = ERR_MALLOC;
        do {
            best_similarity = merge_threhsold;
            best_big_cluster_idx = best_cluster_idx = SIZE_MAX;
            for (i = 0; i < kv_size(*clusters); ++i) {
                cluster_t big_cluster = kv_A(*clusters, i);
                for (j = i + 1; j < kv_size(*clusters) && !done[j]; ++j) {
                    cluster_t* small_cluster = &kv_A(*clusters, j);
                    if (get_shared_count(big_cluster.minimizers.a, big_cluster.minimizers.n, small_cluster) / small_cluster->minimizers.n >= best_similarity) {
                        best_cluster_idx = j;
                        best_big_cluster_idx = i;
                    }
                }
            }
            if (best_cluster_idx != SIZE_MAX) { /* merge cluster at index best_cluster_idx into cluster at index best_big_cluster_idx */
                assert(best_big_cluster_idx != SIZE_MAX);
                if (!err) err = two_way_merge(kv_A(*clusters, best_cluster_idx).minimizers.a, kv_A(*clusters, best_cluster_idx).minimizers.n, &buffer, &kv_A(*clusters, best_big_cluster_idx));
                if (!err) { /* merge ids */
                    len_t big_size, small_size;
                    big_size = kv_A(*clusters, best_big_cluster_idx).ids.n;
                    small_size = kv_A(*clusters, best_cluster_idx).ids.n;
                    kv_reserve(read_id_t, kv_A(*clusters, best_big_cluster_idx).ids, big_size + small_size);
                    memcpy(kv_A(*clusters, best_big_cluster_idx).ids.a + big_size, kv_A(*clusters, best_cluster_idx).ids.a, small_size * sizeof(read_id_t));
                    kv_A(*clusters, best_big_cluster_idx).ids.n = big_size + small_size;
                    done[best_cluster_idx] = TRUE;
                }
            }
        } while (best_cluster_idx != SIZE_MAX);
        assert(!done[0]); /* due to sorting the first cluster is the bigger one and is always active */
        i = seek_done(done, clusters->n, 0);
        j = seek_not_done(done, clusters->n, clusters->n - 1);
        while(i < j) { /* swap to fill hole */
            cluster_t tmp = clusters->a[i];
            clusters->a[i] = clusters->a[j];
            clusters->a[j] = tmp;
            done[i] = FALSE;
            done[j] = TRUE;
            i = seek_done(done, j, i + 1);
            j = seek_not_done(done, clusters->n, j - 1);
        }
        for (j = i; j < clusters->n; ++j) { /* deallocate merged clusters which are now the tail */
            kv_destroy(kv_A(*clusters, j).minimizers);
            kv_destroy(kv_A(*clusters, j).ids); 
        }
        clusters->n = i; /* resize vector to true clusters */
        if (done) free(done);
    }
    kv_destroy(buffer);
    for (i = 0; i < kv_size(*clusters); ++i) kv_destroy(kv_A(*clusters, i).minimizers); /* deallocate all remaining minimizers */
    return err;
}

/*
    FIXME: this function writes a csv file only mapping reads ids to clusters. 
    Extend it to split the original Fastq into the actual clusters or write/re-use a separate tool for that.
*/
int cluster_save(
    clusters_t const *const clusters, 
    char const *const output_filename
) {
    int err;
    FILE *comma_output;
    size_t i, j;
    /* kvec_t(char) filename; */
    assert(clusters);
    err = OK;
    /* 
    kv_init(filename);
    of_len = strlen(output_folder);
    i = of_len - 1;
    while (i < of_len && output_folder[i] == '/') { --of_len; }
    kv_reserve(char, filename, of_len + 17);
    if (!err && !filename.a) err = ERR_MALLOC;
    if (!err && memcpy(filename.a, output_folder, of_len) != filename.a) err = ERR_RUNTIME;
    if (!err && memcpy(filename.a + of_len, "/cluster_ids.csv", 17) != filename.a + of_len) err = ERR_RUNTIME;
    if (!err && (comma_output = fopen(filename.a, "w")) == NULL) err = ERR_FILE;
    */
    if (output_filename) {
        if ((comma_output = fopen(output_filename, "w")) == NULL) err = ERR_FILE;
    } else {
        comma_output = stdout;
    } 
    if (!err) {
        fprintf(comma_output, "read ID,cluster ID\n");
        for (i = 0; i < clusters->n; ++i) {
            for (j = 0; j < clusters->a[i].ids.n; ++j) {
                fprintf(comma_output, "%llu,%lu\n", clusters->a[i].ids.a[j], i);
            }
        }
    }
    if (comma_output) {
        fflush(comma_output);
        if (comma_output != stdout) fclose(comma_output);
    }
    /* kv_destroy(filename); */
    return err;
}
