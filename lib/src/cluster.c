
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
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

#include <assert.h>

#define cmp(x, y) ((x).minimizers.n > (y).minimizers.n)
KSORT_INIT(cluster, cluster_t, cmp)

typedef uint32_t clid_t;
typedef kvec_t(clid_t) clustersidv_t;
KHASHL_MAP_INIT(KH_LOCAL, mm2cluster_t, mm2cls, mm_t, clustersidv_t, kh_hash_uint64, kh_eq_generic)

typedef struct {
    clid_t cluster_id;
    uint32_t count;
} idw_t;

typedef kvec_t(idw_t) clustersidwv_t;
KHASHL_MAP_INIT(KH_LOCAL, mm2cluster_weighted_t, mm2clsw, mm_t, clustersidwv_t, kh_hash_uint64, kh_eq_generic)

typedef kvec_t(size_t) hitv_t;

size_t seek_done(unsigned char const *const flags, const size_t size, const size_t start);
size_t seek_not_done(unsigned char const *const flags, const size_t size, const size_t rstart);
double get_shared_count(mm_t const *const sketch, const size_t sketch_size, cluster_t const *const cluster);
int two_way_merge(mm_t const *const sketch, const size_t sketch_size, mmv_t *const buffer, cluster_t* cluster);
size_t clsid_lower_bound(clid_t a[], size_t n, clid_t threshold_value);
size_t clsid_lower_bound_weighted(idw_t a[], size_t n, clid_t threshold_value);
size_t minimizer_lower_bound(mm_t a[], size_t n, mm_t threshold_value);

int cluster_reads(
    char const *const index_filename, 
    const double similarity_threshold, 
    clusters_t *const clusters
) {
    int err, fd, absent;
    struct stat filestat;
    void *index;
    sketch_metadata_t *len_id;
    mm_t *mm, *mm_itr, minimizer;
    uint64_t nsketches;
    cluster_t empty_cluster;
    size_t i, j, k, best_value, back_idx, insertion_idx;
    clid_t best_cluster_idx;
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
            kv_push(clid_t, kh_val(mm2clusters, map_itr), 0); /* minimizers in first cluster at index 0 */
            ++mm_itr;
        }
        if (!err) kv_push(read_id_t, kv_A(*clusters, 0).ids, len_id->id);
        mm += len_id->size;
        ++len_id;
    }
    for (i = 1; !err && i < nsketches; ++i) {
        if (kv_capacity(hits) == clusters->n) kv_reserve(size_t, hits, 2*clusters->n);
        kv_resize(size_t, hits, clusters->n);
        assert(kv_size(hits) != 0);
        kv_reset(hits, 0);
        mm_itr = mm;
        kv_clear(buffer); /* buffer is to avoid reading multiple times from mapped memory */
        for (j = 0; j < len_id->size; ++j) { /* for each minimizer in read */
            map_itr = mm2cls_get(mm2clusters, *mm_itr);
            if (map_itr < kh_end(mm2clusters)) { /* minimizer already seen before */
                clustersidv_t cluster_ids = kh_val(mm2clusters, map_itr);
                for (k = 0; k < cluster_ids.n; ++k) ++kv_A(hits, kv_A(cluster_ids, k));
            } /* else, if new minimizers not seen before do nothing, add to cluster later */
            kv_push(mm_t, buffer, *mm_itr);
            ++mm_itr;
        }
        best_cluster_idx = 0;
        best_value = kv_A(hits, best_cluster_idx); /* reusing variable for tmp computation */
        for (j = 0; j < kv_size(hits); ++j) { /* equivalent to best_cluster_idx = std::max_element(hits.begin(), hits.end()); */
            if (kv_A(hits, j) > best_value) {
                best_cluster_idx = j;
                best_value = kv_A(hits, j);
            }
        }
        /* fprintf(stderr, "hits.size = %llu, best cluster id = %llu\n", kv_size(hits), best_cluster_idx); */
        if ((double)kv_A(hits, best_cluster_idx) / len_id->size > similarity_threshold) { /* merge read into best cluster */
            kv_push(read_id_t, kv_A(*clusters, best_cluster_idx).ids, len_id->id); /* add read to cluster */
            for (j = 0; j < kv_size(buffer); ++j) { /* for the minimizers in read */
                minimizer = kv_A(buffer, j);
                map_itr = mm2cls_put(mm2clusters, minimizer, &absent);
                if (absent) { /* new minimizer not seen before */
                    kv_init(kh_val(mm2clusters, map_itr));
                    kv_push(clid_t, kh_val(mm2clusters, map_itr), best_cluster_idx);
                } else { 
                    clustersidv_t clsids;
                    clsids = kh_val(mm2clusters, map_itr); /* mm2clusters is a packed struct so pointers to its members can be unaligned */
                    insertion_idx = clsid_lower_bound(clsids.a, clsids.n, best_cluster_idx); /* O(log2(N)) complexity, if lists are short maybe O(N) is better */
                    if (insertion_idx == kv_size(clsids) || clsids.a[insertion_idx] != best_cluster_idx) kv_insert(clid_t, clsids, insertion_idx, best_cluster_idx); /* insert cluster idx to preserve order and set property */
                    kh_val(mm2clusters, map_itr) = clsids; /* push changes due to possible reallocations of clsids's buffer */
                }
                /* add minimizers to clusters by preserving set property (no duplicates) and order */
                mmv_t *cluster_mms = &kv_A(*clusters, best_cluster_idx).minimizers;
                insertion_idx = minimizer_lower_bound(cluster_mms->a, cluster_mms->n, minimizer);
                if (insertion_idx == kv_size(*cluster_mms) || cluster_mms->a[insertion_idx] != minimizer) kv_insert(mm_t, *cluster_mms, insertion_idx, minimizer);
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
                }
                kv_push(clid_t, kh_val(mm2clusters, map_itr), kv_size(*clusters) - 1);
            }
        }
        mm += len_id->size;
        ++len_id;
    }
    
    kv_destroy(empty_cluster.minimizers);
    kv_destroy(empty_cluster.ids);
    kv_destroy(hits);
    kv_destroy(buffer);
    kh_foreach(mm2clusters, map_itr) kv_destroy(kh_val(mm2clusters, map_itr));
    mm2cls_destroy(mm2clusters);
    mm2clusters = NULL;
    if (!err && munmap(index, filestat.st_size) != 0) err = ERR_RUNTIME;
    index = NULL;
    len_id = NULL;
    mm = NULL;
    mm_itr = NULL;
    if (!err && close(fd)) err = ERR_FILE;
    return err;
}

typedef struct {
    clid_t first;
    clid_t second;
} clid_pair_t;
typedef kvec_t(clid_pair_t) clid_pairv_t;
KHASHL_SET_INIT(KH_LOCAL, cls_set_t, clss, clid_t, kh_hash_uint32, kh_eq_generic)
KHASHL_MAP_INIT(KH_LOCAL, cls_map_t, clsm, clid_t, mmv_t, kh_hash_uint32, kh_eq_generic)

int generate_cluster_merging(
    mm2cluster_t *const cluster_map,
    cls_map_t *const cl_set_map
) {
    int err, absent;
    size_t i;
    khint_t map_itr, id_itr;
    clustersidv_t ids;

    assert(cluster_map);
    assert(cl_set_map);
    
    err = OK;
    kh_foreach(cluster_map, map_itr) {
        ids = kh_val(cluster_map, map_itr);
        for (i = 0; i < kv_size(ids); ++i) {
            id_itr = clsm_put(cl_set_map, kv_A(ids, i), &absent);
            if (absent) {
                kv_init(kh_val(cl_set_map, id_itr));
            }
            kv_push(mm_t, kh_val(cl_set_map, id_itr), kh_key(cluster_map, map_itr));
        }
    }
    return err;
}

int detect_overlaps(
    const double merge_threhsold,
    mm2cluster_t *const cluster_map,
    cls_map_t *const cl_set_map,
    clid_pairv_t *const merge_into,
    cls_set_t *const small_hs,
    const hitv_t shared_seed_counts /* pre-allocated by parent function: clusters.size() */
) {
    int err, absent;
    unsigned char contains;
    size_t i, j;
    khint_t set_itr, map_itr, itr, pos;
    size_t max_cluster_id, max_count;

    assert(merge_threhsold);
    assert(cluster_map);
    assert(cl_set_map);
    assert(merge_into);
    assert(small_hs);

    err = OK;
    /* FIXME sort cl_set_map to access shared_seed_counts sequentially */
    kh_foreach(cl_set_map, set_itr) {
        mmv_t hashes = kh_val(cl_set_map, set_itr);
        for (i = 0; i < kv_size(hashes); ++i) {
            map_itr = mm2cls_get(cluster_map, kv_A(hashes, i));
            if (map_itr != kh_end(cluster_map)) {
                clustersidv_t belongs_to = kh_val(cluster_map, map_itr);
                for (j = 0; j < kv_size(belongs_to); ++j) {
                    if (kv_A(belongs_to, j) != kh_key(cl_set_map, set_itr)) {
                        ++kv_A(shared_seed_counts, kv_A(belongs_to, j));
                    }
                }
            }
        }
        if (kv_size(shared_seed_counts)) {
            double shared_perc;
            max_cluster_id = 0;
            max_count = kv_A(shared_seed_counts, max_cluster_id);
            for (i = 0; i < kv_size(shared_seed_counts); ++i) {
                if (kv_A(shared_seed_counts, i) > max_count) {
                    max_cluster_id = i;
                    max_count = kv_A(shared_seed_counts, i);
                }
            }
            shared_perc = (double)kv_size(hashes) / max_count;
            pos = clss_get(small_hs, max_cluster_id);
            if (shared_perc > merge_threhsold && pos == kh_end(small_hs)) {
                clid_pair_t dummy;
                assert(max_cluster_id < UINT32_MAX);
                dummy.first = kh_key(cl_set_map, set_itr);
                itr = clsm_get(cl_set_map, max_cluster_id);
                contains = FALSE;
                for (i = 0; !contains && i < kv_size(*merge_into); ++i) {
                    if (kv_A(*merge_into, i).first == dummy.first && kv_A(*merge_into, i).second == max_cluster_id) {
                        contains = TRUE;
                    }
                }
                if (!contains) {
                    assert(kv_size(hashes) <= kv_size(kh_val(cl_set_map, itr)));
                    if (kv_size(hashes) < kv_size(kh_val(cl_set_map, itr)) || dummy.first < max_cluster_id) {
                        dummy.second = max_cluster_id;
                        kv_push(clid_pair_t, *merge_into, dummy);
                        clss_put(small_hs, dummy.first, &absent);
                    }
                }
            }
        }
        for (i = 0; i < kv_size(shared_seed_counts); ++i) { /* reinit count vector */
            kv_A(shared_seed_counts, i) = 0;
        }
    }
    return err;
}

int apply_cluster_merging(
    clusters_t *const clusters,
    mm2cluster_t *const cluster_map,
    cls_map_t *const cl_set_map,
    clid_pairv_t *const merge_into,
    cls_set_t *const small_hs,
    cls_set_t *const not_large
) {
    int err, absent;
    size_t i, j, k;
    khint_t itr, clsm_map_itr, mm2cls_map_itr;
    kvec_t(unsigned char) available;
    mmv_t mm_vec;
    clustersidv_t cl_vec;;

    assert(clusters);
    assert(cluster_map);
    assert(merge_into);
    assert(small_hs);
    assert(not_large);

    err = OK;
    kv_init(available);
    kv_resize(unsigned char, available, kv_size(*clusters));
    for (i = 0; i < kv_size(available); ++i) kv_A(available, i) = TRUE;
    for (i = 0; i < kv_size(*merge_into); ++i) {
        clid_t target_clid;
        target_clid = kv_A(*merge_into, i).second;
        if (kv_A(available, target_clid)) {
            itr = clss_get(small_hs, target_clid); 
            if (itr == kh_end(small_hs)) {
                clsm_map_itr = clsm_get(cl_set_map, i);
                mm_vec = kh_val(cl_set_map, clsm_map_itr);
                for (j = 0; j < kv_size(mm_vec); ++j) {
                    mm_t mm;
                    mm = kv_A(mm_vec, j);
                    mm2cls_map_itr = mm2cls_get(cluster_map, mm);
                    cl_vec = kh_val(cluster_map, mm2cls_map_itr);
                    absent = TRUE;
                    for (k = 0; absent && k < kv_size(cl_vec); ++k) {
                        if (kv_A(cl_vec, k) == target_clid) {
                            absent = FALSE;
                        }
                    }
                    if (absent) {
                        kv_push(clid_t, cl_vec, target_clid);
                        kh_val(cluster_map, mm2cls_map_itr) = cl_vec;
                    }
                }
            } else {
                clss_put(not_large, target_clid, &absent);
            }
        }
    }
    return err;
}

int cluster_merging(
    const double merge_threhsold,
    clusters_t *const clusters,
    mm2cluster_t *const cluster_map
) {
    int err;
    size_t i, j;
    cls_map_t *cl_set_map;
    clid_pairv_t merge_into;
    cls_set_t *small_hs;
    cls_set_t *not_large;
    khint_t itr;
    hitv_t shared_seed_counts;
    
    assert(clusters);

    err = OK;
    if (merge_threhsold == 0.0) return err;
    kv_init(merge_into);
    cl_set_map = clsm_init();
    small_hs = clss_init();
    not_large = clss_init();
    kv_init(shared_seed_counts);
    kv_resize(size_t, shared_seed_counts, kv_size(*clusters));
    for (i = 0; i < kv_size(shared_seed_counts); ++i) kv_A(shared_seed_counts, i) = 0;
    do {
        kv_clear(merge_into);
        clss_clear(small_hs);
        if (!err) err = generate_cluster_merging(cluster_map, cl_set_map);
        if (!err) err = detect_overlaps(merge_threhsold, cluster_map, cl_set_map, &merge_into, small_hs, shared_seed_counts);
        if (!err) err = apply_cluster_merging(clusters, cluster_map, cl_set_map, &merge_into, small_hs, not_large);
        for (i = 0, j = 0; i < kv_size(merge_into); ++i) {
            itr = clss_get(not_large, kv_A(merge_into, i).second);
            if (itr == kh_end(not_large)) {
                kv_A(merge_into, j) = kv_A(merge_into, i);
                ++j;
            }
        }
        kv_resize(clid_pair_t, merge_into, j);
    } while (!err && !kv_size(merge_into));

    kv_destroy(merge_into);
    clsm_destroy(cl_set_map);
    clss_destroy(small_hs);
    clss_destroy(not_large);
    kv_destroy(shared_seed_counts);
    return err;
}

/*---------------------------- Weighted variant --------------------------------*/

int cluster_reads_weighted(
    char const *const index_filename, 
    const double similarity_threshold, 
    clusters_t *const clusters
) {
    typedef kvec_t(size_t) hitv_t;
    int err, fd, absent;
    struct stat filestat;
    void *index;
    sketch_metadata_t *len_id;
    mm_t *mm, *mm_itr, minimizer;
    uint64_t nsketches;
    cluster_t empty_cluster;
    size_t i, j, k, best_cluster_idx, best_value, back_idx, insertion_idx;
    mmv_t buffer;
    hitv_t hits;
    mm2cluster_weighted_t *mm2clusters;
    khint_t map_itr;
    idw_t dummy;

    assert(index_filename);
    assert(clusters);

    err = OK;
    kv_init(empty_cluster.minimizers);
    kv_init(empty_cluster.ids);
    kv_init(buffer);
    kv_init(hits);
    mm2clusters = mm2clsw_init();

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
            map_itr = mm2clsw_put(mm2clusters, *mm_itr, &absent);
            kv_init(kh_val(mm2clusters, map_itr));
            dummy.cluster_id = 0;
            dummy.count = 1;
            kv_push(idw_t, kh_val(mm2clusters, map_itr), dummy); /* minimizers in first cluster at index 0 */
            ++mm_itr;
        }
        if (!err) kv_push(read_id_t, kv_A(*clusters, 0).ids, len_id->id);
        mm += len_id->size;
        ++len_id;
    }
    for (i = 1; !err && i < nsketches; ++i) {
        size_t read_size_weighted;
        if (kv_capacity(hits) == clusters->n) kv_reserve(size_t, hits, 2*clusters->n);
        kv_resize(size_t, hits, clusters->n);
        assert(kv_size(hits) != 0);
        kv_reset(hits, 0);
        mm_itr = mm;
        kv_clear(buffer); /* buffer is to avoid reading multiple times from mapped memory */
        read_size_weighted = 0;
        for (j = 0; j < len_id->size; ++j) { /* for each minimizer in read */
            map_itr = mm2clsw_get(mm2clusters, *mm_itr);
            if (map_itr < kh_end(mm2clusters)) { /* minimizer already seen before */
                clustersidwv_t cluster_ids = kh_val(mm2clusters, map_itr);
                for (k = 0; k < cluster_ids.n; ++k) {
                    kv_A(hits, kv_A(cluster_ids, k).cluster_id) += kv_A(cluster_ids, k).count;
                    read_size_weighted += kv_A(cluster_ids, k).count;
                }
            } else {/* if new minimizers not seen before do nothing, just add 1 to size */
                ++read_size_weighted;
            }
            kv_push(mm_t, buffer, *mm_itr);
            ++mm_itr;
        }
        best_cluster_idx = 0;
        best_value = kv_A(hits, best_cluster_idx); /* reusing variable for tmp computation */
        for (j = 0; j < kv_size(hits); ++j) { /* equivalent to best_cluster_idx = std::max_element(hits.begin(), hits.end()); */
            if (kv_A(hits, j) > best_value) {
                best_cluster_idx = j;
                best_value = kv_A(hits, j);
            }
        }
        /* fprintf(stderr, "hits.size = %llu, best cluster id = %llu\n", kv_size(hits), best_cluster_idx); */
        if ((double)kv_A(hits, best_cluster_idx) / read_size_weighted > similarity_threshold) { /* merge read into best cluster */
            kv_push(read_id_t, kv_A(*clusters, best_cluster_idx).ids, len_id->id); /* add read to cluster */
            for (j = 0; j < kv_size(buffer); ++j) { /* for the minimizers in read */
                minimizer = kv_A(buffer, j);
                map_itr = mm2clsw_put(mm2clusters, minimizer, &absent);
                dummy.cluster_id = best_cluster_idx;
                if (absent) { /* new minimizer not seen before */
                    kv_init(kh_val(mm2clusters, map_itr));
                    dummy.count = 1;
                    kv_push(idw_t, kh_val(mm2clusters, map_itr), dummy);
                } else { 
                    clustersidwv_t clsids;
                    clsids = kh_val(mm2clusters, map_itr); /* mm2clusters is a packed struct so pointers to its members can be unaligned */
                    insertion_idx = clsid_lower_bound_weighted(clsids.a, clsids.n, best_cluster_idx); /* O(log2(N)) complexity, if lists are short maybe O(N) is better */
                    if (insertion_idx == kv_size(clsids) || clsids.a[insertion_idx].cluster_id != best_cluster_idx) {
                        dummy.count = 1;
                        kv_insert(idw_t, clsids, insertion_idx, dummy); /* insert cluster idx to preserve order and set property */
                    } else {
                        ++kv_A(clsids, insertion_idx).count; /* update count associated to minimizer (which is already in the map) */
                    }
                    kh_val(mm2clusters, map_itr) = clsids; /* push changes due to possible reallocations of clsids's buffer */
                }
                /* add minimizers to clusters by preserving set property (no duplicates) and order */
                mmv_t *cluster_mms = &kv_A(*clusters, best_cluster_idx).minimizers;
                insertion_idx = minimizer_lower_bound(cluster_mms->a, cluster_mms->n, minimizer);
                if (insertion_idx == kv_size(*cluster_mms) || cluster_mms->a[insertion_idx] != minimizer) kv_insert(mm_t, *cluster_mms, insertion_idx, minimizer);
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
                map_itr = mm2clsw_put(mm2clusters, minimizer, &absent);
                if (absent) { /* new minimizers not seen before */
                    kv_init(kh_val(mm2clusters, map_itr));
                }
                dummy.cluster_id = kv_size(*clusters) - 1;
                dummy.count = 1; 
                kv_push(idw_t, kh_val(mm2clusters, map_itr), dummy);
            }
        }
        mm += len_id->size;
        ++len_id;
    }
    
    kv_destroy(empty_cluster.minimizers);
    kv_destroy(empty_cluster.ids);
    kv_destroy(hits);
    kv_destroy(buffer);
    kh_foreach(mm2clusters, map_itr) kv_destroy(kh_val(mm2clusters, map_itr));
    mm2clsw_destroy(mm2clusters);
    mm2clusters = NULL;
    if (!err && munmap(index, filestat.st_size) != 0) err = ERR_RUNTIME;
    index = NULL;
    len_id = NULL;
    mm = NULL;
    mm_itr = NULL;
    if (!err && close(fd)) err = ERR_FILE;
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

int cluster_print(clusters_t const *const clusters) {
    size_t i, j;
    for (i = 0; i < clusters->n; ++i) {
        fprintf(stderr, "cluster[%llu] : [", i);
        if (clusters->a[i].ids.n == 0) return ERR_LOGIC;
        for (j = 0; j < clusters->a[i].ids.n - 1; ++j) {
            fprintf(stderr, "%llu, ", clusters->a[i].ids.a[j]);
        }
        fprintf(stderr, "%llu]\n", clusters->a[i].ids.a[clusters->a[i].ids.n - 1]);
    }
    return OK;
}

/*--------------------------------- private --------------------------------*/

size_t clsid_lower_bound(clid_t a[], size_t n, clid_t threshold_value) {
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

size_t minimizer_lower_bound(mm_t a[], size_t n, mm_t threshold_value) {
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

size_t clsid_lower_bound_weighted(idw_t a[], size_t n, clid_t threshold_value) {
    size_t l, h;
    assert(a);
    l = 0;
    h = n; /* Not n - 1 */
    while (l < h) {
        size_t mid =  l + (h - l) / 2;
        if (threshold_value <= a[mid].cluster_id) {
            h = mid;
        } else {
            l = mid + 1;
        }
    }
    return l;
}

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

/* ----------------------------- slow version ------------------------------ */

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

int cluster_postprocessing(
    const double merge_threhsold,
    clusters_t *const clusters
) {
    int err;
    double best_similarity;
    unsigned char* done;
    size_t i, j, best_cluster_idx, best_big_cluster_idx;
    mmv_t buffer;

    assert(clusters);

    err = OK;
    if (merge_threhsold != 0) {
        kv_init(buffer);
        ks_mergesort_cluster(clusters->n, clusters->a, NULL);
        done = NULL;
        if (!err && (done = calloc(clusters->n, 1)) == NULL) err = ERR_MALLOC;
        if (!err) do {
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
        } while (!err && best_cluster_idx != SIZE_MAX);
        assert(!done[0]); /* due to sorting the first cluster is the bigger one and is always active */
        i = seek_done(done, clusters->n, 0);
        j = seek_not_done(done, clusters->n, clusters->n - 1);
        while(!err && i < j) { /* swap to fill hole */
            cluster_t tmp = clusters->a[i];
            clusters->a[i] = clusters->a[j];
            clusters->a[j] = tmp;
            done[i] = FALSE;
            done[j] = TRUE;
            i = seek_done(done, j, i + 1);
            j = seek_not_done(done, clusters->n, j - 1);
        }
        for (j = i; !err && j < clusters->n; ++j) { /* deallocate merged clusters which are now the tail */
            kv_destroy(kv_A(*clusters, j).minimizers);
            kv_destroy(kv_A(*clusters, j).ids); 
        }
        clusters->n = i; /* resize vector to true clusters */
        if (done) free(done);
        kv_destroy(buffer);
    }
    for (i = 0; i < kv_size(*clusters); ++i) kv_destroy(kv_A(*clusters, i).minimizers); /* deallocate all remaining minimizers */
    return err;
}
