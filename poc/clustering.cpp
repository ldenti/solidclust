#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/mman.h>

#include <algorithm>
#include <iostream>
#include <vector>
#include <unordered_map>

#include "../lib/include/sketch_reads.h"
#include "../lib/include/minimizer.h"

#include <assert.h>

struct cluster_t {
    cluster_t() {}
    std::vector<read_id_t> ids;
    std::vector<mm_t> minimizers;
};

int cluster_save(
    const std::vector<cluster_t>& clusters, 
    char const *const output_filename
);

int cluster_print(const std::vector<cluster_t>& clusters);

int main(int argc, char *argv[])
{
    // argv[1] = index file (containing sketched reads)
    // argv[2] = similarity threshold <double>
    // argv[3] = merge threshold <double>
    // argv[4] = output filename (csv format)
    if (argc != 5) {
        std::cerr << "Options:\n";
        std::cerr << "see inside source code\n";
        return 1;
    }
    std::string index_filename = argv[1];
    double similarity_threshold = std::strtod(argv[2], NULL);
    double merge_threshold = std::strtod(argv[3], NULL);
    std::string output_filename = argv[4];
    std::vector<cluster_t> clusters;
    { // >>> memory management
    int fd = -1;
    struct stat filestat;
    void *index;

    if ((fd = open(index_filename.c_str(), O_RDWR)) < 0) {
        fprintf(stderr, "opening temporary index file failed\n");
        return 2;
    }
    if (fstat(fd, &filestat) != 0) {
        fprintf(stderr, "stat failed\n");
        return 2;
    }
    if ((index = mmap(NULL, filestat.st_size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0)) == MAP_FAILED) return 2;

    std::unordered_map<mm_t, std::vector<uint32_t>> mm2clusters; // map minimizers to list of clusters they appear in so far
    
    std::size_t nsketches = *((uint64_t*)index);
    sketch_metadata_t *len_id = (sketch_metadata_t*)((char*)index + sizeof(nsketches)); /* init iterator to start of metadata */
    mm_t *mm = (mm_t*)((char*)index + sizeof(nsketches) + nsketches * sizeof(sketch_metadata_t)); /* jump to start of sketches */
    /* implement isonclust3 algorithm */
    clusters.push_back(cluster_t()); // kv_push(cluster_t, *clusters, empty_cluster);
    clusters.back().minimizers.reserve(len_id->size); // kv_reserve(mm_t, kv_A(*clusters, 0).minimizers, len_id->size);
    auto mm_itr = mm;
    for (std::size_t i = 0; i < len_id->size; ++i) {
        clusters.back().minimizers.push_back(*mm_itr);
        mm2clusters[*mm_itr] = std::vector<uint32_t>{0};
        ++mm_itr;
    }
    // if (!err && memcpy(kv_A(*clusters, 0).minimizers.a, mm, len_id->size * sizeof(mm_t)) != kv_A(*clusters, 0).minimizers.a) err = ERR_RUNTIME;
    // if (!err) kv_A(*clusters, 0).minimizers.n = len_id->size;
    clusters.back().ids.push_back(len_id->id); // if (!err) kv_push(read_id_t, kv_A(*clusters, 0).ids, len_id->id);
    mm += len_id->size;
    ++len_id;
    std::vector<std::size_t> hits;
    std::vector<mm_t> mm_buffer;
    for (std::size_t i = 1; i < nsketches; ++i) { // for all sketched reads
        hits.resize(clusters.size());
        for (auto& v : hits) v = 0;
        mm_itr = mm;
        mm_buffer.clear();
        for (std::size_t j = 0; j < len_id->size; ++j) { // for the minimizers in read
            auto itr = mm2clusters.find(*mm_itr);
            if (itr == mm2clusters.end()) { // new minimizers not seen before
                // do nothing, add to cluster later
            } else { // minimizer already seen before
                auto& cluster_ids = itr->second;
                for (auto id : cluster_ids) {
                    ++hits[id];
                }
            }
            mm_buffer.push_back(*mm_itr);
            ++mm_itr;
        }
        auto max_itr = std::max_element(hits.begin(), hits.end());
        std::size_t best_cluster_idx = max_itr - hits.begin();
        // for (auto v : hits) std::cerr << v << " ";
        // std::cerr << "\nhits.size() = " << hits.size() << ", best cluster id = " << best_cluster_idx << "\n";
        if (static_cast<double>(*max_itr) / len_id->size > similarity_threshold) { // merge read into best cluster
            clusters.at(best_cluster_idx).ids.push_back(len_id->id); // add read to cluster
            for (auto mm : mm_buffer) { // for the minimizers in read
                auto itr = mm2clusters.find(mm);
                if (itr == mm2clusters.end()) { // new minimizers not seen before
                    mm2clusters[mm] = std::vector<uint32_t>{static_cast<uint32_t>(best_cluster_idx)};
                } else { // minimizer already seen before
                    auto lb = std::lower_bound(itr->second.begin(), itr->second.end(), best_cluster_idx); // O(log2(N)) complexity, if lists are short maybe O(N) is better
                    if (lb == itr->second.end() or *lb != best_cluster_idx) itr->second.insert(lb, best_cluster_idx); // insert cluster idx to preserve order and set property
                }
                // add minimizers to clusters by preserving set property (no duplicates) and order
                auto& cluster_mms = clusters.at(best_cluster_idx).minimizers;
                auto lb = std::lower_bound(cluster_mms.begin(), cluster_mms.end(), mm);
                if (lb == cluster_mms.end() or *lb != mm) cluster_mms.insert(lb, mm);
            }
        } else { // new cluster
            clusters.emplace_back();
            clusters.back().minimizers = mm_buffer;
            clusters.back().ids.push_back(len_id->id);
            for (auto mm : mm_buffer) { // for the minimizers in read
                auto itr = mm2clusters.find(mm);
                if (itr == mm2clusters.end()) { // new minimizers not seen before
                    mm2clusters[mm] = std::vector<uint32_t>{static_cast<uint32_t>(clusters.size() - 1)};
                } else { // minimizer already seen before
                    mm2clusters[mm].push_back(static_cast<uint32_t>(clusters.size() - 1));
                }
            }
        }
        mm += len_id->size;
        ++len_id;
    }
    len_id = NULL;
    mm = NULL;
    if (munmap(index, filestat.st_size) != 0) return 3;
    if (close(fd)) return 3;
    } // <<< memory management

    if (merge_threshold != 0) {
        /* test it iif the previous part is fast enough */
    } 

    if (cluster_save(clusters, output_filename.c_str())) {
        std::cerr << "Error writing output file\n";
        return 5;
    }
    // if (cluster_print(clusters)) {
    //     std::cerr << "Error while printing clusters\n";
    //     return 5;
    // }
    return 0;
}

int cluster_print(const std::vector<cluster_t>& clusters) {
    size_t i, j;
    for (i = 0; i < clusters.size(); ++i) {
        fprintf(stderr, "cluster[%lu] : [", i);
        if (clusters.at(i).ids.size() == 0) return 1;
        for (j = 0; j < clusters.at(i).ids.size() - 1; ++j) {
            fprintf(stderr, "%llu, ", clusters.at(i).ids.at(j));
        }
        fprintf(stderr, "%llu]\n", clusters.at(i).ids.at(clusters.at(i).ids.size() - 1));
    }
    return 0;
}

int cluster_save(
    const std::vector<cluster_t>& clusters, 
    char const *const output_filename
) {
    int err;
    FILE *comma_output;
    size_t i, j;
    err = 0;
    if (output_filename) {
        if ((comma_output = fopen(output_filename, "w")) == NULL) err = 1;
    } else {
        comma_output = stdout;
    } 
    if (!err) {
        fprintf(comma_output, "read ID,cluster ID\n");
        for (i = 0; i < clusters.size(); ++i) {
            for (j = 0; j < clusters.at(i).ids.size(); ++j) {
                fprintf(comma_output, "%llu,%lu\n", clusters.at(i).ids.at(j), i);
            }
        }
    }
    if (comma_output) {
        fflush(comma_output);
        if (comma_output != stdout) fclose(comma_output);
    }
    return err;
}