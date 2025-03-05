#ifndef IOC_KMER_H
#define IOC_KMER_H

#include <stdint.h>

static const uint8_t to_int[128] = {
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 
    5, 0, 5, 1, 5, 5, 5, 2, 5, 5, 5, 5, 5, 5, 5, 5, 
    5, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 
    5, 0, 5, 1, 5, 5, 5, 2, 5, 5, 5, 5, 5, 5, 5, 5, 
    5, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5
};

uint8_t reverse_char(uint8_t c);

void d23(const uint64_t kmer, int k, char *kk);

void d2s(const uint64_t kmer, int k, char *kk);

uint64_t k2d(const char *kmer, uint8_t k);

uint64_t rc(uint64_t kmer, const uint8_t k);

uint64_t lsappend(const uint64_t kmer, const uint64_t c, const uint64_t k);

uint64_t rsprepend(const uint64_t kmer, const uint64_t c, const uint64_t k);

#endif
