#ifndef MINIMIZER_H
#define MINIMIZER_H

#include <stdint.h>
#include <stddef.h>
#include "../../bundled/klib/include/kvec.h"

typedef uint64_t mm_t;
typedef kvec_t(mm_t) mmv_t;

int minimizer_from_string(
    char const *const seq, 
    unsigned char const *const filter, 
    size_t seq_len, 
    const uint8_t k, 
    const uint8_t w, 
    const unsigned char canonical, 
    const uint64_t seed, 
    size_t *mm_count, 
    mmv_t *accumulator
);

#endif /* MINIMIZER_H */