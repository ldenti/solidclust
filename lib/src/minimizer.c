#include <assert.h>
#include "../../bundled/xxhash/include/xxhash.h"
#include "../include/exception.h"
#include "../include/minimizer.h"

const uint8_t seq_nt4_table[256] = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
};

typedef struct {
    uint64_t hash;
    uint64_t itself;
    size_t idx;
} mm_hash_t;

void mm_hash_clear(mm_hash_t *h) {
    assert(h);
    h->hash = UINT64_MAX;
    h->itself = 0;
}

void append_to_accumulator(mm_hash_t const *const mm, size_t idx, unsigned char const *const filter, mmv_t *accumulator) {
    assert(mm);
    assert(accumulator);
    if (filter && !filter[idx]) return;
    kv_push(mm_t, *accumulator, mm->itself);
}

int minimizer_from_string(
    char const *const seq, 
    unsigned char const *const filter, 
    size_t seq_len, 
    const uint8_t k, 
    const uint8_t w, 
    const unsigned char canonical, 
    const uint64_t seed, 
    mmv_t *accumulator
) {
    int err, c;
    unsigned char z, find_brand_new_min;
    size_t i, j, buf_pos, min_pos;
    size_t nbases_since_last_break;
    uint16_t shift;
    uint64_t mask;
    uint64_t km[2];
    mm_hash_t current;
    mm_hash_t *buffer;
    assert(seq);
    assert(accumulator);
    err = OK;
    shift = 2 * (k - 1);
    mask = (1ULL << (2 * k)) - 1;
    km[0] = km[1] = 0;
    nbases_since_last_break = 0;
    find_brand_new_min = FALSE;
    buf_pos = 0;
    min_pos = w;
    z = 0;
    if (!err && !(buffer = malloc(w * sizeof(mm_hash_t)))) err = ERR_MALLOC;
    for (i = 0; !err && i < seq_len; ++i) {
        c = seq_nt4_table[(int)seq[i]];
        if (c < 4) { /* treat a bad quality as a break */
            /* mm_hash_clear(&current); */
            km[0] = (km[0] << 2 | c) & mask;            /* forward k-mer */
            km[1] = (km[1] >> 2) | (3ULL ^ c) << shift; /* reverse k-mer */
            if (canonical && km[0] != km[1]) z = km[0] < km[1] ? 0 : 1;  /* strand, if symmetric k-mer then use previous strand */
            ++nbases_since_last_break;
            if (nbases_since_last_break >= k) {
                current.itself = km[z];
                current.hash = XXH64(&km[z], sizeof(km[z]), seed);
                if (nbases_since_last_break == k + 1) {
                    /* have seen the first window after a break, time to search for the minimum note that the current m-mer is checked by the next if */
                    min_pos = 0;
                    for (j = 0; j < w; ++j) if (buffer[j].hash < buffer[min_pos].hash) min_pos = j;
                }
                if (nbases_since_last_break >= k + 1) {  /* time to update the minimum, if necessary */
                    if (buf_pos % w == min_pos) {  /* old minimum outside window */
                        append_to_accumulator(&buffer[min_pos], i, filter, accumulator);  /* we save the old minimum, length on the right is k by definition */
                        find_brand_new_min = TRUE;
                    } else if (current.hash < buffer[min_pos].hash) {
                        append_to_accumulator(&buffer[min_pos], i, filter, accumulator);
                        min_pos = buf_pos;  /* actual update is outside if */
                    }
                }
                buffer[buf_pos++] = current;
                buf_pos %= w;  /* circular buffer */
                if (find_brand_new_min) {  /* find new minimum if the old one dropped out the window */
                    find_brand_new_min = FALSE;
                    min_pos = buf_pos;
                    for (j = (buf_pos + 1) % w; j < w; ++j) if (buffer[min_pos].hash > buffer[j].hash) min_pos = j;
                    for (j = 0; j <= buf_pos; ++j) if (buffer[min_pos].hash > buffer[j].hash) min_pos = j;
                }
            }
        } else {
            nbases_since_last_break = 0;
            if (min_pos < w) append_to_accumulator(&buffer[min_pos], i, filter, accumulator); /* push current minimum if available */
            min_pos = w;
            buf_pos = 0; /* is it still needed ? */
            /* 
                we always restart at the beginning of the buffer -> this allows to
                use min_pos as the position of the minimizer inside the first k-mer 
            */
        }
    }
    if (nbases_since_last_break == k) {  /* contig.length == 1 */
        for (j = 0; j < w; ++j) {
            if (buffer[j].hash < buffer[min_pos].hash) {
                min_pos = j;
            }
        }
    }
    if (i != 0 && min_pos < w) append_to_accumulator(&buffer[min_pos], i - 1, filter, accumulator);  /* push last minimum if available */
    if (buffer) free(buffer);
    return err;
}