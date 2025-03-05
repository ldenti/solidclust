#ifndef IOC_RADSORT_H
#define IOC_RADSORT_H

/* Code adapted from https://gist.github.com/attractivechaos/2886685 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct {
  uint32_t a;
  uint32_t b;
} pair_t;

#define rstype_t pair_t  /* type of the array */
#define rskey(x) ((x).b) /* specify how to get the integer from rstype_t */

#define RS_MIN_SIZE 64 /* for an array smaller than this, use insertion sort */

typedef struct {
  rstype_t *b, *e; /* begin and end of each bucket */
} rsbucket_t;

/* Insertion sort */
void rs_insertsort(rstype_t *beg, rstype_t *end);

/* Sort between [$beg, $end); take radix from ">>$s&((1<<$n_bits)-1)" */
void rs_sort(rstype_t *beg, rstype_t *end, int n_bits, int s);

void radix_sort(rstype_t *beg, rstype_t *end);

#endif
