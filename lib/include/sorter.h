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
  rstype_t *b;
  rstype_t *e; /* begin and end of each bucket */
} rsbucket_t;

void radix_sort(rstype_t *beg, rstype_t *end);

#endif
