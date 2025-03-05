#include "rsort.h"

void rs_insertsort(rstype_t *beg, rstype_t *end) {
  rstype_t *i;
  for (i = beg + 1; i < end; ++i)
    if (rskey(*i) < rskey(*(i - 1))) {
      rstype_t *j, tmp = *i;
      for (j = i; j > beg && rskey(tmp) < rskey(*(j - 1)); --j)
        *j = *(j - 1);
      *j = tmp;
    }
}

// sort between [$beg, $end); take radix from ">>$s&((1<<$n_bits)-1)"
void rs_sort(rstype_t *beg, rstype_t *end, int n_bits, int s) {
  rstype_t *i;
  int size = 1 << n_bits, m = size - 1;
  rsbucket_t *k, b[size], *be = b + size; // b[] keeps all the buckets

  for (k = b; k != be; ++k)
    k->b = k->e = beg;
  for (i = beg; i != end; ++i)
    ++b[rskey(*i) >> s & m].e;  // count radix
  for (k = b + 1; k != be; ++k) // set start and end of each bucket
    k->e += (k - 1)->e - beg, k->b = (k - 1)->e;
  for (k = b; k != be;) { // in-place classification based on radix
    if (k->b != k->e) {   // the bucket is not full
      rsbucket_t *l;
      if ((l = b + (rskey(*k->b) >> s & m)) != k) { // different bucket
        rstype_t tmp = *k->b, swap;
        do { // swap until we find an element in bucket $k
          swap = tmp;
          tmp = *l->b;
          *l->b++ = swap;
          l = b + (rskey(tmp) >> s & m);
        } while (l != k);
        *k->b++ = tmp; // push the found element to $k
      } else
        ++k->b; // move to the next element in the bucket
    } else
      ++k; // move to the next bucket
  }
  for (b->b = beg, k = b + 1; k != be; ++k)
    k->b = (k - 1)->e; // reset k->b
  if (s) {             // if $s is non-zero, we need to sort buckets
    s = s > n_bits ? s - n_bits : 0;
    for (k = b; k != be; ++k)
      if (k->e - k->b > RS_MIN_SIZE)
        rs_sort(k->b, k->e, n_bits, s);
      else if (k->e - k->b > 1)
        rs_insertsort(k->b, k->e);
  }
}

void radix_sort(rstype_t *beg, rstype_t *end) {
  if (end - beg <= RS_MIN_SIZE)
    rs_insertsort(beg, end);
  else
    rs_sort(beg, end, 8, sizeof(rskey(*beg)) * 8 - 8);
}
