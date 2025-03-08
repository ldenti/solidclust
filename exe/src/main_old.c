#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <zlib.h>

#include "kseq.h"
#include "ksort.h"

#include "kmer.h"
#include "rsort.h"

#ifndef MIN
#define MIN(a, b) ((a) < (b) ? a : b)
#endif

KSEQ_INIT(gzFile, gzread)
KSORT_INIT_GENERIC(uint64_t)

void err(char *s) {
    fprintf(stderr, "%s\n", s);
    exit(1);
}

/* Run-length encoded bit vector storing distance between two ones */
typedef struct {
    uint64_t *p;
    int n;
    int m;
} rle_t;

rle_t *rle_init() {
    rle_t *r = malloc(sizeof(rle_t));
    r->p = calloc(128, sizeof *r->p);
    r->m = 128;
    r->n = 0;
  return r;
}

void rle_insert(rle_t *r, uint64_t rl) {
    if (r->n == r->m) {
        /* fprintf(stderr, "Reallocating rle..\n"); */
        r->p = realloc(r->p, r->m * 2 * sizeof *r->p);
        r->m *= 2;
    }
    r->p[r->n] = rl;
    ++r->n;
}

void rle_clean(rle_t *r) { r->n = 0; }

int rle_write(rle_t *r, FILE *fp) {
    if ((fwrite(&r->n, sizeof r->n, 1, fp)) != 1) err("Error while writing rle (n)..");
    if ((fwrite(r->p, sizeof r->p, r->n, fp)) != r->n) err("Error while writing rle runs..");
    /* add separator */
    if ((fwrite("|", 1, 1, fp)) != 1) err("Error while writing rle separator..");
    return 0;
}

rle_t *rle_load(FILE *fp, rle_t *r) {
    char c;
    assert(fp);
    assert(r);
    if ((fread(&r->n, sizeof r->n, 1, fp)) != 1) err("Error while loading n..");
    if (r->n > r->m) {
        r->p = realloc(r->p, r->n * sizeof(uint64_t));
        r->m = r->n;
    }
    if ((fread(r->p, sizeof r->p, r->n, fp)) != r->n) err("Error while loading runs..");
    if ((fread(&c, 1, 1, fp)) != 1) err("Error while loading sep..");
    if (c != '|') err("Error while loading rle. Separator mismatch (2)..");
    return r;
}

void rle_print(rle_t *r) {
    int i;
    assert(r);
    printf("%d/%d: ", r->n, r->m);
    for (i = 0; i < r->n; ++i) {
        printf("%llu ", r->p[i]);
    }
    printf("\n");
}

/* Returns number of bits set in both vectors */
int rle_intersection(rle_t *r1, rle_t *r2) {
    int s;
    int i1, i2;
    int pone1, pone2;
    assert(r1);
    assert(r2);
    /* rle_print(r1); */
    /* rle_print(r2); */
    s = 0;
    i1 = 0, 
    i2 = 0;
    pone1 = r1->p[i1];
    pone2 = r2->p[i2];

    while (i1 < r1->n && i2 < r2->n) {
        /* printf("Checking %d (%d) vs %d (%d).. ", pone1, i1, pone2, i2); */
        if (pone1 == pone2) { 
            /* printf("Match\n"); */
            ++s;
            ++i1;
            pone1 += r1->p[i1] + 1;
            ++i2;
            pone2 += r2->p[i2] + 1;
        } else if (pone1 < pone2) { 
            /* printf("Moving 1\n"); */
            ++i1;
            pone1 += r1->p[i1] + 1;
        } else { 
            /* printf("Moving 2\n"); */
            ++i2;
            pone2 += r2->p[i2] + 1;
        }
    }
    /* printf("Intersection: %d\n", s); */
    return s;
}

void rle_destroy(rle_t *r) {
    free(r->p);
    free(r);
}

/* kept just for debugging. we dont need it */
typedef struct {
    uint64_t *blocks;
    int n;
    int b;
} bv_t;

void bv_clean(bv_t *bv) {
    int i;
    assert(bv);
    for (i = 0; i < bv->b; ++i) bv->blocks[i] = 0;
}

bv_t *bv_init(int n) {
    bv_t *bv = malloc(sizeof(bv_t));
    bv->n = n;
    bv->b = n / 64;
    bv->b = bv->b > 0 ? bv->b : 1;
    bv->blocks = malloc(bv->b * sizeof(uint64_t));
    bv_clean(bv);
    return bv;
}

void bv_set(bv_t *bv, int p) {
    int b;
    int off;
    assert(bv);

    b = p / 64;
    off = 64 - p % 64 - 1;
    bv->blocks[b] |= ((uint64_t)1 << off);
}

void bv_print(bv_t *bv) {
    int i, b;
    assert(bv);

    for (b = 0; b < bv->b; ++b) {
        for (i = sizeof(uint64_t) * 8 - 1; i >= 0; --i) {
            printf("%llu", (bv->blocks[b] >> i) & 1);
        }
        printf(" ");
    }
    printf("\n");
}

void bv_destroy(bv_t *bv) {
    if (bv && bv->blocks) free(bv->blocks);
    if (bv) free(bv);
}

void rle2bv(rle_t *r, bv_t *b) {
    int i, p;
    assert(r);
    assert(b);

    /* Assuming b to have correct size */
    bv_clean(b);
    p = 0;
    for (i = 0; i < r->n; ++i) {
        p += r->p[i];
        bv_set(b, p);
        ++p;
    }
}

/* List of minimizers */
typedef struct {
    uint64_t *mmers;
    int size;
    int m;
} mmers_t;

mmers_t *ms_init() {
    mmers_t *ms = malloc(sizeof(mmers_t));
    ms->mmers = calloc(128, sizeof *ms->mmers);
    ms->size = 0;
    ms->m = 128;
    return ms;
}

void ms_clean(mmers_t *ms) { 
    assert(ms);
    ms->size = 0; 
}

void ms_destroy(mmers_t *ms) {
    if (ms && ms->mmers) free(ms->mmers);
    if (ms) free(ms);
}

void ms_add(mmers_t *ms, uint64_t x) {
    assert(ms);
    if (ms->size == ms->m) {
        /* fprintf(stderr, "Reallocating mmers..\n"); */
        ms->mmers = realloc(ms->mmers, 2 * ms->m * sizeof *ms->mmers);
        ms->m *= 2;
    }
    ms->mmers[ms->size] = x;
    ++ms->size;
}

double get_score(char *qual, int s, int e) {
    int i;
    double q, score;
    assert(qual);

    q = 0.0;
    score = 1.0;
    for (i = s; i < e; ++i) {
        q = (double)qual[i] - 33.0;
        /**
           ==16787==  Uninitialised value was created by a heap allocation
           ==16787==    at 0x4847C5D: realloc (vg_replace_malloc.c:1801)
           ==16787==    by 0x401E42: kseq_read (main.c:14)
           ==16787==    by 0x4034A2: main (main.c:102)
         **/
        score *= 1.0 - pow(10.0, -q / 10);
    }
    return score;
}

uint64_t get_mmer(uint64_t kmer_d, int w, int m) {
    int i;
    uint64_t mmer, mmer_d, rcmmer_d, cmmer_d, bm;
    mmer = -1;
    mmer_d = 0;   /* mmmer */
    rcmmer_d = 0; /* reverse and complemented mmer */
    cmmer_d = 0;  /* canonical mmer */
    bm = (((uint64_t)1 << (m * 2)) - 1) << (2 * (w - m));
    for (i = 0; i < w - m + 1; ++i) {
        mmer_d = ((kmer_d & bm) >> (2 * (w - m - i)));
        rcmmer_d = rc(mmer_d, m);
        cmmer_d = MIN(mmer_d, rcmmer_d);
        if (cmmer_d < mmer)
            mmer = cmmer_d;
        bm = bm >> 2;
    }
    /* assert(mmer != -1); */
    return mmer;
}

void get_mmers(char *seq, char *qual, int l, int w, int m, mmers_t **mmers) {
    char *kmer;       /* first kmer on sequence (plain) */
    uint64_t kmer_d;  /* kmer */
    uint8_t c;        /* new character to append */
    int p;            /* current position on segment */
    uint64_t mmer;
    assert(seq);
    assert(qual);
    assert(mmers);

    kmer = malloc(sizeof(char) * (w + 1));
    kmer_d = 0;
    p = 0;
    mmer = 0;
    strncpy(kmer, seq, w);
    kmer_d = k2d(kmer, w);
    if (get_score(qual, 0, w) > 0.95) {
        mmer = get_mmer(kmer_d, w, m);
        /* printf("%ld %ld\n", kmer_d, mmer); */
        ms_add(*mmers, mmer);
    }
    for (p = w; p < l; ++p) {
        c = to_int[(int)seq[p]];
        kmer_d = lsappend(kmer_d, c, w);
        if (get_score(qual, p, p + w) > 0.95) {
            mmer = get_mmer(kmer_d, w, m);
            /* printf("%ld %ld\n", kmer_d, mmer); */
            ms_add(*mmers, mmer);
        }
    }
    free(kmer);
}

/* Merge two rle vectors and replace r2 */
void rle_merge(rle_t *r1, rle_t *r2, int m) {
    rle_t *rr;
    int i1, i2;
    int pone1, pone2;
    int pp, rl;
    assert(r1);
    assert(r2);
    /* rle_print(r1); */
    /* rle_print(r2); */
    /* bv_t *bv = bv_init((uint64_t)1 << (2 * m)); */
    /* rle2bv(r1, bv); */
    /* bv_print(bv); */
    /* rle2bv(r2, bv); */
    /* bv_print(bv); */

    rr = rle_init();
    i1 = 0;
    i2 = 0;
    pone1 = r1->p[i1];
    pone2 = r2->p[i2];
    pp = 0;
    while (i1 < r1->n && i2 < r2->n) {
        if (pone1 == pone2) {
            rl = pone1 - pp;
            rle_insert(rr, rl);
            pp = pone1 + 1;
            ++i1;
            pone1 += r1->p[i1] + 1;
            ++i2;
            pone2 += r2->p[i2] + 1;
        } else if (pone1 < pone2) {
            rl = pone1 - pp;
            rle_insert(rr, rl);
            pp = pone1 + 1;
            ++i1;
            pone1 += r1->p[i1] + 1;
        } else {
            rl = pone2 - pp;
            rle_insert(rr, rl);
            pp = pone2 + 1;
            ++i2;
            pone2 += r2->p[i2] + 1;
        }
    }
    /* rle2bv(rr, bv); */
    /* bv_print(bv); */

    if (rr->n > r2->m) {
        r2->p = realloc(r2->p, rr->n * sizeof *r2->p);
    }
    memcpy(r2->p, rr->p, rr->n * sizeof *r2->p);
    r2->m = rr->n;
    r2->n = rr->n;

    rle_destroy(rr);
}

int old_main(int argc, char **argv) {
    int i, r, c;
    char *fastq_file;
    int w = 21, m = 15; /* atoi(argv[2]); */

    fastq_file = argv[1];
    gzFile fp = gzopen(fastq_file, "r");
    kseq_t *seq = kseq_init(fp);
    int l = 0;

    mmers_t *mmers = ms_init(); /* minimizers list for a read */
    rle_t *rle = rle_init();    /* runlength encoded representation of a read */

    /* Store list of [read id in input ordering, number of mmers]*/
    pair_t *mpr = malloc((1 << 24) * sizeof(pair_t));
    int mpr_m = 1 << 24;
    int mpr_n = 0;

    FILE *ofp = fopen("reads.XXX", "wb");

    
    fwrite("IOR|", 1, 4, ofp);

    while ((l = kseq_read(seq)) >= 0) {
        get_mmers(seq->seq.s, seq->qual.s, l, w, m, &mmers);
        int i;

        /* char xx[m + 1]; */
        /* for (i = 0; i < mmers->size; ++i) { */
        /*   d2s(mmers->mmers[i], m, xx); */
        /*   printf("%ld:%s ", mmers->mmers[i], xx); */
        /* } */
        /* printf("\n"); */
        ks_introsort(uint64_t, mmers->size, mmers->mmers);
        /* for (i = 0; i < mmers->size; ++i) { */
        /*   printf("%ld ", mmers->mmers[i]); */
        /* } */
        /* printf("\n"); */
        /* Flag repeated minimizers */
        for (i = 0; i < mmers->size;) {
            int mmer = mmers->mmers[i];
            ++i;
            while (i < mmers->size && mmers->mmers[i] == mmer) {
                mmers->mmers[i] = -1;
                ++i;
            }
        }

        /* for (i = 0; i < mmers->size; ++i) { */
        /*   printf("%ld ", mmers->mmers[i]); */
        /* } */
        /* printf("\n"); */

        /* Inserting minimizers in rle bitvector */
        int nm = 0, rl, p = 0;
        for (i = 0; i < mmers->size; ++i) {
            if (mmers->mmers[i] == -1)
                continue;
            ++nm;
            rl = mmers->mmers[i] - p;
            rle_insert(rle, rl);
            p = mmers->mmers[i] + 1;
        }
        rle_insert(rle, (uint64_t)1 << (2 * m));

        /* Append to file */
        rle_write(rle, ofp);

        /* Store read information (number of) */
        if (mpr_n == mpr_m) {
            mpr = malloc(mpr_m * 2 * sizeof(pair_t));
            mpr_m *= 2;
        }
        mpr[mpr_n] = (pair_t){.a = mpr_n, .b = nm};

        ++mpr_n;

        rle_clean(rle);
        ms_clean(mmers);
    }
    fclose(ofp);
    kseq_destroy(seq);
    gzclose(fp);
    ms_destroy(mmers);

    fprintf(stderr, "We have %d reads\n", mpr_n);

    /* mpr_n is now the total number of reads */

    /* Compute index for fast extraction from binary file */
    int *rofx = malloc(mpr_n * sizeof(int));
    rofx[0] = 0;
    for (i = 1; i < mpr_n; ++i) {
        /* n mmers in a read means (n+1) runs of zeros in rle */
        rofx[i] = rofx[i - 1] + (mpr[i - 1].b + 1);
    }

    /* Sort reads by number of minimizers */
    radix_sort(mpr, mpr + mpr_n);
    /* for (i = 0; i < mpr_n; ++i) */
    /*   printf("%d,%d ", mpr[i].a, mpr[i].b); */
    /* printf("\n"); */

    /* Start clustering */
    int cn = 0; /* number of clusters */
    int cm = 1024;
    int **rpc = malloc(1024 * sizeof(int *));
    int *rpc_m = malloc(1024 * sizeof(int));
    int *rpc_n = malloc(1024 * sizeof(int));
    int *cofx = malloc(1024 * sizeof(int)); /* offsets in cluster file */
    for (i = 0; i < cm; ++i) {
        rpc[i] = malloc(1024 * sizeof(int));
        rpc_m[i] = 1024;
        rpc_n[i] = 0;
        cofx[i] = 0;
    }

    rle_t *cluster_rle = rle_init();
    FILE *ifp = fopen("reads.XXX", "rb");
    int cr;  /* current read in input order */
    int ofx; /* offset in file */
    int rr;
    FILE *cfp = fopen("clusters.XXX", "wb");
    fwrite("IOC|", 1, 4, cfp);

    char sep; /* separator for sanity check */
    for (r = 0; r < mpr_n; ++r) {
        if ((r + 1) % 1000 == 0)
            fprintf(stderr, "Analyzed %d reads..\n", r + 1);
        rr = mpr_n - r - 1;
        cr = mpr[rr].a;
        ofx = 4 * cr + 8 * rofx[cr] + cr;
        fseek(ifp, 4 + ofx - 1, SEEK_SET);
        if ((fread(&sep, 1, 1, ifp)) != 1)
            err("Error while loading sep..");
        if (sep != '|')
            err("Error while loading rle. Separator mismatch (1)..");
        rle_load(ifp, rle);
        /* rle2bv(rle, bv); */
        /* bv_print(bv); */

        int best_c = -1;
        int intersection;
        float jsim;
        float best_sim = 0.0;
        for (c = 0; c < cn; ++c) {
            /* if we are here, we have at least one cluster and we opened the file in reading mode */
            fseek(cfp, cofx[c] - 1, SEEK_SET);
            if ((fread(&sep, 1, 1, cfp)) != 1)
                err("Error while loading cluster sep..");
            if (sep != '|')
                err("Error while loading cluster rle. Separator mismatch (1)..");
            rle_load(cfp, cluster_rle);
            intersection = rle_intersection(rle, cluster_rle);
            jsim = intersection / (float)(rle->n + cluster_rle->n - intersection);
            if (jsim > 0.5 && jsim > best_sim) {
                best_c = c;
                best_sim = jsim;
            }
        }

        if (best_c == -1) {
            if (cn == cm) {
                rpc = realloc(rpc, 2 * cm * sizeof(int *));
                rpc_m = realloc(rpc_m, 2 * cm * sizeof(int *));
                rpc_n = realloc(rpc_n, 2 * cm * sizeof(int));
                cofx = realloc(cofx, 2 * cm * sizeof(int));
                for (i = cn; i < 2 * cm; ++i)
                {
                    rpc[i] = malloc(1024 * sizeof(int));
                    rpc_m[i] = 1024;
                    rpc_n[i] = 0;
                    cofx[i] = 0;
                }
                cm *= 2;
            }
            /* add read id to cluster. we have space since cluster is new */
            rpc[cn][rpc_n[cn]] = cr;
            ++rpc_n[cn];
            /* dump cluster: close, open in append mode, store offset, write, close, open in read mode */
            fclose(cfp);
            cfp = fopen("clusters.XXX", "ab");
            cofx[cn] = ftell(cfp);
            rle_write(rle, cfp);
            fclose(cfp);
            cfp = fopen("clusters.XXX", "rb");
            ++cn;
        } else {
            rle_merge(rle, cluster_rle, m);

            /* add read id to cluster. we may not have space */
            if (rpc_n[best_c] == rpc_m[best_c]) {
                rpc[best_c] = realloc(rpc[best_c], rpc_m[best_c] * 2 * sizeof(int));
                rpc_m[best_c] *= 2;
            }
            rpc[best_c][rpc_n[best_c]] = cr;
            ++rpc_n[best_c];

            fclose(cfp);
            cfp = fopen("clusters.XXX", "ab");
            cofx[best_c] = ftell(cfp);
            rle_write(cluster_rle, cfp);
            fclose(cfp);
            cfp = fopen("clusters.XXX", "rb");
        }
    }
    fprintf(stderr, "Analyzed %d reads..\n", mpr_n);
    fclose(ifp);
    fclose(cfp);

    fprintf(stderr, "We created %d clusters\n", cn);

    /* invert and dump cluster map */
    int *cpr = malloc(mpr_n * sizeof(int));
    for (c = 0; c < cn; ++c) {
        for (r = 0; r < rpc_n[c]; ++r) {
            cpr[rpc[c][r]] = c;
            printf("%d\t%d\n", rpc[c][r], c);
        }
    }
    /* reiterate over reads */
    /* output in batches depending on cluster */
    int flimit = 1024;
    FILE **clusters = malloc(flimit * sizeof(FILE *));
    char fn[128];
    int c_inf = 0, cc, c_sup = 0, batch_size;

    struct stat st = {0};
    if (stat("CLUSTERS", &st) == -1) {
        mkdir("CLUSTERS", 0700);
    }

    while (c_inf < cn) {
        batch_size = c_inf + flimit > cn ? cn - c_inf : flimit;
        c_sup = c_inf + batch_size;
        printf("Dumping clusters [%d:%d)\n", c_inf, c_sup);
        for (i = 0; i < batch_size; ++i) {
            sprintf(fn, "CLUSTERS/c%d.fq", c_inf + i);
            clusters[i] = fopen(fn, "w");
        }

        fp = gzopen(fastq_file, "r");
        seq = kseq_init(fp);
        l = 0;
        FILE *f;
        int r = 0;
        while ((l = kseq_read(seq)) >= 0) {
            cc = cpr[r];
            ++r;
            if (c_inf <= cc && cc < c_sup) {
                f = clusters[cc - c_inf];
                fprintf(f, "@%s\n", seq->name.s);
                fprintf(f, "%s\n", seq->seq.s);
                fprintf(f, "+\n");
                fprintf(f, "%s\n", seq->qual.s);
            }
        }
        kseq_destroy(seq);
        gzclose(fp);

        for (i = 0; i < batch_size; ++i) {
            fclose(clusters[i]);
        }
        c_inf += batch_size;
    }
    free(clusters);

    free(cpr);
    for (i = 0; i < cm; ++i) {
        free(rpc[i]);
    }
    free(rpc);
    free(cofx);
    free(rpc_m);
    free(rpc_n);
    rle_destroy(cluster_rle);
    rle_destroy(rle);
    free(rofx);
    free(mpr);

    return 0;
}
