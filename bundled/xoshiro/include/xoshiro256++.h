#ifndef XOSHIRO256PP_H
#define XOSHIRO256PP_H 

#include <stdint.h>

typedef uint64_t xoshiro256pp_seed_t;
typedef uint64_t xoshiro256pp_t;

typedef struct  {
    uint64_t s[4];
} xoshiro256pp_state_t;

void xoshiro256pp_seed_state(xoshiro256pp_state_t* state, xoshiro256pp_seed_t seed);
uint64_t xoshiro256pp_next(xoshiro256pp_state_t* state);
void xoshiro256pp_jump(xoshiro256pp_state_t* state);
void xoshiro256pp_long_jump(xoshiro256pp_state_t* state);

#endif