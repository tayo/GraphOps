#ifndef PAGERANKCPUCODE_H
#define PAGERANKCPUCODE_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>

#include <MaxSLiCInterface.h>
#include "Pagerank.max"
#include "Maxfiles.h"

//#include "bfs.h"

//uint64_t micro_time(void);

/////////////////////////////////////
//
// Types
//
typedef struct {
  unsigned char c[6];
} /*__attribute__((packed))*/ uint48_t;
typedef struct {
  unsigned char c[6];
} /*__attribute__((packed))*/ float48_t;
typedef struct {
  unsigned char c[3];
} /*__attribute__((packed))*/ uint24_t;

typedef uint48_t  node_t;
typedef uint48_t  edge_t;
typedef uint24_t  prop_uint24_t;
typedef uint48_t  prop_uint48_t;
typedef float48_t prop_f48_t;

void assign24 (uint24_t* ptr, unsigned long long val);
unsigned long long getval24(uint24_t* ptr);
void assign48 (uint48_t* ptr, unsigned long long val);
unsigned long long getval48(uint48_t* ptr);
void assign_f48 (float48_t* ptr, float val);
float getval_f48(float48_t* ptr);

struct timespec ts_diff(struct timespec start, struct timespec end);

void init_prop_uint24_root(prop_uint24_t* prop_ptr, uint64_t N, 
                           uint64_t root, unsigned long long);
void init_prop_uint24(prop_uint24_t* prop_ptr, uint64_t N, unsigned long long);

void init_prop_f48_root(prop_f48_t* prop_ptr, uint64_t N, 
                        uint64_t root, float init_val);
void init_prop_f48(prop_f48_t* prop_ptr, uint64_t N, float init_val);
void init_rep_f48(prop_f48_t* prop_ptr, uint64_t M, float init_val);

void print_graph(node_t *n_ptr, int nb, int eb, int p24b, int p48b,
                 int b_len, const char *title_str);
void cleanup();
#endif
