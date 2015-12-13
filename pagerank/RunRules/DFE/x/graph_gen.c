#include "graph_gen.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

void create_RMAT_graph(node_t* node_ptr, edge_t* edge_ptr,
                       double a, double b, double c, uint32_t N, uint32_t M)
{
  //srand(1204);
  srand(104);

  double d;
  assert(a + b + c < 1);
  d = 1 - (a + b + c);

  uint32_t* src = malloc(M * sizeof(uint32_t));
  uint32_t* dst = malloc(M * sizeof(uint32_t));
  uint32_t* deg = malloc(N * sizeof(uint32_t));
  memset((void*)deg, 0, N*sizeof(uint32_t));

  uint32_t SCALE = (uint32_t) log2((double) N);

  // 1. edge-gen
  for (uint32_t i = 0; i < M; i++) {
    uint32_t u = 1;
    uint32_t v = 1;
    uint32_t step = N/2;
    double av = a;
    double bv = b;
    double cv = c;
    double dv = d;

    double p = drand48();
    if (p < av) { // do nothing
    } else if (p < (av + bv)) {
      v += step;
    } else if (p < (av + bv + cv)) {
      u += step;
    } else {
      v += step;
      u += step;
    }
    for (uint32_t j = 1; j < SCALE; j++) {
      step = step / 2;
      double var = 0.1;
      av *= 0.95 + var * drand48();		// vary abcd by 10%
      bv *= 0.95 + var * drand48();
      cv *= 0.95 + var * drand48();
      dv *= 0.95 + var * drand48();

      double S = av + bv + cv + dv;
      av = av / S;
      bv = bv / S;
      cv = cv / S;
      dv = dv / S;

      // choose partition
      p = drand48();
      if (p < av) { // do nothing
      } else if (p < (av + bv)) {
        v += step;
      } else if (p < (av + bv + cv)) {
        u += step;
      } else {
        v += step;
        u += step;
      }
    }

    src[i] = u - 1;
    dst[i] = v - 1;

    // avoid self edges
    if (src[i] == dst[i]) {
      i = i - 1;
      continue;
    }
  }

  // 2. permute.  skip this step

  // 3. count degree
  for (uint32_t i = 0; i < M; i++) {
    deg[src[i]]++;
  }

  // 4. Write the graph's data structures
  assign48(&node_ptr[0], 0); //node_ptr[0] = 0;
  for(uint32_t i=1; i<=N; i++) {
    assign48(&node_ptr[i], getval48(&node_ptr[i-1]) + deg[i-1]);
  }

  for(uint32_t i=0; i<M; i++) {
    uint32_t u = src[i];
    uint32_t v = dst[i];

    uint32_t pos = deg[u]--;
    //printf("pos[%d]=%d  ", i, pos);
    assert(pos >= 0);
    assign48(&edge_ptr[ getval48(&node_ptr[u]) + pos - 1], v);
  }

  free(src); free(dst); free(deg);
}


/*
 * N: number of nodes
 * M: number of edges
 *    M = N * avg_degree
 */
void create_uniform_random_graph(node_t* node_ptr, edge_t* edge_ptr,
                                 uint32_t N, uint32_t M)
{
  //srand(1204);
  srand(1203);
  //srand(1205);
  //srand(1206);

  uint32_t* src = malloc(M * sizeof(uint32_t));
  uint32_t* dst = malloc(M * sizeof(uint32_t));
  uint32_t* deg = malloc(N * sizeof(uint32_t));
  memset((void*)deg, 0, N*sizeof(uint32_t));

  for(uint32_t i=0; i<M; i++) {
    src[i] = rand() % N;
    dst[i] = rand() % N;

    deg[src[i]]++;
  }

  assign48(&node_ptr[0], 0); //node_ptr[0] = 0;
  for(uint32_t i=1; i<=N; i++) {
    assign48(&node_ptr[i], getval48(&node_ptr[i-1]) + deg[i-1]);
    //node_ptr[i] = node_ptr[i-1] + deg[i-1];
  }

  for(uint32_t i=0; i<M; i++) {
    uint32_t u = src[i];
    uint32_t v = dst[i];

    uint32_t pos = deg[u]--;
    assert(pos > 0);
    assign48(&edge_ptr[ getval48(&node_ptr[u]) + pos - 1], v);
    //edge_ptr[node_ptr[u] + pos - 1] = v;
  }

#ifdef DBG
  printf("node_ptr: %p\nedge_ptr: %p\n", node_ptr, edge_ptr);
  for (int i = 0; i < N; i++) {
    printf("node %d: ", i);
    for (int j = getval48(&node_ptr[i]); j < getval48(&node_ptr[i+1]); j++) {
      printf("%d ", getval48(&edge_ptr[j]));
    }
    printf("  --> degree[%d] is %d\n", i, deg[i]);
  }

  printf("N: %d, M: %d\n", N, M);
  printf("test for uninitialized data..\n");
  for (int i = N-1; i < ((N/48)*48+48); i++) {
    printf("empty node %d: %d\n", i, getval48(&node_ptr[i]));
  }
#endif

  free(src); free(dst); free(deg);
}
