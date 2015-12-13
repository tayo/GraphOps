#ifndef GRAPH_GEN_H
#define GRAPH_GEN_H

#include "PagerankCpuCode.h"

/*
 * N: number of nodes
 * M: number of edges  (M = N * avg_degree)
 */
void create_uniform_random_graph(node_t* node_ptr, edge_t* edge_ptr,
                                 uint32_t N, uint32_t M);

void create_RMAT_graph(node_t* node_ptr, edge_t* edge_ptr,
                       double a, double b, double c, uint32_t N, uint32_t M);
#endif
