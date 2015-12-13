#ifndef EDGE_LIST_INPUT_H
#define EDGE_LIST_INPUT_H

#include "PagerankCpuCode.h"

void read_graph_size(const char* filename);
void read_edge_list(node_t* node_ptr, edge_t* edge_ptr, 
                    uint32_t N, uint32_t M,
                    const char* filename);
#endif
