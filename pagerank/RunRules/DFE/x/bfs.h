#ifndef BFS_H
#define BFS_H

#include "PagerankCpuCode.h"

void sw_bfs(node_t* node_ptr, edge_t* edge_ptr, prop_uint24_t* lvls,
            uint32_t N, uint32_t M, int root);

#endif
