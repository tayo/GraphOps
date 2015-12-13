#include "bfs.h"

#include <stdio.h>

void sw_bfs(node_t* n_ptr, edge_t* e_ptr, prop_uint24_t *lvls,
            uint32_t N, uint32_t M, int root) {

  int burst_len_bytes = 192; //384;  //max_get_burst_size(max_file_t*..)
  int node_bytes = (N+1)*sizeof(node_t);
  int node_bursts = ceil((double)node_bytes/(double)burst_len_bytes);
  int edge_bytes = (M)*sizeof(edge_t);
  int edge_bursts = ceil((double)edge_bytes/(double)burst_len_bytes);
  //int num_props = 4;
  int prop24_bytes = (N)*sizeof(prop_uint24_t);
  int prop24_bursts = ceil((double)prop24_bytes/(double)burst_len_bytes);

  int cur_lvl = 0;
  int edge_count = 0;
  int upd_count = 0;

  int mode = 0;
  int *set0 = malloc(N*sizeof(int)); //initially, current set
  int *set1 = malloc(N*sizeof(int)); //initially, next set
  memset(set0, 0, N*sizeof(int));
  memset(set1, 0, N*sizeof(int));
  set0[root] = 1;
  int *cur_set, *next_set;

  int is_done;
  do {
    is_done = 1;
    if (!mode) {
      cur_set = set0; 
      next_set = set1;
    } else {
      cur_set = set1; 
      next_set = set0;
    }

    //printf("cur_lvl=%d..\n", cur_lvl);
    for (uint32_t i = 0; i < N; i++) {
      if (cur_set[i]) {
        //printf("  node %d: ", i);
        assign24(&lvls[i], cur_lvl);
      }
    }
    for (int i = 0; i < N; i++) {
      if (cur_set[i]) {
        for (uint64_t j = getval48(&n_ptr[i]); j < getval48(&n_ptr[i+1]); j++) {
          int cur_edge = (int)getval48(&e_ptr[j]);
          //printf(" %d", cur_edge);
          edge_count++;
          //printf("node=%d, cur_edge=%d, lvl[cur_edge]=%d\n",
          //       i, cur_edge, (int)(getval24(&lvls[cur_edge])));
          if ((int)(getval24(&lvls[cur_edge])) > 200) {
            //printf("<upd>");
            if (next_set[cur_edge] == 0) {
              next_set[cur_edge] = 1;
              upd_count++;
            }
            is_done = 0;
          }
        }
        //printf("\n");
      }
    }
    //printf("Frontier Size = %d, Edge Count = %d\n", upd_count, edge_count);
    cur_lvl++;
    memset(cur_set, 0, N*sizeof(int));
    mode = 1-mode;
    upd_count = 0;
    edge_count = 0;
    if (cur_lvl == 100) {
      printf("Error in SW BFS..\n");
      free(set0); free(set1);
      break;
    }
  } while (!is_done);

  //if (cur_lvl < 100) printf("\nSuccess..\n");

  free(set0); free(set1);
  print_graph(n_ptr, node_bursts, edge_bursts, prop24_bursts,
              0, burst_len_bytes, "ReferenceBFS");
}
