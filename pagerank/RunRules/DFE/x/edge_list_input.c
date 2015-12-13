#include "edge_list_input.h"

extern int NUM_NODES;
extern int NUM_EDGES;

void read_graph_size(const char* filename) {
  FILE* f = NULL;
  if (NULL == (f = fopen(filename, "r"))) {
    fprintf(stderr, "Error opening file %s for reading.\n", filename);
    cleanup();
    exit(EXIT_FAILURE);
  }

  char line[128];

  // read N, M and allocate space for them
  // N: nodes/vertices      M: edges
  uint32_t N, M;
  while (fgets(line, 128, f)) {
    // skip header comments
    if ((line[0]=='%') | (line[0]=='#')) {
      continue;
    }

    sscanf(line, "%d %d", &N, &M);
    break;
  }
  NUM_NODES = N;
  NUM_EDGES = M;
}
  
void read_edge_list(node_t* node_ptr, edge_t* edge_ptr, 
                    uint32_t N, uint32_t M,
                    const char* filename) {

  FILE* f = NULL;
  if (NULL == (f = fopen(filename, "r"))) {
    fprintf(stderr, "Error opening file %s for reading.\n", filename);
    cleanup();
    exit(EXIT_FAILURE);
  }
  char line[128];
  printf("Reading edge list.. N=%d, M=%d\n", N, M);
  uint32_t* src = malloc(M*sizeof(uint32_t));
  uint32_t* dst = malloc(M*sizeof(uint32_t));
  uint32_t* deg = malloc(N*sizeof(uint32_t));
  memset((void*)deg, 0, N*sizeof(unsigned char));

  uint32_t i=0, s, d;
  uint32_t skip_size=0;
  while (fgets(line, 128, f)) {
    //printf("i=%d  ", i);
    if ((line[0]=='%')|(line[0]=='#')) continue;
    // skip the graph size line
    if (skip_size == 0) {
      skip_size=1;
      continue;
    }
    sscanf(line, "%d %d", &s, &d);
    /*
    // debug
    if (i > 16518940) {
      printf(" i=%d, s=%d, d=%d, deg[%d]=%d\n", i, s, d, s, deg[s]);
    }
    if (i < 20) {
      printf(" i=%d, s=%d, d=%d, deg[%d]=%d\n", i, s, d, s, deg[s]);
    }
    */
    src[i] = s;
    dst[i] = d;
    deg[s]++;
    i++;
  }

  printf("done reading. about to assign to node_ptr..\n");
  assign48(&node_ptr[0], 0);
  for(uint32_t i=1; i<=N; i++) {
    assign48(&node_ptr[i], getval48(&node_ptr[i-1]) + deg[i-1]);
  }

  for (i=0; i < M; i++) {
    uint32_t u = src[i];
    uint32_t v = dst[i];
    uint32_t pos = deg[u]--;
    /*
    // debug
    if (i < 20) {
      printf("    i=%d, u=%d, v=%d, pos=%d\n", i, u, v, pos);
    }
    */
    /*
    if ((pos <= 0)|(i<20)) {
      printf("pos failure: i=%d, u=%d, v=%d, pos=%d\n", i, u, v, pos);
    }
    */
    assert(pos > 0);
    assign48(&edge_ptr[ getval48(&node_ptr[u]) + pos - 1], v);
  }

  free(src); free(dst); free(deg);
  fclose(f);
}
