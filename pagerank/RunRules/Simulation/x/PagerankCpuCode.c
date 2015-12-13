#include "PagerankCpuCode.h"

#define AVG_DEG   8 //4
//#define NUM_NODES (512*1024) //(128*1024) 
//#define NUM_EDGES (NUM_NODES*AVG_DEG)
int NUM_NODES = (512*1024); //(128*1024) 
int NUM_EDGES; // (NUM_NODES*AVG_DEG)

#include "bfs.h"
#include "graph_gen.h"
#include "edge_list_input.h"
#include <omp.h>

// Globals
void *data_in, *data_out;

void assign24 (uint24_t* ptr, unsigned long long val) {
  unsigned long long tmp = 0xff;
  ptr->c[0] = val & tmp;
  ptr->c[1] = (val & (tmp<< 8)) >>  8;
  ptr->c[2] = (val & (tmp<<16)) >> 16;
}
unsigned long long getval24(uint24_t* ptr) {
  unsigned long long retval = 0;
  retval += (unsigned long long)(ptr->c[0]);
  retval += (unsigned long long)(ptr->c[1]) <<  8;
  retval += (unsigned long long)(ptr->c[2]) << 16;
  return retval;
}

void assign48 (uint48_t* ptr, unsigned long long val) {
  unsigned long long tmp = 0xff;
  ptr->c[0] = val & tmp;
  ptr->c[1] = (val & (tmp<< 8)) >>  8;
  ptr->c[2] = (val & (tmp<<16)) >> 16;
  ptr->c[3] = (val & (tmp<<24)) >> 24;
  ptr->c[4] = (val & (tmp<<32)) >> 32;
  ptr->c[5] = (val & (tmp<<40)) >> 40;
}
unsigned long long getval48(uint48_t* ptr) {
  unsigned long long retval = 0;
  retval += (unsigned long long)(ptr->c[0]);
  retval += (unsigned long long)(ptr->c[1]) <<  8;
  retval += (unsigned long long)(ptr->c[2]) << 16;
  retval += (unsigned long long)(ptr->c[3]) << 24;
  retval += (unsigned long long)(ptr->c[4]) << 32;
  retval += (unsigned long long)(ptr->c[5]) << 40;
  return retval;
}

void assign_f48(float48_t* ptr, float val) {
  uint32_t tmp = 0xff;
  uint32_t *p = (uint32_t*)(&val); //needed to do bit-wise operations
  ptr->c[0] =  *p & tmp;
  ptr->c[1] = (*p & (tmp<< 8)) >>  8;
  ptr->c[2] = (*p & (tmp<<16)) >> 16;
  ptr->c[3] = (*p & (tmp<<24)) >> 24;
  //ptr->c[4] = 0;
  //ptr->c[5] = 0;
}
float getval_f48(float48_t* ptr) {
  uint32_t retval_int = 0; //need 32-bit zero quantity
  retval_int |=  (ptr->c[0]);
  retval_int |= ((ptr->c[1]) << 8);
  retval_int |= ((ptr->c[2]) << 16);
  retval_int |= ((ptr->c[3]) << 24);
  float* p = (float*)(&retval_int);
  float retval = *p;
  return retval;
}


/////////////////////////////////////
//
// Helper Functions
//

void cleanup() {
  free(data_in);
  free(data_out);
}

/*
uint64_t micro_time(void) {
	struct timeval t;
	//struct timezone z;
	//gettimeofday(&t, &z);
	gettimeofday(&t, NULL);
	return t.tv_sec * 1000000 + t.tv_usec;
}
*/
struct timespec ts_diff(struct timespec start, struct timespec end) {
  struct timespec temp;
  if ((end.tv_nsec-start.tv_nsec)<0) {
    temp.tv_sec = end.tv_sec-start.tv_sec-1;
    temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
  } else {
    temp.tv_sec = end.tv_sec-start.tv_sec;
    temp.tv_nsec = end.tv_nsec-start.tv_nsec;
  }
  return temp;
}


void init_prop_uint24_root(prop_uint24_t* prop_ptr, uint64_t N, 
                           uint64_t root, unsigned long long init_val) {
  for (uint64_t i = 0; i < N; i++) {
    if (i == root)
      assign24(&prop_ptr[i], 0);   //set root to zero
    else
      assign24(&prop_ptr[i], init_val); //set other nodes to 'init_val'
  }
}
void init_prop_uint24(prop_uint24_t* prop_ptr, uint64_t N, 
                      unsigned long long init_val) {
  for (uint64_t i = 0; i < N; i++) {
    assign24(&prop_ptr[i], init_val); 
  }
}

void init_prop_f48_root(prop_f48_t* prop_ptr, uint64_t N, 
                        uint64_t root, float init_val) {
  for (uint64_t i = 0; i < N; i++) {
    if (i == root)
      assign_f48(&prop_ptr[i], 1.0); //set root to 1.0
    else
      assign_f48(&prop_ptr[i], init_val); //set other nodes to 'init_val'
  }
}
void init_prop_f48(prop_f48_t* prop_ptr, uint64_t N, float init_val) {
  for (uint64_t i = 0; i < N; i++) {
    assign_f48(&prop_ptr[i], init_val);
  }
}
void init_rep_f48(prop_f48_t* rep_ptr, uint64_t M, float init_val) {
  for (uint64_t i = 0; i < M; i++) {
    assign_f48(&rep_ptr[i], init_val);
  }
}

// Go through float48 prop array and put degree of node in most significant 16 bits
void prepend_deg(prop_f48_t* prop_ptr, node_t* n_ptr) {
  int i = 0;
  uint16_t tmp = 0xff;

  for (i = 0; i < NUM_NODES; i++) {
    uint16_t deg = (uint16_t)(getval48(&n_ptr[i+1]) - getval48(&n_ptr[i]));
    float48_t* ptr = (float48_t*)(&prop_ptr[i]);
    ptr->c[4] = (deg & tmp);
    ptr->c[5] = (deg & (tmp<<8)) >> 8;
  }
}

// Go through float48 rep array and put degree of node in most significant 16 bits
void prepend_deg_rep(prop_f48_t* rep_ptr, node_t* n_ptr, edge_t* e_ptr) {
  int i = 0;
  uint16_t tmp = 0xff;
  for (i = 0; i < NUM_NODES; i++) {
    for (uint64_t j = getval48(&n_ptr[i]); j < getval48(&n_ptr[i+1]); j++) {
      uint64_t node = getval48(&e_ptr[j]);
      uint16_t deg = (uint16_t)(getval48(&n_ptr[node+1]) - getval48(&n_ptr[node]));
      float48_t* ptr = (float48_t*)(&rep_ptr[j]);
      ptr->c[4] = (deg & tmp);
      ptr->c[5] = (deg & (tmp<<8)) >> 8;
    }
  }
}

// retrieve the most significant 16-bits as a 16-bit integer
uint16_t getdeg16(float48_t* ptr) {
  uint16_t retval = 0;
  retval |= (ptr->c[4]);
  retval |= ( (ptr->c[5]) << 8);
  return retval;
}

// Divide every 32-bit float in rep array by its 16-bit integer (degree)
void div_by_deg_rep(prop_f48_t* rep_ptr, node_t* n_ptr) {
  for (int i = 0; i < NUM_NODES; i++) {
    for (uint64_t j = getval48(&n_ptr[i]); j < getval48(&n_ptr[i+1]); j++) {
      float48_t* ptr = (float48_t*)(&rep_ptr[j]);
      float f        = getval_f48(ptr);
      uint16_t deg   = getdeg16(ptr);
      if (deg == 0) deg = (uint16_t)1;  // prevent division by zero
      assign_f48(ptr, f/((float)deg) );
    }
  }
}

// pagerank: accumulate diff between prop array and rep array
float acc_diff(node_t *n_ptr, edge_t* e_ptr, 
               prop_f48_t* prop_ptr, prop_f48_t* rep_ptr) {

  // This is an O(1) lookup.  This data structure can be created off-line.
  // use this data structure to determine where in edge list we can find a rep
  // array entry for each entry
  int *pos_map = malloc(NUM_NODES*sizeof(int));
  for (int i = 0; i < NUM_NODES; i++) 
    pos_map[i]=-1;
  for (int j = 0; j < NUM_EDGES; j++) {  // populate the mappings
    int dst_node = getval48(&e_ptr[j]);
    if (pos_map[dst_node] == -1) pos_map[dst_node] = j;
  }

  float diff = 0.0;
  for (int i = 0; i < NUM_NODES; i++) {
    float48_t* ptr = (float48_t*)(&prop_ptr[i]);
    float prop_val = getval_f48(ptr);
    float rep_val = 0.0;
    uint16_t deg = 0;
    // need to find location of a rep array entry for this node
    int rep_pos = pos_map[i];
    float48_t* r_ptr = (float48_t*)(&rep_ptr[rep_pos]);
    rep_val = getval_f48(r_ptr);
    deg = getdeg16(r_ptr);

    /*
    // brute force search takes a VERY long time for large graphs
    for (uint64_t j = 0; j < NUM_EDGES; j++) {
      if (i == getval48(&e_ptr[j])) {
        float48_t* r_ptr = (float48_t*)(&rep_ptr[j]);
        rep_val = getval_f48(r_ptr);
        deg = getdeg16(r_ptr);
        break;
      }
    }
    */
//if (i<10)  /**************DEBUG********************/
//printf("@ i=%d, deg=%d, prop_val=%f, rep_val=%f, rep_val*deg=%f, diff=%f..\n", 
//   i, deg, prop_val, rep_val, rep_val*(float)deg, diff);

    // rep_val has been divided by deg. Account for this when calculating the diff
    diff += fabsf(prop_val - ((float)deg*rep_val) );
  }
  free(pos_map);
  return diff;
}

/* this version assumes bi-directional graph
// pagerank: accumulate diff between prop array and rep array
float acc_diff(node_t *n_ptr, edge_t* e_ptr, 
               prop_f48_t* prop_ptr, prop_f48_t* rep_ptr) {
  float diff = 0.0;
  for (int i = 0; i < NUM_NODES; i++) {
    float48_t* ptr = (float48_t*)(&prop_ptr[i]);
    float prop_val = getval_f48(ptr);
    // need to find location of a rep array entry for this node
    
    float rep_val = 0.0;
    uint64_t edge_idx = getval48(&n_ptr[i]); // take first nbr
    uint64_t nbr_node = getval48(&e_ptr[edge_idx]); // look at that nbr's nbrs
    for (uint64_t j = getval48(&n_ptr[nbr_node]); 
                  j < getval48(&n_ptr[nbr_node+1]); j++) {
      if (i == getval48(&e_ptr[j])) {
        float48_t* r_ptr = (float48_t*)(&rep_ptr[j]);
        rep_val = getval_f48(r_ptr);
        break;
      }
    }
    diff += abs(prop_val - rep_val);
  }
  return diff;
}
*/

/*
// Optimized scatter
// 1. Make shadow structures (node list, edge list) on host that use proper 4-byte
//    aligned addresses
//    uint32_t* node_list, edge_list;
// 2. Perform parallel scatter using pragma omp
      #pragma omp parallel
      for (int i = 0; i < NUM_NODES; i++) [
        for (int j = node_list[i]; j < node_list[j]; j++) {
          uint64_t n = edge_list[j];
          rep[j] = deg[n]; //could do several updates do different replicated vectors here
*/
uint32_t* node_list; 
uint32_t* edge_list; 
float*    prop_ar; 
float*    rep_ar;
void prepare_scat(node_t *n_ptr, edge_t *e_ptr, 
                 prop_f48_t* prop_ptr, prop_f48_t* rep_ptr) {
  node_list = (uint32_t*)malloc(NUM_NODES*sizeof(uint32_t));
  prop_ar   = (float*)malloc(NUM_NODES*sizeof(float));
  edge_list = (uint32_t*)malloc(NUM_EDGES*sizeof(uint32_t));
  rep_ar    = (float*)malloc(NUM_EDGES*sizeof(float));

  for (int i = 0; i < NUM_NODES; i++) {
    node_list[i] = getval48(&n_ptr[i]);
    prop_ar[i]   = getval_f48(&prop_ptr[i]);
  }
  for (int j = 0; j < NUM_EDGES; j++) {
    edge_list[j] = getval48(&e_ptr[j]);
    //rep_ar[j]    = getval_f48(&rep_ptr[j]);
  }
}
void test_scat() { // prints a few values to ensure scatter is correct
  for (int i = 0; i < 32; i++) {
    printf("node %d: [ ", i);
    for (int j = node_list[i]; j < node_list[i+1]; j++) {
      printf("%d(%f)  ", edge_list[j], rep_ar[j]);
    }
    printf("] pagerank=%f\n", prop_ar[i]);
  }
}
void scatter_opt() {
  int i=0, j=0;
  omp_set_num_threads(16);
  int chunk = NUM_NODES/12;
#pragma omp parallel private(i,j) shared(node_list,edge_list,rep_ar,prop_ar)
  {
#pragma omp for nowait //schedule(static)
    for (i = 0; i < NUM_NODES; i++) {
      for (j = node_list[i]; j < node_list[i+1]; j++) {
        //uint32_t dn = edge_list[j];
        //rep_ar[j] = prop_ar[dn];
        rep_ar[j] = prop_ar[edge_list[j]];
      }
    }
  }
}
void scatter_opt2() {
  int j=0;
  omp_set_num_threads(12);
  int chunk = NUM_NODES/12;
#pragma omp parallel private(j) shared(node_list,edge_list,rep_ar,prop_ar)
  {
#pragma omp for nowait //schedule(static)
    for (j = 0; j < NUM_EDGES; j++) {
      //uint32_t dn = edge_list[j];
      //rep_ar[j] = prop_ar[dn];
      rep_ar[j] = prop_ar[edge_list[j]];
    }
  }
}
void end_scat() {
  free(node_list);  free(edge_list);
  free(prop_ar);    free(rep_ar);
}




// scatter values from prop array to rep array
void rep_scatter(node_t *n_ptr, edge_t *e_ptr, 
                 prop_f48_t* prop_ptr, prop_f48_t* rep_ptr) {
  for (int i = 0; i < NUM_NODES; i++) {
    for (uint64_t j = getval48(&n_ptr[i]); j < getval48(&n_ptr[i+1]); j++) {
      uint64_t dn = getval48(&e_ptr[j]);
      float48_t* p_ptr = (float48_t*)(&prop_ptr[dn]);
      float prop_val = getval_f48(p_ptr);
      float48_t* r_ptr = (float48_t*)(&rep_ptr[j]);
      assign_f48(r_ptr, prop_val);
    }
  }
}



void print_graph(node_t *n_ptr, int nb, int eb, int p24b, int p48b, int b_len,
                 const char *title_str) {
  FILE *f;
  int len = strlen(title_str);
  char *newstr = malloc(len+5);
  strcpy(newstr, title_str);
  strcat(newstr, ".txt");
  if (NULL == (f = fopen(newstr, "w"))) {
    fprintf(stderr, "Error opening file %s for writing.\n", newstr);
    cleanup();
    exit(EXIT_FAILURE);
  }

  edge_t *e_ptr =  (edge_t*)n_ptr +  (nb*b_len)/sizeof(*e_ptr);

  prop_uint24_t *p_ptr0 = (prop_uint24_t*)e_ptr + ((eb*b_len)/sizeof(prop_uint24_t));
  prop_f48_t    *p_ptr1 = (prop_f48_t*)p_ptr0 + ((p24b*b_len)/sizeof(prop_f48_t));
  prop_f48_t    *p_ptr2 = p_ptr1 + ((p48b*b_len)/sizeof(prop_f48_t));

  fprintf(f, title_str);
  fprintf(f, "\n\n-------------------------------------\n\n");
  for (int i = 0; i < NUM_NODES; i++) {
    fprintf(f, "node %d:  [ ", i);
    for (uint64_t j = getval48(&n_ptr[i]); j < getval48(&n_ptr[i+1]); j++) {
      fprintf(f, "%lld(%d|%f) ", getval48(&e_ptr[j]), 
                                 getdeg16(&p_ptr2[j]),
                                 getval_f48(&p_ptr2[j]) );
    }
    fprintf(f, "]  page_rank=%f  \n", getval_f48(&p_ptr1[i]));
  }
  fprintf(f,"\n\n-------------------------------------\n\n");
  fclose(f);
  free(newstr);
}



int main(int argc, char **argv) {

  // assign NUM_NODES and NUM_EDGES
  NUM_EDGES = NUM_NODES*AVG_DEG;
  if (argc == 2) 
    read_graph_size(argv[1]);

  int burst_len_bytes = 192; //384;  //max_get_burst_size(max_file_t*..)
  int node_bytes = (NUM_NODES+1)*sizeof(node_t);
  int node_bursts = ceil((double)node_bytes/(double)burst_len_bytes);
  int edge_bytes = (NUM_EDGES)*sizeof(edge_t);
  int edge_bursts = ceil((double)edge_bytes/(double)burst_len_bytes);

  // Application-specific: Betweeness Centrality
  // Node properties: levels , page_rank
  // [ lvls ] 
  // lvls: 24-bit uints    others: 48-bit floats (highest 16 bits are zero)
  // 4 property arrays needed
  int num_prop24s = 1;
  int num_prop48s = 1; // page_rank property: 48-bit float
  int num_rep48s  = 1; // replicated arrays (num_edges elements, instead of N)
  int prop24_bytes = (NUM_NODES)*sizeof(prop_uint24_t);
  int prop48_bytes = (NUM_NODES)*sizeof(prop_f48_t);
  int rep48_bytes  = (NUM_EDGES)*sizeof(prop_f48_t);
  int prop24_bursts = ceil((double)prop24_bytes/(double)burst_len_bytes);
  int prop48_bursts = ceil((double)prop48_bytes/(double)burst_len_bytes);
  int rep48_bursts  = ceil((double)rep48_bytes/(double)burst_len_bytes);
  int prop_bursts = num_prop24s*prop24_bursts + 
                    num_prop48s*prop48_bursts +
                    num_rep48s *rep48_bursts;

  int size_bytes = (node_bursts + 
                    edge_bursts + 
                    prop_bursts) * burst_len_bytes; 
  printf("N = %d, M = %d\n", NUM_NODES, NUM_EDGES);
  //printf("\nsize_bytes=%d , n_bursts=%d , e_bursts=%d , p_bursts=%d\n", 
  //       size_bytes, node_bursts, edge_bursts, prop_bursts);

	data_in  = malloc(size_bytes);
	data_out = malloc(size_bytes);
	if(!data_in || !data_out) {
		fprintf(stderr, "Failed to allocate memory for data I/O.\n");
		return 1;
	}
  memset(data_out, 0, size_bytes);
  node_t *node_ptr = (node_t*)data_in;
  edge_t *edge_ptr = (edge_t*)data_in + ((node_bursts*burst_len_bytes)/sizeof(edge_t));
  //printf("node_ptr: %p\nedge_ptr: %p\n", node_ptr, edge_ptr);


  if (argc == 2) {
    // Read graph in from edge list
    read_edge_list(node_ptr, edge_ptr, NUM_NODES, NUM_EDGES, argv[1]);
  }
  else {
    // Generate graph
    create_uniform_random_graph(node_ptr, edge_ptr, NUM_NODES, NUM_EDGES);

    //create_RMAT_graph(node_ptr, edge_ptr, 0.4, 0.25, 0.25, NUM_NODES, NUM_EDGES);
      // 4M node, 32M edges: pagerank was very slow (~7s runtime per iteration)
      // and did not converge
    //create_RMAT_graph(node_ptr, edge_ptr, 0.6, 0.20, 0.15, 
    //                  NUM_NODES, NUM_EDGES);
  }


  // Initialize Property arrays
  //prop24
  prop_uint24_t *prop_ptr0 = (prop_uint24_t*)edge_ptr + 
                        ((edge_bursts*burst_len_bytes)/sizeof(prop_uint24_t));
  //prop48
  prop_f48_t *prop_ptr1 = (prop_f48_t*)prop_ptr0 +  // pointer to property array
                        ((prop24_bursts*burst_len_bytes)/sizeof(prop_f48_t));
  //rep48
  prop_f48_t *prop_ptr2 = prop_ptr1 + // pointer to the replicated array
                        ((prop48_bursts*burst_len_bytes)/sizeof(prop_f48_t));

  // page_rank initial value: 1/N
  float init_val = (float)1/(float)NUM_NODES;
  init_prop_f48(prop_ptr1, NUM_NODES, init_val); //page_rank property array
    prepend_deg(prop_ptr1, node_ptr);  //prepend degree in high 16 bits
  init_rep_f48(prop_ptr2, NUM_EDGES, init_val); //page_rank replicated array
    prepend_deg_rep(prop_ptr2, node_ptr, edge_ptr); //<--problem here
    div_by_deg_rep(prop_ptr2, node_ptr);

  printf("Finished generating input graph...\n");fflush(stdout);
  //print_graph(node_ptr, node_bursts, edge_bursts, 
  //            prop24_bursts, prop48_bursts, burst_len_bytes,
  //            "InputGraph");


  // Initialize timing data structures
  struct timespec ta_start, ta_end, ta_d;
  struct timespec tb_start, tb_end, tb_d;
  uint64_t elapsed_nsec;

  printf("Initializing maxfile...");fflush(stdout);
  max_file_t *maxfile;
  if (NULL == (maxfile = Pagerank_init())) {
    fprintf(stderr, "Problem initializing maxfile.\n"); fflush(stderr);
    exit(EXIT_FAILURE);
  }
  printf("done\n");fflush(stdout);

  printf("Creating engine, loading maxfile...");fflush(stdout);
  max_engine_t *eng;
  if (NULL == (eng = max_load(maxfile, "*"))) {
    fprintf(stderr, "Problem creating engine.\n"); fflush(stderr);
    max_file_free(maxfile);
    exit(EXIT_FAILURE);
  }
  printf("done\n\n");fflush(stdout);

  printf("Preparing params for loading data into DFE memory...");fflush(stdout);
  int num_elems_32 = size_bytes / sizeof(uint32_t);
    Pagerank_writeLMem_actions_t wr_actions;
    wr_actions.param_size = num_elems_32;
    wr_actions.param_start = 0;
    wr_actions.instream_fromcpu = (uint32_t*)data_in;
  printf("done\n");fflush(stdout);

  //PageRank parameters
  // e:   acceptable error
  // d:   damping factor (kernel parameter)
  // max: maximum number of iterations allowed for convergence
  float d = 0.85;
  float e = 0.001;
  int max = 5;
  float term = (1.0 - d)/(float)NUM_NODES;

  //Debug
  uint64_t numBrsts0;
  uint64_t numBrsts1;
  uint64_t numBrsts2;
  uint64_t numBrsts3;
  uint64_t numBrsts4;
  printf("Preparing params for running engine on FPGA...");fflush(stdout);
    Pagerank_actions_t grph_actions;
    grph_actions.param_nodeAddr  = 0;
    // address of normal page_rank array
    grph_actions.param_propAddr  = node_bursts + edge_bursts +  
                                   prop24_bursts;   // lvl array
    // address of replicated page_rank array
    grph_actions.param_repAddr   = node_bursts + edge_bursts +  
                                   prop24_bursts +  // lvl array
                                   prop48_bursts;   // BCRep array
    grph_actions.param_NumNodes  = NUM_NODES; 
    grph_actions.param_StartCnt  = 100; 
    grph_actions.param_prTerm    = term;  // pagerank term
    grph_actions.param_d         = d;     // damping factor: d
    //cycle at which to send stop interrupt
    grph_actions.param_StopCnt   = (512*1024*1024);
    //number of cycles for done signal stability
    grph_actions.param_uDVal     = 20;

    //Scalar Outputs
    grph_actions.outscalar_MemUnit0_numBrsts    = &numBrsts0;
    grph_actions.outscalar_MemUnit1_numBrsts    = &numBrsts1;
    grph_actions.outscalar_MemUnit2_numBrsts    = &numBrsts2;
    grph_actions.outscalar_MemUnit3_numBrsts    = &numBrsts3;
    grph_actions.outscalar_MemUnit4_numBrsts    = &numBrsts4;

  max_config_set_string(MAX_CONFIG_DEBUG_DIRECTORY,  "./");
  printf("done\n");fflush(stdout);

  printf("Preparing params for reading DFE memory...");fflush(stdout);
    Pagerank_readLMem_actions_t rd_actions;
    rd_actions.param_size = num_elems_32;
    rd_actions.param_start = 0;
    rd_actions.outstream_tocpu = (uint32_t*)data_out;
  printf("done\n");fflush(stdout);

  // global timers
  uint64_t wr_time=0, run_time=0, rd_time=0, diff_time=0, scatter_time=0;
  uint64_t opt_scat_time=0;

  int iter = 0; // iteration number
  float diff = 0.0;

  clock_gettime(CLOCK_REALTIME, &tb_start);
  do {
    diff = 0.0;
    printf("\nBeginning iteration %d..\n", iter);

    //   write to fpga
    printf("Loading data into DFE memory...");fflush(stdout);
      clock_gettime(CLOCK_REALTIME, &ta_start);
      Pagerank_writeLMem_run(eng, &wr_actions);
      clock_gettime(CLOCK_REALTIME, &ta_end);
    ta_d = ts_diff(ta_start, ta_end);
    elapsed_nsec = ta_d.tv_sec*(1E9) + ta_d.tv_nsec;
    printf("done..%ld:%ld -> %ldns\n",
           ta_d.tv_sec, ta_d.tv_nsec, elapsed_nsec);fflush(stdout);
    wr_time += elapsed_nsec;

    //   call computation
    printf("----\nRunning engine...");fflush(stdout);
      clock_gettime(CLOCK_REALTIME, &ta_start);
      Pagerank_run( eng , &grph_actions ); // Advanced Static Interface
      clock_gettime(CLOCK_REALTIME, &ta_end);
    printf("done\n");fflush(stdout);
    printf("End: %ld:%ld\n", ta_end.tv_sec, ta_end.tv_nsec);
    ta_d = ts_diff(ta_start, ta_end);
    elapsed_nsec = ta_d.tv_sec*(1E9) + ta_d.tv_nsec;
    printf("Process time: %ld:%ld -> %ldns\n----\n",
           ta_d.tv_sec, ta_d.tv_nsec, elapsed_nsec);
    run_time += elapsed_nsec;

    // read from fpga
    printf("Reading DFE memory...");fflush(stdout);
      clock_gettime(CLOCK_REALTIME, &ta_start);
      Pagerank_readLMem_run(eng, &rd_actions);
      clock_gettime(CLOCK_REALTIME, &ta_end);
    ta_d = ts_diff(ta_start, ta_end);
    elapsed_nsec = ta_d.tv_sec*(1E9) + ta_d.tv_nsec;
    printf("done..%ld:%ld -> %ldns\n",
           ta_d.tv_sec, ta_d.tv_nsec, elapsed_nsec);fflush(stdout);
    rd_time += elapsed_nsec;

    // copy output data to input data
    memcpy((void*)data_in, (void*)data_out, size_bytes);

    // accumulate the diff: old values are still in rep array
    printf("Accumulating diff...");fflush(stdout);
      clock_gettime(CLOCK_REALTIME, &ta_start);
      diff = acc_diff(node_ptr, edge_ptr, prop_ptr1, prop_ptr2);
      clock_gettime(CLOCK_REALTIME, &ta_end);
    ta_d = ts_diff(ta_start, ta_end);
    elapsed_nsec = ta_d.tv_sec*(1E9) + ta_d.tv_nsec;
    printf("done..%ld:%ld -> %ldns\n",
           ta_d.tv_sec, ta_d.tv_nsec, elapsed_nsec);fflush(stdout);
    diff_time += elapsed_nsec;

    // scatter: from prop array to rep array
    printf("Scattering from prop array to rep array and dividing...");fflush(stdout);
      clock_gettime(CLOCK_REALTIME, &ta_start);
      rep_scatter(node_ptr, edge_ptr, prop_ptr1, prop_ptr2);
      div_by_deg_rep(prop_ptr2, node_ptr);  // rep arrays store pr/deg
      clock_gettime(CLOCK_REALTIME, &ta_end);
    ta_d = ts_diff(ta_start, ta_end);
    elapsed_nsec = ta_d.tv_sec*(1E9) + ta_d.tv_nsec;
    printf("done..%ld:%ld -> %ldns\n",
           ta_d.tv_sec, ta_d.tv_nsec, elapsed_nsec);fflush(stdout);
    scatter_time += elapsed_nsec;

    // Optimized scatter function
    printf("Preparing optimized scatter...");fflush(stdout);
    prepare_scat(node_ptr, edge_ptr, prop_ptr1, prop_ptr2);
    printf("done...\n");fflush(stdout);
    //test_scat();
    printf("Optimized scatter for time...");fflush(stdout);
      clock_gettime(CLOCK_REALTIME, &ta_start);
    scatter_opt2();
    //scatter_opt();
      clock_gettime(CLOCK_REALTIME, &ta_end);
    //test_scat();
    end_scat();
    ta_d = ts_diff(ta_start, ta_end);
    elapsed_nsec = ta_d.tv_sec*(1E9) + ta_d.tv_nsec;
    printf("done..%ld:%ld -> %ldns\n",
           ta_d.tv_sec, ta_d.tv_nsec, elapsed_nsec);fflush(stdout);
    opt_scat_time += elapsed_nsec;


    printf("Ending iteration %d..diff=%f\n", iter, diff);
    iter++;
  } while ((diff > e) && (iter < max));

  clock_gettime(CLOCK_REALTIME, &tb_end);
  tb_d = ts_diff(tb_start, tb_end);
  elapsed_nsec = tb_d.tv_sec*(1E9) + tb_d.tv_nsec;
  printf("\n----\nOverall time: %ld:%ld -> %ldns\n\n",
         tb_d.tv_sec, tb_d.tv_nsec, elapsed_nsec);
  printf("run_time = %f ms\n", (float)run_time/(float)1000000);
  printf("diff_time = %f ms\n", (float)diff_time/(float)1000000);
  printf("scatter_time = %f ms\n", (float)scatter_time/(float)1000000);
  printf("opt_scatter_time = %f ms\n", (float)opt_scat_time/(float)1000000);
  printf("\nwrite_to_fpga = %f ms\n", (float)wr_time/(float)1000000);
  printf("read_from_fpga = %f ms\n", (float)rd_time/(float)1000000);

  printf("MemUnits:   num0Brsts:%d, numBrsts1:%d, numBrsts2:%d, numBrsts3:%d, numBrsts4:%d\n",
         numBrsts0, numBrsts1, numBrsts2, numBrsts3, numBrsts4);

  printf("\nUnloading engine...");fflush(stdout);
  max_unload(eng);
  printf("done\n");fflush(stdout);

  print_graph((node_t*)data_in, node_bursts, 
              edge_bursts, prop24_bursts, prop48_bursts,
              burst_len_bytes, "OutputGraph");

  cleanup();

  return 0;
}
