/*
 *  Main header file for percolation code.
 */

/*
 *  System size L
 */

//#define L 432

/*
 *  Use 1D decomposition over NPROC processes across first dimension
 *  For an LxL simulation, the local arrays are of size MxN
 */

#define NPROC 3

//#define M L/NPROC
//#define N L

/*
 *  Prototypes for supplied functions
 */

/*
 *  Visualisation
 */

void percwrite(char *percfile, int **map, int ncluster,int L);
void percwritedynamic(char *percfile, int **map, int l, int ncluster);
void initialize_oldmap(int **old, int M,int L);
void initialize_map(int **map, int L, double r, double rho, int *nhole);
void compute_newmap(int **old, int **new, int newval, int oldval, int M, int N, int *nchangelocal);
void test_percolation(int **map, int L, int *perc);

/*
 *  Random numbers
 */

void rinit(int ijkl);
float uni(void);
