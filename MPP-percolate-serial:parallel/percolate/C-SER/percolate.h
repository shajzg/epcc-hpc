/*
 *  Main header file for percolation exercise.
 */

/*
 *  System size L
 */

#define L 768

/*
 *  Although overall system is square, i.e. size L x L, we will define
 *  different variables for the first and second dimensions. This is
 *  because, in the parallel code, the local arrays will not be
 *  square. For example, using a simple 1D decomposition over NPROC
 *  processes, then M = L/NPROC and N = L
 */

#define M L
#define N L

/*
 *  Prototypes for supplied functions
 */

/*
 *  Visualisation
 */

void percwrite(char *mapfile, int map[L][L], int ncluster);
void percwritedynamic(char *mapfile, int **map, int l, int ncluster);

/*
 *  Random numbers
 */

void rinit(int ijkl);
float uni(void);
