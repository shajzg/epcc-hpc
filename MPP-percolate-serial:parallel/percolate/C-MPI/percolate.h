/*
 *  Main header file for percolation code.
 */

/*
 *  System size L
 */

#define L 768

/*
 *  Use 1D decomposition over NPROC processes across first dimension
 *  For an LxL simulation, the local arrays are of size MxN
 */

#define NPROC 4

#define M L/NPROC
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
