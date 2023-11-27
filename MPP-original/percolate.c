#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <unistd.h>
#include "arralloc.h"
#include "percolate.h"

int main(int argc, char *argv[])
{
  /*
   * Initialize basic parameters to be used.
   */

  int **map;
  int L = 432;
  int seed = 6543;
  int i, j, nhole, step, oldval, newval;
  int nchangelocal, nchange, printfreq;
  int ilft, irt, perc;
  int size, rank;
  int tag = 1;
  int direction, disp;
  int ndims = 2;
  int dims[ndims];
  int period[ndims];
  int reorder = 0;
  int left, right, up, down;
  int height, width;
  int count, blocklength;
  int coords[2];
  int M, N;
  int opt;
  double rho, r;

  /*
   * Initialize MPI routines.
   */

  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);
  MPI_Request sendrequest[size];
  MPI_Request newrequest, request[4];
  MPI_Datatype vector, halo_j, halo_i;
  MPI_Status status, stats[4];
  MPI_Comm comm2d;

  /*
   * Initialize parameter for topology. Set index j to be periodic.
   */

  direction = 1;
  disp = 1;
  dims[0] = 0;
  dims[1] = 0;
  period[0] = 0;
  period[1] = 1;

  /*
   * Create dims and Cartesian Topology.
   */

  MPI_Dims_create(size, ndims, dims);
  MPI_Cart_create(comm, ndims, dims, period, reorder, &comm2d);
  MPI_Cart_shift(comm2d, direction, disp, &left, &right);
  MPI_Cart_shift(comm2d, 0, disp, &up, &down);

  /*
   * Initalize seed and mapsize using arguments from input.
   */

  while ((opt = getopt(argc, argv, ":s:l:")) != -1)
  {
    switch (opt)
    {
    case 'l':
      L = atoi(optarg);
      break;
    case 's':
      seed = atoi(optarg);
      break;
    }
  }

  if(L < 0){
    printf("Map size should be no smaller than zero");
    MPI_Finalize();
    return 0;
  }
  if(seed < 0){
    printf("Seed should be no smaller than zero");
    MPI_Finalize();
    return 0;
  }

  /*
   * Initialize map.
   */

  map = (int **)arralloc(sizeof(int), 2, L + size, L + size);

  /*
   * Initialize parameters for derived data type vector and submit it to MPI.
   * The count and blocklength are set to be ceiling of integer part to solve 
   * uneven distribution.
   */

  height = dims[0];
  width = dims[1];
  count = (int)ceil((double)L / height);
  blocklength = (int)ceil((double)L / width);
  M = count;
  N = blocklength;
  MPI_Type_vector(count, blocklength, L + size, MPI_INT, &vector);
  MPI_Type_commit(&vector);

  MPI_Type_vector(1, blocklength, N + 2, MPI_INT, &halo_j);
  MPI_Type_commit(&halo_j);
  MPI_Type_vector(count, 1, N + 2, MPI_INT, &halo_i);
  MPI_Type_commit(&halo_i);

  /*
   * Initialize old map. The 2-d array new is used as a temporary map to update old map.
   */

  int **old, **new, **smallmap;
  new = (int **)arralloc(sizeof(int), 2, M + 2, N + 2);
  smallmap = (int **)arralloc(sizeof(int), 2, M + 2, L + size);
  old = (int **)arralloc(sizeof(int), 2, M + 2, N + 2);

  initialize_oldmap(old, M, N);

  printfreq = 100;

  if (rank == 0)
  {
    printf("percolate: running on %d process(es)\n", size);

    /*
      *  Set most important value: the rock density rho (between 0 and 1)
      */

    rho = 0.411;

    printf("percolate: L = %d, rho = %f, seed = %d\n",
           L, rho, seed);

    rinit(seed);

    nhole = 0;

    /*
      *  Initialise map with density rho. Zero indicates rock, a positive
      *  value indicates a hole. For the algorithm to work, all the holes
      *  must be initialised with a unique integer
      */

    initialize_map(map, L, r, rho, &nhole);

    printf("percolate: rho = %f, actual density = %f\n",
           rho, 1.0 - ((double)nhole) / ((double)L * L));
  }

  /*
    * Controller process distribute map to other processes. 
    * Each process has an old map which is a rectangular part of map.
    * Coordinate of the first point of rectangle is calculated according
    * to the coordinate of the process in the topology.
    */

  if (rank == 0)
  {
    for (int k = 0; k < size; k++)
    {
      MPI_Cart_coords(comm2d, k, 2, coords);

      int x = coords[0] * count;
      int y = coords[1] * blocklength;

      MPI_Bsend(&(map[x][y]), 1, vector, k, 0, comm);
    }
  }

  /*
    * Each process receives data from process 0 using old map.
    */

  MPI_Recv(&smallmap[1][1], 1, vector, 0, 0, comm, &status);

  for (int i = 1; i <= M; i++)
  {
    for (int j = 1; j <= N; j++)
    {
      old[i][j] = smallmap[i][j];
    }
  }

  /*
    * Initialize step and number of changes.
    */

  step = 1;
  nchange = 1;

  /*
    * Start timing. Before timing, add barrier to all processes.
    */

  MPI_Barrier(MPI_COMM_WORLD);
  double time1 = MPI_Wtime();

  /*
    * Start to execute the core part of the code, which includes swapping halos
    * and calculating new maps. The while loop stops when there is no change any more.
    */

  MPI_Cart_coords(comm2d, rank, 2, coords);

  int x = coords[0] * count;
  int y = coords[1] * blocklength;
  int row = M, col = N;

  /*
    * Ensure size of old map plus its coordinate in map can not surpass mapsize.
    * This step is to make sure that halos that we would send are correct.
    */

  if (x + M > L)
  {
    row = L - x;
  }

  if (y + N > L)
  {
    col = L - y;
  }

  while (nchange != 0)
  {

    /*
      * Send halos to the processes up, down, left and right respectively.
      */

    MPI_Issend(&old[1][col], 1, halo_i, right, tag, comm, &request[0]);
    MPI_Issend(&old[1][1], 1, halo_i, left, tag, comm, &request[1]);
    MPI_Issend(&old[1][1], 1, halo_j, up, tag, comm, &request[2]);
    MPI_Issend(&old[row][1], 1, halo_j, down, tag, comm, &request[3]);

    /*
      * Receive halos from processes up, down, left and right respectively.
      */

    MPI_Recv(&old[0][1], 1, halo_j, up, tag, comm, &stats[0]);
    MPI_Recv(&old[row + 1][1], 1, halo_j, down, tag, comm, &stats[1]);
    MPI_Recv(&old[1][0], 1, halo_i, left, tag, comm, &stats[2]);
    MPI_Recv(&old[1][col + 1], 1, halo_i, right, tag, comm, &stats[3]);

    /*
      * Initialize number of local changes.
      */

    nchangelocal = 0;

    /*
      * Start computation.
      */

    compute_newmap(old, new, newval, oldval, M, N, &nchangelocal);

    /*
       *  Compute global number of changes on every rank
       */

    MPI_Allreduce(&nchangelocal, &nchange, 1, MPI_INT, MPI_SUM, comm);

    /*
       *  Report progress every now and then
       */

    if (step % printfreq == 0)
    {
      if (rank == 0)
      {
        printf("percolate: number of changes on step %d is %d\n",
               step, nchange);
      }
    }

    /*
       * Wait for nonblocking send to finish.
       */

    MPI_Waitall(4, request, stats);

    /*
       *  Copy back in preparation for next step, omitting halos
       */

    for (i = 1; i <= M; i++)
    {
      for (j = 1; j <= N; j++)
      {
        old[i][j] = new[i][j];
      }
    }

    step++;
  }

  /*
   * Finish timing.
   */

  MPI_Barrier(MPI_COMM_WORLD);
  double time2 = MPI_Wtime();

  if (rank == 0)
  {
    printf(" time per step is %g\n", (time2 - time1) / step);
  }

  for (int i = 1; i <= M; i++)
  {
    for (int j = 1; j <= N; j++)
    {
      smallmap[i][j] = old[i][j];
    }
  }

  /*
   *  Now each processs sends old map back to map
   */

  MPI_Issend(&smallmap[1][1], 1, vector, 0, 0, comm, &newrequest);

  if (rank == 0)
  {

    for (int i = 0; i < size; i++)
    {
      MPI_Cart_coords(comm2d, i, 2, coords);

      int x = coords[0] * count;
      int y = coords[1] * blocklength;

      MPI_Recv(&map[x][y], 1, vector, i, 0, comm, &status);
    }
  }

  /*
   *  Test to see if percolation occurred by looking for positive numbers
   *  that appear on both the top and bottom edges
   */

  if (rank == 0)
  {

    perc = 0;

    test_percolation(map, L, &perc);

    if (perc != 0)
    {
      printf("percolate: cluster DOES percolate\n");
    }
    else
    {
      printf("percolate: cluster DOES NOT percolate\n");
    }

    /*
       *  Write the map to the file "map.pgm", displaying the two
       *  largest clusters. If the last argument here was 3, it would
       *  display the three largest clusters etc. The picture looks
       *  cleanest with only a single cluster, but multiple clusters
       *  are useful for debugging.
       */

    percwrite("map.pgm", map, 8, L);
  }

  MPI_Finalize();

  return 0;
}

void initialize_oldmap(int **old, int M, int N)
{
  for (int i = 0; i < M + 2; i++)
  {
    for (int j = 0; j < N + 2; j++)
    {
      old[i][j] = 0;
    }
  }
}

void initialize_map(int **map, int L, double r, double rho, int *nhole)
{
  for (int i = 0; i < L; i++)
  {
    for (int j = 0; j < L; j++)
    {
      r = uni();

      if (r < rho)
      {
        map[i][j] = 0;
      }
      else
      {
        *nhole = *nhole + 1;
        map[i][j] = *nhole;
      }
    }
  }
}

void compute_newmap(int **old, int **new, int newval, int oldval, int M, int N, int *nchangelocal)
{
  for (int i = 1; i <= M; i++)
  {
    for (int j = 1; j <= N; j++)
    {
      oldval = old[i][j];
      newval = oldval;

      /*                                                                                                                                                                                           
               * Set new[i][j] to be the maximum value of old[i][j]                                                                                                                                       
               * and its four nearest neighbours                                                                                                                                                         
               */

      if (oldval != 0)
      {
        if (old[i - 1][j] > newval)
          newval = old[i - 1][j];
        if (old[i + 1][j] > newval)
          newval = old[i + 1][j];
        if (old[i][j - 1] > newval)
          newval = old[i][j - 1];
        if (old[i][j + 1] > newval)
          newval = old[i][j + 1];

        if (newval != oldval)
        {
          *nchangelocal += 1;
        }
      }

      new[i][j] = newval;
    }
  }
}
void test_percolation(int **map, int L, int *perc)
{
  for (int ilft = 0; ilft < L; ilft++)
  {
    if (map[0][ilft] > 0)
    {
      for (int irt = 0; irt < L; irt++)
      {
        if (map[0][ilft] == map[L - 1][irt])
        {
          *perc = 1;
        }
      }
    }
  }
}
