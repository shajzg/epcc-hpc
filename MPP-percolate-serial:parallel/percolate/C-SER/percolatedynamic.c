#include <stdio.h>
#include <stdlib.h>

#include "percolate.h"
#include "arralloc.h"

/*
 * Simple serial program to test for percolation of a cluster.
 *
 * This version illustrates the use of dynamic array allocation via the
 * arralloc routine.
 *
 * In a real program, the array dimensions would be variables
 * (e.g. taken from the command line arguments). For simplicty, Here
 * we just use constant values as in the static allocation version.
 */

int main(int argc, char *argv[])
{
  /*
   *  Define the main arrays for the simulation
   */

  // int old[M+2][N+2], new[M+2][N+2];

  int **old, **new;

  /*
   *  Additional array WITHOUT halos for initialisation and IO. This
   *  is of size LxL because, even in our parallel program, we do
   *  these two steps in serial
   */

  // int map[L][L];

  int **map;

  /*
   *  Variables that define the simulation
   */

  int seed;
  double rho;

  /*
   *  Local variables
   */

  int i, j, nhole, step, maxstep, oldval, newval, nchange, printfreq;
  int jleft, jright, perc;
  double r;

  if (argc != 2)
    {
      printf("Usage: percolate <seed>\n");
      return 0;
    }

  /*
   *  Update for a fixed number of steps, periodically print progress
   */

  maxstep = 5*L;
  printfreq = 100;

  /*
   *  Set most important value: the rock density rho (between 0 and 1)
   */

  rho = 0.408;

  /*
   *  Set the randum number seed and initialise the generator
   */

  seed = atoi(argv[1]);

  printf("percolate: L = %d, rho = %f, seed = %d, maxstep = %d\n",
	 L, rho, seed, maxstep);

  /*
   * Allocate 2D LxL integer arrays dynamically
   */

  old = (int **) arralloc(sizeof(int), 2, M+2, N+2);
  new = (int **) arralloc(sizeof(int), 2, M+2, N+2);

  map = (int **) arralloc(sizeof(int), 2, L, L);

  if (NULL == old || NULL == new || NULL == map)
    {
      printf("percolate: array allocation failed\n");
      return 1;
    }

  rinit(seed);

  /*
   *  Initialise map with density rho. Zero indicates rock, a positive
   *  value indicates a hole. For the algorithm to work, all the holes
   *  must be initialised with a unique integer.
   */

  nhole = 0;

  for (i=0; i < L; i++)
    {
      for (j=0; j < L; j++)
	{
	  r=uni();
	  
	  if(r < rho)
	    {
	      map[i][j] = 0;
	    }
	  else
	    {
	      nhole++;
	      map[i][j] = nhole;
	    }
	}
    }

  printf("percolate: rho = %f, actual density = %f\n",
	 rho, 1.0 - ((double) nhole)/((double) L*L) );

  /*
   * Initialise the old array: copy the LxL array map to the centre of
   * old, and set the halo values to zero.
   */

  for (i=1; i <= M; i++)
    {
      for (j=1; j <= N; j++)
	{
	  old[i][j] = map[i-1][j-1];
	}
    }

  for (i=0; i <= M+1; i++)  // zero the bottom and top halos
    {
      old[i][0]   = 0;
      old[i][N+1] = 0;
    }

  for (j=0; j <= N+1; j++)  // zero the left and right halos
    {
      old[0][j]   = 0;
      old[M+1][j] = 0;
    }

  step = 1;
  nchange = 1;

  while (step <= maxstep)
    {
      nchange = 0;

      for (i=1; i<=M; i++)
	{
	  for (j=1; j<=N; j++)
	    {
	      oldval = old[i][j];
	      newval = oldval;

	      /*
	       * Set new[i][j] to be the maximum value of old[i][j]
	       * and its four nearest neighbours
	       */

	      if (oldval != 0)
		{
		  if (old[i-1][j] > newval) newval = old[i-1][j];
		  if (old[i+1][j] > newval) newval = old[i+1][j];
		  if (old[i][j-1] > newval) newval = old[i][j-1];
		  if (old[i][j+1] > newval) newval = old[i][j+1];

		  if (newval != oldval)
		    {
		      ++nchange;
		    }
		}

	      new[i][j] = newval;
	    }
	}

      /*
       *  Print progress every now and then
       */

      if (step % printfreq == 0)
	{
	  printf("percolate: changes on step %d is %d\n",
		 step, nchange);
	}

      /*
       *  Copy back in preparation for next step, omitting halos
       */

      for (i=1; i<=M; i++)
	{
	  for (j=1; j<=N; j++)
	    {
	      old[i][j] = new[i][j];
	    }
	}

      step++;
    }

  /*
   *  We set a maximum number of steps to ensure the algorithm always
   *  terminates. However, if we hit this limit before the algorithm
   *  has finished then there must have been a problem (e.g. the value
   *  of maxstep is too small)
   */

  if (nchange != 0)
    {
      printf("percolate: WARNING max steps = %d reached but nchange != 0\n",
	     maxstep);
    }

  /*
   *  Copy the centre of old, excluding the halos, into map
   */
  
  for (i=1; i<=M; i++)
    {
      for (j=1; j<=N; j++)
	{
	  map[i-1][j-1] = old[i][j];
	}
    }

  /*
   *  Test to see if percolation occurred by looking for positive numbers
   *  that appear on both the left and right-hand sides
   */

  perc = 0;

  for (jleft=0; jleft < L; jleft++)
    {
      if (map[0][jleft] > 0)
	{
	  for (jright=0; jright < L; jright++)
	    {
	      if (map[L-1][jright] == map[0][jleft])
		{
		  perc = 1;
		}
	    }
	}
    }

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

  percwritedynamic("map.pgm", map, L, 2);

  /*
   * Free the arrays
   */

  free(old);
  free(new);
  free(map);

  return 0;
}
