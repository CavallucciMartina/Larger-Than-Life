/****************************************************************************
 *
 * omp-ltl-serial.c - Larger than Life serial version ( GoL parameters)
 *
 * Written by Martina Cavallucci <martina.cavallucci(at)studio.unibo.it>
 * --------------------------------------------------------------------------
 *
 * Serial computation of the Larger Than Life
 * Compile with:
 *
 * gcc -fopenmp -std=c99 -Wall -Wpedantic ltl-serial.c -o ltl-serial
 *
 * (note that this is a serial version, so -fopenmp is actually
 * ignored, but will be required for a parallel version).
 *
 * Run with:
 * (Test it with the Game of Life parameter )
 * ./ltl-serial (R =) 1 (B1 = B2 =) 3 (D1 =) 3  (D2 = )4 nsteps input_file output_file
 *
 ****************************************************************************/
 #include "hpc.h"
 #include <stdio.h>
 #include <stdlib.h>
 #include <string.h>
 #include <assert.h>
 #include <ctype.h> /* for isdigit */

 typedef unsigned char cell_t;

 /* The struct of the bmap_t with a size and a point of unsigned char */
 typedef struct {
     int n;
     cell_t *bmap;
 } bmap_t;

 /* Returns a pointer to the cell of coordinates (i,j) in the bitmap grid
    with the HALO included.
  */

 cell_t *IDX(cell_t *grid, int n, int i, int j, int r)
 {
     //return bmap + i * n + j;
 	return grid + ((i)*(n +(2*r))+(j));
 }

/* Fill the ghost cells of |grid| in order to have cyclic boundary conditions*/

void fill_ghost( cell_t *grid, int n, int r)
{
	 int ISTART = r;
   int JSTART = r;
	 int IEND= n + r - 1;
   int JEND = n + r - 1;
    int i, j;
    /* copy top and bottom */
    for (j=JSTART; j<JEND; j++) {
      for( int k = 1; k < r + 1; k++){
        *IDX(grid, n, IEND + k, j, r) = *IDX(grid, n, ISTART + k - 1, j, r);
        *IDX(grid, n, ISTART - k, j, r) = *IDX(grid, n, IEND - k + 1, j, r);
      }
    }
    /* copy left and right */
    for (i=ISTART; i<IEND; i++) {
      for(int k=1; k< r + 1; k++){
        *IDX(grid, n, i, JEND + k, r) = *IDX(grid, n, i, JSTART + k - 1, r);
          *IDX(grid, n, i, JSTART - k, r) = *IDX(grid, n, i, JEND + k - 1, r);
    }
  }
    /* copy corners */
      for (i = 0; i < r ; i++) {
        for (j = 0; j < r; j++) {
              // *IDX(grid, n, i ,j , r) = // 0 0 // 0 1  // 1 0 // 1  1
               *IDX(grid, n, i ,j , r) = *IDX(grid, n, i + r + 1, j + r + 1 , r);
              /// *IDX(grid, n, i + r + n ,j + r + n, r) = //5 5 // 6 6 // 5  6 // 6 5
              *IDX(grid, n, i + r + n ,j + r + n, r) = *IDX(grid, n,i + r, j + r, r);
               //*IDX(grid, n, i ,j + r + n, r) = // 0 5 // 0 6 // 1  5  // 1  6
                *IDX(grid, n, i ,j + r + n, r) = *IDX(grid, n, i + r + 1, j + r  , r);
            //   *IDX(grid, n, i + r + n ,j, r) = // 5 0 // 6 0 // 5  1 // 6  1
                *IDX(grid, n, i + r + n ,j, r) = *IDX(grid, n,i + r , j + r + 1 , r);
      }
    }
}

  /**
  * Write the content of the bmap_t structure pointed to by ltl to the
  * file f in PBM format and and allocates space
  *for the ghost cell that will be assigned.
  * The caller is responsible for passing a
  * pointer f to a file opened for writing
  */
 void write_ltl( bmap_t* grid, FILE *f , int  r )
 {
     int i, j;
     const int n = grid->n;
     fprintf(f, "P1\n");
     fprintf(f, "# produced by ltl\n");
     fprintf(f, "%d %d\n", n, n);
     for (i= r ; i< n + r; i++) {
        for (j= r ; j< n + r; j++) {
             fprintf(f, "%d ", *IDX(grid->bmap, n, i, j, r));
         }
         fprintf(f, "\n");
     }
 }
 /*Count how many neighbors are living in the range of i j element */
 int count_neighbors(cell_t* cur,int n, int i, int j, int r )
 {
   int nbors = 0;
  for(int k = 1 ; k < r + 1 ; k++){
      for(int f= 1 ; f < r + 1 ; f++){
          nbors = nbors +
                 *IDX(cur,n, i-k,j-f, r) + *IDX(cur, n, i-k,j, r) + *IDX(cur, n,i-k,j+f, r) +
                 *IDX(cur,n,i  ,j-f,r) +                            *IDX(cur,n,i  ,j+f, r) +
                 *IDX(cur,n,i + k,j - f,r) + *IDX(cur,n,i+k,j,r) + *IDX(cur,n,i+k,j+f,r);
          }
}

	return nbors;
}
 /*Compute of the Larger than life*/
 void compute_ltl( cell_t *cur, cell_t *next, int nc, int r, int b1, int b2, int d1, int d2)
 {
    int i, j, element, c = 0;
    const int n = nc;
   	 for (i = r; i < n + r ; i++)
   	  {
       for (j = r; j < n + r; j++)
        {
           	c = count_neighbors(cur, n, i, j , r);
            element = *IDX(cur, n, i, j, r);

            if( !element && c >= b1 && c <= b2){ // if it can relive
               *IDX(next, n, i, j, r) = 1; //Set it as live
            }else if( element && c + 1 >= d1 && c + 1 <= d2) // if the cell remaining live
            {
            	*IDX(next, n, i, j,r ) = 1;
            }else{
              *IDX(next, n, i, j, r) = 0; // set it as died
            }
       }
     }
 	}

 /**
  * Read a PBM file from file f. The caller is responsible for passing
  * a pointer f to a file opened for reading. This function is not very
  * robust; it may fail on perfectly legal PBM images, but should work
  * for the images produced by gen-input.c. Also, it should work with
  * PBM images produced by Gimp (you must save them in "ASCII format"
  * when prompted).
  */
 void read_ltl( bmap_t *ltl, FILE* f, int r)
 {

    char buf[2048];
    char *s;
    int n, i, j;
    int width, height;

     /* Get the file type (must be "P1") */
    s = fgets(buf, sizeof(buf), f);
    if (0 != strcmp(s, "P1\n")) {
        fprintf(stderr, "FATAL: Unsupported file type \"%s\"\n", buf);
        exit(-1);
    }
    /* Get any comment and ignore it; does not work if there are
       leading spaces in the comment line */
    do {
        s = fgets(buf, sizeof(buf), f);
    } while (s[0] == '#');
    /* Get width, height; since we are assuming square images, we
       reject the input if width != height. */
    sscanf(s, "%d %d", &width, &height);
    if ( width != height ) {
        fprintf(stderr, "FATAL: image width (%d) and height (%d) must be equal\n", width, height);
        exit(-1);
    }
    ltl->n = n = width;
    const int ISTART = r;
  	const int IEND   = n + r;
  	const int JSTART = r;
  	const int JEND   = n + r ;
    int ng = n + (2*r);
    ltl->bmap = (cell_t*)malloc( ng * ng * sizeof(cell_t));
    /* scan bitmap; each pixel is represented by a single numeric
       character ('0' or '1'); spaces and other separators are ignored
       (Gimp produces PBM files with no spaces between digits) */
    for (i = ISTART; i < IEND ; i++) {
         for (j = JSTART; j < JEND ; j++) {
            int val;
            do {
                val = fgetc(f);
                if ( EOF == val ) {
                    fprintf(stderr, "FATAL: error reading input\n");
                    exit(-1);
                }
            } while ( !isdigit(val) );
            *IDX(ltl->bmap, n, i, j, r) = (val - '0');
        }
    }
 }

 int main( int argc, char* argv[] )
 {
     int R, B1, B2, D1, D2, nsteps,s;
     const char *infile, *outfile;
     FILE *in, *out;
     bmap_t cur, next;
     double tstart, tend;

     if ( argc != 9 ) {
         fprintf(stderr, "Usage: %s R B1 B2 D1 D2 nsteps infile outfile\n", argv[0]);
         return -1;
     }
     R = atoi(argv[1]);
     B1 = atoi(argv[2]);
     B2 = atoi(argv[3]);
     D1 = atoi(argv[4]);
     D2 = atoi(argv[5]);
     nsteps = atoi(argv[6]);
     infile = argv[7];
     outfile = argv[8];

     assert(  R <= 8  );
     assert(  0 <= B1 );
     assert( B1 <= B2 );
     assert(  1 <= D1 );
     assert( D1 <= D2 );

     in = fopen(infile, "r");
     if (in == NULL) {
         fprintf(stderr, "FATAL: can not open \"%s\" for reading\n", infile);
         exit(-1);
     }

     read_ltl(&cur, in, R);
     fclose(in);

     fprintf(stderr, "Size of input image: %d x %d\n", cur.n, cur.n);
     fprintf(stderr, "Model parameters: R=%d B1=%d B2=%d D1=%d D2=%d nsteps=%d\n",
             R, B1, B2, D1, D2, nsteps);

     out = fopen(outfile, "w");
     if ( out == NULL ) {
         fprintf(stderr, "FATAL: can not open \"%s\" for writing", outfile);
         exit(-1);
     }
      const int n = cur.n;
       int ng = n + (2*R);
        next.n = n;
     next.bmap = (cell_t*)malloc( ng * ng * sizeof(cell_t));

  /*   const int ISTART = R;
    const int IEND   = next.n + R - 1;
    const int JSTART = R;
    const int JEND   = next.n + R - 1;
     for (int i = ISTART; i < IEND + 1 ; i++) {
          for (int j = JSTART; j < JEND + 1 ; j++) {

          printf("%d",*IDX(next.bmap, next.n, i, j, R));
         }
       }*/
     tstart = hpc_gettime();

    for (s = 0; s < nsteps; s++) {
    	bmap_t tmp;
    	fill_ghost(cur.bmap, n , R);
     	compute_ltl(cur.bmap, next.bmap, n, R, B1, B2, D1, D2);
        tmp = cur;
        cur = next;
        next = tmp;

    }


    tend = hpc_gettime();
    fprintf(stderr, "Execution time %f\n", tend - tstart);

    write_ltl(&cur, out, R);

    fclose(out);

    free(cur.bmap);
    free(next.bmap);

     return 0;
 }
