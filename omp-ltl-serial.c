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

 /* This struct is defined here as an example; it is possible to modify
    this definition, or remove it altogether if you prefer to pass
    around the pointer to the bitmap directly. */
 typedef struct {
     int n;
     cell_t *bmap;
 } bmap_t;

 /* Returns a pointer to the cell of coordinates (i,j) in the bitmap
    bmap */
    /*Ritorna il valore di una cella date le coordinate */
 cell_t *IDX(cell_t *bmap, int n, int i, int j)
 {
     return bmap + i*n + j;
 }
  /**
  * Write the content of the bmap_t structure pointed to by ltl to the
  * file f in PBM format. The caller is responsible for passing a
  * pointer f to a file opened for writing
  */
  /*Scrivere le matrici di char con la mappa ltl*/
 void write_ltl( bmap_t* ltl, FILE *f )
 {
     int i, j;
     const int n = ltl->n;
     fprintf(f, "P1\n");
     fprintf(f, "# produced by ltl\n");
     fprintf(f, "%d %d\n", n, n);
     for (i=0; i<n; i++) {
         for (j=0; j<n; j++) {
             fprintf(f, "%d ", *IDX(ltl->bmap, n, i, j));
         }
         fprintf(f, "\n");
     }
 }
 /*Count how many neighbors are living in the range */
 int CountLivingNeighbors(cell_t* cur,int n, int i, int j, int r )
 {
	int count = 0;
	int x = i - r ;
	int y = j -r;
	if( j < r){
		y = 1;
	}
	if( i < r ){
		x = 1;
	}
	int z = (r*2) + 1;
	for (int d = 0 ; d < z ; d++ ){
		for( int k = 0; k < z ; k ++){

			if (*IDX(cur, n, x + d, y + k) == 1){
				count ++;

			}
		}
	}
	return count;
}
 /*Compute of the Larger than life*/
 void compute_ltl( cell_t *cur, cell_t *next, int nc, int r, int b1, int b2, int d1, int d2)
 {
    int i, j;
    const int n = nc;
   	 for (i = 0; i < n; i++)
   	  {
       for (j = 0; j < n; j++)
        {
          if(*IDX(cur, n, i, j) == 0)//if the cell is died
          {
           	int c = CountLivingNeighbors(cur, n, i, j , r);
            if( b1 <= c && c <= b2){ // if it can relive
               *IDX(next, n, i, j) = 1; //Set it as live
            }else // if the cell remaining died
            {
            	*IDX(next, n, i, j) = *IDX(cur, n, i, j) ; 
            }
          }
          if(*IDX(cur, n, i, j) == 1)//if is living
          {
            int c = CountLivingNeighbors(cur,n, i, j , r) + 1;
            if( d1 <= c && c <= d2){ // if it has to remain live
            	*IDX(next, n, i, j) = *IDX(cur, n, i, j) ; //Remaining live
            }else
            {
                *IDX(next, n, i, j) = 0; // set it as died
            }
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
 void read_ltl( bmap_t *ltl, FILE* f )
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
    ltl->bmap = (cell_t*)malloc( n * n * sizeof(cell_t));
    /* scan bitmap; each pixel is represented by a single numeric
       character ('0' or '1'); spaces and other separators are ignored
       (Gimp produces PBM files with no spaces between digits) */
    for (i = 0; i < n; i++) {
         for (j = 0; j < n; j++) {
            int val;
            do {
                val = fgetc(f);
                if ( EOF == val ) {
                    fprintf(stderr, "FATAL: error reading input\n");
                    exit(-1);
                }
            } while ( !isdigit(val) );
            *IDX(ltl->bmap, n, i, j) = (val - '0');
        }
    }
 }

 int main( int argc, char* argv[] )
 {
     int R, B1, B2, D1, D2, nsteps,s;
     const char *infile, *outfile;
     FILE *in, *out;
     bmap_t cur, next, tmp;
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
     read_ltl(&cur, in);
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
     next.n = n;
     next.bmap = (cell_t*)malloc( n * n * sizeof(cell_t));
     tstart = hpc_gettime();
    for (s = 0; s < nsteps; s++) {
     	compute_ltl(cur.bmap,next.bmap, n, R,B1,B2,D1,D2);
        tmp = cur;
        cur = next;
        next = tmp;
    }
    tend = hpc_gettime();
    fprintf(stderr, "Execution time %f\n", tend - tstart);

    write_ltl(&cur, out);

    fclose(out);

    free(cur.bmap);
    free(next.bmap);

     return 0;
 }
