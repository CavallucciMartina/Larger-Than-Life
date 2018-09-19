/****************************************************************************
 *
 * cuda-ltl.cu - Larger than Life with CUDA using shared memory
 *
 * Written by Martina Cavallucci <martina.cavallucci(at)studio.unibo.it>
 * --------------------------------------------------------------------------
 *
 * This version of the Larger than Life uses shared memory and should
 * work correctly with any domain size n.
 * Compile with:
 *
 * nvcc cuda-ltl.cu -o cuda-ltl
 *
 * Run with:
 * (Test it with the Game of Life parameter )
 * ./cuda-ltl (R =) 1 (B1 = B2 =) 3 (D1 =) 3  (D2 = )4 nsteps input_file output_file
 *
 ****************************************************************************/
 #include "hpc.h"
 #include <stdio.h>
 #include <stdlib.h>
 #include <string.h>
 #include <assert.h>
 #include <ctype.h> /* for isdigit */

#define BLKSIZE 32
/*we assume the presence of 8 rows / columns of ghost cells per side,
 then using only those that serve*/
#define HALO 8
/* We use 1D blocks to copy ghost cells; in this case we can use up to
   1024 threads per block (for the GPUs on the lab machine) */

#define BLKSIZE_GHOST 1024

typedef unsigned char cell_t;

/* The struct of the bmap_t with a size and a point of unsigned char */

typedef struct {
  int n;
  cell_t *bmap;
} bmap_t;

/* The following function makes indexing of the two-dimensional CA
    grid easier. Instead of writing, e.g., grid[i][j] (which you can
    not do anyway, since the CA grids are passed around as pointers to
    linear blocks of data), you write IDX(grid, n, i, j) to get a
    pointer to grid[i][j]. This function assumes that the size of the
    CA grid is (n+2)*(n+2), where the first and last rows/columns are
    ghost cells.

    Note the use of both the __device__ and __host__ qualifiers: this
    function can be called both from host and device code. */

__device__ __host__ cell_t *IDX(cell_t *grid, int n, int i, int j)
{
   	return (grid + i*(n + 2 * HALO)+j);
}

/* Fill the ghost cells of |grid| in order to have cyclic boundary conditions*/
__global__ void copy_top_bottom(cell_t *grid, int n)
{
  const int end = HALO + n - 1;
  const int j = HALO + threadIdx.x + blockIdx.x * blockDim.x;
  const int k =  threadIdx.y + blockIdx.y * blockDim.y + 1;
  /* Copy top and bottom */

  if( k < HALO + 1){
    if ( j < end + 1) {
      *IDX(grid, n, end + k, j) = *IDX(grid, n, HALO + k - 1, j);
      *IDX(grid, n, HALO - k, j) = *IDX(grid, n, end - k + 1, j);
      }
    }
}
__global__ void copy_left_right(cell_t *grid, int n )
{

  const int end = HALO + n - 1;
  const int i = HALO + threadIdx.y + blockIdx.y * blockDim.y;
  const int k =  threadIdx.y + blockIdx.y * blockDim.y + 1;
  /* Copy left and right */

  if( k < HALO + 1){
    if ( i < end + 1 ) {
      *IDX(grid, n, i, end + k) = *IDX(grid, n, i, HALO + k - 1);
      *IDX(grid, n, i, HALO - k) = *IDX(grid, n, i, end + k - 1);
    }
  }
}
__global__ void copy_corners(cell_t *grid, int n)
{
  const int i =  threadIdx.y + blockIdx.y * blockDim.y;
  const int j = threadIdx.x + blockIdx.x * blockDim.x;

  /* Copy corners*/
  if( i < HALO){
    if (j < HALO){
      *IDX(grid, n, i ,j) = *IDX(grid, n, i + HALO + 1, j + HALO + 1 );
      *IDX(grid, n, i + HALO + n ,j + HALO + n) = *IDX(grid, n, i + HALO, j + HALO);
      *IDX(grid, n, i ,j + HALO + n) = *IDX(grid, n, i + HALO + 1, j + HALO );
      *IDX(grid, n, i + HALO + n ,j) = *IDX(grid, n, i + HALO , j + HALO + 1 );
    }
  }
}
/**
  * Write the content of the bmap_t structure pointed to by ltl to the
  * file f in PBM format and and allocates space
  * for the ghost cell that will be assigned.
  * The caller is responsible for passing a
  * pointer f to a file opened for writing
*/
void write_ltl( bmap_t* grid, FILE *f , int r )
{
  const int n = grid->n;
  fprintf(f, "P1\n");
  fprintf(f, "# produced by ltl\n");
  fprintf(f, "%d %d\n", n, n);
  for (int i = r ; i < n + r; i++) {
    for (int j = r ; j < n + r; j++) {
         fprintf(f, "%d ", *IDX(grid->bmap, n, i, j));
     }
         fprintf(f, "\n");
     }
 }

 /*Compute of the Larger than life*/
__global__ void compute_ltl( cell_t *cur, cell_t *next, int n, int r, int b1, int b2, int d1, int d2)
{
   /*we assume the presence of 8 rows / columns of ghost cells per side,
    then using only those that serve*/

   __shared__ cell_t buf[BLKSIZE+2*HALO][BLKSIZE+2*HALO];

   /* "global" indexes */
   const int gi = HALO + threadIdx.y + blockIdx.y * blockDim.y;
   const int gj = HALO + threadIdx.x + blockIdx.x * blockDim.x;
   /* "local" indexes */
   const int li = HALO + threadIdx.y;
   const int lj = HALO + threadIdx.x;
   int nbors = 0;

    /*Copy elements from global memory to local memory of block*/
    if ( gi < n + 2 * HALO && gj < n + 2 * HALO ) {
        buf[li][lj] = *IDX(cur, n, gi, gj);

        if (li < 2 * HALO) { /* left-right */
            buf[li - HALO   ][lj] = *IDX(cur, n, gi - HALO, gj);
            buf[li + BLKSIZE][lj] = (gi + BLKSIZE < n + 2 * HALO ? *IDX(cur, n, gi + BLKSIZE, gj) : 0);
        }
        if (lj < 2 * HALO) { /* top-bottom */
            buf[li][lj - HALO   ] = *IDX(cur, n, gi, gj - HALO);
            buf[li][lj + BLKSIZE] = (gj + BLKSIZE < n + 2 * HALO ? *IDX(cur, n, gi, gj + BLKSIZE) : 0);
        }
        if (li < 2 * HALO && lj < 2 * HALO) { /* corners */
          buf[li - HALO][lj - HALO   ] = *IDX(cur, n, gi - HALO, gj - HALO);
          buf[li - HALO][lj + BLKSIZE] = (gj + BLKSIZE < n + 2 * HALO ? *IDX(cur, n, gi - HALO, gj + BLKSIZE) : 0);
          buf[li + BLKSIZE][lj - HALO] = (gi + BLKSIZE < n + 2 * HALO ? *IDX(cur, n, gi + BLKSIZE, gj - HALO) : 0);
          buf[li + BLKSIZE][lj + BLKSIZE] = (gi + BLKSIZE < n + 2 * HALO && gj + BLKSIZE < n + 2 * HALO ? *IDX(cur, n, gi + BLKSIZE, gj + BLKSIZE) : 0);
        }
    }
    __syncthreads();

    const int globali = r + threadIdx.y + blockIdx.y * blockDim.y;
    const int globalj = r + threadIdx.x + blockIdx.x * blockDim.x;
    const int localy = r + threadIdx.y;
    const int localx = r + threadIdx.x;
    int i,j;
    for(i = localy - r ; i < localy + r ; i++){
      for(j = localx - r; j < localx + r; j++){
          nbors = nbors +
          buf[i][j] ;
        }
    }

    if( !buf[localx][localy] && nbors >= b1 && nbors <= b2){ // if it can relive
      *IDX(next, n, globali, globalj) = 1; //Set it as live

    }else if(buf[localx][localy] && nbors + 1 >= d1 && nbors + 1 <= d2) // if the cell remaining live
    {
      *IDX(next, n, globali, globalj) = 1;// set it as live
    }else{
      *IDX(next, n, globali, globalj) = 0; // set it as died
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
    int ng = n + (2*HALO);
    ltl->bmap = (cell_t*)malloc( ng * ng * sizeof(cell_t));
    /* scan bitmap; each pixel is represented by a single numeric
       character ('0' or '1'); spaces and other separators are ignored
       (Gimp produces PBM files with no spaces between digits) */
    for (i = HALO; i < n + HALO ; i++) {
         for (j = HALO; j < n + HALO ; j++) {
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
     bmap_t cur;
     cell_t *d_cur, *d_next,*d_tmp;
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
     dim3 cpyBlock(BLKSIZE_GHOST,HALO);
     dim3 cpyGrid((n + BLKSIZE-1)/BLKSIZE, (n + BLKSIZE-1)/BLKSIZE,1);
     dim3 stepBlock(BLKSIZE,BLKSIZE);
     dim3 stepGrid((n + BLKSIZE-1)/BLKSIZE, (n + BLKSIZE-1)/BLKSIZE);

     const size_t size = (n+2*HALO)*(n+2*HALO)*sizeof(*(cur.bmap));
    /* Allocate space for device copy of cur and next grids */
    cudaMalloc((void**)&d_cur, size);
    cudaMalloc((void**)&d_next, size);

    /* Copy initial grid to d_cur */
	  cudaMemcpy(d_cur,cur.bmap, size, cudaMemcpyHostToDevice);
    tstart = hpc_gettime();

    for (s = 0; s < nsteps; s++) {
       copy_top_bottom<<<cpyGrid,cpyBlock>>>(d_cur, n);
       copy_left_right<<<cpyGrid,cpyBlock>>>(d_cur, n);
       copy_corners<<<cpyGrid,cpyBlock>>>(d_cur, n);
       compute_ltl<<<stepGrid,stepBlock>>>(d_cur, d_next, n, R, B1, B2, D1, D2);
       d_tmp = d_cur;
       d_cur = d_next;
       d_next = d_tmp;
    }
    cudaDeviceSynchronize();
    tend = hpc_gettime();
    fprintf(stderr, "Execution time %f\n", tend - tstart);
    cudaMemcpy(cur.bmap, d_cur, size, cudaMemcpyDeviceToHost);
    write_ltl(&cur, out, R);
    fclose(out);
    free(cur.bmap);
    cudaFree(d_cur);
    cudaFree(d_next);
     return 0;
 }
