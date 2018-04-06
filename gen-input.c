/****************************************************************************
 *
 * gen-input.c - Input generator for LtL
 *
 * Written in 2017 by Moreno Marzolla <moreno.marzolla(at)unibo.it>
 *
 * To the extent possible under law, the author(s) have dedicated all 
 * copyright and related and neighboring rights to this software to the 
 * public domain worldwide. This software is distributed without any warranty.
 *
 * You should have received a copy of the CC0 Public Domain Dedication
 * along with this software. If not, see 
 * <http://creativecommons.org/publicdomain/zero/1.0/>. 
 *
 * --------------------------------------------------------------------------
 *
 * Compile with:
 * gcc -std=c99 -Wall -Wpedantic gen-input.c -o gen-input
 *
 * Run with:
 * ./gen-input size rho > file.pbm
 *
 * Where:
 * - size is the size (number of rows/columns) of the image (integer)
 * - rho is the density of black dots in the image (real, 0 < rho < 1)
 *
 * Example:
 * ./gen-input 512 0.3 > in.pbm
 *
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>

int main( int argc, char* argv[] )
{    
    int i, j, n;
    float rho;

    if ( argc != 3 ) {
        fprintf(stderr, "Usage: %s size rho > outfile\n", argv[0]);
        return -1;
    }

    n = atoi(argv[1]);
    rho = atof(argv[2]);

    printf("P1\n");
    printf("# produced by gen-input.c %d %f\n", n, rho);
    printf("%d %d\n", n, n);
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            printf("%d ", (((float)rand())/RAND_MAX < rho) );
        }
        printf("\n");
    }
    return 0;
}
