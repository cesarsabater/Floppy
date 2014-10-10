#include <stdio.h>

void lin_solve_opt( int N, int b, float **x, float **x0, float a, float c)
{
	int i, j, k;
	printf("lin_solve version %d\n", versionXXXX);
	#pragma scop
	for ( k=0 ; k<20 ; k++ ) {
		for ( i=1 ; i<=N ; i++ ) { 
			for ( j=1 ; j<=N ; j++ ) {
				  x[i][j] = (x0[i][j] + a*(x[i-1][j]+x[i+1][j]+x[i][j-1]+x[i][j+1]))/c;
			}
		}
		//set_bnd ( N, b, x );
	}
	#pragma endscop
}

