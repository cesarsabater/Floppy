#include <stdio.h>
#include <omp.h>

#define SWAP(x0,x) {float ** tmp=x0;x0=x;x=tmp;}
#define FOR_EACH_CELL for ( i=1 ; i<=N ; i++ ) { for ( j=1 ; j<=N ; j++ ) {
#define END_FOR }}

// grid
extern int **grid;
extern int slot_size;
extern int iter_from_level(int);
// display
extern double time1, time2;

// generator
extern void lin_solve_safe(int N, int b, float **x, float **x0, float a, float c);

void add_source ( int N, float **x, float **s, float dt )
{
	int i, j, size=(N+2);
	for ( i=0 ; i< size ; i++ ) for (j=0; j < size; j++) 
		x[i][j] += dt*s[i][j];
}

void set_bnd ( int N, int b, float * x )
{
	/*
	WARNING: This may be an incorrect way to simulate how a border behaves
	int i;

	for ( i=1 ; i<=N ; i++ ) {
		x[IX(0  ,i)] = b==1 ? -x[IX(1,i)] : x[IX(1,i)];
	 	x[IX(N+1,i)] = b==1 ? -x[IX(N,i)] : x[IX(N,i)]; 
		x[IX(i,0  )] = b==2 ? -x[IX(i,1)] : x[IX(i,1)]; 
	* x[IX(i,N+1)] = b==2 ? -x[IX(i,N)] : x[IX(i,N)]; 
	}
	x[IX(0  ,0  )] = 0.5f*(x[IX(1,0  )]+x[IX(0  ,1)]);
	x[IX(0  ,N+1)] = 0.5f*(x[IX(1,N+1)]+x[IX(0  ,N)]);
	x[IX(N+1,0  )] = 0.5f*(x[IX(N,0  )]+x[IX(N+1,1)]);
	x[IX(N+1,N+1)] = 0.5f*(x[IX(N,N+1)]+x[IX(N+1,N)]);
	*/
}

void lin_solve_complex( int N, int b, float **x, float **x0, float a, float c)
{
	int i, j, k;
	//printf("lin_solve version %d\n", 0);
	for ( k=0 ; k<20 ; k++ ) {
		for ( i=1 ; i<=N ; i++ ) {
			for ( j=1 ; j<=N ; j++ ) {
						x[i][j] = (x0[i][j] + a*(x[i-1][j]+x[i+1][j]+x[i][j-1]+x[i][j+1]))/c;
			}
		//set_bnd ( N, b, x );
		}
	}
}
void lin_solve_original( int N, int b, float **x, float **x0, float a, float c)
{
	int i, j, k;
	int gi, gj;
	//printf("lin_solve version %d\n", 0);
	for ( k=0 ; k<20 ; k++ ) {
		for ( i=1 ; i<=N ; i++ ) {
			gi = (i-1)/slot_size; 
			for ( j=1 ; j<=N ; j++ ) {
					gj = (j-1)/slot_size;
					if (k < iter_from_level(grid[gi][gj])) { 
						x[i][j] = (x0[i][j] + a*(x[i-1][j]+x[i+1][j]+x[i][j-1]+x[i][j+1]))/c;
					} 
			}
		}
		//set_bnd ( N, b, x );
	}
}


void diffuse ( int N, int b, float ** x, float ** x0, float diff, float dt)
{
	float a=dt*diff*N*N;
	lin_solve_safe( N, b, x, x0, a, 1+4*a);
}

void advect ( int N, int b, float **d, float **d0, float **u, float **v, float dt )
{
	int i, j, i0, j0, i1, j1;
	float x, y, s0, t0, s1, t1, dt0;

	dt0 = dt*N;
	for ( i=1 ; i<=N ; i++ ) { for ( j=1 ; j<=N ; j++ ) {
		x = i-dt0*u[i][j]; y = j-dt0*v[i][j];
		if (x<0.5f) x=0.5f; if (x>N+0.5f) x=N+0.5f; i0=(int)x; i1=i0+1;
		if (y<0.5f) y=0.5f; if (y>N+0.5f) y=N+0.5f; j0=(int)y; j1=j0+1;
		s1 = x-i0; s0 = 1-s1; t1 = y-j0; t0 = 1-t1;
		d[i][j] = s0*(t0*d0[i0][j0]+t1*d0[i0][j1])+
					 s1*(t0*d0[i1][j0]+t1*d0[i1][j1]);
	} }
	// set_bnd ( N, b, d );
}

void project ( int N, float **u, float **v, float **p, float **div)
{
	int i, j;

	for ( i=1 ; i<=N ; i++ ) { for ( j=1 ; j<=N ; j++ ) {
		div[i][j] = -0.5f*(u[i+1][j]-u[i-1][j]+v[i][j+1]-v[i][j-1])/N;
		p[i][j] = 0;
	} }	
	//set_bnd ( N, 0, div ); set_bnd ( N, 0, p );

	lin_solve_complex(N, 0, p, div, 1, 4);

	for ( i=1 ; i<=N ; i++ ) { for ( j=1 ; j<=N ; j++ ) {
		u[i][j] -= 0.5f*N*(p[i+1][j]-p[i-1][j]);
		v[i][j] -= 0.5f*N*(p[i][j+1]-p[i][j-1]);
	} }
	//set_bnd ( N, 1, u ); set_bnd ( N, 2, v );
}

void dens_step(int N, float **x, float **x0, float **u, float **v, float diff, float dt)
{
	//time1 = omp_get_wtime();
	add_source ( N, x, x0, dt );
	SWAP ( x0, x ); diffuse ( N, 0, x, x0, diff, dt);
	SWAP ( x0, x ); advect ( N, 0, x, x0, u, v, dt);
	//time2 = omp_get_wtime();
	//printf( "Number of seconds: %f\n", time2 - time1 );
}

void vel_step ( int N, float **u, float **v, float **u0, float **v0, float visc, float dt)
{
	add_source ( N, u, u0, dt ); add_source ( N, v, v0, dt );
	SWAP(u0, u); diffuse(N, 1, u, u0, visc, dt);
	SWAP(v0, v); diffuse(N, 2, v, v0, visc, dt);
	//time1 = omp_get_wtime();
	project(N, u, v, u0, v0);
	//time2 = omp_get_wtime();
	//printf( "Number of seconds: %f\n", time2 - time1 );
	SWAP ( u0, u ); SWAP ( v0, v );
	//time1 = omp_get_wtime();
	advect(N, 1, u, u0, u0, v0, dt ); advect(N, 2, v, v0, u0, v0, dt );
	//time2 = omp_get_wtime();
	project ( N, u, v, u0, v0);
	//printf( "Number of seconds: %f\n", time2 - time1 );
}


