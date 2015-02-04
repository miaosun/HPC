/*
 * Jacobi_parallel.c
 *
 *  Created on: Feb 2, 2015
 *      Author: miaosun
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "mpi.h"

#define n 1000                         /* number of grid points */
#define TOLERANCE 0.00000001
#define PI 3.14159265

static double rms_norm(double **a);
static double jacobi_parallel(double **p, int rank, int n_proc, int ite);
static double analytical_norm(int iter);

int main(int argc, char **argv)
{
	int rank, n_proc, n_grid, i, j, ite;
	double time1, time2;
	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &n_proc);

	double **p;
	p = ( double** )malloc( n*sizeof(double *));
	for(i=0; i<n; i++)
		p[i] = (double *) malloc(n*sizeof(double));

	for(i=0; i<n; i++)                 /* clear the array  */
	{
		for (j=0; j<n; j++)
			p[i][j] = 0;
	}

	for(i=1; i<n-1; i++)       // boundary condition for x=0
	{
		p[i][0] = sin(PI*i/n) * sin(PI*i/n);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	time1 = MPI_Wtime();

	double old_norm = 0.0;
	double new_norm = 0.0;
	for(ite = 1; ite < 10000; )
	{
		old_norm = jacobi_parallel(p, rank, n_proc, ite);
		ite++;
		new_norm = jacobi_parallel(p, rank, n_proc, ite);

		if(rank == n_proc-1)
		{
			if(fabs(new_norm - old_norm) < TOLERANCE)
			{
				printf("\nConverged after %d iterations\n", ite);

				double analytical = analytical_norm(ite);
				printf("\nError: %f\n", fabs(analytical-new_norm));
				break;
			}
		}
	}

	time2 = MPI_Wtime();

	if(rank == n_proc-1)
	{
		printf("Time spent: %gs\n\n", time2-time1);
	}

	free(p);
	exit(0);
	MPI_Finalize();
}

//double rms_norm (double a[n][n] )
double rms_norm(double **a)
{
	int i, j;
	double v = 0.0;

	for ( j = 0; j < n; j++ )
	{
		for ( i = 0; i < n; i++ )
		{
			v = v + a[i][j] * a[i][j];
		}
	}
	v = sqrt ( v / ( double ) ( n * n )  );

	return v;
}

double analytical_norm(int iter)
{
	int ii, i, j, k;
	double **a;
	a = ( double** )malloc( n*sizeof(double *));
	for(ii=0; ii<n; ii++)
		a[ii] = (double *) malloc(n*sizeof(double));

	double h = 1.0 / (n-1);
	double res = 0.0;
	long double b = 0.0;
	long double c = 0.0;
	for(i=0; i<n; i++)
	{
		for(j=0; j<n; j++)
		{
			double x = i*h;
			double y = j*h;
			for(k=1; k<=iter; k++)
			{
				if(k != 2)
				{
					b = 4*(cos(k*PI)-1)*(1/tanh(k*PI))*sin(k*PI*y)*(tanh(k*PI)*cosh(k*PI*x)-sinh(k*PI*x));
					c = (PI*(pow(k,3)-4*k));

					res += b/c;
				}
			}
			a[i][j] = res;
		}
	}

	return rms_norm(a);
}

double jacobi_parallel(double **p, int rank, int n_proc, int ite)
{
	int i, j, n_grid, s_index, count;

	MPI_Status status;

	n_grid = (n-2)/n_proc + 2;
	s_index = (n-2)/n_proc;

	count = n;

	int downtag = 1;
	int uptag = 2;

	// each processor updating values inside its range [exclusive]
	for(i = rank*(s_index)+1; i<=rank*(s_index)+s_index; i++)
	{
		for(j=1; j<n-1; j++)
		{
			p[i][j] = 0.25*(p[i+1][j]+p[i-1][j]+p[i][j+1]+p[i][j-1]);
		}
	}

	if(rank == 0)
	{
		MPI_Sendrecv(&p[s_index][0], count, MPI_DOUBLE, 1, downtag, &p[s_index+1][0], count, MPI_DOUBLE, 1, uptag, MPI_COMM_WORLD, &status);
	}
	else if(rank == n_proc-1)
	{
		MPI_Sendrecv(&p[n-n_grid+1][0], count, MPI_DOUBLE, rank-1, uptag, &p[n-n_grid][0], count, MPI_DOUBLE, rank-1, downtag, MPI_COMM_WORLD, &status);
	}
	else
	{
		MPI_Sendrecv(&p[rank*s_index+1][0], count, MPI_DOUBLE, rank-1, uptag, &p[rank*s_index][0], count, MPI_DOUBLE, rank-1, downtag, MPI_COMM_WORLD, &status);

		MPI_Sendrecv(&p[rank*s_index+s_index][0], count, MPI_DOUBLE, rank+1, downtag, &p[rank*s_index+s_index+1][0], count, MPI_DOUBLE, rank+1, uptag, MPI_COMM_WORLD, &status);
	}

	return rms_norm(p);
}
