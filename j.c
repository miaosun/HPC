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

#define n 20                         /* number of grid points */
#define TOLERANCE 0.00000001
#define PI 3.14159265

static double rms_norm(double *a);
static double jacobi_parallel(double *p, int rank, int n_proc, int ite);

int main(int argc, char **argv)
{
	int rank, n_proc, n_grid, i, j, ite;
	double time1, time2;
	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &n_proc);

	//double p[n][n];
	double *p;
	p = malloc(n*n*sizeof(double));

	for(i=0; i<n; i++)                 /* clear the array  */
	{
		for (j=0; j<n; j++)
			p[i*n+j] = 0;
	}

	for(i=1; i<n-1; i++)       // boundary condition for x=0
	{
		p[i*n] = sin(PI*i/n) * sin(PI*i/n);
	}

	time1 = MPI_Wtime();

	double old_norm = 0.0;
	double new_norm = 0.0;
	for(ite = 1; ite < 10000; )
	{
		old_norm = jacobi_parallel(p, rank, n_proc, ite);
		printf("\nOld_norm: %f\n", old_norm);
		ite++;
		new_norm = jacobi_parallel(p, rank, n_proc, ite);
		printf("\nNew_norm: %f\n", old_norm);

		if(fabs(new_norm - old_norm) < TOLERANCE)
		{
			printf("\nConverged after %d iterations\n", ite);
			break;
		}
	}

	time2 = MPI_Wtime();

	if(rank == 0)
	{
		printf("No. iterations: %d\nTime spent: %gs\n", ite, time2-time1);
	}
	/*
	if(rank == 0)
	{
		n_grid = (n-2) / n_proc + 2;
		count = n-2;

		sendBuf = (double *) malloc(n*sizeof(double));
		for(i=0; i<n; i++)
		{
			sendBuf[i] = p[n-2][i];
		}
		//for(j=1; j<n_proc; j++)
		//	MPI_Send(&count, 1, MPI_INT, j, 1, MPI_COMM_WORLD);
	}
	//else
	//	MPI_Recv(&count, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);

	MPI_Bcast(&count, 1, MPI_INT, 0, MPI_COMM_WORLD);
	rcvBuf = (double *) malloc(count*sizeof(double));

	if(rank != 0 && rank != n_proc-1)
	{
		MPI_SendRecv(sendBuf, count, MPI_DOUBLE, rank+1, 1, rcvBuf, count, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, &status);

	}

	if(rank == 0)
	{

	}

	if(rank == n_proc-1)
	{

	}
	 */

	MPI_Finalize();
}

//double rms_norm (double a[n][n] )
double rms_norm(double *a)
{
	int i;
	int j;
	double v;

	v = 0.0;

	for ( j = 0; j < n; j++ )
	{
		for ( i = 0; i < n; i++ )
		{
			v = v + a[i*n+j] * a[i*n+j];
		}
	}
	v = sqrt ( v / ( double ) ( n * n )  );

	return v;
}

//double jacobi_parallel(double p[n][n], int rank, int n_proc, int ite)
double jacobi_parallel(double *p, int rank, int n_proc, int ite)
{
	int i, j, n_grid, s_index, count;
	double *sendBuf, *rcvBuf, *sendBuf1, *sendBuf2, *rcvBuf1, *rcvBuf2;
	//double sendBuf[n], rcvBuf[n], sendBuf1[n], sendBuf2[n], rcvBuf1[n], rcvBuf2[n];
	MPI_Status status;

	n_grid = (n-2)/n_proc + 2;
	s_index = (n-2)/n_proc;

	count = n-2;

	sendBuf = (double *) malloc(count*sizeof(double));
	rcvBuf = (double *) malloc(count*sizeof(double));
	sendBuf1 = (double *) malloc(count*sizeof(double));
	sendBuf2 = (double *) malloc(count*sizeof(double));
	rcvBuf1 = (double *) malloc(count*sizeof(double));
	rcvBuf2 = (double *) malloc(count*sizeof(double));


	if(rank == 0)
	{
		for(j=1; j<n-1; j++)
		{
			sendBuf[j*n-1] = p[n_grid*n-1+j];
		}
		MPI_Sendrecv(sendBuf, count, MPI_DOUBLE, 1, 1, rcvBuf, count, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &status);
	}
	else if(rank == n_proc-1)
	{
		for(j=1; j<n-1; j++)
		{
			sendBuf[j*n-1] = p[(n-n_grid)*n+1+j];
		}
		MPI_Sendrecv(sendBuf, count, MPI_DOUBLE, rank-1, 1, rcvBuf, count, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, &status);
	}
	else
	{
		//for(i = rank*(s_index); i<=rank*(s_index)+(s_index+1); i++)
		//for(int i=0; i<n_grid; i++)
		//{

		//}
		for(j=1; j<n-1; j++)
		{
			sendBuf1[j*n-1] = p[rank*s_index*n+1+i];
			MPI_Sendrecv(sendBuf1, count, MPI_DOUBLE, rank-1, 1, rcvBuf1, count, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, &status);
			sendBuf2[j*n-1] = p[rank*s_index*n+s_index+1+i];
			MPI_Sendrecv(sendBuf2, count, MPI_DOUBLE, rank+1, 1, rcvBuf2, count, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, &status);
		}
	}

	for(i = rank*(s_index); i<=rank*(s_index)+(s_index+1); i++)
	{
		for(j=1; j<n-1; j++)
		{
			p[i*n+j] = 0.25*(p[i*n+1+j]+p[i*n-1+j]+p[i*n+j+1]+p[i*n+j-1]);
		}
	}

	return rms_norm(p);
}
