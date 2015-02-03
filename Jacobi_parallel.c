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

static double rms_norm(double **a);
static double jacobi_parallel(double **p, int rank, int n_proc, int ite);

int main(int argc, char **argv)
{
	int rank, n_proc, n_grid, i, j, ite;
	double time1, time2;
	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &n_proc);

	//double p[n][n];
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

	time1 = MPI_Wtime();

	if(rank==0)
	{
		for(i=0; i<n; i++)
		{
			for(j=0; j<n; j++)
				printf("%f ", p[i][j]);
			printf("\n");
		}
	}

	double old_norm = 0.0;
	double new_norm = 0.0;
	for(ite = 1; ite < 10000; )
	{
		printf("\nbefore jacobi iteration\n");
		old_norm = jacobi_parallel(p, rank, n_proc, ite);
		printf("\nafter jacobi iteration\n");
		ite++;
		new_norm = jacobi_parallel(p, rank, n_proc, ite);

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

	free(p);
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

//double jacobi_parallel(double p[n][n], int rank, int n_proc, int ite)
double jacobi_parallel(double **p, int rank, int n_proc, int ite)
{
	int i, j, n_grid, s_index, count;
	double *sendBuf, *rcvBuf, *sendBuf1, *sendBuf2, *rcvBuf1, *rcvBuf2, *sendBuf3, *rcvBuf3;
	//double sendBuf[n], rcvBuf[n], sendBuf1[n], sendBuf2[n], rcvBuf1[n], rcvBuf2[n];
	MPI_Status status;

	n_grid = (n-2)/n_proc + 2;
	s_index = (n-2)/n_proc;

	count = n-2;

/*
	rcvBuf = (double *) malloc(count*sizeof(double));


	rcvBuf1 = (double *) malloc(count*sizeof(double));
	rcvBuf2 = (double *) malloc(count*sizeof(double));
	rcvBuf3 = (double *) malloc(count*sizeof(double));
*/
	int downtag = 1;
	int uptag = 2;

	if(rank == 0)
	{
		sendBuf = (double *) malloc(count*sizeof(double));
		for(j=1; j<n-1; j++)
		{
			sendBuf[j-1] = p[s_index][j];
		}
		//MPI_Sendrecv(sendBuf, count, MPI_DOUBLE, 1, 1, rcvBuf, count, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &status);
		printf("\nrank 0, before MPI_Send\n");
		MPI_Send(sendBuf, count, MPI_DOUBLE, 1, downtag, MPI_COMM_WORLD);
		free(sendBuf);
		printf("\nrank 0, before MPI_Recv\n");
		//MPI_Recv(rcvBuf, count, MPI_DOUBLE, 1, uptag, MPI_COMM_WORLD, &status);
		MPI_Recv(&p[s_index+1][0], count, MPI_DOUBLE, 1, uptag, MPI_COMM_WORLD, &status);
		printf("\nrank 0, after MPI_Recv\n");
	}
	else if(rank == n_proc-1)
	{
		sendBuf3 = (double *) malloc(count*sizeof(double));
		for(j=1; j<n-1; j++)
		{
			sendBuf3[j-1] = p[n-s_index-1][j];
		}
		//MPI_Sendrecv(sendBuf, count, MPI_DOUBLE, rank-1, 1, rcvBuf, count, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, &status);
		MPI_Send(sendBuf3, count, MPI_DOUBLE, rank-1, uptag, MPI_COMM_WORLD);
		free(sendBuf3);
		printf("\nrank last, before MPI_Recv\n");
		//MPI_Recv(rcvBuf3, count, MPI_DOUBLE, rank-1, downtag, MPI_COMM_WORLD, &status);
		MPI_Recv(&p[n-n_grid][0], count, MPI_DOUBLE, rank-1, downtag, MPI_COMM_WORLD, &status);
		printf("\nrank last, after MPI_Recv\n");
	}
	else
	{
		//for(i = rank*(s_index); i<=rank*(s_index)+(s_index+1); i++)
		//for(int i=0; i<n_grid; i++)
		//{

		//}
		sendBuf1 = (double *) malloc(count*sizeof(double));
		sendBuf2 = (double *) malloc(count*sizeof(double));
		for(j=1; j<n-1; j++)
		{
			sendBuf1[j-1] = p[rank*s_index+1][i];
			//MPI_Sendrecv(sendBuf1, count, MPI_DOUBLE, rank-1, 1, rcvBuf1, count, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, &status);
			MPI_Send(sendBuf1, count, MPI_DOUBLE, rank-1, uptag, MPI_COMM_WORLD);
			free(sendBuf1);
			printf("\nrank %d, before MPI_Recv\n", rank);

			//MPI_Recv(rcvBuf1, count, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, &status);
			MPI_Recv(&p[rank*s_index][0], count, MPI_DOUBLE, rank-1, downtag, MPI_COMM_WORLD, &status);
			printf("\nrank %d, after MPI_Recv\n", rank);


			sendBuf2[j-1] = p[rank*s_index+s_index][i];
			//MPI_Sendrecv(sendBuf2, count, MPI_DOUBLE, rank+1, 1, rcvBuf2, count, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, &status);
			MPI_Send(sendBuf2, count, MPI_DOUBLE, rank-1, downtag, MPI_COMM_WORLD);
			free(sendBuf2);
			printf("\nrank %d, before MPI_Recv\n", rank);
			//MPI_Recv(rcvBuf2, count, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, &status);
			MPI_Recv(&p[rank*s_index+s_index+1][0], count, MPI_DOUBLE, rank-1, uptag, MPI_COMM_WORLD, &status);
			printf("\nrank %d, after MPI_Recv\n", rank);
		}
	}
/*
	for(i = rank*(s_index); i<=rank*(s_index)+(s_index+1); i++)
	{
		for(j=1; j<n-1; j++)
		{
			p[i][j] = 0.25*(p[i+1][j]+p[i-1][j]+p[i][j+1]+p[i][j-1]);
		}
	}
*/
	return rms_norm(p);
}
