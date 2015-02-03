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

#define n 50                         /* number of grid points */
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

	double old_norm = 0.0;
	double new_norm = 0.0;
	for(ite = 1; ite < 10000; )
	{
		//printf("\nbefore jacobi iteration\n");
		old_norm = jacobi_parallel(p, rank, n_proc, ite);
		ite++;
		new_norm = jacobi_parallel(p, rank, n_proc, ite);

		if(fabs(new_norm - old_norm) < TOLERANCE)
		{
			printf("\nConverged after %d iterations\n", ite);
			break;
		}
	}

	time2 = MPI_Wtime();

	//if(rank == 0)
	//{
	printf("No. iterations: %d\nTime spent: %gs\n", ite, time2-time1);

	//}

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

//double jacobi_parallel(double p[n][n], int rank, int n_proc, int ite)
double jacobi_parallel(double **p, int rank, int n_proc, int ite)
{
	int i, j, n_grid, s_index, count;
	double *sendBuf, *rcvBuf, *sendBuf1, *sendBuf2, *rcvBuf1, *rcvBuf2, *sendBuf3, *rcvBuf3;
	//double sendBuf[n], rcvBuf[n], sendBuf1[n], sendBuf2[n], rcvBuf1[n], rcvBuf2[n];
	MPI_Status status;

	n_grid = (n-2)/n_proc + 2;
	s_index = (n-2)/n_proc;

	count = n;

	int downtag = 1;
	int uptag = 2;

	if(rank == 0)
	{/*
		for(i = 1; i<=s_index; i++)
		{
			for(j=1; j<n-1; j++)
			{
				p[i][j] = 0.25*(p[i+1][j]+p[i-1][j]+p[i][j+1]+p[i][j-1]);
			}
		}*/

		sendBuf = (double *) malloc(count*sizeof(double));
		for(j=0; j<n; j++)
		{
			sendBuf[j] = p[s_index][j];
		}
		//MPI_Sendrecv(sendBuf, count, MPI_DOUBLE, 1, 1, rcvBuf, count, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &status);

		//MPI_Send(sendBuf, count, MPI_DOUBLE, 1, downtag, MPI_COMM_WORLD);
		//free(sendBuf);
//		MPI_Send(&p[s_index][0], count, MPI_DOUBLE, 1, downtag, MPI_COMM_WORLD);

		//MPI_Recv(rcvBuf, count, MPI_DOUBLE, 1, uptag, MPI_COMM_WORLD, &status);
//		MPI_Recv(&p[s_index+1][0], count, MPI_DOUBLE, 1, uptag, MPI_COMM_WORLD, &status);
//		MPI_Sendrecv(&p[s_index][0], count, MPI_DOUBLE, 1, downtag, &p[s_index+1][0], count, MPI_DOUBLE, 1, uptag, MPI_COMM_WORLD, &status);

	}
	else if(rank == n_proc-1)
	{/*
		for(i = rank*(s_index)+1; i<=rank*(s_index)+s_index; i++)
		{
			for(j=1; j<n-1; j++)
			{
				p[i][j] = 0.25*(p[i+1][j]+p[i-1][j]+p[i][j+1]+p[i][j-1]);
			}
		}*/

		sendBuf3 = (double *) malloc(count*sizeof(double));
		for(j=0; j<n; j++)
		{
			sendBuf3[j] = p[n-s_index-1][j];
		}
		//MPI_Sendrecv(sendBuf, count, MPI_DOUBLE, rank-1, 1, rcvBuf, count, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, &status);
		//MPI_Send(sendBuf3, count, MPI_DOUBLE, rank-1, uptag, MPI_COMM_WORLD);
		//free(sendBuf3);
//		MPI_Send(&p[n-s_index-1][0], count, MPI_DOUBLE, rank-1, uptag, MPI_COMM_WORLD);

		//MPI_Recv(rcvBuf3, count, MPI_DOUBLE, rank-1, downtag, MPI_COMM_WORLD, &status);
//		MPI_Recv(&p[n-n_grid][0], count, MPI_DOUBLE, rank-1, downtag, MPI_COMM_WORLD, &status);
//		MPI_Sendrecv(&p[n-s_index-1][0], count, MPI_DOUBLE, rank-1, uptag, &p[n-n_grid][0], count, MPI_DOUBLE, rank-1, downtag, MPI_COMM_WORLD, &status);

	}
	else
	{
		sendBuf1 = (double *) malloc(count*sizeof(double));
		sendBuf2 = (double *) malloc(count*sizeof(double));
		for(j=0; j<n; j++)
		{
			sendBuf1[j] = p[rank*s_index+1][j];
			//MPI_Sendrecv(sendBuf1, count, MPI_DOUBLE, rank-1, 1, rcvBuf1, count, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, &status);
			//MPI_Send(sendBuf1, count, MPI_DOUBLE, rank-1, uptag, MPI_COMM_WORLD);
			//free(sendBuf1);
	//		MPI_Send(&p[rank*s_index+1][0], count, MPI_DOUBLE, rank-1, uptag, MPI_COMM_WORLD);


			//MPI_Recv(rcvBuf1, count, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, &status);
	//		MPI_Recv(&p[rank*s_index][0], count, MPI_DOUBLE, rank-1, downtag, MPI_COMM_WORLD, &status);

//			MPI_Sendrecv(&p[rank*s_index+1][0], count, MPI_DOUBLE, rank-1, uptag, &p[rank*s_index][0], count, MPI_DOUBLE, rank-1, downtag, MPI_COMM_WORLD, &status);

			sendBuf2[j] = p[rank*s_index+s_index][j];
			//MPI_Sendrecv(sendBuf2, count, MPI_DOUBLE, rank+1, 1, rcvBuf2, count, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, &status);
			//MPI_Send(sendBuf2, count, MPI_DOUBLE, rank-1, downtag, MPI_COMM_WORLD);
			free(sendBuf2);
	//		MPI_Send(&p[rank*s_index+s_index][0], count, MPI_DOUBLE, rank-1, downtag, MPI_COMM_WORLD);

			//MPI_Recv(rcvBuf2, count, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, &status);
	//		MPI_Recv(&p[rank*s_index+s_index+1][0], count, MPI_DOUBLE, rank-1, uptag, MPI_COMM_WORLD, &status);
//			MPI_Sendrecv(&p[rank*s_index+s_index][0], count, MPI_DOUBLE, rank+1, downtag, &p[rank*s_index+s_index+1][0], count, MPI_DOUBLE, rank+1, uptag, MPI_COMM_WORLD, &status);

		}
	}


	for(i = rank*(s_index)+1; i<=rank*(s_index)+s_index; i++)
	{
		for(j=1; j<n-1; j++)
		{
			p[i][j] = 0.25*(p[i+1][j]+p[i-1][j]+p[i][j+1]+p[i][j-1]);
		}
	}

	return rms_norm(p);
}
