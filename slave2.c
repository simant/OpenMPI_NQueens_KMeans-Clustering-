/* 
 * File:   main.c
 * Author: simant
 *
 * Created on April 25, 2012, 10:57 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>

#define true 1;
#define false 0;

/*
 * 
 */
int count = 0;

struct nQueens {
    int n;
    int col[20];
    int numFit, numUnFit, numSolution;
};

int Fitness(struct nQueens *q, int i) {
    // Check if queen in row k threatens queen in row i
    int k;
    for (k = 0; k < i; k++)
        if (q->col[i] == q->col[k] || abs(q->col[i] - q->col[k]) == i - k)
            return false; // does threaten, so not Fit
    return true; // no threats, so Fit
}


// Display each solution as it's found, and statistics

void OutputSolution(struct nQueens *q, int numprocs, char *processor_name) {
    int i;
    //printf("Process  on %s out of %d\n", processor_name, numprocs);

    printf("%d\t\t%d\t\t%d\t\t", q->numSolution, q->numFit, q->numUnFit);
    for (i = 0; i < q->n; i++) {
        printf("(%d,%d)", i + 1, q->col[i] + 1);
        //fprintf(fp,"%d\t%d\t",i+1,q->col[i]);
        count++;
    }
    printf("\n");
}

void queens(struct nQueens *nq, int i, int numprocs, char *processor_name) {

    int j;
    //nq.col[i]=7;
    //printf("%d\n", nq.numFit);
    if (Fitness(nq, i - 1)) { // Continue only if columns 0,...,i-1 are Fit.
        nq->numFit++;
        if (i == nq->n) { // Have a complete solution.
            nq->numSolution++;
            //OutputSolution(nq, numprocs, processor_name);
        } else {
            for (j = 0; j < nq->n; j++) { // place queen in
                nq->col[i] = j; // row i, column j
                queens(nq, i + 1, numprocs, processor_name); // and continue to next row
            }
        }
    } else {
        nq->numUnFit++;
    }

}

//Now display results

void finish(struct nQueens *nq, int rank, int numprocs, char *processor_name, clock_t *start) {

    double time_used;
    clock_t end = clock();
    time_used = (((double) end - (double) *start)) / ((double) CLOCKS_PER_SEC);
    printf("Process %d out of %d on machine %s gives %d solutions in %f seconds\n", rank, numprocs, processor_name, nq->numSolution, time_used);
}

int main(int argc, char** argv) {
    struct nQueens nq;
    int N, i;
    int numprocs, rank, namelen, err;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Status status;
    clock_t start;
    double cpu_time_used;
    start = clock();
    //start1 = MPI_Wtime();
    err = MPI_Init(&argc, &argv); /* Initialize MPI */
    if (err != MPI_SUCCESS) {
        printf("MPI_init failed!\n");
        exit(1);
    }

    err = MPI_Comm_size(MPI_COMM_WORLD, &numprocs); /* Get nr of tasks */
    if (err != MPI_SUCCESS) {
        printf("MPI_Comm_size failed!\n");
        exit(1);
    }

    err = MPI_Comm_rank(MPI_COMM_WORLD, &rank); /* Get id of this process */
    if (err != MPI_SUCCESS) {
        printf("MPI_Comm_rank failed!\n");
        exit(1);
    }

    MPI_Get_processor_name(processor_name, &namelen); /*Get name of the machine on which process is running*/


    N = atoi(argv[1]);

    if (N > 3) {
        nq.n = N;
        nq.numUnFit = 0;
        nq.numFit = 0;
        nq.numSolution = 0;
        nq.col[0] = rank;
        printf("Process:%d out of %d processes started on %s\n", rank, numprocs, processor_name);
        queens(&nq, 1, numprocs, processor_name);
        finish(&nq, rank, numprocs, processor_name, &start);
        MPI_Finalize();
    }
    return (EXIT_SUCCESS);

}
