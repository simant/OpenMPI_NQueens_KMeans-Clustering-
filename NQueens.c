/* 
 * File:   main.c
 * Author: Simant
 *
 * Created on April 19, 2012, 12:11 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define true 1;
#define false 0;

/*
 * 
 */

struct nQueens {
    int n;
    int col[20];
    int numPromising, numNonPromising, numSolution;
};

int promising(struct nQueens *q, int i) {
    // Check if queen in row k threatens queen in row i
    int k;
    for (k = 0; k < i; k++)
        if (q->col[i] == q->col[k] || abs(q->col[i] - q->col[k]) == i - k)
            return false; // does threaten, so not promising
    return true; // no threats, so promising
}


// Display each solution as it's found, and statistics

void OutputSolution(struct nQueens *q) {
    int i;
    printf("%d\t%d\t%d\t", q->numSolution, q->numPromising, q->numNonPromising);
    for (i = 0; i < q->n; i++)
        printf("(%d,%d)", i + 1, q->col[i] + 1);
    printf("\n");
}

void queens(struct nQueens *nq, int i) {

    int j;
    //nq.col[i]=7;
    //printf("%d\n", nq.numPromising);
    if (promising(nq, i - 1)) { // Continue only if columns 0,...,i-1 are promising.
        nq->numPromising++;
        if (i == nq->n) { // Have a complete solution.
            nq->numSolution++;
            OutputSolution(nq);
        } else {
            for (j = 0; j < nq->n; j++) { // place queen in
                 nq->col[i] = j; // row i, column j
                queens(nq, i + 1); // and continue to next row
            }
        }
    } else {
        nq->numNonPromising++;
    }

}

/*
struct nQueens start(struct nQueens nq, int N) {
    //initializing solution set
    nq.n = N;
    nq.numNonPromising = 0;
    nq.numPromising = 0;
    nq.numSolution = 0;

}
 */

void finish(struct nQueens *nq) {
    printf("No. of solutions =%d\n", nq->numSolution);
}

int main(int argc, char** argv) {
    struct nQueens nq;
    int N;
    clock_t start, end;
    double cpu_time_used;
    printf("Enter Number of Queens for which solution is required:");
    scanf("%d", &N);
    start= clock();
    if (N > 0) {
        printf("#Solutions\t #P\t #NonP\t Coordinates\n");
        //nq.n = N;
        nq.n = N;
        nq.numNonPromising = 0;
        nq.numPromising = 0;
        nq.numSolution = 0;
       //nq.col[0] = 3;
        //nq = start(nq, N);
        queens(&nq, 0);
        finish(&nq);
    }
    end = clock();
    cpu_time_used = ((double) (clock() - start))/CLOCKS_PER_SEC;
    printf("Total CPU time taken to execute is %f\n",cpu_time_used);
    return (EXIT_SUCCESS);
}
