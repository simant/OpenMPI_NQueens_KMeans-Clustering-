/* 
 * File:   main.c
 * Author: Simant
 *
 * Created on April 16, 2012, 10:18 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

#define true 1;
#define false 0;

extern int count;

int menu() {
    int ch;
    printf("Please select one of the options:\n");
    printf("1:Display the solutions(warning may take more time)\n2:Display only No. of solutions\n");
    scanf("%d", &ch);
    return ch;

}

/*this  program compiles the slave despatches the compiled code to all the slave nodes
 * and execute them on each node with initializations
 */

int main(int argc, char** argv) {
    int proc, i = 0, j = 0, num, num_ip = 0;
    char ipaddress[2][15], buf[2];
    char N[1];
    int ch;
    char str1[50], str2[20], str3[50], str4[6];

    FILE *fp;
    //clock_t start,end;
    double start, end;
    double totaltime;

    //reading ip addresses of slave nodes from hostfile 
    fp = fopen("hostfile", "r");
    if (!fp) {
        printf("Problem in opening the host file");
        exit(0);
    }
    rewind(fp);
    // while(!fp) num_ip++;fscanf(fp, "%s\n", ipaddress[0]);printf("%s\t %d\n",ipaddress[0],num_ip);
    while (!feof(fp)) {
        fscanf(fp, "%s\n", ipaddress[i]);
        printf("%s\n", ipaddress[i]);
        i++;

    }
    //i=i+1;
    fclose(fp);
    ch = menu();
    if (ch == 1) {
        strcpy(str1, "scp slave mpi@");
        strcpy(str2, ":/home/mpi ");
        strcpy(str3, "mpiexec --hostfile hostfile -n ");
        strcpy(str4, " slave");

        //compile the slave program
        if (system("mpicc slave.c -o slave")) {
            printf("compilation of slave1 program failed\n");
            exit(0);
        }
    }
    else if (ch == 2) {
        strcpy(str1, "scp slave2 mpi@");
        strcpy(str2, ":/home/mpi ");
        strcpy(str3, "mpiexec --hostfile hostfile -n ");
        strcpy(str4, " slave2");

        //compile the slave program
        if (system("mpicc slave2.c -o slave2")) {
            printf("compilation of slave2S program failed\n");
            exit(0);
        }
    }
    else {
        printf("Invalid choice\n");
        menu();
    }

    //dispatching object files to all slave nodes
    for (j = 1; j < i; j++) {
        strcat(str1, ipaddress[j]);
        strcat(str1, str2);
        printf("\n%s\n", str1);
        int e = system(str1);
        if (e == -1) {
            printf("Can't dispatch slave program to the slave nodes\n");
            exit(0);
        }
    }


    // printf("\n%s\n",str1);
    printf("Enter number of queens for which solution is required:");
    scanf("%s", N);
    num = atoi(N);

    //validation of user input
    if (num <= 3) {
        printf("There is no solution for number of queens less than 4\n");
        exit(1);
    } else {
        if (num > 14) { //for more no. of queens more time is taken
            printf("Warning: The execution will take some time to execute\n");
        }

        strcat(str3, N);
        strcat(str3, str4);
        strcat(str3, " ");
        strcat(str3, N);
        printf("\n%s\n", str3);
        start = MPI_Wtime();
        //start = clock();
        while (system(str3)) {
            printf("Slave programs can't be executed on slave nodes\n");
            //exit(0);	
        }
        end = clock();
        totaltime = (MPI_Wtime() - start);
        printf("Total time used in performing execution of %d processes is %f seconds\n", num, totaltime);
    }
    return (EXIT_SUCCESS);
}
