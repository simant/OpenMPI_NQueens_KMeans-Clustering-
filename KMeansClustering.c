/*program to find k-means cluster using parallel computing and openmpi*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>    
#include <sys/types.h>  
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>     
#include <assert.h>
#define MAX_CHAR_PER_LINE 100
#include <mpi.h>
int      _debug;


int     mpi_kmeans(float**, int, int, int, float, int*, float**, MPI_Comm);
float** mpi_read(int, char*, int*, int*, MPI_Comm);
int     mpi_write(int, char*, int, int, int, float**, int*, int, MPI_Comm);


/*---< file_read() >---------------------------------------------------------*/
float** file_read(int   isBinaryFile,  /* flag: 0 or 1 */
                  char *filename,      /* input file name */
                  int  *numObjs,       /* no. data objects (local) */
                  int  *numCoords)     /* no. coordinates */
{
    float **objects;
    int     i, j, len;
    ssize_t numBytesRead;

    if (isBinaryFile) {  /* input file is in raw binary format -------------*/
        int infile;
        if ((infile = open(filename, O_RDONLY, "0600")) == -1) {
            fprintf(stderr, "Error: no such file (%s)\n", filename);
            return NULL;
        }
        numBytesRead = read(infile, numObjs,    sizeof(int));
        assert(numBytesRead == sizeof(int));
        numBytesRead = read(infile, numCoords, sizeof(int));
        assert(numBytesRead == sizeof(int));
        if (_debug) {
            printf("File %s numObjs   = %d\n",filename,*numObjs);
            printf("File %s numCoords = %d\n",filename,*numCoords);
        }

        /* allocate space for objects[][] and read all objects */
        len = (*numObjs) * (*numCoords);
        objects    = (float**)malloc((*numObjs) * sizeof(float*));
        assert(objects != NULL);
        objects[0] = (float*) malloc(len * sizeof(float));
        assert(objects[0] != NULL);
        for (i=1; i<(*numObjs); i++)
            objects[i] = objects[i-1] + (*numCoords);

        numBytesRead = read(infile, objects[0], len*sizeof(float));
        assert(numBytesRead == len*sizeof(float));

        close(infile);
    }
    else {  /* input file is in ASCII format -------------------------------*/
        FILE *infile;
        char *line, *ret;
        int   lineLen;

        if ((infile = fopen(filename, "r")) == NULL) {
            fprintf(stderr, "Error: no such file (%s)\n", filename);
            return NULL;
        }

        /* first find the number of objects */
        lineLen = MAX_CHAR_PER_LINE;
        line = (char*) malloc(lineLen);
        assert(line != NULL);

        (*numObjs) = 0;
        while (fgets(line, lineLen, infile) != NULL) {
            /* check each line to find the max line length */
            while (strlen(line) == lineLen-1) {
                /* this line read is not complete */
                len = strlen(line);
                fseek(infile, -len, SEEK_CUR);

                /* increase lineLen */
                lineLen += MAX_CHAR_PER_LINE;
                line = (char*) realloc(line, lineLen);
                assert(line != NULL);

                ret = fgets(line, lineLen, infile);
                assert(ret != NULL);
            }

            if (strtok(line, " \t\n") != 0)
                (*numObjs)++;
        }
        rewind(infile);
        if (_debug) printf("lineLen = %d\n",lineLen);

        /* find the no. objects of each object */
        (*numCoords) = 0;
        while (fgets(line, lineLen, infile) != NULL) {
            if (strtok(line, " \t\n") != 0) {
                /* ignore the id (first coordiinate): numCoords = 1; */
                while (strtok(NULL, " ,\t\n") != NULL) (*numCoords)++;
                break; /* this makes read from 1st object */
            }
        }
        rewind(infile);
        if (_debug) {
            printf("File %s numObjs   = %d\n",filename,*numObjs);
            printf("File %s numCoords = %d\n",filename,*numCoords);
        }

        /* allocate space for objects[][] and read all objects */
        len = (*numObjs) * (*numCoords);
        objects    = (float**)malloc((*numObjs) * sizeof(float*));
        assert(objects != NULL);
        objects[0] = (float*) malloc(len * sizeof(float));
        assert(objects[0] != NULL);
        for (i=1; i<(*numObjs); i++)
            objects[i] = objects[i-1] + (*numCoords);

        i = 0;
        /* read all objects */
        while (fgets(line, lineLen, infile) != NULL) {
            if (strtok(line, " \t\n") == NULL) continue;
            for (j=0; j<(*numCoords); j++)
                objects[i][j] = atof(strtok(NULL, " ,\t\n"));
            i++;
        }

        fclose(infile);
        free(line);
    }

    return objects;
}


/*----< euclid_dist_2() >----------------------------------------------------*/
/* square of Euclid distance between two multi-dimensional points            */
__inline static
float euclid_dist_2(int    numdims,  /* no. dimensions */
                    float *coord1,   /* [numdims] */
                    float *coord2)   /* [numdims] */
{
    int i;
    float ans=0.0;

    for (i=0; i<numdims; i++)
        ans += (coord1[i]-coord2[i]) * (coord1[i]-coord2[i]);

    return(ans);
}

/*----< find_nearest_cluster() >---------------------------------------------*/
__inline static
int find_nearest_cluster(int     numClusters, /* no. clusters */
                         int     numCoords,   /* no. coordinates */
                         float  *object,      /* [numCoords] */
                         float **clusters)    /* [numClusters][numCoords] */
{
    int   index, i;
    float dist, min_dist;

    /* find the cluster id that has min distance to object */
    index    = 0;
    min_dist = euclid_dist_2(numCoords, object, clusters[0]);

    for (i=1; i<numClusters; i++) {
        dist = euclid_dist_2(numCoords, object, clusters[i]);
        /* no need square root */
        if (dist < min_dist) { /* find the min and its array index */
            min_dist = dist;
            index    = i;
        }
    }
    return(index);
}

/*----< mpi_kmeans() >-------------------------------------------------------*/
int mpi_kmeans(float    **objects,     /* in: [numObjs][numCoords] */
               int        numCoords,   /* no. coordinates */
               int        numObjs,     /* no. objects */
               int        numClusters, /* no. clusters */
               float      threshold,   /* % objects change membership */
               int       *membership,  /* out: [numObjs] */
               float    **clusters,    /* out: [numClusters][numCoords] */
               MPI_Comm   comm)        /* MPI communicator */
{
    int      i, j, rank, index, loop=0, total_numObjs;
    int     *newClusterSize; /* [numClusters]: no. objects assigned in each
                                new cluster */
    int     *clusterSize;    /* [numClusters]: temp buffer for Allreduce */
    float    delta;          /* % of objects change their clusters */
    float    delta_tmp;
    float  **newClusters;    /* [numClusters][numCoords] */
    extern int _debug;

    if (_debug) MPI_Comm_rank(comm, &rank);

    /* initialize membership[] */
    for (i=0; i<numObjs; i++) membership[i] = -1;

    /* need to initialize newClusterSize and newClusters[0] to all 0 */
    newClusterSize = (int*) calloc(numClusters, sizeof(int));
    assert(newClusterSize != NULL);
    clusterSize    = (int*) calloc(numClusters, sizeof(int));
    assert(clusterSize != NULL);

    newClusters    = (float**) malloc(numClusters *            sizeof(float*));
    assert(newClusters != NULL);
    newClusters[0] = (float*)  calloc(numClusters * numCoords, sizeof(float));
    assert(newClusters[0] != NULL);
    for (i=1; i<numClusters; i++)
        newClusters[i] = newClusters[i-1] + numCoords;

    MPI_Allreduce(&numObjs, &total_numObjs, 1, MPI_INT, MPI_SUM, comm);
    if (_debug) printf("%2d: numObjs=%d total_numObjs=%d numClusters=%d numCoords=%d\n",rank,numObjs,total_numObjs,numClusters,numCoords);

    do {
        double curT = MPI_Wtime();
        delta = 0.0;
        for (i=0; i<numObjs; i++) {
            /* find the array index of nestest cluster center */
            index = find_nearest_cluster(numClusters, numCoords, objects[i],
                                         clusters);

            /* if membership changes, increase delta by 1 */
            if (membership[i] != index) delta += 1.0;

            /* assign the membership to object i */
            membership[i] = index;

            /* update new cluster centers : sum of objects located within */
            newClusterSize[index]++;
            for (j=0; j<numCoords; j++)
                newClusters[index][j] += objects[i][j];
        }

        /* sum all data objects in newClusters */
        MPI_Allreduce(newClusters[0], clusters[0], numClusters*numCoords,
                      MPI_FLOAT, MPI_SUM, comm);
        MPI_Allreduce(newClusterSize, clusterSize, numClusters, MPI_INT,
                      MPI_SUM, comm);

        /* average the sum and replace old cluster centers with newClusters */
        for (i=0; i<numClusters; i++) {
            for (j=0; j<numCoords; j++) {
                if (clusterSize[i] > 1)
                    clusters[i][j] /= clusterSize[i];
                newClusters[i][j] = 0.0;   /* set back to 0 */
            }
            newClusterSize[i] = 0;   /* set back to 0 */
        }
            
        MPI_Allreduce(&delta, &delta_tmp, 1, MPI_FLOAT, MPI_SUM, comm);
        delta = delta_tmp / total_numObjs;

        if (_debug) {
            double maxTime;
            curT = MPI_Wtime() - curT;
            MPI_Reduce(&curT, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
            if (rank == 0) printf("%2d: loop=%d time=%f sec\n",rank,loop,curT);
        }
    } while (delta > threshold && loop++ < 500);

    if (_debug && rank == 0) printf("%2d: delta=%f threshold=%f loop=%d\n",rank,delta,threshold,loop);

    free(newClusters[0]);
    free(newClusters);
    free(newClusterSize);
    free(clusterSize);

    return 1;
}



float** mpi_read(int       isBinaryFile,  /* flag: 0 or 1 */
                 char     *filename,      /* input file name */
                 int      *numObjs,       /* no. data objects (local) */
                 int      *numCoords,     /* no. coordinates */
                 MPI_Comm  comm)
{
    float    **objects;
    int        i, j, len, divd, rem;
    int        rank, nproc;
    MPI_Status status;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nproc);

    if (isBinaryFile) {  /* using MPI-IO to read file concurrently */
        int            err;
        MPI_Offset     disp;
        MPI_Datatype   filetype;
        MPI_File       fh;

        err = MPI_File_open(comm, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
        if (err != MPI_SUCCESS) {
            char errstr[MPI_MAX_ERROR_STRING];
            int  errlen;
            MPI_Error_string(err, errstr, &errlen);
            printf("Error at opening file %s (%s)\n",filename,errstr);
            MPI_Finalize();
            exit(1);
        }

        /* read numObjs & numCoords from the 1st 2 integers */
        MPI_File_read(fh, numObjs,   1, MPI_INT, &status);
        MPI_File_read(fh, numCoords, 1, MPI_INT, &status);

        if (*numObjs <= 0 || *numCoords <= 0) {
            printf("Error: file format (%s)\n",filename);
            MPI_Finalize();
            exit(1);
        }

        divd = (*numObjs) / nproc;
        rem  = (*numObjs) % nproc;
        len  = (rank < rem) ? rank*(divd+1) : rank*divd + rem;
        disp = 2 * sizeof(int) + len * (*numCoords) * sizeof(float);

        /* adjust numObjs to be local size */
        (*numObjs) = (rank < rem) ? divd+1 : divd;

        /* allocate space for data points */
        objects    = (float**)malloc((*numObjs)              * sizeof(float*));
        assert(objects != NULL);
        objects[0] = (float*) malloc((*numObjs)*(*numCoords) * sizeof(float));
        assert(objects[0] != NULL);
        for (i=1; i<(*numObjs); i++)
            objects[i] = objects[i-1] + (*numCoords);

        /* define a file type for file view */
        MPI_Type_contiguous((*numObjs), MPI_FLOAT, &filetype);
        MPI_Type_commit(&filetype);

        MPI_File_set_view(fh, disp, MPI_FLOAT, filetype, "native",
                          MPI_INFO_NULL);
        MPI_File_read_all(fh, objects[0], (*numObjs)*(*numCoords),
                          MPI_FLOAT, &status);
        MPI_Type_free(&filetype);
        MPI_File_close(&fh);
    }
    else { /* ASCII format: let proc 0 read and distribute to others */
        if (rank == 0) {
            objects = file_read(0, filename, numObjs, numCoords);
            if (objects == NULL) *numObjs = -1;
        }

        /* broadcast global numObjs and numCoords to the rest proc */
        MPI_Bcast(numObjs,   1, MPI_INT, 0, comm);
        MPI_Bcast(numCoords, 1, MPI_INT, 0, comm);

        if (*numObjs == -1) {
            MPI_Finalize();
            exit(1);
        }

        divd = (*numObjs) / nproc;
        rem  = (*numObjs) % nproc;

        if (rank == 0) {
            int index = (rem > 0) ? divd+1 : divd;

            /* index is the numObjs partitioned locally in proc 0 */
            (*numObjs) = index;

            /* distribute objects[] to other processes */
            for (i=1; i<nproc; i++) {
                int msg_size = (i < rem) ? (divd+1) : divd;
                MPI_Send(objects[index], msg_size*(*numCoords), MPI_FLOAT,
                         i, i, comm);
                index += msg_size;
            }

            /* reduce the objects[] to local size */
            objects[0] = realloc(objects[0],
                                 (*numObjs)*(*numCoords)*sizeof(float));
            assert(objects[0] != NULL);
            objects    = realloc(objects, (*numObjs)*sizeof(float*));
            assert(objects != NULL);
        }
        else {
            /*  local numObjs */
            (*numObjs) = (rank < rem) ? divd+1 : divd;

            /* allocate space for data points */
            objects    = (float**)malloc((*numObjs)            *sizeof(float*));
            assert(objects != NULL);
            objects[0] = (float*) malloc((*numObjs)*(*numCoords)*sizeof(float));
            assert(objects[0] != NULL);
            for (i=1; i<(*numObjs); i++)
                objects[i] = objects[i-1] + (*numCoords);

            MPI_Recv(objects[0], (*numObjs)*(*numCoords), MPI_FLOAT, 0,
                     rank, comm, &status);
        }
    }

    return objects;
}


/*---< mpi_write() >---------------------------------------------------------*/
int mpi_write(int        isOutFileBinary, /* flag: 0 or 1 */
              char      *filename,     /* input file name */
              int        numClusters,  /* no. clusters */
              int        numObjs,      /* no. data objects */
              int        numCoords,    /* no. coordinates (local) */
              float    **clusters,     /* [numClusters][numCoords] centers */
              int       *membership,   /* [numObjs] */
              int        totalNumObjs, /* total no. data objects */
              MPI_Comm   comm)
{
    int        divd, rem, len, err;
    int        i, j, k, rank, nproc;
    char       outFileName[1024], fs_type[32], str[32], *delim;
    MPI_File   fh;
    MPI_Status status;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nproc);

    delim = strchr(filename, ':');
    if (delim != NULL) {
        strncpy(fs_type, filename, delim-filename);
        fs_type[delim-filename] = '\0';
        /* real file name starts from delim+1 */
        delim++;
    }
    else
        delim = filename;

    /* output: the coordinates of the cluster centres ----------------------*/
    /* only proc 0 do this, because clusters[] are the same across all proc */
    if (rank == 0) {
        printf("Writing coordinates of K=%d cluster centers to file \"%s.cluster_centres\"\n",
               numClusters, delim);
        sprintf(outFileName, "%s.cluster_centres", filename);
        err = MPI_File_open(MPI_COMM_SELF, outFileName, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
        if (err != MPI_SUCCESS) {
            char errstr[MPI_MAX_ERROR_STRING];
            int  errlen;
            MPI_Error_string(err, errstr, &errlen);
            printf("Error at opening file %s (%s)\n", outFileName,errstr);
            MPI_Finalize();
            exit(1);
        }
        if (isOutFileBinary) {
            MPI_File_write(fh, &numClusters, 1, MPI_INT, &status);
            MPI_File_write(fh, &numCoords,   1, MPI_INT, &status);
            MPI_File_write(fh, clusters[0], numClusters*numCoords, MPI_FLOAT, &status);
        }
        else {
            for (i=0; i<numClusters; i++) {
                sprintf(str, "%d ", i);
                MPI_File_write(fh, str, strlen(str), MPI_CHAR, &status);
                for (j=0; j<numCoords; j++) {
                    sprintf(str, "%f ", clusters[i][j]);
                    MPI_File_write(fh, str, strlen(str), MPI_CHAR, &status);
                }
                MPI_File_write(fh, "\n", 1, MPI_CHAR, &status);
            }
        }
        MPI_File_close(&fh);
    }

    /* output: the closest cluster centre to each of the data points --------*/
    if (rank == 0)
        printf("Writing membership of N=%d data objects to file \"%s.membership\"\n",
               totalNumObjs, delim);

    if (isOutFileBinary) {
        sprintf(outFileName, "%s.membership", filename);
        err = MPI_File_open(comm, outFileName, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
        if (err != MPI_SUCCESS) {
            char errstr[MPI_MAX_ERROR_STRING];
            int  errlen;
            MPI_Error_string(err, errstr, &errlen);
            printf("Error at opening file %s (%s)\n", outFileName,errstr);
            MPI_Finalize();
            exit(1);
        }
        /* write numObjs from the 1st integer */
        if (rank == 0)
            MPI_File_write(fh, &totalNumObjs, 1, MPI_INT, &status);

        divd = totalNumObjs / nproc;
        rem  = totalNumObjs % nproc;
        len  = (rank < rem) ? rank*(divd+1) : rank*divd + rem;

        MPI_Offset disp = (len + 1) * sizeof(int);
        MPI_File_seek(fh, disp, MPI_SEEK_SET);
        MPI_File_write(fh, membership, numObjs, MPI_INT, &status);
        MPI_File_close(&fh);
    }
    else {
        if (rank == 0) { /* gather membership[] from all processes ----------*/
            int divd = totalNumObjs / nproc;
            int rem  = totalNumObjs % nproc;

            sprintf(outFileName, "%s.membership", filename);
            err = MPI_File_open(MPI_COMM_SELF, outFileName, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
            if (err != MPI_SUCCESS) {
                char errstr[MPI_MAX_ERROR_STRING];
                int  errlen;
                MPI_Error_string(err, errstr, &errlen);
                printf("Error at opening file %s (%s)\n", outFileName,errstr);
                MPI_Finalize();
                exit(1);
            }

            /* first, print out local membership[] */
            for (j=0; j<numObjs; j++) {
                sprintf(str, "%d %d\n", j, membership[j]);
                MPI_File_write(fh, str, strlen(str), MPI_CHAR, &status);
            }

            k = numObjs;
            for (i=1; i<nproc; i++) {
                numObjs = (i < rem) ? divd+1 : divd;
                MPI_Recv(membership, numObjs, MPI_INT, i, i, comm, &status);

                for (j=0; j<numObjs; j++) {
                    sprintf(str, "%d %d\n", k++, membership[j]);
                    MPI_File_write(fh, str, strlen(str), MPI_CHAR, &status);
                }
            }
            MPI_File_close(&fh);
        }
        else {
            MPI_Send(membership, numObjs, MPI_INT, 0, rank, comm);
        }
    }
    return 1;
}
/*---< usage() >------------------------------------------------------------*/
static void usage(char *argv0, float threshold) {
    char *help =
        "Usage: %s [switches] -i filename -n num_clusters\n"
        "       -i filename    : file containing data to be clustered\n"
        "       -b             : input file is in binary format (default no)\n"
        "       -r             : output file in binary format (default no)\n"
        "       -n num_clusters: number of clusters (K must > 1)\n"
        "       -t threshold   : threshold value (default %.4f)\n"
        "       -o             : output timing results (default no)\n"
        "       -d             : enable debug mode\n";
    fprintf(stderr, help, argv0, threshold);
}

/*---< main() >-------------------------------------------------------------*/
int main(int argc, char **argv) {
           int     opt;
    extern char   *optarg;
    extern int     optind;
           int     i, j;
           int     isInFileBinary, isOutFileBinary;
           int     is_output_timing, is_print_usage;

           int     numClusters, numCoords, numObjs, totalNumObjs;
           int    *membership;    /* [numObjs] */
           char   *filename;
           float **objects;       /* [numObjs][numCoords] data objects */
           float **clusters;      /* [numClusters][numCoords] cluster center */
           float   threshold;
           double  timing, io_timing, clustering_timing;

           int        rank, nproc, mpi_namelen;
           char       mpi_name[MPI_MAX_PROCESSOR_NAME];
           MPI_Status status;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Get_processor_name(mpi_name,&mpi_namelen);

    /* some default values */
    _debug           = 0;
    threshold        = 0.001;
    numClusters      = 0;
    isInFileBinary   = 0;
    isOutFileBinary  = 0;
    is_output_timing = 0;
    is_print_usage   = 0;
    filename         = NULL;

    while ( (opt=getopt(argc,argv,"p:i:n:t:abdorh"))!= EOF) {
        switch (opt) {
            case 'i': filename=optarg;
                      break;
            case 'b': isInFileBinary = 1;
                      break;
            case 'r': isOutFileBinary = 1;
                      break;
            case 't': threshold=atof(optarg);
                      break;
            case 'n': numClusters = atoi(optarg);
                      break;
            case 'o': is_output_timing = 1;
                      break;
            case 'd': _debug = 1;
                      break;
            case 'h': is_print_usage = 1;
                      break;
            default: is_print_usage = 1;
                      break;
        }
    }

    if (filename == 0 || numClusters <= 1 || is_print_usage == 1) {
        if (rank == 0) usage(argv[0], threshold);
        MPI_Finalize();
        exit(1);
    }

    if (_debug) printf("Proc %d of %d running on %s\n", rank, nproc, mpi_name);

    MPI_Barrier(MPI_COMM_WORLD);
    io_timing = MPI_Wtime();

    /* read data points from file ------------------------------------------*/
    objects = mpi_read(isInFileBinary, filename, &numObjs, &numCoords,
                       MPI_COMM_WORLD);

    if (_debug) { /* print the first 4 objects' coordinates */
        int num = (numObjs < 4) ? numObjs : 4;
        for (i=0; i<num; i++) {
            char strline[1024], strfloat[16];
            sprintf(strline,"%d: objects[%d]= ",rank,i);
            for (j=0; j<numCoords; j++) {
                sprintf(strfloat,"%10f",objects[i][j]);
                strcat(strline, strfloat);
            }
            strcat(strline, "\n");
            printf("%s",strline);
        }
    }

    timing            = MPI_Wtime();
    io_timing         = timing - io_timing;
    clustering_timing = timing;

    /* allocate a 2D space for clusters[] (coordinates of cluster centers)
       this array should be the same across all processes                  */
    clusters    = (float**) malloc(numClusters *             sizeof(float*));
    assert(clusters != NULL);
    clusters[0] = (float*)  malloc(numClusters * numCoords * sizeof(float));
    assert(clusters[0] != NULL);
    for (i=1; i<numClusters; i++)
        clusters[i] = clusters[i-1] + numCoords;

    MPI_Allreduce(&numObjs, &totalNumObjs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    /* pick first numClusters elements in feature[] as initial cluster centers*/
    if (rank == 0) {
        for (i=0; i<numClusters; i++)
            for (j=0; j<numCoords; j++)
                clusters[i][j] = objects[i][j];
    }
    MPI_Bcast(clusters[0], numClusters*numCoords, MPI_FLOAT, 0, MPI_COMM_WORLD);

    /* membership: the cluster id for each data object */
    membership = (int*) malloc(numObjs * sizeof(int));
    assert(membership != NULL);

    /* start the core computation -------------------------------------------*/
    mpi_kmeans(objects, numCoords, numObjs, numClusters, threshold, membership,
               clusters, MPI_COMM_WORLD);

    free(objects[0]);
    free(objects);

    timing            = MPI_Wtime();
    clustering_timing = timing - clustering_timing;

    /* output: the coordinates of the cluster centres ----------------------*/
    mpi_write(isOutFileBinary, filename, numClusters, numObjs, numCoords,
              clusters, membership, totalNumObjs, MPI_COMM_WORLD);

    free(membership);
    free(clusters[0]);
    free(clusters);

    /*---- output performance numbers ---------------------------------------*/
    if (is_output_timing) {
        double max_io_timing, max_clustering_timing;

        io_timing += MPI_Wtime() - timing;

        /* get the max timing measured among all processes */
        MPI_Reduce(&io_timing, &max_io_timing, 1, MPI_DOUBLE,
                   MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&clustering_timing, &max_clustering_timing, 1, MPI_DOUBLE,
                   MPI_MAX, 0, MPI_COMM_WORLD);

        if (rank == 0) {
            printf("\nPerforming **** Simple Kmeans  (MPI) ****\n");
            printf("Num of processes = %d\n", nproc);
            printf("Input file:        %s\n", filename);
            printf("numObjs          = %d\n", totalNumObjs);
            printf("numCoords        = %d\n", numCoords);
            printf("numClusters      = %d\n", numClusters);
            printf("threshold        = %.4f\n", threshold);

            printf("I/O time           = %10.4f sec\n", max_io_timing);
            printf("Computation timing = %10.4f sec\n", max_clustering_timing);
        }
    }

    MPI_Finalize();
    return(0);
}

