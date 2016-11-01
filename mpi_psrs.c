//
//  mpi_psrs.c
//  psrs
//
//  Created by Margaret Guo on 10/25/16.
//  Copyright © 2016 Margaret Guo. All rights reserved.
//

/**************************
 Algorithm of Parallel sorting by regular sampling:
 Start: Set p processors and generate random n keys.
 Phase 0: Divide n keys into p parts and assign them to each threads.
 Phase 1: Using quicksort to sort the list of keys(its size is n/p) locally and generate p pivots.
 Phase 2: Merge p sets of p pivots together using merge sort with processor 0, then generate p-1 pivots .
 Phase 3: Divide list of data into p parts in each processors.
 Phase 4: All 'ith' sublist are gathered and then merged by merge sort.
 Phase 5: The root processor gathers data from all processors and assembles the sorted list of n keys.
 *************************/

#include <mpi.h>
#include <stdio.h>
#include <sys/time.h>//gettimeofday()
#include <stdlib.h>
#include <string.h>
//#include "mpi_psrs.h"

static double getCurTime(void){
  struct timeval tim;
  gettimeofday(&tim, NULL);
  return tim.tv_sec+(tim.tv_usec/1000000.0);
}
///https://github.com/finlir/sample_sort/blob/master/sort.c

int cmpfunc (const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}////https://www.tutorialspoint.com/c_standard_library/c_function_qsort.htm

void multimerge(int **partitions, int *mergedarray, int *recv_leng, int length, int size){
  int min = partitions[0][0];
  int k = 0;
  int minp = 0;
  int *base;
  int *count;
  
  base = malloc(sizeof(int)*size);
  count = malloc(sizeof(int)*size);

  for (int i = 0; i < size; ++i){
    base[i] = partitions[i][0];
    count[i] = 0;
  }
  while (k < length){
    for (int i = 0; i < size; ++i){
      if (base[i] <= min && count[i] < recv_leng[i]){
        min = base[i];
        minp = i;
      }
      mergedarray[k++] = min;
      count[minp]++;
    }
    if (count[minp] < recv_leng[minp]){
      base[minp] = partitions[minp][count[minp]];
      min = base[minp];
    }
    else{
      for (int i = 0; i < size; ++i){
        if (count[i] < recv_leng[i]){
          min = base[i];
          break;
        }
      }
    }
  }
}

int checktest(int *sorted, int length){
  for (int i = 0; i < length-1; ++i){
    if (sorted[i] > sorted[i+1]){
      return 0;
    }
  }
  return 1;
}

int main(int argc, char** argv) {
  
  int seed;
  int dataSize;
  int *sub_data;
  int *dataLeng;
  int *dataStart;
  int *data;
  int *pivots1;
  int *pivots2;
  int *pivots3;
  double startTime, endTime;
  int startT1, endT1;
    
  // Initialize the MPI environment. The two arguments to MPI Init are not
  // currently used by MPI implementations, but are there in case future
  // implementations might need the arguments.
  //MPI_Init(&argc, &argv);
  MPI_Init(NULL, NULL);
  seed = 10000000;
  dataSize = 8000000;

  for (int i = 0; i < argc; ++i){
    if (strcmp(argv[i],"-size") == 0){
      dataSize = atoi(argv[++i]);//set data size
    }
    if (strcmp(argv[i],"-seed") == 0){
      seed = atoi(argv[++i]);//set seed
    }
  }

  // Get the number of processes
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // Get the rank of the process
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Get the name of the processor
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);

  if (rank == 0){
    srandom(seed);
    data = malloc(sizeof(int)*dataSize);
    dataLeng = malloc(sizeof(int)*size);
    dataStart = malloc(sizeof(int)*size);
    for (int i = 0; i < dataSize; ++i){
      data[i] = random()%dataSize;
    }
    //int a = checktest(data, dataSize);
    //printf("%d\n", a);
    startTime = getCurTime();
  }
  for (int i = 0; i < size; ++i){
    dataLeng[i] = dataSize/size;
    dataStart[i] = i*dataSize/size;
  }
  dataLeng[size-1] += dataSize%size;


  /********* PHASE 1 *********/
  sub_data = malloc(sizeof(int)*dataLeng[rank]);
  pivots1 = malloc(sizeof(int)*size);

  startT1 = getCurTime();
  MPI_Scatter(data, dataLeng[0], MPI_INT, sub_data,
            dataLeng[rank], MPI_INT, 0, MPI_COMM_WORLD);
  qsort(data, dataLeng[rank], sizeof(int), cmpfunc); 
  for (int i = 0; i < size; ++i){
    pivots1[i] = sub_data[i*dataLeng[rank]/size];
  }
  endT1 = getCurTime();
  

  /********* PHASE 2 *********/
  if (rank == 0) {
    pivots2 = malloc(sizeof(float) * size * size);
  }
  MPI_Gather(pivots1, size, MPI_INT, pivots2, size, MPI_INT, 0,
            MPI_COMM_WORLD);
  if (rank == 0){
    qsort(pivots2, size*size, sizeof(int), cmpfunc);
    pivots3 = malloc(sizeof(int)*(size-1));
    for (int i = 0; i < size-1; ++i){
      pivots3[i] = pivots2[size*(i+1) + size/2 - 1];
    }
    MPI_Bcast(pivots3, size-1, MPI_INT, 0, MPI_COMM_WORLD);
  }

  /********* PHASE 3 *********/
  int *sub_array_start;
  int *sub_array_length;
  sub_array_start = malloc(sizeof(int)*size);
  sub_array_length = malloc(sizeof(int)*size);
  int index = 0;
  for (int i = 0; i < rank-1; ++i){
    sub_array_start[i] = index;
    sub_array_length[i] = 0;
    while (index < dataSize && data[index] < pivots3[i]){
        sub_array_length[i]++;
        index++;
    }
  }
  sub_array_start[size-1] = index;
  sub_array_length[size-1] = dataLeng[rank] - index;

  /********* PHASE 4 *********/
  int *recv;
  int *recv_leng;
  int *recv_start;
  int **partitions;
  int *p;
  int *sorted;
  int slength = 0;
  int *all;
  recv = malloc(sizeof(int)*dataSize);
  recv_leng = malloc(sizeof(int)*size);
  recv_start = malloc(sizeof(int)*size);
  
  MPI_Gather(&sub_array_start[rank], 1, MPI_INT, recv_start, 1, MPI_INT, rank,
            MPI_COMM_WORLD);//gather the start of all sub arrays
  MPI_Gather(&sub_array_length[rank], 1, MPI_INT, recv_leng, 1, MPI_INT, rank,
            MPI_COMM_WORLD);//gather the length of all sub arrays

  MPI_Gather(&sub_data[recv_start[rank]], recv_leng[rank], MPI_INT, recv, 
            recv_leng[rank], MPI_INT, rank, MPI_COMM_WORLD);// gather all data

  p = recv;
  partitions = malloc(sizeof(int *)*size);
  sorted = malloc(sizeof(int)*dataSize);
  for (int i = 0; i < size; ++i){
    partitions[i] = p + recv_leng[i];    
    p += recv_leng[i];
    slength += recv_leng[i];
  }
  multimerge(partitions, sorted, recv_leng, slength, size);

  //to test if the sorted array is in order
  //int a = checktest(sorted, slength);
  //printf("%d\n", a);

  all = malloc(sizeof(int)*dataSize);  
  MPI_Gather(sorted, slength, MPI_INT, all, slength, MPI_INT, 0, MPI_COMM_WORLD);
  if (rank == 0){
    endTime = getCurTime();
    printf("%lf %lf\n", startTime, endTime);
  }
  //to test if the sorted array is in order
  //int a = checktest(sorted, slength);
  //printf("%d\n", a);

  // Finalize the MPI environment. No more MPI calls can be made after this
  MPI_Finalize();
  return 0;
}
