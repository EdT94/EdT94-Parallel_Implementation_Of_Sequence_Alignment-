#include <mpi.h>
#include <stdio.h>
#include <omp.h>
#include <string.h>
#include <stdlib.h>
#include "myProto.h"



int main(int argc, char *argv[])
{
	float weights[4];
	int numOfSeq2, numOfProcs, rank;
	char seq1[SEQUENCE_1_MAX_SIZE], **allSeq2;
	double starttime, endtime;
	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numOfProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if(rank == 0)
		starttime = MPI_Wtime();
	if (numOfProcs != 2)
	{
		printf("Run the example with 2 processes only\n");
		MPI_Abort(MPI_COMM_WORLD, __LINE__);
	}
	
	if(rank == 0)
	{
		readDataFromFile(READ_FILE, weights, seq1, &numOfSeq2, &allSeq2);
		if(numOfSeq2 != 4)
		{
			printf("Run the example with 4 seq2 sequences only\n");
			MPI_Abort(MPI_COMM_WORLD, __LINE__);
		}
	}
			
    	sendAndRecieveData(seq1, &allSeq2, &numOfSeq2, weights, rank, &status, sizeof(weights)/sizeof(weights[0]));
    		
    	computeAlignmentAndSaveToFile(seq1, allSeq2, numOfSeq2 , weights, rank, &status);
    	
	
	if(rank == 0)
	{
		endtime = MPI_Wtime();
		printf("\ntotal calculation time: %.2f seconds\n", endtime-starttime);
	}
		
	MPI_Finalize();
	
	return 0;
}

