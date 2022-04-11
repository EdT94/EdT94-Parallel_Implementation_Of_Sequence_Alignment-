#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "myProto.h"



/* Reading data from the input file by process 0 */
void readDataFromFile(const char *fileName, float * weights, char * seq1, int * numOfSeq2, char *** allSeq2)
{
	FILE* fp;
	
	char buffer1[SEQUENCE_1_MAX_SIZE];
	char buffer2[SEQUENCE_2_MAX_SIZE];
	
	if ((fp = fopen(fileName, "r")) == 0)
	{
		printf("cannot open file %s for reading\n", fileName);
		return;
	}
	//weights[0] = starsWeight, weights[1] = colonsWeight, weights[2] = pointsWeight, weights[3] = spacesWeight
	fscanf(fp, "%f %f %f %f", &weights[0], &weights[1], &weights[2], &weights[3]);
	fscanf(fp, "%s",buffer1);
	
	strcpy(seq1,buffer1);
	
	fscanf(fp, "%d", numOfSeq2);
	
	*allSeq2 = (char**) malloc(sizeof(char*) * (*numOfSeq2));
	if (*allSeq2 == NULL)
	{
		printf("Failed to allocate memory for '*allSeq2' in 'readDataFromFile' function\n");
		exit(0);
	}
	
	for(int i = 0 ; i < *numOfSeq2 ; i++)
	{
		fscanf(fp, "%s", buffer2);
		(*allSeq2)[i] = (char*) malloc(strlen(buffer2));
		if ((*allSeq2)[i] == NULL)
		{
			printf("Failed to allocate memory for '(*allSeq2)[%d]' in 'readDataFromFile' function\n", i);
			exit(0);
		}
		strcpy((*allSeq2)[i], buffer2);
	}
	
	
	fclose(fp);

}



/* This function performed by each process. Each of the processes prepares its own data and sends them for calculation. */
void  computeAlignmentAndSaveToFile(char * seq1, char **allSeq2,int numOfSeq2, float * weights, int rank, MPI_Status * status)
{
	int seq1Length, seq2Length , offsetOfMaxScore[numOfSeq2], limit, nthElement[numOfSeq2], kthElement[numOfSeq2], index; 
	float score[numOfSeq2];
	double starttime, endtime;
	int * offsets;	//*(offsets+i) = score for this offset
	
	seq1Length = strlen(seq1);
	
	if(rank == 0)
		limit = numOfSeq2/2; //numOfseq2 of 0 is 4	
	else	
		limit = numOfSeq2; //numOfseq2 of 1 is 2		
	for(index = 0 ; index < limit ;index++)
	{
		score[index] = 0.0;
		nthElement[index] = kthElement[index] = offsetOfMaxScore[index] = 0;
		seq2Length = strlen(allSeq2[index]);
		starttime = MPI_Wtime();
	
		//for each offset compute the best mutant sequence 
		computeBestMutantSequence(seq1,allSeq2[index], seq1Length, seq2Length, weights, rank, &nthElement[index], &kthElement[index], &offsetOfMaxScore[index], &score[index]);
		
		
		endtime   = MPI_Wtime();
		
      		printf("\ncomputation by process %d for seq2 number %d took %.2f seconds:\noffset %d produces maximum score of %.2f with MS(%d, %d)\n", rank, index+1, endtime-starttime, offsetOfMaxScore[index], score[index], nthElement[index], kthElement[index]);	
	}
	
	collectData(limit, rank, offsetOfMaxScore, nthElement, kthElement, score, status, numOfSeq2);
	
	if(rank == 0)
	{
		for(int i = 0 ; i < numOfSeq2 ; i++)//free allSeq2 and write data to file
			free(allSeq2[i]);
		writeToFile(WRITE_FILE, offsetOfMaxScore, nthElement, kthElement, score, numOfSeq2);
	}					
}



/* Function that collects calculated data by process 0 */
void collectData(int limit, int rank, int * offsetOfMaxScore, int * nthElement, int * kthElement, float * score, MPI_Status * status,int numOfSeq2)
{
	for(int index = 0 ; index < limit ;index++)
	{
		if(rank == 0)
		{
			//recieves data from process 1 and stores it in the second half of it's array
			MPI_Recv(&offsetOfMaxScore[numOfSeq2/2+index], 1, MPI_INT, 1, 0, MPI_COMM_WORLD,status);
			MPI_Recv(&nthElement[numOfSeq2/2+index], 1, MPI_INT, 1, 0, MPI_COMM_WORLD,status);
			MPI_Recv(&kthElement[numOfSeq2/2+index], 1, MPI_INT, 1, 0, MPI_COMM_WORLD,status);
			MPI_Recv(&score[numOfSeq2/2+index], 1, MPI_FLOAT, 1, 0, MPI_COMM_WORLD,status);
			
		}
		else
		{	
			//sends data to process 0
			MPI_Send(&offsetOfMaxScore[index], 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
			MPI_Send(&nthElement[index],1,MPI_INT,0, 0, MPI_COMM_WORLD);
			MPI_Send(&kthElement[index],1,MPI_INT,0, 0, MPI_COMM_WORLD);
			MPI_Send(&score[index],1,MPI_FLOAT,0, 0, MPI_COMM_WORLD);
		}
	
	}
}



/* A function that gets the strings and their length, calculates the mutant sequence for each offset it's score with the help of cuda and omp. Each thread calculates different mutant sequence as the schedule is dynamic */
void computeBestMutantSequence(char * seq1, char * seq2, int seq1Length, int seq2Length, float * weights, int rank, int * nthElement, int * kthElement, int * offsetOfMaxScore, float * score)
{
	
	char copyOfSeq2[strlen(seq2)], stringAfterDeletedElement[strlen(seq2)];
	strcpy(copyOfSeq2,seq2);
	float * mutantsScore;
	float maxPossibleScore = (strlen(seq2) - 2) * weights[0];
	int n, k, localOffsetOfMaxScore, maxScoreFound = 0;
	for( n = 0 ;n < strlen(seq2)-1;n++)
	{
		#pragma omp parallel for private(mutantsScore, stringAfterDeletedElement, copyOfSeq2) schedule(dynamic)
		for( k = n + 1; k < strlen(seq2) ; k++)
		{
			if(maxScoreFound!=1)
			{
				strcpy(stringAfterDeletedElement, copyOfSeq2 + n + 1);
				copyOfSeq2[n] = '\0';
				strcat(copyOfSeq2,stringAfterDeletedElement);
			
				//at this point copyOfSeq2 is after nth deleted element
				strcpy(stringAfterDeletedElement, copyOfSeq2 + k);
				copyOfSeq2[k-1] = '\0';
				strcat(copyOfSeq2,stringAfterDeletedElement);
			
				//for current sequence with deleted chars find the offset that produces best score using cuda
				computeOnGPU(&mutantsScore, seq1,copyOfSeq2, seq1Length, sizeof(copyOfSeq2)/sizeof(copyOfSeq2[0])-2, weights, WEIGHTS_ARRAY_SIZE, &localOffsetOfMaxScore);
			
				#pragma omp critical
				{
					if( *score < *(mutantsScore + localOffsetOfMaxScore))
					{
						*score = *(mutantsScore + localOffsetOfMaxScore);
						*nthElement = n+1;
						*kthElement = k+1;
						*offsetOfMaxScore = localOffsetOfMaxScore ;
					}
					strcpy(copyOfSeq2,seq2);
					free(mutantsScore);
					
					/* debugging with omp */
					//printf("checking MS(%d, %d) by thread num %d, score:%.2f  maxPossibleScore:%.2f. Best so far:MS(%d, %d), offset :%d\n ", n, k, omp_get_thread_num(), *score,maxPossibleScore, *nthElement, *kthElement, *offsetOfMaxScore);
					
					/* debugging */
					//printf("checking MS(%d, %d), score:%.2f  maxPossibleScore:%.2f. Best so far:MS(%d, %d), offset :%d\n ", n, k, *score,maxPossibleScore, *nthElement, *kthElement, *offsetOfMaxScore);
					
				}
			
				#pragma omp nowait
				{
					if(*score >= maxPossibleScore)
						maxScoreFound = 1;
					
						
				}			
			}
		}
	}
}



/* Process 0 sends part of the data to process 1 */
void sendAndRecieveData(char * seq1, char ***allSeq2,int * numOfSeq2,float * weights, int rank, MPI_Status *status, int lengthOfWeightsArray)
{

	if(rank ==0)
	{
		MPI_Send(seq1, strlen(seq1) + 1 , MPI_CHAR, 1, 0, MPI_COMM_WORLD);
		int numOfSeq2ToSlave ;
		if(*numOfSeq2 % 2 == 0 )
			numOfSeq2ToSlave = *numOfSeq2 / 2;
		else 
			numOfSeq2ToSlave = *numOfSeq2 / 2 + 1;
			
		MPI_Send(&numOfSeq2ToSlave, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
		
		//send second half of the seq2 sequences
		for( int i = *numOfSeq2/2 ; i < *numOfSeq2 ; i++)
			MPI_Send((*allSeq2)[i],  strlen((*allSeq2)[i]) + 1 , MPI_CHAR, 1, 0, MPI_COMM_WORLD);
		
		//send weights
		for(int i = 0 ; i < lengthOfWeightsArray; i++)
			MPI_Send(&weights[i], 1, MPI_FLOAT, 1, 0, MPI_COMM_WORLD);
	}
	
	else
	{
	
		MPI_Recv(seq1, SEQUENCE_1_MAX_SIZE + 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD,status);
		MPI_Recv(numOfSeq2, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,status);
		char **seq2FromMaster = (char**) malloc(sizeof(char*) * (*numOfSeq2));
		if (*seq2FromMaster == NULL)
		{
			printf("Failed to allocate memory for '*seq2FromMaster' in 'sendAndRecieveData' function\n");
			exit(0);
		}
		char buffer[SEQUENCE_2_MAX_SIZE+1];
		
		//recieve second half of the seq2 sequences
		for( int i = 0 ; i < *numOfSeq2 ; i++)
		{
			MPI_Recv(&buffer, SEQUENCE_2_MAX_SIZE+1, MPI_CHAR, 0, 0, MPI_COMM_WORLD,status);
			seq2FromMaster[i] = (char*) malloc(strlen(buffer));
			strcpy(seq2FromMaster[i], buffer);	
		}
		
		*allSeq2 = seq2FromMaster;
		
		//recieve weights
		for(int i = 0 ; i < lengthOfWeightsArray ; i++)
			MPI_Recv(&weights[i], 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD,status);
			
	}
}



/* Writing the results to the output file */
void writeToFile(const char *fileName, int *offsetOfMaxScore, int *nthElementh, int *kthElementh, float * scores, int numOfSeq2)
{
	FILE* fp;
	if ((fp = fopen(fileName, "w")) == 0)
	{
		printf("cannot open file %s for writing\n", fileName);
		return;
	}
	for(int i = 0 ; i < numOfSeq2 ; i++)
	{
		fprintf(fp, "offset %d produces maximum score of %.2f for seq2 number %d with MS(%d, %d)\n",offsetOfMaxScore[i], scores[i], i+1, nthElementh[i], kthElementh[i]);
	}
	
	
	
	
	fclose(fp);
}

