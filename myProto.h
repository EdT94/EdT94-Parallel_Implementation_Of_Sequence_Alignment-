#pragma once

#define SEQUENCE_1_MAX_SIZE 5000
#define SEQUENCE_2_MAX_SIZE 3000
#define WEIGHTS_ARRAY_SIZE 4
#define READ_FILE "Input.txt"
#define WRITE_FILE "output.txt"




void readDataFromFile(const char *fileName, float * weights, char * seq1, int * numOfSeq2, char *** allSeq2);


void  computeAlignmentAndSaveToFile(char * seq1, char **allSeq2,int numOfSeq2, float * weights, int rank, MPI_Status * status);


void collectData(int limit, int rank, int * offsetOfMaxScore, int * nthElement, int * kthElement, float * score, MPI_Status * status, int numOfSeq2);


void computeBestMutantSequence(char * seq1, char * seq2, int seq1Length, int seq2Length, float * weights, int rank, int * nthElementh, int * kthElement,int * offsetOfMaxScore, float * score);


void sendAndRecieveData(char * seq1, char ***allSeq2,int * numOfSeq2, float * weights, int rank, MPI_Status *status, int lengthOfWeightsArray);


void writeToFile(const char *fileName, int *offsetOfMaxScore, int *nthElementh, int *kthElementh, float * scores, int numOfSeq2);


int computeOnGPU(float ** mutantsScore, char * seq1, char * seq2, int seq1Length, int seq2Length, float * weights,int weightsSize, int * offsetOfMaxScore);





