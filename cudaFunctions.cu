#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <omp.h>


/* This function runs on the GPU, gets two characters and checks if they are both in the 'conservativeGroupLength' group, if so, returns 1, otherwise 0 */
__device__ int inConservative(const char seq1Char, const char seq2Char)
{
	int char1InGroup = 0 ;
	int char2InGroup = 0 ;
	const char * conservativeGroups[9] = { "NDEQ", "NEQK", "STA", "MILV", "QHRK", "NHQK", "FYW", "HY", "MILF" };

 	int conservativeGroupLength = sizeof(conservativeGroups)/sizeof(conservativeGroups[0]);
	//#pragma omp parallel for
	for(int i = 0 ; i < conservativeGroupLength ; i++)
	{
		const char * group = conservativeGroups[i];
		for(int j = 0 ; j < sizeof(group)/sizeof(char) ; j++)
		{
			if(seq1Char == group[j])
				char1InGroup = 1;
			else if(seq2Char == group[j])
				char2InGroup = 1;
					
		}
		if(char1InGroup == 1 && char2InGroup == 1)
			return 1;
		else
		{
			char1InGroup = 0 ;
			char2InGroup = 0 ;
		}
	}
	
	return 0;
	
}



/* This function runs on the GPU, gets two characters and checks if they are both in the 'semiConservativeGroups' group, if so, returns 1, otherwise 0 */
__device__ int inSemiConservative(const char seq1Char,const char seq2Char)
{
	int char1InGroup = 0 ;
	int char2InGroup = 0 ;
	const char * semiConservativeGroups[11] = { "SAG", "ATV", "CSA", "SGND", "STPA", "STNK", "NEQHRK", "NDEQHK", "SNDEQK", "HFY", "FVLIM" };

	int semiConservativeGroupLength = sizeof(semiConservativeGroups)/sizeof(semiConservativeGroups[0]);
	//#pragma omp parallel for
	for(int i = 0 ; i < semiConservativeGroupLength ; i++)
	{
		const char * group = semiConservativeGroups[i];
		for(int j = 0 ; j < sizeof(group)/sizeof(char) ; j++)
		{
			if(seq1Char == group[j])
				char1InGroup = 1;
			else if(seq2Char == group[j])
				char2InGroup = 1;
					
		}
		if(char1InGroup == 1 && char2InGroup == 1)
			return 1;
		else
		{
			char1InGroup = 0 ;
			char2InGroup = 0 ;
		}
	}
	
	return 0;
}



/* A function that gets an array and its length and returns the index with the largest value */
int findIndexOfMaxScore(float * offset, int length)
{
	int indexOfMax = 0;
	
	for( int i = 1 ; i < length ; i ++ )
	{
		if(offset[i] > offset[indexOfMax])
			indexOfMax = i ;
		
	}
	
	return indexOfMax;
}



/* A function that runs on the GPU. Each thread whose value is less than the 'maxOffset' participates in the calculation. Each thread is responsible for specific offset, comparing one letter from seq2 to seq1. Each thread runs 'seq2Length' times and adds to it's own offset the specific weight value for each comparison between letters */
__global__  void computeScoresOnGPU(float * offset, char * seq1, char * seq2, int seq2Length, float * weights, int maxOffset)
{
	int threadNum = blockDim.x * blockIdx.x + threadIdx.x;
	offset[threadNum] = 0;
	
	if(threadNum <= maxOffset)
	{
		//#pragma omp parallel for
		for(int i = 0 ; i < seq2Length ; i++)
		{
			if(seq1[threadNum+i] == seq2[i]) 
				offset[threadNum] += weights[0];
			
				
			else if(inConservative(seq1[threadNum+i], seq2[i]) == 1)	
				offset[threadNum] += weights[1];
			
				
			else if (inSemiConservative(seq1[threadNum+i], seq2[i]) == 1) 
				offset[threadNum] += weights[2];
				 
			else 
				offset[threadNum] += weights[3];
			
		}
	}
}



/* A function that alocates memory on GPU and copies the data to GPU */
void * cudaInit(void * array, size_t size)
{
	// Error code to check return values for CUDA calls
	cudaError_t err = cudaSuccess;
	
	void *dev_data;
	 
      
	// Allocate memory on GPU to copy the data from the host
	err = cudaMalloc((void **)&dev_data, size);
	if (err != cudaSuccess)
	{
		printf("Failed to allocate device memory on GPU - %s\n", cudaGetErrorString(err));
        	exit(EXIT_FAILURE);
   	}
    
   	 // Copy data from host to the GPU memory
   	 err = cudaMemcpy(dev_data, array, size, cudaMemcpyHostToDevice);
    	if (err != cudaSuccess)
    	{
       	 printf( "Failed to Copy data from host to the GPU memory - %s\n", cudaGetErrorString(err));
       	 exit(EXIT_FAILURE);
   	}

	return dev_data;
}



/* Free allocated memory on GPU */
void freeCuda(void * data)
{
	if (cudaFree(data) != cudaSuccess )
	{
		printf("Failed to free device data - %s\n", cudaGetErrorString(cudaGetLastError()));
		exit(EXIT_FAILURE);
	}
}



/* A function that prepares the data for calculation with cuda, and copies it back to the host*/
int computeOnGPU(float ** mutantsScore, char * seq1, char * seq2, int seq1Length, int seq2Length, float * weights,int weightsSize, int * offsetOfMaxScore)
{
	char * dev_seq1, * dev_seq2;
	float * dev_weights, * dev_offset, * mutantsTemp;
	int maxOffset = (seq1Length-seq2Length);
	int threadsPerBlock = 256 ;
	int blocksPerGrid = (maxOffset + threadsPerBlock - 1 ) / threadsPerBlock;
		
	mutantsTemp = (float*)malloc(sizeof(float)*maxOffset);
	if (mutantsTemp == NULL)
	{
		printf("Failed to allocate memory for 'mutantsTemp' in 'computeOnGPU' function\n");
		exit(0);
	}

	size_t seq1LengthCuda = seq1Length * sizeof(char);
	size_t seq2LengthCuda = seq2Length * sizeof(char);
	size_t sizeOfWeightsArrayCuda = weightsSize * sizeof(float);
	size_t sizeOfMutantsArrayCuda = maxOffset * sizeof(float);


	dev_seq1 = (char *)cudaInit((void *)seq1, seq1LengthCuda);
	dev_seq2 = (char *)cudaInit((void *)seq2, seq2LengthCuda);
	dev_weights = (float *)cudaInit((void*)weights, sizeOfWeightsArrayCuda);
	

	if (cudaMalloc((void **)&dev_offset, sizeOfMutantsArrayCuda) != cudaSuccess)
   	{
		printf("Failed to allocate device memory on GPU - %s\n", cudaGetErrorString(cudaGetLastError()));
		exit(EXIT_FAILURE);
    	}
    	
	computeScoresOnGPU<<<blocksPerGrid, threadsPerBlock>>>(dev_offset, dev_seq1, dev_seq2, seq2Length, dev_weights, maxOffset);

	if (cudaGetLastError() != cudaSuccess)
	{
		printf("Failed to launch kernel -  %s\n", cudaGetErrorString(cudaGetLastError()));
		exit(EXIT_FAILURE);
	}
		
	if (cudaMemcpy(mutantsTemp, dev_offset, sizeOfMutantsArrayCuda, cudaMemcpyDeviceToHost) != cudaSuccess)
	{
       	printf( "Failed to Copy data from GPU to host memory - %s\n", 				cudaGetErrorString(cudaGetLastError()));
       	exit(EXIT_FAILURE);
    	}
		 
	*offsetOfMaxScore = findIndexOfMaxScore(mutantsTemp, maxOffset);
	*mutantsScore = mutantsTemp;	
			
	freeCuda(dev_seq1);
	freeCuda(dev_seq2);
	freeCuda(dev_weights);
	freeCuda(dev_offset);
	freeCuda(dev_signs);

	return 0;
	
}










