
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <zlib.h>

int hNumberOfReads = 0;                 // holds number of reads that are processed
int hReadsWritten = 0;                  // number of reads written to output file
char input_filename[400]="/home/linux/cuda-workspace/PreprocessReads/test.txt"; // input file name
char output_filename[400]="/mnt/sdc/Sequences/out.txt";                         // output file name
char GZIP_PIPE[410] = "gzip -> ";       // option -gz will fallback to uncompressed if error on pipe opening
char PIGZ_PIPE[410] = "pigz -> ";       // option -pigz for multithreaded compression will fallback to gz and then to uncompressed
gzFile inputFile;
gzFile templateFile;						// holds a second pointer in input read file for templating output with read name sequence and quality
FILE *outputFile;
bool overwrite_outputfile = false;      // if set true the output file will be overwritten (cmd line option -f)
bool write_gz_output = false;           // specify .gz file compressed output (cmd line option -gz)
bool write_pigz_output = false;         // specify pigz multithreaded compression for output (-pigz option)

bool AdaptorTrimming = true;			// algorithm selectors
bool QualityTrimming = true;
bool gRNATrimming = false;
int gRNA_start;							// holds gRNA position relative to adaptor
int gRNA_length;
bool verbose = false;					// if true print details of processing to screen
bool FastQ_WARNING = false;				// set to true is a FastQ file warning has been issued

short MAX_READ_LENGTH = 203;            // maximal length of reads in input file
short MIN_OUT_READ_LEN = 26;            // minimum length of reads after trimming to write to output file

short MIN_BASE_QUAL = 20;				// parameter for minimum quality for good base
short MIN_QUAL_AFTER_BAD = 25;			// parameter for average base quality after bad base
short PHRED = 33;			        	// subtract from qualities for numerical value of solexa
bool auto_phred_detect = false;         // if true try to detect PHRED format in FastQ input file

int MAX_READS_ON_GPU = 10000000;        // maximum number of reads that fit simultaneously on GPU

#define MAX_ADAPTOR_LENGTH 100          // maximum length of adaptor sequence for trimming, can contain N for any base
#define MAX_FASTQ_BUF_LINE 500          // structure of buffer for read from FastQ file

char Adaptor[MAX_ADAPTOR_LENGTH] = "ATTAGCTAGCTAGATGATCGACTAGCT"; // adaptor on CPU side from parsing command line
short mismatches = 2;					// holds maximum number of mismatches on cpu side

int numGPUs = 0;                        // holds number of CUDA devices in the system [1 TITAN to max 2 TITAN Z]
#define MAX_GPU_NUM 4
int device[MAX_GPU_NUM];                // device numbers to use
cudaStream_t upload_stream[MAX_GPU_NUM];    // cuda streams per device for overlapped asynchronous operation
cudaStream_t compute_stream[MAX_GPU_NUM];
cudaStream_t download_stream[MAX_GPU_NUM];

char *dSequenceArray[MAX_GPU_NUM];	    // device pointer to data array, will hold hReadSequence and hQualstr
char *dQualArray[MAX_GPU_NUM];			// device pointer to data array, will hold hReadSequence and hQualstr
short *dIntArray[MAX_GPU_NUM][2];       // two device pointers to integer array, will hold hReadLen
cudaEvent_t QUAL_COPY_DONE[MAX_GPU_NUM];   // events for device host synchronization
cudaEvent_t SEQ_COPY_DONE[MAX_GPU_NUM];
cudaEvent_t TrimmQ_DONE[MAX_GPU_NUM];
cudaEvent_t TrimmA_DONE[MAX_GPU_NUM];
cudaEvent_t DTH_COPY_DONE[MAX_GPU_NUM];

short *hReadLen[MAX_GPU_NUM][2];        // two buffer to hold read length on host
char *hReadSequence[MAX_GPU_NUM][2];       // buffer holding sequences of reads on host
char *hReadQualstr[MAX_GPU_NUM];     // two buffer holding qualities of reads on host

struct read
{
    char name[MAX_FASTQ_BUF_LINE];
    char sequence[MAX_FASTQ_BUF_LINE];
    char descr[MAX_FASTQ_BUF_LINE];
    char qualstr[MAX_FASTQ_BUF_LINE];
} Read;                                 // buffer one read from FastQ file i 4 lines



// macro checks return value of the CUDA runtime call and exits the application if the call failed.

#define CUDA_CHECK_RETURN(value) {											\
	cudaError_t _m_cudaStat = value;										\
	if (_m_cudaStat != cudaSuccess) {										\
		fprintf(stderr, "Error %s at line %d in file %s\n",					\
				cudaGetErrorString(_m_cudaStat), __LINE__, __FILE__);		\
		exit(1);															\
	} }


// For trimming adaptor sequences store adaptor and length as well as quality cutoffs in device constant memory

__constant__ __device__ char ADAPTOR[MAX_ADAPTOR_LENGTH];     // adaptor sequence can contain N for any match
__constant__ __device__ short ADAPTOR_LEN;                    // length of adaptor sequence in bases
__constant__ __device__ short MAX_MM=2;                       // maximum mismatches allowed not including Ns
__constant__ __device__ short QUAL_PADDING=33;				  // hold PHRED base for qualities
__constant__ __device__ short MIN_QUAL=20;					  // hold minimal quality that is considered good base
__constant__ __device__ short MIN_AVERAGE_QUAL_AFTER_BAD=25;  // minimum quality average of remaining read after bad base

void SetAdaptorSequence(const char *AdaptorSequence, const short mismatches, const short PHRED, const short MIN_BASE_QUAL, const short MIN_QUAL_AFTER_BAD)
{																// define the adaptor, mismatches and quality parameters before trimming
	short AdLen = strnlen(AdaptorSequence, MAX_ADAPTOR_LENGTH); // get length of adaptor sequence

	for(int dev=0; dev < numGPUs; dev++)
	{
		cudaSetDevice(dev);
		if(verbose) printf("[%d] Setting Adaptor sequence %s, length %d on GPU \n", dev, AdaptorSequence, AdLen);
		CUDA_CHECK_RETURN(cudaMemcpyToSymbol(ADAPTOR, AdaptorSequence, AdLen));           // copy to device
		CUDA_CHECK_RETURN(cudaMemcpyToSymbol(ADAPTOR_LEN, &AdLen, sizeof(short)));        // constant memory
		if(verbose) printf("[%d] Setting mismatches %d, PHRED+%d on GPU \n", dev, mismatches, PHRED);
		CUDA_CHECK_RETURN(cudaMemcpyToSymbol(MAX_MM, &mismatches, sizeof(short)));
		CUDA_CHECK_RETURN(cudaMemcpyToSymbol(QUAL_PADDING, &PHRED, sizeof(short)));
		if(verbose) printf("[%d] Setting minimum base quality %d, average read quality after bad %d on GPU \n", dev, MIN_BASE_QUAL, MIN_QUAL_AFTER_BAD);
		CUDA_CHECK_RETURN(cudaMemcpyToSymbol(MIN_QUAL, &MIN_BASE_QUAL, sizeof(short)));
		CUDA_CHECK_RETURN(cudaMemcpyToSymbol(MIN_AVERAGE_QUAL_AFTER_BAD, &MIN_QUAL_AFTER_BAD, sizeof(short)));
	}
}


__global__ void Find_gRNA(char *dDataArray, short *read_length, int reads, short MAX_READ_LENGTH)
{
  for(long idx = blockIdx.x * blockDim.x + threadIdx.x; idx < reads; idx += blockDim.x * gridDim.x)
  {
      if (idx < reads)
      {
          char *sequence = &dDataArray[idx*MAX_READ_LENGTH];
          short rl=read_length[idx];
          short n;
          short p=0;
          short mm=0;
          read_length[idx]=0;
          while (p <= rl-ADAPTOR_LEN)
          {
              mm=0;
              for (n=0; n<ADAPTOR_LEN; n++)
              {
                  if (ADAPTOR[n]!='N' && sequence[p+n]!=ADAPTOR[n])
                  {
                     mm++;
                     if (mm > MAX_MM) break;
                  }
              }
              if (mm <= MAX_MM)
              {
                  read_length[idx]=p;
                  break;
              }
              p++;
          }
          if (mm <= MAX_MM ) continue;
          p=rl-ADAPTOR_LEN+1;
          if (rl >= ADAPTOR_LEN-1)
          {
              mm=0;
              for (n=0; n<ADAPTOR_LEN-1; n++)
              {
                  if (ADAPTOR[n]!='N' && sequence[p+n]!=ADAPTOR[n])
                  {
                     mm++;
                     if (mm >= MAX_MM) break;
                  }
              }
              if (mm < MAX_MM)
              {
                  read_length[idx]=p;
                  continue;
              }
          }
          p++;
          if (rl >= ADAPTOR_LEN-2)
          {
              mm=0;
              for (n=0; n<ADAPTOR_LEN-2; n++)
              {
                  if (ADAPTOR[n]!='N' && sequence[p+n]!=ADAPTOR[n])
                  {
                     mm++;
                     break;
                  }
              }
              if (MAX_MM > (mm+1))
              {
                  read_length[idx]=p;
              }
          }
      }
  }
}

__global__ void TrimmAdaptors(char *dDataArray, short *read_length, int reads, short MAX_READ_LENGTH)
{
  for(long idx = blockIdx.x * blockDim.x + threadIdx.x; idx < reads; idx += blockDim.x * gridDim.x)
  {
      if (idx < reads)
      {
          char *sequence = &dDataArray[idx*MAX_READ_LENGTH];
          short rl=read_length[idx];
          short n;
          short p=0;
          short mm=0;
          while (p <= rl-ADAPTOR_LEN)
          {
              mm=0;
              for (n=0; n<ADAPTOR_LEN; n++)
              {
                  if (ADAPTOR[n]!='N' && sequence[p+n]!=ADAPTOR[n])
                  {
                     mm++;
                     if (mm > MAX_MM) break;
                  }
              }
              if (mm <= MAX_MM)
              {
                  read_length[idx]=p;
                  break;
              }
              p++;
          }
          if (mm <= MAX_MM ) continue;
          p=rl-ADAPTOR_LEN+1;
          if (rl >= ADAPTOR_LEN-1)
          {
              mm=0;
              for (n=0; n<ADAPTOR_LEN-1; n++)
              {
                  if (ADAPTOR[n]!='N' && sequence[p+n]!=ADAPTOR[n])
                  {
                     mm++;
                     if (mm >= MAX_MM) break;
                  }
              }
              if (mm < MAX_MM)
              {
                  read_length[idx]=p;
                  continue;
              }
          }
          p++;
          if (rl >= ADAPTOR_LEN-2)
          {
              mm=0;
              for (n=0; n<ADAPTOR_LEN-2; n++)
              {
                  if (ADAPTOR[n]!='N' && sequence[p+n]!=ADAPTOR[n])
                  {
                     mm++;
                     break;
                  }
              }
              if (MAX_MM > (mm+1))
              {
                  read_length[idx]=p;
              }
          }
      }
  }
}

__global__ void TrimmQuals(char *quals, short *read_length, int reads, short MAX_READ_LENGTH)
{
  for(long idx = blockIdx.x * blockDim.x + threadIdx.x; idx < reads; idx += blockDim.x * gridDim.x)
  {
      if (idx < reads)
      {
          char *q = &quals[idx*MAX_READ_LENGTH];
          short rl=read_length[idx];
          short i;
          for(i=0; i<rl; i++)
          {
              if ((q[i]-QUAL_PADDING)<=MIN_QUAL)
              {
                 short s=0;
                 short r=i+1;
                 while (r<rl)
                 {
                    s=s+q[r]-QUAL_PADDING;
                    r++;
                 }
                 if(s<(MIN_AVERAGE_QUAL_AFTER_BAD*(rl-i+1)))
                 {
                    read_length[idx]=i;
                    break;
                 }
              }
          }

      }
  }
}


void InitGPUs(void)
{
	if(verbose) printf("\n");
	for (int dev=0; dev<numGPUs; dev++)                         // initialize devices and allocate cuda resources
    {
    	if(verbose) printf("[%d] Initializing GPU\n", dev);
    	CUDA_CHECK_RETURN(cudaSetDevice(device[dev]));
    	CUDA_CHECK_RETURN(cudaMalloc((void **)&dSequenceArray[dev], MAX_READS_ON_GPU * MAX_READ_LENGTH * sizeof(char)));    // one data array
    	CUDA_CHECK_RETURN(cudaMalloc((void **)&dQualArray[dev], MAX_READS_ON_GPU * MAX_READ_LENGTH * sizeof(char)));    // one data array
    	CUDA_CHECK_RETURN(cudaMalloc((void **)&dIntArray[dev][0], MAX_READS_ON_GPU * sizeof(short)));  // two dIntArrays
    	CUDA_CHECK_RETURN(cudaMalloc((void **)&dIntArray[dev][1], MAX_READS_ON_GPU * sizeof(short)));

    	CUDA_CHECK_RETURN(cudaStreamCreate(&upload_stream[dev]));                                   // create streams on devices
    	CUDA_CHECK_RETURN(cudaStreamCreate(&compute_stream[dev]));
    	CUDA_CHECK_RETURN(cudaStreamCreate(&download_stream[dev]));

    	CUDA_CHECK_RETURN(cudaEventCreateWithFlags(&QUAL_COPY_DONE[dev], cudaEventDisableTiming || cudaEventBlockingSync));  // events for device host synchronization
    	CUDA_CHECK_RETURN(cudaEventCreateWithFlags(&SEQ_COPY_DONE[dev], cudaEventDisableTiming));
    	CUDA_CHECK_RETURN(cudaEventCreateWithFlags(&TrimmQ_DONE[dev], cudaEventDisableTiming));
    	CUDA_CHECK_RETURN(cudaEventCreateWithFlags(&TrimmA_DONE[dev], cudaEventDisableTiming));
    	CUDA_CHECK_RETURN(cudaEventCreateWithFlags(&DTH_COPY_DONE[dev], cudaEventDisableTiming || cudaEventBlockingSync));

        if(verbose) printf("[%d] Allocating pinned CPU memory for reads\n", dev);
        CUDA_CHECK_RETURN(cudaMallocHost((void **) &hReadLen[dev][0], MAX_READS_ON_GPU * sizeof(short)));   // allocate two buffers for reads length
        CUDA_CHECK_RETURN(cudaMallocHost((void **) &hReadLen[dev][1], MAX_READS_ON_GPU * sizeof(short)));

        CUDA_CHECK_RETURN(cudaMallocHost((void **) &hReadSequence[dev][0], MAX_READS_ON_GPU * MAX_READ_LENGTH * sizeof(char)));   // two buffers for reads sequence
        CUDA_CHECK_RETURN(cudaMallocHost((void **) &hReadSequence[dev][1], MAX_READS_ON_GPU * MAX_READ_LENGTH * sizeof(char)));

        CUDA_CHECK_RETURN(cudaMallocHost((void **) &hReadQualstr[dev], MAX_READS_ON_GPU * MAX_READ_LENGTH * sizeof(char))); // one buffer for qualities

    }
}

void FreeGPUs(void)
{
    for (int dev=0; dev<numGPUs; dev++)                         // release cuda resources from all devices
    {
    	CUDA_CHECK_RETURN(cudaSetDevice(dev));

    	CUDA_CHECK_RETURN(cudaEventDestroy(QUAL_COPY_DONE[dev]));
    	CUDA_CHECK_RETURN(cudaEventDestroy(SEQ_COPY_DONE[dev]));
    	CUDA_CHECK_RETURN(cudaEventDestroy(TrimmQ_DONE[dev]));
    	CUDA_CHECK_RETURN(cudaEventDestroy(TrimmA_DONE[dev]));
    	CUDA_CHECK_RETURN(cudaEventDestroy(DTH_COPY_DONE[dev]));
    	CUDA_CHECK_RETURN(cudaStreamDestroy(upload_stream[dev]));
    	CUDA_CHECK_RETURN(cudaStreamDestroy(compute_stream[dev]));
    	CUDA_CHECK_RETURN(cudaStreamDestroy(download_stream[dev]));

    	CUDA_CHECK_RETURN(cudaFree(dIntArray[dev][0]));
    	CUDA_CHECK_RETURN(cudaFree(dIntArray[dev][1]));
    	CUDA_CHECK_RETURN(cudaFree(dSequenceArray[dev]));
    	CUDA_CHECK_RETURN(cudaFree(dQualArray[dev]));

        cudaFreeHost(hReadLen[dev][0]);                        // free all page locked host memory
        cudaFreeHost(hReadLen[dev][1]);
        cudaFreeHost(hReadSequence[dev][0]);
        cudaFreeHost(hReadSequence[dev][1]);
        cudaFreeHost(hReadQualstr[dev]);
    }

}

int FindGPUs(void)
{
    int deviceCount;
    numGPUs = 0;
    CUDA_CHECK_RETURN(cudaGetDeviceCount(&deviceCount)); // search cuda devices and identify compute 3.5 with 6 GB
    if (deviceCount == 0) return(0);

    for (int dev=0; dev<deviceCount; dev++)
    {
    	cudaSetDevice(dev);
    	cudaDeviceProp deviceProp;
    	cudaGetDeviceProperties(&deviceProp, dev);

    	int driverVersion, runtimeVersion;
    	cudaDriverGetVersion(&driverVersion);
    	cudaRuntimeGetVersion(&runtimeVersion);

        if(((deviceProp.major == 3 && deviceProp.minor >= 5) || (deviceProp.major > 3)) && deviceProp.totalGlobalMem > 6144000000)
        {
            if(verbose)
            {
            	printf("\nDevice %d: %s\n", dev, deviceProp.name);
        	    printf("  CUDA Driver / Runtime Version:     %d.%d / %d.%d\n", driverVersion/1000, (driverVersion%100)/10, runtimeVersion/1000, (runtimeVersion%100)/10);
        	    printf("  CUDA Compute Capability:           %d.%d\n", deviceProp.major, deviceProp.minor);
                    printf("  PCI bus ID:                        %d:%d\n", deviceProp.pciBusID, deviceProp.pciDeviceID);
        	    printf("  Total amount of global memory:     %.0f MBytes (%llu bytes)\n", (float)deviceProp.totalGlobalMem/1048576.0f, (unsigned long long) deviceProp.totalGlobalMem);
        	    printf("  Number of Multiprocessors:         %d\n", deviceProp.multiProcessorCount);
            }
        	device[numGPUs] = dev;
        	numGPUs++;
        	if (numGPUs > MAX_GPU_NUM) break;
        }
    }
    return(numGPUs);
}

void OpenFiles()
{
    // handle compression using zlib gzFile functions for input and pipe gzip for output

    if(!overwrite_outputfile)       // test if output file exists if option -f is not specified
    {
        outputFile = fopen(output_filename, "r");
    	if(outputFile!=NULL)
    	{
    		printf("ERROR output file exists already. To overwrite run with option -f\n");
    		exit(255);
    	}  // do not attempt close file as we would have already exited if opening was successful !
    }
	inputFile = gzopen(input_filename, "r");
    if(inputFile == NULL)
    {
    	printf("ERROR opening input file: %s\n", input_filename);
    	exit(255);
    }
	templateFile = gzopen(input_filename, "r");  // holds a second access to input read file for rereading positions for output
    if(templateFile == NULL)
    {
    	printf("ERROR opening second handle to input file.\n");
    	exit(255);
    }
    if(write_pigz_output)
    {
    	outputFile = popen(strcat(PIGZ_PIPE, output_filename), "w");
        if(outputFile == NULL)
        {
        	printf("Option -pigz: could not open pipe ... falling back on libz.\n\tMake sure pigz is installed with: sudo apt-get install pigz\n");
        	write_gz_output = true;
        	write_pigz_output = false;
        }
    }
	if(write_gz_output && !write_pigz_output)
	{
		outputFile = popen(strcat(GZIP_PIPE, output_filename), "w");
        if(outputFile == NULL)
        {
        	printf("Option -gz: could not open pipe ... falling back to uncompressed output.\n");
        	write_gz_output = false;
        }
	}
	if(!write_gz_output && !write_pigz_output) outputFile = fopen(output_filename, "w");
    if(outputFile == NULL)
    {
    	printf("ERROR opening output file: %s\n", output_filename);
    	exit(255);
    }
}

int ReadBlockOfReads(int max_num_to_read, short *hReadLen, char* hReadSequence, char *hReadQualstr)
{
	int n = 0;                                                                    // returns the number of reads read
	short length;

    while(gzgets(inputFile, Read.name, MAX_FASTQ_BUF_LINE) != NULL)
    {
    	if(Read.name[0] == 13) break;				// end if CR <13> or LF <10> are encountered
    	if(Read.name[0] == 10) break;
    	if(gzgets(inputFile, Read.sequence, MAX_FASTQ_BUF_LINE)==NULL)
        {
            printf("ERROR parsing input!\n");
            exit(255);
        }
        if(gzgets(inputFile, Read.descr, MAX_FASTQ_BUF_LINE)==NULL)
        {
            printf("ERROR parsing input!\n");
            exit(255);
        }
        if(gzgets(inputFile, Read.qualstr, MAX_FASTQ_BUF_LINE)==NULL)
        {
            printf("ERROR parsing input!\n");
            exit(255);
        }

        if(!FastQ_WARNING)                         // write warning if FastQ file contains non-standard lines
        {
            if(Read.name[0]!='@' || Read.descr[0]!='+')
            {
            	printf("Warning: Encountered unusual line identifiers in FastQ file format.\n");
            	FastQ_WARNING = true;
            }
        }
        length = strlen(Read.sequence)-1;
        if(length > MAX_READ_LENGTH) length = MAX_READ_LENGTH;

//        printf("%d, ", length);

        memcpy(&hReadSequence[n*MAX_READ_LENGTH], Read.sequence, length);
        memcpy(&hReadQualstr[n*MAX_READ_LENGTH], Read.qualstr, length);
        hReadLen[n] = length;

        n++;
        if(n == max_num_to_read) break;

    }
    hNumberOfReads = hNumberOfReads + n;                           // update global variable to total number of reads
    return(n);
}

int WriteBlockOfReads(int num_of_reads_to_write, const short *hReadLen, short min_out_read_length)
{
    int n=0;
    int w=0;
    short length;
    char str_buffer[MAX_FASTQ_BUF_LINE];

    if(num_of_reads_to_write == 0) return(0);
    if(verbose) printf("I/O > num_of_reads_to_write = %d\n", num_of_reads_to_write);
    while(gzgets(templateFile, Read.name, MAX_FASTQ_BUF_LINE) != NULL)
    {
        if(gzgets(templateFile, Read.sequence, MAX_FASTQ_BUF_LINE)==NULL)
        {
            printf("ERROR parsing input file for writing trimmed output!\n");
            exit(255);
        }
        if(gzgets(templateFile, Read.descr, MAX_FASTQ_BUF_LINE)==NULL)
        {
            printf("ERROR parsing input file for writing trimmed output!\n");
            exit(255);
        }
        if(gzgets(templateFile, Read.qualstr, MAX_FASTQ_BUF_LINE)==NULL)
        {
            printf("ERROR parsing input file for writing trimmed output!\n");
            exit(255);
        }

        length = hReadLen[n];

        if(gRNATrimming)                         // isolate the gRNA from adaptor match
        {
        	if(length>0)
        	{
        		strncpy(str_buffer, &Read.sequence[length+gRNA_start], gRNA_length);
        		strncpy(&Read.sequence[0], str_buffer, gRNA_length);
        		Read.sequence[gRNA_length]=10;
        		Read.sequence[gRNA_length+1]=0;
        		strncpy(str_buffer, &Read.qualstr[length+gRNA_start], gRNA_length);
        		strncpy(&Read.qualstr[0], str_buffer, gRNA_length);
        		Read.qualstr[gRNA_length]=10;
        		Read.qualstr[gRNA_length+1]=0;

                fprintf(outputFile,"%s", Read.name);
                fprintf(outputFile,"%s", Read.sequence);
                fprintf(outputFile,"%s", Read.descr);
                fprintf(outputFile,"%s", Read.qualstr);

                w++;
        	}
        }
        else                                     // normal quality and adaptor trimming
        {
            if(length >= min_out_read_length)
            {
                Read.sequence[length]=10;        // write only reads that are longer than minimum to output file
                Read.sequence[length+1]=0;
                Read.qualstr[length]=10;
                Read.qualstr[length+1]=0;

                fprintf(outputFile,"%s", Read.name);
                fprintf(outputFile,"%s", Read.sequence);
                fprintf(outputFile,"%s", Read.descr);
                fprintf(outputFile,"%s", Read.qualstr);

                w++;
             }
        }

        n++;
        if(n == num_of_reads_to_write) break;
    }
    hReadsWritten = hReadsWritten + w;
    return(n);
}

void CloseFiles(void)
{
    gzclose(inputFile);
    gzclose(templateFile);
    if(write_gz_output || write_pigz_output) pclose(outputFile); else fclose(outputFile);
}

bool Detect_FastQ_Format_PHRED(void)
{
	int n = 0;                                                                    // returns the number of reads read
	short length;
	uint min_qual = 65;
	uint max_qual = 65;
	bool deciphered = false;

	inputFile = gzopen(input_filename, "r");
    if(inputFile == NULL)
    {
    	printf("ERROR opening input file %s\n", input_filename);
    	exit(255);
    }

    while(gzgets(inputFile, Read.name, MAX_FASTQ_BUF_LINE) != NULL)
    {
    	if(Read.name[0] == 13) break;				// end if CR <13> or LF <10> are encountered
    	if(Read.name[0] == 10) break;
    	if(gzgets(inputFile, Read.sequence, MAX_FASTQ_BUF_LINE)==NULL)
        {
            printf("ERROR parsing input!\n");
            exit(255);
        }
        if(gzgets(inputFile, Read.descr, MAX_FASTQ_BUF_LINE)==NULL)
        {
            printf("ERROR parsing input!\n");
            exit(255);
        }
        if(gzgets(inputFile, Read.qualstr, MAX_FASTQ_BUF_LINE)==NULL)
        {
            printf("ERROR parsing input!\n");
            exit(255);
        }

        if(!FastQ_WARNING)                         // write warning if FastQ file contains non-standard lines
        {
            if(Read.name[0]!='@' || Read.descr[0]!='+')
            {
            	printf("Warning: Encountered unusual line identifiers in FastQ file format.");
            	FastQ_WARNING = true;
            }
        }
        length = strlen(Read.sequence)-1;

        for(n=0; n<length;n++)
        {
        	if(min_qual > (uint)Read.qualstr[n])
        	{
        		min_qual = (uint) Read.qualstr[n];
        	}
        	else
        	{
            	if(max_qual < (uint)Read.qualstr[n]) max_qual = (uint) Read.qualstr[n];

        	}
        }
        if(min_qual < 59)
        {
        	PHRED = 33;                                    // FastQ Sanger and Illumina 1.8+ format
        	deciphered = true;
        	break;
        }
        if(max_qual > 80)
        {
        	PHRED = 64;                                    // FastQ Illumina 1.3+ format
        	deciphered = true;
        	break;
        }
    }

    gzclose(inputFile);
    return(deciphered);
}

int PreprocessReadsOnGPU(void)
{
	bool DONE = false;
	bool ALL_DONE = false;

	int ReadsProcessed = 0;
    int ReadsToProcess[MAX_GPU_NUM];
    int ReadsSentToGPU[MAX_GPU_NUM];
    int dev_side = 0;                             // left or right side of device memories [0] or [1]
    int dev_next_side = 1;
    int side_flip;

    if(verbose) printf("\nPreprocessing reads on GPU\n");


    for(int dev=0; dev<numGPUs; dev++)                         // copy data and run kernels initial round to fill CUDA queue and not waste time !
    {
    	ReadsToProcess[dev] = 0;
    	ReadsSentToGPU[dev] = 0;
    	if(!DONE)
        {
        	ReadsToProcess[dev] = ReadBlockOfReads(MAX_READS_ON_GPU, hReadLen[dev][dev_side], hReadSequence[dev][dev_side], hReadQualstr[dev]);

            if(ReadsToProcess[dev] < MAX_READS_ON_GPU)
            {
            	DONE = true;
            }
            if(ReadsToProcess[dev] == 0) continue;

    	    CUDA_CHECK_RETURN(cudaSetDevice(device[dev]));

    	    if(QualityTrimming) CUDA_CHECK_RETURN(cudaMemcpyAsync(dQualArray[dev], hReadQualstr[dev], ReadsToProcess[dev] * MAX_READ_LENGTH * sizeof(char), cudaMemcpyHostToDevice, upload_stream[dev]));
    	    CUDA_CHECK_RETURN(cudaMemcpyAsync(dIntArray[dev][dev_side], hReadLen[dev][dev_side], ReadsToProcess[dev] * sizeof(short), cudaMemcpyHostToDevice, upload_stream[dev]));
    	    CUDA_CHECK_RETURN(cudaEventRecord(QUAL_COPY_DONE[dev], upload_stream[dev]));

    	    CUDA_CHECK_RETURN(cudaStreamWaitEvent(compute_stream[dev], QUAL_COPY_DONE[dev], 0));
    	    if(QualityTrimming)
    	    {
    	    	if(verbose) printf("[%d] TrimmQuals\n", dev);
    	    	TrimmQuals<<<4096, 256, 0, compute_stream[dev]>>>(dQualArray[dev], dIntArray[dev][dev_side], ReadsToProcess[dev], MAX_READ_LENGTH);
    	    	CUDA_CHECK_RETURN(cudaGetLastError());
    	    }
    	    CUDA_CHECK_RETURN(cudaEventRecord(TrimmQ_DONE[dev], compute_stream[dev]));

    	    if(AdaptorTrimming || gRNATrimming) CUDA_CHECK_RETURN(cudaMemcpyAsync(dSequenceArray[dev], hReadSequence[dev][dev_side], ReadsToProcess[dev] * MAX_READ_LENGTH * sizeof(char), cudaMemcpyHostToDevice, upload_stream[dev]));
    	    CUDA_CHECK_RETURN(cudaEventRecord(SEQ_COPY_DONE[dev], upload_stream[dev]));

    	    CUDA_CHECK_RETURN(cudaStreamWaitEvent(compute_stream[dev], SEQ_COPY_DONE[dev], 0));
    	    if(gRNATrimming)
    	    {
    	    	if(verbose) printf("[%d] Find_gRNA\n", dev);
    	        Find_gRNA<<<4096, 256, 0, compute_stream[dev]>>>(dSequenceArray[dev], dIntArray[dev][dev_side], ReadsToProcess[dev], MAX_READ_LENGTH);
    	        CUDA_CHECK_RETURN(cudaGetLastError());
    	    }
    	    else
    	    {
    	    	if(AdaptorTrimming)
    	    	{
    	    		if(verbose) printf("[%d] TrimmAdaptors\n", dev);
    	    		TrimmAdaptors<<<4096, 256, 0, compute_stream[dev]>>>(dSequenceArray[dev], dIntArray[dev][dev_side], ReadsToProcess[dev], MAX_READ_LENGTH);
    	    		CUDA_CHECK_RETURN(cudaGetLastError());
    	    	}
    	    }
    	    CUDA_CHECK_RETURN(cudaEventRecord(TrimmA_DONE[dev], compute_stream[dev]));

    	    CUDA_CHECK_RETURN(cudaStreamWaitEvent(download_stream[dev], TrimmA_DONE[dev], 0));
    	    CUDA_CHECK_RETURN(cudaMemcpyAsync(hReadLen[dev][dev_side], dIntArray[dev][dev_side], ReadsToProcess[dev] * sizeof(short), cudaMemcpyDeviceToHost, download_stream[dev]));
    	    CUDA_CHECK_RETURN(cudaEventRecord(DTH_COPY_DONE[dev], download_stream[dev]));
    	    if(verbose) printf("[%d] Memcpy DtoH %d reads\n", dev, ReadsToProcess[dev]);

    	    ReadsSentToGPU[dev] = ReadsToProcess[dev];
    	    ReadsToProcess[dev] = 0;
       }
    }

    while(!ALL_DONE)
    {
    	ALL_DONE = true;                                            // setup for test if should stop if there is no more data on the GPU

    	if(!DONE)
        {
            for(int dev=0; dev<numGPUs; dev++)                         // read a block of data, copy QUALS to GPU and run QUAL kernel
        	{
            	ReadsToProcess[dev] = 0;

        	    CUDA_CHECK_RETURN(cudaSetDevice(device[dev]));

            	if(ReadsSentToGPU[dev]>0)
            	{
            		CUDA_CHECK_RETURN(cudaEventSynchronize(TrimmQ_DONE[dev]));
            	}

            	ReadsToProcess[dev] = ReadBlockOfReads(MAX_READS_ON_GPU, hReadLen[dev][dev_next_side], hReadSequence[dev][dev_next_side], hReadQualstr[dev]);
                if(ReadsToProcess[dev] < MAX_READS_ON_GPU)
                {
                	DONE = true;
                }
                if(ReadsToProcess[dev] == 0) break;

                ALL_DONE = false;                                      // there is more work on the GPUs that requires further processing

                if(QualityTrimming) CUDA_CHECK_RETURN(cudaMemcpyAsync(dQualArray[dev], hReadQualstr[dev], ReadsToProcess[dev] * MAX_READ_LENGTH * sizeof(char), cudaMemcpyHostToDevice, upload_stream[dev]));
        	    CUDA_CHECK_RETURN(cudaMemcpyAsync(dIntArray[dev][dev_next_side], hReadLen[dev][dev_next_side], ReadsToProcess[dev] * sizeof(short), cudaMemcpyHostToDevice, upload_stream[dev]));
        	    CUDA_CHECK_RETURN(cudaEventRecord(QUAL_COPY_DONE[dev], upload_stream[dev]));

        	    CUDA_CHECK_RETURN(cudaStreamWaitEvent(compute_stream[dev], QUAL_COPY_DONE[dev], 0));
        	    if(QualityTrimming)
        	    {
        	    	if(verbose) printf("[%d] TrimmQuals\n", dev);
        	    	TrimmQuals<<<4096, 256, 0, compute_stream[dev]>>>(dQualArray[dev], dIntArray[dev][dev_next_side], ReadsToProcess[dev], MAX_READ_LENGTH);
        	    	CUDA_CHECK_RETURN(cudaGetLastError());
        	    }
        	}
        }

        for(int dev=0; dev<numGPUs; dev++)                         // copy sequence and run Adaptor trim kernel
    	{
        	if(ReadsToProcess[dev]>0)
        	{
        	    CUDA_CHECK_RETURN(cudaSetDevice(device[dev]));

        	    if(ReadsSentToGPU[dev]>0)
        	    {
            	    CUDA_CHECK_RETURN(cudaEventSynchronize(TrimmA_DONE[dev]));
        	    }

        	    if(AdaptorTrimming || gRNATrimming) CUDA_CHECK_RETURN(cudaMemcpyAsync(dSequenceArray[dev], hReadSequence[dev][dev_next_side], ReadsToProcess[dev] * MAX_READ_LENGTH * sizeof(char), cudaMemcpyHostToDevice, upload_stream[dev]));
        	    CUDA_CHECK_RETURN(cudaEventRecord(SEQ_COPY_DONE[dev], upload_stream[dev]));

        	    CUDA_CHECK_RETURN(cudaStreamWaitEvent(compute_stream[dev], SEQ_COPY_DONE[dev], 0));
        	    if(gRNATrimming)
        	    {
        	    	if(verbose) printf("[%d] Find_gRNA\n", dev);
        	    	Find_gRNA<<<4096, 256, 0, compute_stream[dev]>>>(dSequenceArray[dev], dIntArray[dev][dev_next_side], ReadsToProcess[dev], MAX_READ_LENGTH);
        	    	CUDA_CHECK_RETURN(cudaGetLastError());
        	    }
        	    else
        	    {
        	    	if(AdaptorTrimming)
        	    	{
        	    		if(verbose) printf("[%d] TrimmAdaptors\n", dev);
        	    		TrimmAdaptors<<<4096, 256, 0, compute_stream[dev]>>>(dSequenceArray[dev], dIntArray[dev][dev_next_side], ReadsToProcess[dev], MAX_READ_LENGTH);
        	    		CUDA_CHECK_RETURN(cudaGetLastError());
        	    	}
        	    }
        	    CUDA_CHECK_RETURN(cudaEventRecord(TrimmA_DONE[dev], compute_stream[dev]));
        	}
    	}

        for(int dev=0; dev<numGPUs; dev++)                          // copy data from devices and write output file
    	{
            if(verbose) printf("[%d] > ReadsSentToGPU = %d\n", dev, ReadsSentToGPU[dev]);

        	if(ReadsSentToGPU[dev]>0)
        	{
        	    CUDA_CHECK_RETURN(cudaSetDevice(device[dev]));
        		CUDA_CHECK_RETURN(cudaEventSynchronize(DTH_COPY_DONE[dev]));
        		if(verbose) printf("[%d] Writing hReadLen[%d][%d]\n", dev, dev, dev_side);
        	    WriteBlockOfReads(ReadsSentToGPU[dev], hReadLen[dev][dev_side], MIN_OUT_READ_LEN);    // write output block
        	    ReadsProcessed += ReadsSentToGPU[dev];             // count how many reads were considered for writing - actual output has minimal length requirement
        	}

        	if(ReadsToProcess[dev]>0)
        	{
        	    CUDA_CHECK_RETURN(cudaStreamWaitEvent(download_stream[dev], TrimmA_DONE[dev], 0));
        	    CUDA_CHECK_RETURN(cudaMemcpyAsync(hReadLen[dev][dev_next_side], dIntArray[dev][dev_next_side], ReadsToProcess[dev] * sizeof(short), cudaMemcpyDeviceToHost, download_stream[dev]));
        	    CUDA_CHECK_RETURN(cudaEventRecord(DTH_COPY_DONE[dev], download_stream[dev]));
        	    if(verbose) printf("[%d] Memcpy DtoH %d reads\n", dev, ReadsToProcess[dev]);
        	}
    	    ReadsSentToGPU[dev] = ReadsToProcess[dev];
    	    ReadsToProcess[dev] = 0;
    	}
        side_flip = dev_side;                                      // dev_side alternate between 0 and 1 to use memory efficiently
        dev_side = dev_next_side;
        dev_next_side = side_flip;
    }
    return(ReadsProcessed);
}

void usage_help(void)
{
	printf("\nProgram Options:\n\n");
	printf("\t-i, -I ....... specify the path for the input file [required]\n");
	printf("\t-o, -O ....... specify the path for the output file [required]\n");
	printf("\t-len, -Len ... specify the minimum length for reads written to the output file\n");
	printf("\t-f ........... over write output file if it exists\n");
	printf("\t-gz .......... write output to file in compressed format (.gz)\n");
	printf("\t-pigz ........ write output to file in compressed format (.gz) using multithreaded compression\n");
	printf("\t-v, -V ....... verbose, prints details on processing to screen\n");
	printf("\t-a, -A ....... adaptor sequence, can contain Ns that will match any base in the read\n");
	printf("\t-mm, -MM ..... number of mismatches allowed in adaptor search (not including Ns)\n");
	printf("\t-Phred ....... sets the base of the qualities PHRED values [default 33], use AUTO for autodetection\n");
	printf("\t-Qmin ........ sets the minimum acceptable quality score for a good base [default 20]\n");
	printf("\t-Qualab ...... sets the minimum average base quality score for remainder of the read after a bad base [default 25]\n");
	printf("\t-dev, -Dev ... specify CUDA device to use (can be used repeatedly for selecting multiple GPUs,\n");
	printf("\t               if not specified an automatic search for GPUs with compute capability 3.5 and 6 GB memory is performed)\n");
	printf("\nAlgorithm selection:\n\n");
	printf("\t/TA .......... select adaptor trimming only (can be also used before option /gRNA to avoid quality trimming)\n");
	printf("\t/TQ .......... select quality trimming only\n");
	printf("\t/TAQ ......... select quality and adaptor trimming [default]\n");
	printf("\t/gRNA ........ select gRNA recovery, followed by start relative to start of adaptor sequence and length of gRNA to be retrieved\n");
	printf("\t               ( /gRNA 12 20 will retrieve 20 nucleotide read sequence starting 12 positions after the 5' end match of adaptor)\n");
	printf("\nProgrammer options variables:\n");
	printf("\t-G ........... set MAX_READS_ON_GPU variable\n");
	printf("\t-R ........... set MAX_READ_LENGTH variable\n");
	printf("\nexample usage:\n");
	printf("\tPreprocessReads -dev 0 -dev 1 -i infile.fastq -o outfilename.fastq -len 26 -a ATTAGATCGATGCTA -mm 2\n\n");
	exit(3);
}

void invalid_command_line_option(const char *errorstr)
{
	printf("\nERROR %s\n", errorstr);
	printf("\nUsage ----------------------------------------------\n");
	usage_help();
}

void parseCommandLineArguments(int argc, char *argv[])
{
	bool inputfile_set = false;
	bool outputfile_set = false;
	bool adaptor_set = false;
	bool mismatch_set = false;
	bool outlen_set = false;
	for(int n=1; n<argc; n++)
	{
		if(argv[n][0] == '-')
		{
			char a = argv[n][1];
			switch (a)
			{
				case 'I':
				case 'i':
					if(inputfile_set) invalid_command_line_option("option -i: input file name already specified.");
					if(++n<argc)
					{
						if(strlen(argv[n])>400) invalid_command_line_option("option -i: input file path length exceeds limit.");
						strcpy(input_filename, argv[n]);
						inputfile_set=true;
						break;
					}
					else invalid_command_line_option("option -i: INPUT file name missing.");
					break;
				case 'O':
				case 'o':
					if(outputfile_set) invalid_command_line_option("option -o: output file name already specified.");
					if(++n<argc)
					{
						if(strlen(argv[n])>400) invalid_command_line_option("option -o: output file path length exceeds limit.");
						strcpy(output_filename, argv[n]);
						outputfile_set=true;
						break;
					}
					else invalid_command_line_option("option -o: OUTPUT file name missing.");
					break;
				case 'l':
				case 'L':
					if(argv[n][2]!='e' || argv[n][3]!='n') invalid_command_line_option("invalid command line option encountered.");
					if(outlen_set) invalid_command_line_option("option -len: minimum read length to include in output file is already specified.");
					if(++n<argc)
					{
						int ol=atoi(argv[n]);
						if(ol<0 || ol>MAX_READ_LENGTH) invalid_command_line_option("option -len: encountered invalid minimum read length to include in output file.");
						MIN_OUT_READ_LEN = (short) ol;
						outlen_set=true;
						break;
					}
					else invalid_command_line_option("option -len: expected minimum output read length not found.");
					break;

				case 'A':
				case 'a':
					if(adaptor_set) invalid_command_line_option("option -a: adaptor sequence already specified.");
					if(++n<argc)
					{
						if(strlen(argv[n])>MAX_ADAPTOR_LENGTH) invalid_command_line_option("option -a: adaptor sequence length exceeds limit.");
						if(strlen(argv[n])<3) invalid_command_line_option("option -a: adaptor sequence length too short.");
						strcpy(Adaptor, argv[n]);
						adaptor_set=true;
						break;
					}
					else invalid_command_line_option("option -a: expected adaptor sequence.");
					break;
				case 'M':
				case 'm':
					if(argv[n][2]!='m' && argv[n][2]!='M') invalid_command_line_option("invalid command line option encountered.");
					if(mismatch_set) invalid_command_line_option("option -mm: number of mismatches is already specified.");
					if(++n<argc)
					{
						int mm=atoi(argv[n]);
						if(mm>20) invalid_command_line_option("option -mm: number of mismatches exceeds limit.");
						mismatches = (short) mm;
						mismatch_set=true;
						break;
					}
					else invalid_command_line_option("option -mm: expected number of mismatches not found.");
					break;
				case 'd':
				case 'D':
					if(argv[n][2]!='e' || argv[n][3]!='v') invalid_command_line_option("invalid command line option encountered.");
					if(numGPUs>=MAX_GPU_NUM) invalid_command_line_option("option -dev: number of CUDA devices exceeds maximum.");
					if(++n<argc)
					{
						int d=atoi(argv[n]);
						if(d<0 || d>20) invalid_command_line_option("option -dev: invalid device number.");
						device[numGPUs] = d;
						numGPUs++;
						break;
					}
					else invalid_command_line_option("option -dev: expected CUDA device number not found.");
					break;
				case 'f':
					overwrite_outputfile = true;                                        // -f
					break;
				case 'g':                                                               // -gz
					if((argv[n][2]=='z') && (strlen(argv[n])==3))
					{
						write_gz_output = true;
						break;
					}
					else invalid_command_line_option("invalid command line option encountered.");
					break;
				case 'p':
				case 'P':
					if((strcmp(argv[n], "-Pigz")==0) || (strcmp(argv[n], "-pigz")==0))  // -pigz
					{
						write_pigz_output = true;
						break;
					}
					else                                                                // -phred
					{
						if(argv[n][2]!='h' || argv[n][3]!='r' || argv[n][4]!='e' || argv[n][5]!='d') invalid_command_line_option("invalid command line option encountered.");
						if(++n<argc)
						{
							if((strcmp(argv[n], "AUTO")==0) || (strcmp(argv[n], "auto")==0))
							{
								auto_phred_detect = true;
								break;
							}
							else
							{
								int d=atoi(argv[n]);
								if(d<0 || d>100) invalid_command_line_option("option -Phred: invalid phred base.");
								PHRED = d;
								break;
							}
						}
						else invalid_command_line_option("option -phred: expected number or AUTO not found.");
					}
					break;
				case 'Q':
				case 'q':
					if((strcmp(argv[n], "-Qmin")==0) || (strcmp(argv[n], "-qmin")==0))
					{
						if(++n<argc)
						{
							int qm=atoi(argv[n]);
							if(qm<10 || qm>64) invalid_command_line_option("option -qmin: outside 10 to 64 limit.");
							MIN_BASE_QUAL = (short) qm;
							break;
						}
						else invalid_command_line_option("option -qmin: expected number not found.");
					}
					else
					{
						if((strcmp(argv[n], "-Qualab")==0) || (strcmp(argv[n], "-qualab")==0))
						{
							if(++n<argc)
							{
								int qm=atoi(argv[n]);
								if(qm<10 || qm>64) invalid_command_line_option("option -qualab: outside 10 to 64 limit.");
								MIN_QUAL_AFTER_BAD = (short) qm;
								break;
							}
							else invalid_command_line_option("option -qmin: expected number not found.");
						}
					}
					break;
				case 'v':
				case 'V':
					verbose = true;
					break;
				case 'H':
				case 'h':
				case '?':
					printf("\nHELP ----------------------------------------------\n");
					usage_help();
					break;
				case 'G':                                                           // programmer option to set MAX_NUM_READS_ON_GPU parameter
					if(++n<argc)
					{
						int mrg=atoi(argv[n]);
						if(mrg<=0) invalid_command_line_option("invalid command line option encountered.");
						MAX_READS_ON_GPU = mrg;
						break;
					}
					else invalid_command_line_option("invalid command line option encountered.");
					break;
				case 'R':                                                           // programmer option to set MAX_READ_LENGTH parameter
					if(++n<argc)
					{
						int mrl=atoi(argv[n]);
						if(mrl<=0 || mrl>10000) invalid_command_line_option("invalid command line option encountered.");
						MAX_READ_LENGTH = (short) mrl;
						break;
					}
					else invalid_command_line_option("invalid command line option encountered.");
					break;

				default:
					invalid_command_line_option("invalid command line option encountered.");
			}
		}
		else
		{
			if(argv[n][0] == '/')
			{
				char a = argv[n][1];
				switch (a)
				{
					case 'T':
						AdaptorTrimming = false;
						QualityTrimming = false;
						gRNATrimming = false;
						if(strlen(argv[n]) == 3)
						{
							if(argv[n][2]=='A')
							{
								AdaptorTrimming = true;
								break;
							}
							else
							{
								if(argv[n][2]=='Q')
								{
									QualityTrimming = true;
									break;
								}
								else invalid_command_line_option("option /T: algorithm selection is invalid.");
							}
						}
						else
						{
							if(strlen(argv[n])==3 && argv[n][2]=='A' && argv[n][3]=='Q')
							{
								AdaptorTrimming = true;
								QualityTrimming = true;
								break;
							}
							else invalid_command_line_option("option /T: algorithm selection is invalid.");
						}
						break;
					case 'g':
						if(strlen(argv[n]) == 5 && argv[n][2]=='R' && argv[n][3]=='N' && argv[n][4]=='A')
						{
							gRNATrimming = true;
							if(++n < argc)
							{
								gRNA_start = atoi(argv[n]);
							}
							else invalid_command_line_option("option /gRNA: expected start position relative to adaptor 5' not found.");
							if(++n < argc)
							{
								gRNA_length = atoi(argv[n]);
							}
							else invalid_command_line_option("option /gRNA: expected gRNA length not found.");
							break;
						}
						else invalid_command_line_option("option /g: algorithm selection is invalid.");
						break;
				}
			}
			else invalid_command_line_option("parsing command line options.");
		}
	}
	if(!inputfile_set || !outputfile_set || ((AdaptorTrimming || gRNATrimming) && !adaptor_set)) invalid_command_line_option("missing required command line options.");
	if(strcmp(input_filename, output_filename)==0)
	{
		printf("ERROR: encountered identical input and output file paths!\n");         // input file cannot be identical to output file
		exit(253);
	}
	printf("Setting parameters:\n");
	printf("\tInput file name  : %s\n", input_filename);
	printf("\tOutput file name : %s\n", output_filename);
	if(AdaptorTrimming || gRNATrimming) printf("\tAdaptor sequence : %s\n", Adaptor);
	if(verbose)
	{
		if(AdaptorTrimming || gRNATrimming) printf("\tMismatches = %d\n", mismatches);
		printf("\tTrimming: ");
		if(QualityTrimming) printf("Quality ");
		if(AdaptorTrimming) printf("Adaptors ");
		if(gRNATrimming) printf("Find gRNA");
		printf("\n");
		if(QualityTrimming)
		{
			printf("\tMinimum accepted base quality = %d; average quality of read after bad base = %d\n", MIN_BASE_QUAL, MIN_QUAL_AFTER_BAD);
		}
		printf("\tFastQ Quality Encoding: ");
		if(!auto_phred_detect) printf("PHRED+%d\n", PHRED); else printf("<AUTO DETECTION>\n");
		printf("\tMAX_READ_LENGTH  : %d\n", MAX_READ_LENGTH);
		printf("\tMinimum read length for output = %d\n", MIN_OUT_READ_LEN);
		if(overwrite_outputfile) printf("\tForce overwrite output file if it exists.\n");
	}
	if(numGPUs>0)
	{
		printf("\nUsing CUDA devices: ");
		for(int n=0; n<numGPUs; n++) printf("[%d]", device[n]);
		printf("\n");
	}
}

int main(int argc, char *argv[])
{
	time_t start_t = time(NULL);                        // for timing the process
	time_t end_t;
    int reads_processed = 0;

    parseCommandLineArguments(argc, argv);                        // get parameters from the command line
    if(numGPUs == 0)
    {
        if(FindGPUs() == 0)                             // check for GPUs
        {
            printf("\nERROR: No suitable CUDA devices found ... exiting\n");
            exit(1);
        }
        else
        {
            if(numGPUs == 1)
            	printf("\nOne CUDA device with at least compute capability 3.5 and 6 GB memory found.\n");
            else
                printf("\n%d CUDA devices with at least compute capability 3.5 and 6 GB memory found.\n", numGPUs);
        }
    }

    if(auto_phred_detect)
    {
    	if(Detect_FastQ_Format_PHRED())
    	{
    		if(verbose) printf("<AUTO DETECTION> FastQ Quality Encoding: PHRED+%d\n", PHRED);
    	}
    	else
    	{
    		printf("\nERROR: Unable to detect FastQ format version.\n");
    		exit(32);
    	}
    }

    OpenFiles();                                    // access files
    InitGPUs();
    SetAdaptorSequence(Adaptor, mismatches, PHRED, MIN_BASE_QUAL, MIN_QUAL_AFTER_BAD);

    reads_processed = PreprocessReadsOnGPU();

    printf("\nReads processed: %d\n", reads_processed);
    printf("Reads in input file %s : %d\n",input_filename, hNumberOfReads);
    printf("Reads written to output file %s : %d\n", output_filename, hReadsWritten);
    if(FastQ_WARNING) printf("Encountered potential issues with the FastQ file format in %s\n", input_filename);

    CloseFiles();
    FreeGPUs();

    end_t = time(NULL);
    printf("\nTotal time elapsed: %.0lf seconds\n", difftime(end_t, start_t));
//    printf("\n<Press RTN key to exit>\n");
//    getchar();

    return(0);
}
