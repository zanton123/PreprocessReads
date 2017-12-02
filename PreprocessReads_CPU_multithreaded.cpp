//#define __CODEBLOCKS__
    // uncomment for alternate CODE BLOCKS MinGW code path to work around popen() and strnlen()
#ifndef __CODEBLOCKS__
char VERSION[] = "0.1.cpu";
#else
char VERSION[] = "0.1.cpu.CODEBLOCKS";
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef __CODEBLOCKS__
size_t strnlen(const char *str, size_t max)
{
    const char *end =(char *) memchr(str, 0, max);
    return end ? (size_t) (end-str) : max;
}
#endif

#include <time.h>
#include <iostream>       // std::cout
#include <thread>         // std::thread
#include <mutex>          // std::mutex, std::unique_lock
#include <stdexcept>      // std::logic_error

#include <zlib.h>

int hNumberOfReads = 0;   // holds number of reads that are processed
int hReadsWritten = 0;    // number of reads written to output file
int ReadStartWriteNext = 0;// used to preserve order of reads in the output file when running many threads
std::mutex ReadFile_mutex;
std::mutex WriteFile_mutex;
char input_filename[400]="/home/linux/cuda-workspace/PreprocessReads/test.txt"; // input file name
char output_filename[400]="/mnt/sdc/Sequences/out.txt";                         // output file name
char GZIP_PIPE[410] = "gzip -> ";       // option -gz will fallback to uncompressed if error on pipe opening
char PIGZ_PIPE[410] = "pigz -> ";       // option -pigz for multithreaded compression will fallback to gz and then to uncompressed
gzFile inputFile; 						// input FastQ file with read name sequence and quality
FILE *outputFile;						// output FastQ file
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

#define MAX_READ_LENGTH 256             // maximal length of reads in input file
int MIN_OUT_READ_LEN = 26;              // minimum length of reads after trimming to write to output file

int MIN_BASE_QUAL = 20;			    	// parameter for minimum quality for good base
int MIN_QUAL_AFTER_BAD = 25;			// parameter for average base quality after bad base
int PHRED = 33;			        	    // subtract from qualities for numerical value of solexa
bool auto_phred_detect = false;         // if true try to detect PHRED format in FastQ input file

#define MAX_ADAPTOR_LENGTH 100          // maximum length of adaptor sequence for trimming, can contain N for any base
#define MAX_FASTQ_BUF_LINE 500          // structure of buffer for read from FastQ file
#define MAX_READ_NAME_LENGTH 256		// maximal length of the name of reads to be stored in output file
char Adaptor[MAX_ADAPTOR_LENGTH] = "ATCTCGTATGCCGTCTTCTGCTT"; // adaptor on CPU side from parsing command line
short mismatches = 2;					// holds maximum number of mismatches on cpu side

#define MAX_THREAD_NUM 64				// for AMD Threadripper and Intel XEON CPUs

int num_Threads=8;
int MAX_READS_PER_THREAD = 10000;       // maximum number of reads that are processed simultaneously per thread

struct Read_t {
	int ReadLen;
	char ReadSequence[MAX_READ_LENGTH];
	char ReadQualstr[MAX_READ_LENGTH];
	char ReadName[MAX_READ_NAME_LENGTH];
};

struct worker_data 
{
	bool AdaptorTrimming;			// algorithm selectors
	bool QualityTrimming;
	bool gRNATrimming;
	bool verbose;
	
	int MAX_READS_PER_THREAD;
	Read_t *Reads;
	int ReadsProcessed;

	char ADAPTOR[MAX_ADAPTOR_LENGTH];
	int ADAPTOR_LEN;
	int MAX_MM;
	
	int QUAL_PADDING;
	int MIN_QUAL;
	int MIN_QUAL_AFTER_BAD;
	
	int MIN_OUT_READ_LEN;
	int ErrNo;
};

std::thread tid[MAX_THREAD_NUM];
worker_data data[MAX_THREAD_NUM];

struct read
{
    char name[MAX_FASTQ_BUF_LINE];
    char sequence[MAX_FASTQ_BUF_LINE];
    char descr[MAX_FASTQ_BUF_LINE];
    char qualstr[MAX_FASTQ_BUF_LINE];
};                                 // buffer one read from FastQ file i 4 lines

void Find_gRNA(Read_t *Reads, int ReadsToProcess, char *ADAPTOR, int ADAPTOR_LEN, int MAX_MM)
{
  for(long idx = 0; idx < ReadsToProcess; idx++)
  {
      if (idx < ReadsToProcess)
      {
          char *sequence = Reads[idx].ReadSequence; 
          int rl=Reads[idx].ReadLen;
          int n=0;
          int p=0;
          int mm=0;
          Reads[idx].ReadLen=0;
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
                  Reads[idx].ReadLen=p;
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
                  Reads[idx].ReadLen=p;
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
                  Reads[idx].ReadLen=p;
              }
          }
      }
  }
}

void TrimmAdaptors(Read_t *Reads, int ReadsToProcess, char *ADAPTOR, int ADAPTOR_LEN, int MAX_MM)
{
  for(long idx=0; idx < ReadsToProcess; idx++)
  {
      if (idx < ReadsToProcess)
      {
          char *sequence = Reads[idx].ReadSequence;
          int rl=Reads[idx].ReadLen;
          int n=0;
          int p=0;
          int mm=0;
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
                  Reads[idx].ReadLen=p;
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
                  Reads[idx].ReadLen=p;
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
                  Reads[idx].ReadLen=p;
              }
          }
      }
  }
}

void TrimmQuals(Read_t *Reads, int ReadsToProcess, int QUAL_PADDING, int MIN_QUAL, int MIN_QUAL_AFTER_BAD)
{
  for(int idx=0; idx < ReadsToProcess; idx++)
  {
      if (idx < ReadsToProcess)
      {
          char *q = Reads[idx].ReadQualstr;
          int rl=Reads[idx].ReadLen;
          for(int i=0; i<rl; i++)
          {
              if ((q[i]-QUAL_PADDING)<=MIN_QUAL)
              {
                 int s=0;
                 int r=i+1;
                 while (r<rl)
                 {
                    s=s+q[r]-QUAL_PADDING;
                    r++;
                 }
                 if(s<(MIN_QUAL_AFTER_BAD*(rl-i+1)))
                 {
//                  printf("MIN_QUAL_AFTER_BAD %d, QUAL_PADDING %d, MIN_QUAL %d, idx %d, Reads[idx].ReadLen %d, q[i] %d s %d i %d rl %d %s", MIN_QUAL_AFTER_BAD, QUAL_PADDING, MIN_QUAL, idx, Reads[idx].ReadLen, q[i], s, i, rl, Reads[idx].ReadName);
                    Reads[idx].ReadLen=i;
                    break;
                 }
              }
          }

      }
  }
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
#ifndef __CODEBLOCKS__
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
#else  // __CODEBLOCKS__ code path
    if(write_gz_output || write_pigz_output)
        outputFile = (FILE *) gzopen(output_filename, "w");
    else
        outputFile = fopen(output_filename, "w");
#endif
    if(outputFile == NULL)
    {
    	printf("ERROR opening output file: %s\n", output_filename);
    	exit(255);
    }
}

int ReadBlockOfReads(int max_num_to_read, Read_t *Reads)
{
	int n = 0;                                                                    // returns the number of reads read
	int length;
	struct read Read;

    while(gzgets(inputFile, Read.name, MAX_FASTQ_BUF_LINE) != NULL)
    {
    	if(Read.name[0] == 13) break;				// end if CR <13> or LF <10> are encountered
    	if(Read.name[0] == 10) break;
    	if(gzgets(inputFile, Read.sequence, MAX_FASTQ_BUF_LINE) == NULL)
        {
            printf("ERROR parsing input!\n");
            exit(255);
        }
        if(gzgets(inputFile, Read.descr, MAX_FASTQ_BUF_LINE) == NULL)
        {
            printf("ERROR parsing input!\n");
            exit(255);
        }
        if(gzgets(inputFile, Read.qualstr, MAX_FASTQ_BUF_LINE) == NULL)
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

        memcpy(Reads[n].ReadSequence, Read.sequence, length);
        memcpy(Reads[n].ReadQualstr, Read.qualstr, length);
        Reads[n].ReadLen = length;
        length=strlen(Read.name)-1;
        if(length > MAX_READ_NAME_LENGTH-1) length = MAX_READ_NAME_LENGTH-1;
        memcpy(Reads[n].ReadName, Read.name, length);
        Reads[n].ReadName[length]=0;

        n++;
        if(n >= max_num_to_read) break;

    }
    hNumberOfReads = hNumberOfReads + n;                           // update global variable to total number of reads
    return(n);
}

int WriteBlockOfReads(int num_of_reads_to_write, const Read_t *Reads, int min_out_read_length)
{
    int n=0;
    int w=0;
    int length;
    char str_buffer[MAX_FASTQ_BUF_LINE];
    struct read Read;									// temporary storage for read to be written

    if(num_of_reads_to_write == 0) return(0);
    if(verbose) printf("I/O > num_of_reads_to_write = %d\n", num_of_reads_to_write);
    Read.descr[0]='+'; Read.descr[1]=10; Read.descr[2]=0;      // generate generic read descriptor line
    while(n < num_of_reads_to_write)
    {
        length=strnlen(Reads[n].ReadName, MAX_READ_NAME_LENGTH); // regenerate writeable read name
        memcpy(Read.name, Reads[n].ReadName, length);
		Read.name[length]=10;
		Read.name[length+1]=0;

        length = Reads[n].ReadLen;

        if(gRNATrimming)                        // isolate the gRNA from adaptor match
        {
        	if(length>0)
        	{
        		memcpy(str_buffer, &Reads[n].ReadSequence[length+gRNA_start], gRNA_length);
        		memcpy(Read.sequence, str_buffer, gRNA_length);
        		Read.sequence[gRNA_length]=10;
        		Read.sequence[gRNA_length+1]=0;

        		memcpy(str_buffer, &Reads[n].ReadQualstr[length+gRNA_start], gRNA_length);
        		memcpy(Read.qualstr, str_buffer, gRNA_length);
        		Read.qualstr[gRNA_length]=10;
        		Read.qualstr[gRNA_length+1]=0;
#ifndef __CODEBLOCKS__
                fprintf(outputFile,"%s", Read.name);
                fprintf(outputFile,"%s", Read.sequence);
                fprintf(outputFile,"%s", Read.descr);
                fprintf(outputFile,"%s", Read.qualstr);
#else
                if(write_gz_output || write_pigz_output)
                {
                    gzprintf((gzFile) outputFile,"%s", Read.name);
                    gzprintf((gzFile) outputFile,"%s", Read.sequence);
                    gzprintf((gzFile) outputFile,"%s", Read.descr);
                    gzprintf((gzFile) outputFile,"%s", Read.qualstr);
                }
                else
                {
                    fprintf(outputFile,"%s", Read.name);
                    fprintf(outputFile,"%s", Read.sequence);
                    fprintf(outputFile,"%s", Read.descr);
                    fprintf(outputFile,"%s", Read.qualstr);
                }
#endif
                w++;
        	}
        }
        else                                     // normal quality and adaptor trimming
        {
            if(length >= min_out_read_length)
            {
                memcpy(Read.sequence, Reads[n].ReadSequence, length);
                memcpy(Read.qualstr, Reads[n].ReadQualstr, length);

                Read.sequence[length]=10;        // write only reads that are longer than minimum to output file
                Read.sequence[length+1]=0;

                Read.qualstr[length]=10;
                Read.qualstr[length+1]=0;
#ifndef __CODEBLOCKS__
                fprintf(outputFile,"%s", Read.name);
                fprintf(outputFile,"%s", Read.sequence);
                fprintf(outputFile,"%s", Read.descr);
                fprintf(outputFile,"%s", Read.qualstr);
#else
                if(write_gz_output || write_pigz_output)
                {
                    gzprintf((gzFile) outputFile,"%s", Read.name);
                    gzprintf((gzFile) outputFile,"%s", Read.sequence);
                    gzprintf((gzFile) outputFile,"%s", Read.descr);
                    gzprintf((gzFile) outputFile,"%s", Read.qualstr);
                }
                else
                {
                    fprintf(outputFile,"%s", Read.name);
                    fprintf(outputFile,"%s", Read.sequence);
                    fprintf(outputFile,"%s", Read.descr);
                    fprintf(outputFile,"%s", Read.qualstr);
                }
#endif
                w++;
             }
//           else printf("Read %s (%d) [%d] not written!\n", Read.sequence, length, min_out_read_length);
        }
        n++;
    }
    hReadsWritten = hReadsWritten + w;
    return(n);
}

void CloseFiles(void)
{
    gzclose(inputFile);
#ifndef __CODEBLOCKS__
    if(write_gz_output || write_pigz_output) pclose(outputFile); else fclose(outputFile);
#else   // __CODEBLOCKS__ code path
     if(write_gz_output || write_pigz_output) gzclose((gzFile) outputFile); else fclose(outputFile);
#endif
}

bool Detect_FastQ_Format_PHRED(void)
{
	int n = 0;                                                                    // returns the number of reads read
	int length;
	unsigned int min_qual = 65;
	unsigned int max_qual = 65;
	bool deciphered = false;
	struct read Read;

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
        	if(min_qual > (unsigned int)Read.qualstr[n])
        	{
        		min_qual = (unsigned int) Read.qualstr[n];
        	}
        	else
        	{
            	if(max_qual < (unsigned int)Read.qualstr[n]) max_qual = (unsigned int) Read.qualstr[n];

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

void WorkerThread(int tID)
{
	bool DONE=false;
	bool wait_for_write;
	int ReadsToProcess=0;
	int ReadStartProcess=0;
	data[tID].ReadsProcessed=0;
	data[tID].ErrNo=0;
	std::unique_lock<std::mutex> read_lock(ReadFile_mutex);
	read_lock.unlock();
	std::unique_lock<std::mutex> write_lock(WriteFile_mutex);
	write_lock.unlock();
	while(!DONE)
	{
		read_lock.lock();  // read a block of reads in
		ReadStartProcess=hNumberOfReads;  // save the total number of reads read before to know when to write processed reads
		ReadsToProcess=ReadBlockOfReads(data[tID].MAX_READS_PER_THREAD, data[tID].Reads); // read a new block of reads
		read_lock.unlock();
		if(ReadsToProcess==0) DONE=true;
		else
		{
			// trimm the reads
			if(data[tID].QualityTrimming) TrimmQuals(data[tID].Reads, ReadsToProcess, data[tID].QUAL_PADDING, data[tID].MIN_QUAL, data[tID].MIN_QUAL_AFTER_BAD);
			if(data[tID].gRNATrimming) Find_gRNA(data[tID].Reads, ReadsToProcess, data[tID].ADAPTOR, data[tID].ADAPTOR_LEN, data[tID].MAX_MM);
			else
				if(data[tID].AdaptorTrimming) TrimmAdaptors(data[tID].Reads, ReadsToProcess, data[tID].ADAPTOR, data[tID].ADAPTOR_LEN, data[tID].MAX_MM);

			wait_for_write=true;	// write reads to file in the correct order based on ReadStartProcess
			while(wait_for_write)
			{
				write_lock.lock();
				if(ReadStartWriteNext==ReadStartProcess)    // unlock and retry for read start for block
				{
					// write a block of reads
					if (data[tID].verbose) printf("TR[%d] ", tID);
					int reads_written=WriteBlockOfReads(ReadsToProcess, data[tID].Reads, data[tID].MIN_OUT_READ_LEN);
					ReadStartWriteNext+=ReadsToProcess;      // update for start of next block of reads
					data[tID].ReadsProcessed+=reads_written;				 // update data record
					if(reads_written!=ReadsToProcess) data[tID].ErrNo=-1;  // error flag if number written is incorrect
					wait_for_write=false;
				}
				write_lock.unlock();
				if(wait_for_write) std::this_thread::yield();
			}
		}
	}
}

int PreprocessReadsWithThreads(void)
{
    if(verbose) printf("\nPreprocessing reads multithreaded on CPU\n");
	int ReadsProcessed = 0;
																  // initialize threads and allocate memory
	if(num_Threads>MAX_THREAD_NUM)
	{
		printf("Number of threads (%d) greater than maximal number of threads that can be processed (%d), Exiting ...\n", num_Threads, MAX_THREAD_NUM);
		exit(255);
	}
	for (int tr=0; tr<num_Threads; tr++)
    {
		if(verbose) printf("Thread[%d] Allocating CPU memory for reads\n", tr);
		data[tr].Reads = (Read_t *) malloc( MAX_READS_PER_THREAD * sizeof(Read_t) ); // buffer array for reads
		if(data[tr].Reads==NULL)
		{
			printf("ERROR - Memory allocation failed.\n");
			if(tr<1)
			{
				printf("CRITICAL ERROR: No thread running after memory allocation failure!\nExiting ...\n");
				exit(1);
			}
			else
			{
				num_Threads=tr;
				if(verbose) printf("Running with %d threads.\n", num_Threads);
				continue;
			}
		}
		else
		{
			data[tr].MAX_READS_PER_THREAD=MAX_READS_PER_THREAD;
			data[tr].AdaptorTrimming = AdaptorTrimming;
			data[tr].QualityTrimming = QualityTrimming;
			data[tr].gRNATrimming = gRNATrimming;
			data[tr].verbose = verbose;

			strncpy(data[tr].ADAPTOR, Adaptor, MAX_ADAPTOR_LENGTH);
			data[tr].ADAPTOR_LEN=strnlen(Adaptor, MAX_ADAPTOR_LENGTH);
			data[tr].MAX_MM=mismatches;
			
			data[tr].QUAL_PADDING=PHRED;
			data[tr].MIN_QUAL=MIN_BASE_QUAL;
			data[tr].MIN_QUAL_AFTER_BAD=MIN_QUAL_AFTER_BAD;
			
			data[tr].MIN_OUT_READ_LEN=MIN_OUT_READ_LEN;
			data[tr].ErrNo=0;
		}
	}
	for (int tr=0; tr<num_Threads; tr++) tid[tr]=std::thread(WorkerThread, tr);    // launch threads
	for (int tr=0; tr<num_Threads; tr++) tid[tr].join();		         		   // wait on ALL threads to complete

	for (int tr=0; tr<num_Threads; tr++)						                   // collect data and cleanup
	{	
		
//		printf("idx %d, max_per_th %d, ad:%s (%d), mm=%d, Phred=%d, minq=%d, qab=%d, min_len_out=%d, ErrNo:%d\n", tr, data[tr].MAX_READS_PER_THREAD, data[tr].ADAPTOR, data[tr].ADAPTOR_LEN, data[tr].MAX_MM, data[tr].QUAL_PADDING, data[tr].MIN_QUAL, data[tr].MIN_QUAL_AFTER_BAD, data[tr].MIN_OUT_READ_LEN, data[tr].ErrNo);
		
		ReadsProcessed+=data[tr].ReadsProcessed;
		if(data[tr].Reads!=NULL) free(data[tr].Reads);
		if(data[tr].ErrNo!=0) printf("\nWARNING: Thread[%d] carries ErrorNo %d\n", tr, data[tr].ErrNo);
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
	printf("\t-Threads ..... specify the number of CPU threads to use (default 8)\n");
	printf("\nAlgorithm selection:\n\n");
	printf("\t/TA .......... select adaptor trimming only (can be also used before option /gRNA to avoid quality trimming)\n");
	printf("\t/TQ .......... select quality trimming only\n");
	printf("\t/TAQ ......... select quality and adaptor trimming [default]\n");
	printf("\t/gRNA ........ select gRNA recovery, followed by start relative to start of adaptor sequence and length of gRNA to be retrieved\n");
	printf("\t               ( /gRNA 12 20 will retrieve 20 nucleotide read sequence starting 12 positions after the 5' end match of adaptor)\n");
	printf("\nProgrammer options variables:\n");
	printf("\t-G ........... set MAX_READS_PER_THREAD variable\n");
	printf("\t-R ........... set MAX_READ_LENGTH variable, this option is unused on CPU and set to %d \n", MAX_READ_LENGTH);
	printf("\nexample usage:\n");
	printf("\tPreprocessReads -threads 4 -i infile.fastq -o outfilename.fastq -len 26 -a ATTAGATCGATGCTA -mm 2\n\n");
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
	int arg_len;
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
				case 't':   // -Threads -threads option is only available on CPU version of PreprocessReads
				case 'T':   // default is 8 Threads if not specified to keep compatibility to GPU version, see also -Dev -dev
					arg_len=strlen(argv[n]);
					if(arg_len!=8) invalid_command_line_option("invalid command line option encountered.");
					if(strncmp(&argv[n][2],"hreads", arg_len-2)) invalid_command_line_option("invalid command line option encountered.");
					if(++n<argc)
					{
						int t=atoi(argv[n]);
						if(t<1 || t>MAX_THREAD_NUM) invalid_command_line_option("option -Threads: invalid thread number.");
						num_Threads=t;
						break;
					}
					else invalid_command_line_option("option -Threads: expected thread number not specified.");
					break;
				case 'd':  // -Dev -dev optiom is parsed but not used for compatibility with GPU version of PreprocessReads
				case 'D':
					if(argv[n][2]!='e' || argv[n][3]!='v') invalid_command_line_option("invalid command line option encountered.");
// CPU version:  	if(numGPUs>=MAX_GPU_NUM) invalid_command_line_option("option -dev: number of CUDA devices exceeds maximum.");
					if(++n<argc)
					{
						int d=atoi(argv[n]);
						if(d<0 || d>20) invalid_command_line_option("option -dev: invalid device number.");
// CPU version: 		device[numGPUs] = d;
//						numGPUs++;
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
								if(d<1 || d>100) invalid_command_line_option("option -Phred: invalid phred base.");
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
				case 'G':                                                           // programmer option to set MAX_NUM_READS_PER_THREAD parameter
					if(++n<argc)
					{
						int mrt=atoi(argv[n]);
						if(mrt <= 0) invalid_command_line_option("invalid command line option encountered.");
						MAX_READS_PER_THREAD=mrt;
						break;
					}
					else invalid_command_line_option("invalid command line option encountered.");
					break;
				case 'R':                                                           // programmer option to set MAX_READ_LENGTH parameter
					if(++n<argc)
					{
						int mrl=atoi(argv[n]);
						if(mrl<=0 || mrl>10000) invalid_command_line_option("invalid command line option encountered.");
// CPU version ignore	MAX_READ_LENGTH = (short) mrl;
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
	printf("\nUsing %d CPU threads.\n", num_Threads);
}

int main(int argc, char *argv[])
{
	printf("PreprocessReads multithreaded CPU version %s\n", VERSION);
	time_t start_t = time(NULL);                        // for timing the process
	time_t end_t;
    int reads_processed = 0;

    parseCommandLineArguments(argc, argv);                        // get parameters from the command line
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

    OpenFiles();
                                        // access files
	if(verbose)
	{
		int AdLen = strnlen(Adaptor, MAX_ADAPTOR_LENGTH); // get length of adaptor sequence
		printf("Using Adaptor sequence %s, length %d \n", Adaptor, AdLen);
		printf("\tmismatches %d, PHRED+%d on GPU \n", mismatches, PHRED);
		printf("\tminimum base quality %d, average read quality after bad %d \n", MIN_BASE_QUAL, MIN_QUAL_AFTER_BAD);
	}
	
    reads_processed = PreprocessReadsWithThreads();

    CloseFiles();

    printf("\nReads processed: %d\n", reads_processed);
    printf("Reads in input file %s : %d\n",input_filename, hNumberOfReads);
    printf("Reads written to output file %s : %d\n", output_filename, hReadsWritten);
    if(FastQ_WARNING) printf("\nWARNING: Encountered potential issues with the FastQ file format in %s\n", input_filename);

    end_t = time(NULL);
    printf("\nTotal time elapsed: %.0lf seconds\n", difftime(end_t, start_t));
//    printf("\n<Press RTN key to exit>\n");
//    getchar();

    return(0);
}
