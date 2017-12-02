#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#define MAX_FASTQ_BUF_LINE 500          // structure of buffer for read from FastQ file

void Usage(void) {
	printf("\nUSAGE: generate_reads <NUMBER> <LENGTH> <FILENAME>\n");
	printf("\tNUMBER ..... specifies the number of reads that will be generated.\n");
	printf("\tLENGTH ..... specified the length of the reads that will be generated.\n");
	printf("\tFILENAME ... file name for output. For filenames endining in .gz zlib coompression will be used.\n");
	printf("\n\t-f ......... force overwrite output file if it exists.\n");
	exit(1);
}

int numReadVarieties=4;

char *ReadName[4]={
	"@ReadVariety0::\n",
	"@ReadVariety1::\n",
	"@ReadVariety2::\n",
	"@ReadVariety3::\n"};
	
char *ReadSequence[4]={
	"ATCAGCATCGCACTCAGCATCAGACATCACGCGGGGCA\n",
	"GGGACTGACTAGCTACGGGGGGGGGCATATTTTTTTTTTTTTTTTTTTTT\n",
	"GGTGGCGTTGGCTAGAGGGAAAAAATCCCATCCCAACCACCGTAGCTGAGGAGGGACCG\n",
	"TTGATGCACACACAGACTACGGGCATCGACATGCATCAGCTACGACATA\n" };
	
char *ReadQualstr[4]={
	"A.da.aca.a,a;.ac,a.c;a.adjadjhfufHHUHh\n",
	"ADADJCNAHCYHDJAAHUDHU_-ADA_D_A_CLA::acpaadkklaAAAA\n",
	"gggcagaidakkcaidadakkcKAHDACNIKAJAKJKA.,;;;a.;;,..,,,,lllko\n",
	"ATGGAAAACBACbabandm,adladlkdadokldaldkhfkjfjakdlq\n" };

char ReadDescr[]="+\n";

int main(int argc, char *argv[]) {
	bool gz=false;
	printf("Arg num: %d\n", argc);
	if(argc<4) Usage();
	for(int i=0; i<argc; i++) {
		printf("argv[%d]: %s\n", i, argv[i]);
	}
	int numReads=atoi(argv[1]);
	if(numReads==0 || numReads>100000000) {
		printf("Number of reads should be between 1 and 100.000.000\n");
		Usage();
	}
	int ReadLen=atoi(argv[2]);
	if(ReadLen<1 || ReadLen>256) {
		printf("Length of reads should be between 1 and 256\n");
		Usage();
	}
	FILE *outputfile=fopen(argv[3],"r");
	if(outputfile!=NULL) {
		fclose(outputfile);
		if( argc!=5) {
			printf("Output file exists, use option -f to overwrite.\n");
			Usage();
		}
		else if(argc==5) {
			if(strcmp(argv[4], "-f")!=0) {
				printf("Output file exists, use option -f to overwrite.\n");
				Usage();
			}
		}
	}

	int l=strlen(argv[3]);
	gzFile gzoutput;
	if(strcmp(&argv[3][l-3], ".gz")==0) {
		gzoutput=gzopen(argv[3],"w");
		if (!gzoutput) {
			printf("ERROR: Failed to open output file: %s", argv[3]);
			exit(255);
		}
		gz=true;
	}
	else {
		outputfile=fopen(argv[3],"w");
		if(!outputfile) {
			printf("ERROR: Failed to open output file: %s", argv[3]);
			exit(255);
		}
		gz=false;
	}
	int variety=0;
	for(int r=0; r<numReads; r++){
		if(gz) {
			gzprintf(gzoutput,"%s", ReadName[variety]);
			gzprintf(gzoutput,"%s", ReadSequence[variety]);
			gzprintf(gzoutput,"%s", ReadDescr);
			gzprintf(gzoutput,"%s", ReadQualstr[variety]);
		}
		else {
			fprintf(outputfile,"%s", ReadName[variety]);
			fprintf(outputfile,"%s", ReadSequence[variety]);
			fprintf(outputfile,"%s", ReadDescr);
			fprintf(outputfile,"%s", ReadQualstr[variety]);
		}
		variety++;
		if(variety>=numReadVarieties) variety=0;
	}
	if(gz) gzclose(gzoutput); else fclose(outputfile);		
}
		
