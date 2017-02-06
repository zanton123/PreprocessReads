# PreprocessReads
A CUDA based high throughput NGS short read preprocessing tool. Massive parallel sequencing generates million of short sequence reads that have to be processed and aligned to reference genomes. Ever greater amount of data creates a need for efficient preprocessign tools to maximize the use of the information and to avoid alignment and interpretation mistakes to the best possible extent. Although read preprocessing is not complex it has been less developed than NGS read aligners. Owing to its position as the initial step of analysis pipelines, preprocessing offers considerable opportunities for reducing noise from sequencing data and for controlling data formats and the analysis strategy. In addition *throughtful* preprocessing can contribute to efforts for futher reducing alignment errors and runtime even in already optimized pipelines.

PreprocessReads assists with trimming reads by reported base quality, trimming of adaptor sequences on the 3'-end, and the recovery of short insert sequences such as gRNA sequences from CRISPR/Cas library screens. The tool draws GPU acceleration and can process 100 milion 100 nucleotide Solexa reads in approximately 5 minutes on a NVIDIA Keppler or Maxwell accelerator.

PreprocessReads might be of interest for educational purposes as it adopts a simple and extensible source format. At this time the source does not integrate with NVIDIA's nvBIO library but is constructed as an easy to understand single CUDA C++ source file that draws on the CUDA runtime API for parallel processing of data.


**Usage: PreprocessReads [algorithm selection] -i [fastq or fastq.gz input file] -o [output file] [options]**

Program Options:

	-i, -I ....... specify the path for the input file [required]
	-o, -O ....... specify the path for the output file [required]
	-len, -Len ... specify the minimum length for reads written to the output file
	-f ........... over write output file if it exists
	-gz .......... write output to file in compressed format (.gz)
	-pigz ........ write output to file in compressed format (.gz) using multithreaded compression
	-v, -V ....... verbose, prints details on processing to screen
	-a, -A ....... adaptor sequence, can contain Ns that will match any base in the read
	-mm, -MM ..... number of mismatches allowed in adaptor search (not including Ns)
	-Phred ....... sets the base of the qualities PHRED values [default 33], use AUTO for autodetection
	-Qmin ........ sets the minimum acceptable quality score for a good base [default 20]
	-Qualab ...... sets the minimum average base quality score for remainder of the read after a bad base [default 25]
	-dev, -Dev ... specify CUDA device to use (can be used repeatedly for selecting multiple GPUs,
	               if not specified an automatic search for GPUs with compute capability 3.5 and 6 GB memory is performed)

Algorithm selection:

	/TA .......... select adaptor trimming only (can be also used BEFORE option /gRNA to avoid quality trimming)
	/TQ .......... select quality trimming only
	/TAQ ......... select quality and adaptor trimming [default]
	/gRNA ........ select gRNA recovery, followed by start relative to start of adaptor sequence and length of gRNA to be retrieved
	               ( /gRNA 12 20 will retrieve 20 nucleotide read sequence starting 12 positions after the 5' end match of adaptor)

Programmer options variables:

	-G ........... set MAX_READS_ON_GPU variable
	-R ........... set MAX_READ_LENGTH variable

example usage:

	PreprocessReads -dev 0 -dev 1 -i infile.fastq -o outfilename.fastq -len 26 -a ATCTCGTATGCCGTCTTCTGCTT -mm 2



###How to build PreprocessReads:

PreprocessReads depends on the core cuda and c libraries as well as libz (**install zlib1g-dev**). To write compressed output files gzip should be in the path, else we fall back to uncompressed output. For using the multithreaded compression option **-pigz** you need also the package pigz installed and in the path. The source contains the eclipse project files for building with the NVIDIA CUDA IDE on Linux. If you would want to build from the terminal the following commands can be used:

```
cd
git clone https://github.com/zanton123/PreprocessReads.git
cd PreprocessReads
nvcc -lz -o PreprocessReads PreprocessReads.cu
```

This will build an executable PreprocessReads for your the CUDA card on your machine. If your development and production machines are different consult the CUDA toolkit documentation on the *-arch* and *-gencode* parameters. Also note that the source code presently expects at _**minimum 6 gigabyte GPU memory and compute capability of 3.5**_. PreprocessReads has been tested on the following NVIDIA GPUs:`
```
GTX Titan
GTX Titan Z
GTX Titan X (Maxwell)
```
