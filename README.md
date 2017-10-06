# SavvySuite
Suite of tools for analysing off-target reads to find CNVs, homozygous regions, and shared haplotypes

## Compiling
This code requires the htsjdk library and the JAMA matrix maths library. The easiest way to get everything required is to download the GATK Jar. Compiling the code is then done by:
```
javac -cp GenomeAnalysisTK.jar:. *.java
```
For all subsequent examples, it is assumed that both the GenomeAnalysisTK.jar and the compiled *.class files are in the classpath. If you have a large server, it is also sensible to add the -XX:ParallelGCThreads=2 -XX:ConcGCThreads=2 options to java, to prevent it creating too many garbage collection threads.
## Usage
This suite contains three separate tools for analysing off-target reads.
### SavvyCNV
This software analyses the read depth of off-target reads to detect CNVs. It requires a reasonable number of samples sequenced using the same method. (Don't mix samples sequenced using different methods - it won't work well.) The sample data must be provided in aligned BAM files. First, each BAM file must be converted to a coverage summary file, with the following command:
```
java -Xmx1g CoverageBinner sample.bam >sample.coverageBinner
```
This step can be performed in parallel on each sample, and will produce a file approximately 7MB in size. The operation requires very little RAM.

To perform the analysis, the following command should be used:
```
java -Xmx30g SavvyCNV -d <size> *.coverageBinner >cnv_list.csv 2>log_messages.txt
```
The <size> parameter is the size of the chunks that the genome is split into. If you have targeted sequencing with three million reads and about 50% off target reads, then a chunk size of 200,000 is appropriate. It is sensible to process male and female samples separately if CNVs in the X/Y chromosomes are to be detected.
In addition, the following arguments can be provided:
+ -trans <transition probability> - This is the transition probability to use for the Viterbi algorithm. The default is 0.00001. To increase the sensitivity and false positive rate, increase this parameter.
+ -cutoff <noise cutoff> - This is the noise threshold above which a chunk of the genome will be excluded from analysis. The default of 0.25 is probably best in most situations.
+ -g - Switches on the generation of graphs for all samples that have a detected CNV. The generated graphs will be placed in the same directory as the *.coverageBinner file.
+ -a - The same as ```-g```, but produces a graph for every sample.
+ -mosaic - Switches the software into mosaic mode. Normally, the state probability calculations assume that the relative dosage is either <=0.5, 1, or >=1.5. Dosage levels must cross the mid-point between 1 and 0.5 or 1.5 before they become evidence of a CNV. This increases sensitivity and specificity at the cost of being able to detect mosaic CNVs. With this switch, mosaics can be detected. The size parameter will need to be increased, and small CNVs will not be detected as effectively.

The output cnv_list.csv contains a tab-separated list of detected CNVs. The columns in the output are:
1. Chromosome
2. CNV start position
3. CNV end position
4. Deletion/duplication
5. Number of genome chunks used as evidence for the CNV
6. Width of CNV in chunks (including gaps containing noisy chunks, such as centromeres)
7. The phred score of the CNV
8. The phred score divided by the width of the CNV in chunks. Most valid CNVs have a value more than ten in this column
9. The relative dosage in the CNV, so 0.5 for a heterozygous deletion etc.
10. The filename containing the sample summary.

The log_messages.txt file contains log messages, and also a summary of each sample. Each sample has a line containing the following columns:
1. Noisyness of the sample before excluding known CNVs
2. Noisyness of the sample after excluding known CNVs (the CNV calling is performed twice, as the calling depends on the noisyness of the sample). This figure should be below 0.2 for good results.
3. Number of deletions found in the sample
4. Number of duplications found in the sample
5. The filename containing the sample summary.

###SavvyHomozygosity
This software analyses the off-target reads to determine homozygous regions of the genome for a single sample. It can be run as follows:
```
java -Xmx5g SavvyHomozygosity whole_genome.vcf sample.bam >sample.bed
```
The whole_genome.vcf file must contain whole-genome sequencing data with many (a few hundred would be appropriate) samples, in order to calculate population linkage disequilibrium. The sample.bam contains the sample targeted sequencing data. The software produces a BED file containing homozygous regions.
