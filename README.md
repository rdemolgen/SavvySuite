# SavvySuite
Suite of tools for analysing off-target reads to find CNVs, homozygous regions, and shared haplotypes.

This software was written by Matthew Wakeling at the University of Exeter, and was presented at the 2017 ASHG meeting in Orlando, Florida.

To cite this software, please use the following references. For SavvyCNV, use:
> https://www.biorxiv.org/content/10.1101/617605v1

For all other parts of SavvySuite, cite the ASHG conference presentation:
> Wakeling, MN, De Franco E, Hattersley AT, Ellard S. Making the most of targeted sequencing: detecting CNVs and homozygous regions using offtarget reads with SavvyCNV. 67th Annual Meeting of the American Society of Human Genetics. Orlando, FL; 17–21 October 2017.

and if desired, cite this web site too:
> Wakeling MN. SavvySuite. 2018. https://github.com/rdemolgen/SavvySuite.

## Running Java
This code requires the htsjdk library and the JAMA matrix maths library. The easiest way to get everything required is to download the GATK Jar. All operations require this GATK jar and the SavvySuite java to be in the Java classpath, in order for Java to find it.

This can be done in two ways. The first option is to set the CLASSPATH environment variable:
```
export CLASSPATH=/path/to/GenomeAnalysisTK.jar:/path/to/SavvySuite/directory
```
The ":" character separates the two parts of this path, to specify that code can be found in the two places. The second option is to specify the "-cp" option every time you run java, like this:
```
java -cp /path/to/GenomeAnalysisTK.jar:/path/to/SavvySuite/directory blah blah blah
```
For all subsequent code fragments, where "java" or "javac" is specified, it is assumed that the classpath is correctly configured as specified above, either by adding the "-cp" option or using the CLASSPATH environment variable.

If you have a large server, it is also sensible to add the -XX:ParallelGCThreads=2 -XX:ConcGCThreads=2 options to java, to prevent it creating too many garbage collection threads. This is a small performance enhancement, and is done like this:
```
java -XX:ParallelGCThreads=2 -XX:ConcGCThreads=2 blah blah blah
```

## Compiling
Compiling the code is then done by:
```
javac *.java
```
in the SavvySuite directory.

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
java -Xmx30g SavvyCNV -d (size) *.coverageBinner >cnv_list.csv 2>log_messages.txt
```
The (size) parameter is the size of the chunks that the genome is split into. If you have targeted sequencing with three million reads and about 50% off target reads, then a chunk size of 200,000 is appropriate. It is sensible to process male and female samples separately if CNVs in the X/Y chromosomes are to be detected. If you have a lot of samples and SavvyCNV takes a long time, then please use the SelectControlSamples software described below.
In addition, the following arguments can be provided:
+ -trans (transition probability) - This is the transition probability to use for the Viterbi algorithm. The default is 0.00001. To increase the sensitivity and false positive rate, increase this parameter.
+ -cutoff (noise cutoff) - This is the noise threshold above which a chunk of the genome will be excluded from the calculation of how noisy a sample is. The default of 0.25 is probably best in most situations.
+ -g - Switches on the generation of graphs for all samples that have a detected CNV. The generated graphs will be placed in the same directory as the *.coverageBinner file.
+ -a - The same as ```-g```, but produces a graph for every sample.
+ -cytoBands (file) - The location of a cytoBands file for the genome reference, for example downloadable from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz - this will be plotted in the graph behind the data.
+ -mosaic - Switches the software into mosaic mode. Normally, the state probability calculations assume that the relative dosage is either <=0.5, 1, or >=1.5. Dosage levels must cross the mid-point between 1 and 0.5 or 1.5 before they become evidence of a CNV. This increases sensitivity and specificity at the cost of being able to detect mosaic CNVs. With this switch, mosaics can be detected. The size parameter will need to be increased, and small CNVs will not be detected as effectively.
+ -sv (number) - This changes the number of singular vectors that are removed for noise reduction. The default is 5. This must be less than the number of samples.
+ -minReads (number) - This sets the minimum number of reads that a genome chunk must have on average across the samples in order to be analysed. The default is 20.
+ -minProb (number) - This sets the greatest (lowest) probability that a single chunk can contribute to a CNV. The number is a phred score. This is the largest quality score that a very small CNV can have.
+ -chr (chromosome name) - This limits the analysis to just one chromosome - all reads outside that chromosome will be discarded.
+ -case - All samples listed after this option (until a -control option) are marked as case samples, and CNV calling will be performed on them. This is the default.
+ -control - All samples listed after this option (until a -case option) are marked as control samples, and CNV calling will not be performed on them.
+ -data - An additional file will be created for all case samples, which contains the raw data that SavvyCNV used to call CNVs. The output file is named after the input CoverageBinner file, with ".<bin_size>.data" appended. The file contains 8 columns, which are:
  1. Chromosome
  2. Bin start position
  3. Bin end position
  4. Noise-corrected normalised read depth
  5. Estimated error in normalised read depth
  6. Uncorrected normalised read depth
  7. Phred score probability of a deletion in this bin
  8. Phred score probability of a duplication in this bin

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
10. The filename of the input CoverageBinner file for this sample

The log_messages.txt file contains log messages, and also a summary of each sample. Each sample has a line containing the following columns:
1. Noisyness of the sample before excluding known CNVs
2. Noisyness of the sample after excluding known CNVs (the CNV calling is performed twice, as the calling depends on the noisyness of the sample). This figure should be below 0.2 for good results. If too many of the samples have a high noisyness, then increase the chunk size parameter (-d <size>).
3. Number of deletions found in the sample
4. Number of duplications found in the sample
5. The number of reads in the sample.
6. The filename containing the sample summary.

### SelectControlSamples
This software selects a subset of samples from a larger control pool that match a given sample or set of samples best. This is useful if you have a lot of control samples, but SavvyCNV takes a long time to run with all of them. SavvyCNV uses Singular Vector Decomposition (SVD), which takes time approximately proportional to the number of samples to the power of 2.4, so if for example 150 samples takes 30 minutes, then 300 samples will take a bit more than two and a half hours. SelectControlSamples prepares a statistical summary of a large collection of samples, using SVD on a subset of those samples to identify the read depth patterns, then extending it to the rest of the samples. This summary can then be interrogated to identify the most similar set of samples to a sample that you wish to call CNVs on.

The software has two modes of operation. The first mode takes a list of CoverageBinner files and produces a summary file, having performed SVD. This operation takes considerable time, but only needs to be performed once. The second mode of operation takes the summary file and a list of CoverageBinner files, and identifies the CoverageBinner files that were used to create the summary that are most similar. This operation is fast. Samples from the summary file that are also specified in this operation are automatically excluded from the output, as it is not sensible to use a sample as its own control when searching for CNVs.

The software should be run as follows:
```
java -Xmx30g SelectControlSamples -d (size) *.coverageBinner >summary_file
```
to create the summary file. The (size) parameter is the size of the chunks that the genome is split into, in the same way as with SavvyCNV. The following additional options are available:
+ -minReads (number) - This sets the minimum number of reads that a genome chunk must have on average across the samples in order to be analysed. The default is 20.
+ -subset (number) - This sets the number of samples that the SVD will be performed on, with a default of 50. Set this according to how much time you want the software to take - 150 is not unreasonable. The software will take time proportional to this parameter to the power of 2.4, plus time proportional to the total number of samples on the command line.
+ -chr (chromosome name) - This limits the analysis to just one chromosome - all reads outside that chromosome will be discarded.
+ -cross - Instead of producing a summary file, the software will produce a full table of the difference metrics between all samples. The output file can be visualised in gnuplot by running "plot 'file' with image".
+ -svs - Instead of producing a summary file the software will print the singular vectors of all the samples. One row is printed for each sample, and the multiple values from the singular vector (as many as specified in the -subset option) are in multiple columns, ordered from the largest singular value downwards.

To use the summary file and produce a list of best-matching samples, the software should be run as follows:
```
java -Xmx30g SelectControlSamples -summary summary_file *.coverageBinner
```
You should only specify a small number of samples (or preferably one) on this command, as the samples from the summary file are selected based on the sum of the differences to the samples on the command line. The selected samples will match a single specified sample better than they will match multiple specified samples. The block size, minReads, and chromosome options are saved in the summary file so they do not need to be specified here. The following additional options are available:
+ -subset (number) - This sets the number of samples that this software will print out.
+ -stats - Without this option, the output of this command is a list of CoverageBinner files that went into the summary file. With this option, a second column separated by a tab character is added to the output, which is the sum distance to the specified samples.
+ -cross - Instead of selecting matching samples, all other arguments are ignored, and the software will produce a full table of the difference metrics between all samples in the summary file.
+ -svs - Instead of selecting matching samples, all other arguments are ignored, and the software will print the singular vectors of all the samples in the summary file. One row is printed for each sample, and the multiple values from the singular vector (as many as specified in the -subset option) are in multiple columns, ordered from the largest singular value downwards.

SavvyCnv can therefore be run on a sample by running:
```
java -Xmx30g SavvyCNV -d (size) case_sample.coverageBinner -control `java -Xmx30g SelectControlSamples -summary summary_file` >cnv_list.csv 2>log_messages.txt
```

### SavvyCNVJointCaller
This software performs joint calling of CNVs in multiple samples. It is not suited to calling CNVs in large numbers of samples, but is intended to be used with multiple members of the same family. The algorithm favours CNVs starting and ending in the same location in multiple samples. This should reduce the incidence of (for example) false positive de novo CNVs being detected when CNVs are called independently in two parents and a child, if the CNV detected in the child is falsely detected as slightly larger than the CNV inherited from a parent.

This software uses the output from SavvyCNV when it is given the "-data" option. SavvyCNV creates these files named after the input CoverageBinner files, with ".<bin_size>.data" appended.

To perform the analysis, the following command should be used:
```
java -Xmx30g SavvyCNVJointCaller *.coverageBinner.<bin_size>.data >cnv_list.csv 2>log_messages.txt
```
Extra options are:
+ -trans (transition probability) - This is the transition probability to use for the Viterbi algorithm. The default is 0.00001. To increase the sensitivity and false positive rate, increase this parameter.
+ -mosaic - Switches the software into mosaic mode. Normally, the state probability calculations assume that the relative dosage is either <=0.5, 1, or >=1.5. Dosage levels must cross the mid-point between 1 and 0.5 or 1.5 before they become evidence of a CNV. This increases sensitivity and specificity at the cost of being able to detect mosaic CNVs. With this switch, mosaics can be detected. The size parameter will need to be increased, and small CNVs will not be detected as effectively.
+ -minProb (number) - This sets the greatest (lowest) probability that a single chunk can contribute to a CNV. The number is a phred score. This is the largest quality score that a very small CNV can have.

The output cnv_list.csv contains a tab-separated list of detected CNVs. The columns in the output are:
1. Chromosome
2. CNV start position
3. CNV end position
4. Deletion/duplication/normal
5. Number of genome chunks used as evidence for the CNV
6. The phred score for a deletion for this sample (positive means a deletion is likely, negative means it is unlikely)
7. The phred score for a duplication for this sample
8. The phred score divided by the width of the CNV in chunks. Most valid CNVs have a value more than ten in this column
9. The relative dosage for this sample in this region, so 0.5 for a heterozygous deletion etc.
10. The filename of the input data file for this sample

### PrepareLinkageData
This software is used to pre-process linkage disequilibrium data, to get it into a format that can be read quickly by SavvyHomozygosity and SavvySharedHaplotypes. It requires a vcf file containing whole genome genotype data for many samples (a few hundred would be appropriate). It can be run as follows:
```
java -Xmx5g PrepareLinkageData whole_genome.vcf >linkage_data
```
If you have enough genomes sequenced from a similar ethnicity to the samples analysed using targeted or exome sequencing, then use the above. If not, then you can download a pre-prepared data file from https://mega.nz/#!lppUyCpZ!LM2CrnCUZFl3O-WeARQtWdMPVxnV93vQ4f0ZDdSUaak which was converted from the 1000 genomes data.

### SavvyHomozygosity
This software analyses the off-target reads to determine homozygous regions of the genome for a single sample. It can be run as follows:
```
java -Xmx5g SavvyHomozygosity linkage_data sample.bam >sample.bed
```
The sample.bam contains the sample targeted sequencing data. The software produces a BED file containing homozygous regions. The analysis requires at least a million off-target reads, and preferably more, so a sample with three million reads and about 50% off-target is about right. The software searches for read pairs, so the amount of detail that can be extracted (and consequently the time taken) is proportional to the square of the number of off-target reads. Samples with too few reads do not produce a valid result.

### SavvySharedHaplotypes
This software analyses the off-target reads from two samples to determine areas with homozygous shared haplotypes, using the same methodology as SavvyHomozygosity. It can be run as follows:
```
java -Xmx5g SavvySharedHomozygosity linkage_data sample1.bam sample2.bam >shared.bed
```

### SavvyVcfHomozygosity
This software analyses a VCF file to determine homozygous regions of the genome. It can be run as follows:
```
java -Xmx5g SavvyVcfHomozygosity linkage_data input.vcf (options) (samplenames) >sample.bed
```
If a single samplename is used with no options, then the output will be a BED file containing homozygous regions for the sample. If multiple samplenames are used, then only regions that are homozygous in all samples are output. Adding the "-p" option alters this behaviour, so that regions where all the named samples share a haplotype are output.

The linkage_data file is used to filter the variants in the VCF file to only allow variants that are close to a Hardy-Weinberg equilibrium, to reduce the number of false heterozygous variants in homozygous regions. The variants are then analysed using a multi-resolution hidden markov model to find regions of at least 16 variants where the proportion of heterozygous variants is no more than 10%.
