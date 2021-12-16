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
+ -headers - Makes SavvyCNV output a header at the top of the CNV list it generates.

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

#### Controlling sensitivity and specificity
The sensitivity and specificity of SavvyCNV can be controlled by adjusting the command line arguments. Sensitivity and specificity can be traded off against each other - that is, you can adjust the software to be highly sensitive with the consequence of having low specificity, or to be highly specific with the consequence of being less sensitive. The command line arguments that allow this to be done are as follows:
1. -d (size) - the chunk size parameter primarily determines how much of the genome is analysed, and how noisy the resulting data is. If a very large chunk size is chosen (for instance 2Mbp) then each chunk will have a large amount of data available, which decreases the noise present in each chunk. The trade-off is that CNVs that are smaller than the chunk size (or more realistically a little larger) cannot easily be detected, so a large chunk size favours specificity over sensitivity. A smaller chunk size (for instance 200kbp) will have less data in each chunk, increasing the noise in each chunk, but increasing the chances of smaller CNVs being detected. If the chunk size is decreased further, then some chunks will have insufficient data to analyse, and are rejected from the analysis. This can be useful if the data is from targeted or exome sequencing, as targeted regions of the genome will have a higher read depth, and chunks covering these targeted regions will remain in the analysis. For exome or targeted sequencing with a mean targeted read depth of 100 or more, a chunk size of 200bp or 400bp will be appropriate, and will find on-target CNVs only - off-target regions will be rejected from the analysis for insufficient data. The noise level for each sample is reported in the statistics output from SavvyCNV. For reliable detection of CNVs, the noise level should be less than around 0.15 - this ensures that the majority of chunks have a read depth anomaly close to zero enough to distinguish a true deletion or duplication. If SavvyCNV is run in mosaic mode, then the noise level needs to be less, which can be achieved by increasing the chunk size. The noise level should be less than the difference in dosage that you would like to detect.
2. -minReads (number) - the minimum number of reads determines which chunks are excluded from the analysis. The default is 20, as that appears to be the number of reads where a deletion/duplication can be distinguished from the noise. The theoretical best possible noise level for a chunk is (according to the Poisson distribution) 1/sqrt(read_count). For 20 reads this is 0.22, so the majority of chunks will have a normalised read depth between 0.78 and 1.22. SavvyCNV models the error in each chunk in each sample using the overall noise level of each sample combined with the actual spread of normalised read depth for all the samples for that chunk, and this error is often slightly larger than the theoretical minimum from the Poisson distribution. Increasing this parameter from the default will cause SavvyCNV to reject more chunks that have insufficient data, increasing the specificity, with the consequence that CNVs covered by the rejected chunks cannot be detected.
3. -trans (probability) - This is the main control of sensitivity and specificity. SavvyCNV uses an HMM to detect CNVs, and the transition probability is the modelled probability that one part of the genome transitions over to a CNV or back to normal dosage. An HMM finds the arrangement with the greatest probability, so chunks with a sufficient read depth anomaly must be found to override the transition probability twice - once to enter the CNV and again to exit the CNV. If the transition probability is set very low (for instance 0.0000000001) then transitions are discouraged and significant evidence must be found in the data to force a CNV detection, therefore sensitivity is reduced and specificity is increased. A higher transition probability (for instance 0.1) allows transitions to occur with little evidence, increasing sensitivity and reducing specificity. It should also be noted that a high transition probability increases the chances that a CNV will be split into multiple parts by a chunk that erroneously has a normal read depth.
4. -minProb (phred probability) - This limits how much evidence a single chunk can contribute towards the analysis. If this is set higher than the square of the transition probability, then it will become impossible for a CNV to be detected from the evidence in a single chunk.

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

### ValidateGenomeCnvsUsingAlleleBalance
This software analyses a VCF file, and determines whether the read balance of the variants present agrees with a list of CNVs provided to it. For example, you should expect to not find any heterozygous variants inside a deletion, and you should expect to see a 1:2 read balance for variants inside a duplication. It can be run as follows:
```
java -Xmx5g ValidateGenomeCnvsUsingAlleleBalance [-columns chr start end sample] cnv_list linkage_data vcf_files...
```
The options are:
+ -columns is an optional argument. By default, the software will read the cnv_list file, and take the chromosome, start, end, and sample name from columns 1, 2, 3, and 4. The -columns option allows these values to be taken from a different column instead.
+ The cnv_list file is a text file containing tab-separated values. The software needs at least the chromosome, start, end, and sample name for each CNV that is to be validated. If a sample name is "all" in the list, then the likelihood of that CNV existing is tested in all the samples in all the provided VCF files. The file should have a one-line header.
+ The linkage_data file is used to filter the variants in the VCF file to only allow common variants that have a good allele balance in most samples. The keyword "none" can be used instead, and the software will use all available variants without filtering.
+ The list of vcf files will be iterated through to find data for all the samples.
So for instance, if the cnv_list file contains the following (with tabs separating columns):
```
chr start end sample
6 144200000 144400000 all
```
then you can run:
```
java -Xmx5g ValidateGenomeCnvsUsingAlleleBalance cnv_list none vcf_files...
```
and it will evaluate the likelihood of a 6q24 CNV existing in all the samples found in all the VCF files, from the allele balance of the variants.

The output contains multiple rows, one for each match between a cnv_list line and a sample in a VCF file. Each line starts with a copy of the line from the cnv_list, followed by these extra columns:
1. Sample - the sample name that the data is taken from. This will usually be the same as the data from the cnv_list earlier in the line, but it is useful when the cnv_list specifies "all" samples.
2. Variants - the total number of variants found in the area of the alleged CNV.
3. Heterozygous - the number of heterozygous variants found in the area of the alleged CNV.
4. Reads - the total number of reads found in the heterozygous reads.
5. 25% - the signal strength for the heterozygous variants being in a 1:3 allele ratio.
6. 33% - the signal strength for the heterozygous variants being in a 1:2 allele ratio.
7. 50% - the signal strength for the heterozygous variants being in a 1:1 allele ratio.
8. HetProportion - the proportion of the variants that are heterozygous.
9. PeakPosition - the allele ratio with the highest signal strength. This is the best estimate of the allele ratio of variants in the area of the alleged CNV.
10. PeakHeight - the highest signal strength value.
A region may or may not have sufficient data to evaluate a CNV. This is the purpose of the Variants, Heterozygous, and Reads columns. A region with sufficient data should have at least several heterozygous variants, and ideally at least 1000 reads.

The data can show the following kinds of results:
1. Homozygous area or a deletion (or a normal region on the X/Y chromosome in males). This is indicated by the HetProportion column having a very low value. It is not possible to distinguish between a homozygous region and a deleted (monosomy) region using this tool. For this, use SavvyCNV or another read-depth CNV detection method instead.
2. Normal (disomy) region (or duplication on the X/Y chromosome in males). This is indicated by the HetProportion value being around 0.6 (could vary from 0.3 to 0.8), the 50% value being higher than the 33% or 25% value, and the PeakPosition value being close to 0.5.
3. Duplicated (trisomy) region. This is indicated by the 33% value being higher than the 50% or 25% value, and the PeakPosition value being close to 0.33. The HetProportion value may vary from 0.3 to 0.9.
4. Double-duplicated region. This is indicated by the 25% value being higher than the 50% or 33% value, and the PeakPosition value being close to 0.25. The HetProportion value may vary from 0.3 to 0.9.
5. Contaminated sample. This is indicated by a very high HetProportion, and by a 25% value higher than the 50% or 33% value. The main distinguishing feature of sample contamination is that the whole genome will have these values, rather than just a small alleged CNV. Contamination is better detected by VerifyBamID and analysed using SavvyContaminationFinder instead of this software.
6. Mosaic UPD. This is indicated by a very low 50% value. At different levels of mosaicism, this can be indistinguishable from a trisomy or double-duplicated region using this tool. The best way to distinguish this is to look at the read depth using SavvyCNV or another read-depth CNV detection method.

### SavvyContaminationFinder
This software analyses a VCF file containing variants from multiple samples, and determines the level of cross-contamination between samples. This can be used to identify the source of contamination of a sample in order to analyse the root cause of the contamination. See also SavvyContaminationRepairer below, which can correct variants that are incorrectly called because of the cross-contamination.

The software uses the allele depth (AD) information in the VCF file, so multi-sample VCF files produced by GATK GenotypeGVCFs or CombineVariants are suitable. It uses the allele depths for the sample being analysed and the genotypes (derived from allele depths) for the remaining samples in the file to build a model of the mixture of sample contributions to the sample being analysed. It then finds the mixture that best matches the allele depths, using a minimisation algorithm based on preconditioned conjugate gradients. The variants are then corrected based on the calculated contribution from the uncontaminated source sample, and the minimisation repeated, until it converges on a solution. The arguments for the software are:
1. The VCF file containing variants from multiple samples.
2. The name of the sample to analyse, or "all" to analyse all of the samples in the file.
3. (Optional, default: no limit) The maximum mean read depth of a variant to be included in the analysis. For whole genome sequencing, set this to approximately double the mean read depth of your sequencing, and this will filter out variants in collapsed repeat regions and other non-diploid regions. For targeted/exome sequencing, set this to a very high number, as targeted regions may have a high read depth.
4. (Optional, default: 1) The number of CPU threads to process the data with. This is only useful if argument 2 is "all".

The output contains multiple rows, describing the content of each sample in the VCF file that was analysed. The columns are:
1. The sample being described.
2. The sample contributing reads to the sequence data.
3. The proportion of the sample being described that originated from the sample contributing reads.

Given uncontaminated samples, you should expect each described sample to have a near-zero proportion value for every contributing sample, except of course itself, which should have a proportion value near 1. For a contaminated sample, the contaminating sample will have a proportion value representing the amount of that sample that has been mixed into the contaminated sample, and the contribution from the contaminated sample itself will be less than 1 accordingly.

### SavvyContaminationRepairer
This software corrects the variants in a VCF file taking account of cross-contamination as detected by SavvyContaminationFinder. The genotypes (GT), genotype quality (GQ) and phred-scale likelihood (PL) fields for each sample are updated with new values. The old GT and PL fields are copied to new OLDGT and OLDPL fields. The software reads the contamination values from the standard input, in the same format as produced by SavvyContaminationFinder. In a contaminated sample, the GQ and PL values will naturally be lower than an uncontaminated sample, showing the lower confidence in the original genotype. The software takes one or two arguments. They are:
1. The VCF file containing the contaminated sample data.
2. (Optional) The output VCF file. If this argument is not provided, then the VCF file will be output on the standard output, and an index cannot be created.

A typical workflow for processing contaminated samples is to run the following:
```
java SavvyContaminationFinder contaminated.vcf all >contamination.txt
java SavvyContaminationRepairer contaminated.vcf decontaminated.vcf <contamination.txt
```
This will produce the contamination.txt file which describes the contamination source for each sample, which should be checked to make sure the contamination originates from a sample that is present in the VCF file. It will then produce decontaminated.vcf which has most of the genotypes corrected.

Contamination tends to cause homozygous reference (i.e. no variant present) locations to have a false positive heterozygous variant call, if the contaminant has that variant. It also tends to cause homozygous variant locations to be mis-called as heterozygous, if the contaminant is not also homozygous variant. If a homozygous variant is rare, then this is increasingly likely, meaning that contamination causes the most likely disease-causing variants to be most likely to be mis-called. SavvyContaminationRepairer is effective at correcting both of these faults, correcting 98.5 to 99% of variants in GIAB not-difficult regions at contamination levels between 10% and 20% in tests. Heterozygous variants are less likely to be mis-called due to contamination.
                                                                                       
### CoverageOffTarget
This software analyses a set of CoverageBinner files, and determines how many off-target reads there are in them. Firstly, the software scans through the files, and identifies all the 200bp genome chunks that have fewer than a threshold number of reads on average (default 5). Then it prints statistics for each input file. The arguments are:
1. -threshold <number> to set the threshold number of reads. For each 200bp chunk of the genome, the average number of reads (across all the input samples) is calculated, and a chunk is counted as on-target if it is above this threshold, and off-target if it is below. The default is 5, and it can be a floating point number.
2. -readLength <number> to tell the software the read length of the samples. This enables the software to calculate the mean read depth in off-target regions.

All remaining arguments are taken to be input files to process. The software will first print the number of chunks that are on-target, and the proportion of the genome that is, and then for each input sample the following columns:
1. The name of the CoverageBinner file that contains the data.
2. The total number of reads in the file.
3. The number of reads that are in on-target chunks.
4. The proportion of reads that are in on-target chunks.
5. The number of reads that are off-target.
6. The mean read depth in off-target regions. If the -readLength argument is not supplied, then this column is not produced.

### FindLargeInsertSizes
This software searches for read pairs in a BAM file that have an abnormal insert size. If multiple abnormal read pairs in the same location agree, this may indicate the presence of a structural variant, such as a CNV, inversion, or translocation. This detection method is reasonably good at detecting small CNVs, although it will not be able to detect CNVs for which there are no abnormal read pairs. This can occur if targeted sequence data is analysed (and therefore the CNV boundaries are not covered), or in WGS data if the CNV boundaries are in two regions of the genome with similar sequences. Small CNVs are hard to detect using read depth methods or variant allele fraction methods, and this read pair method will provide a better sensitivity and specificity. For large CNVs, a read-depth method or variant allele fraction method is able to collect much more evidence than boundary-based methods such as a read pair method like this, or split read methods. Therefore, this software produces many false positives of large CNVs and translocations due to incorrect mapping of reads in repetitive regions in the genome, and these should be interpreted with caution.

This software analyses multiple BAM files. Events will be detected using the first BAM file specified on the command line, and then the statistics for that event are calculated for each of the BAM files, to work out whether the same event is present in each sample. This allows for a more effective method of determining whether events are de novo in a sample for instance, if the proband is specified as the first sample on the command line, and the parents are specified in addition.

The command line arguments are:
1. -limit (chromosome) (start) (end) - This will limit the analysis to just the specified region of the genome.
2. -stddev - This will add additional columns to the output, containing information on the standard deviations of the read pair insert sizes.
3. -maxInsert (number) - This changes the maximum insert size of a read pair for the read pair to be deemed "normal". Read pairs with an insert size greater than this are analysed. The default is 1000bp, but if you know the distribution of insert sizes for your sequence data, you can reduce this to match the very top of the distribution.
4. -minMq (number) - This sets the minimum mapping quality of reads. All reads with a mapping quality below this number are discarded. The default is 20.
5. All other arguments are treated as BAM file names. Events are detected in the first BAM file and confirmed in all BAM files.

The software produces one line of output for each event found. Each line is tab-separated with the following columns:
1. The chromosome:position location of the event as it was found.
2. "F" or "R" to indicate whether the reads participating in the event at this location are paired in the forward or reverse direction.
3. The chromosome:position location of the read pairs of the reads in the event - this is an approximation of the position of "other end" of the event.
4. "F" or "R" to indicate whether the read pairs at the other end of the event are paired in the forward or reverse direction.
5. The number of read pairs that agree with the event in the first BAM file.
6. The distance between the two ends of the event, of "?" if the two ends are in different chromosomes.
7. A description of the event. The possible values are:
  a. Deletion_left_edge
  b. Deletion_right_edge
  c. Duplication_left_edge
  d. Duplication_right_edge
  e. Inversion_inner_left
  f. Inversion_inner_right
  g. Inversion_outer_left
  h. Inversion_outer_right
  i. Translocation
8. For each BAM file on the command line, two columns are added:
  a. The number of read pairs found that agree with the event
  b. The total number of reads in the vicinity of the event
