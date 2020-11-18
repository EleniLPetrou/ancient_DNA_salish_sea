# Download and index Atlantic herring genome

I used the chromosome-level assembly of the Atlantic herring genome (Pettersson et al. 2019) as a reference genome for this project. Here is the NCBI link to the genome : https://www.ncbi.nlm.nih.gov/genome/genomes/15477
Assembly ID : GCA_900700415.1

After I downloaded the genome from NCBI, I indexed it using *bowtie2*. Indexing compresses the size of the file and makes queries fast.

## Index Atlantic herring genome

``` bash

BASEDIR=/media/ubuntu/Herring_aDNA/atlantic_herring_genome # Path to directory with genome.
REF=GCA_900700415.1_Ch_v2.0.2_genomic.fna # Genome fasta file name
BASENAME=GCA_900700415 # Write bt2 data to files with this dir/basename


bowtie2-build $BASEDIR'/'$REF $BASEDIR'/'$BASENAME

```


# Align ancient samples to the genome

## Explanation of terms:

*bowtie2 -q -x <bt2-idx> -U <r> -S <sam>*
  
-q query input files are in fastq format

-x <bt2-idx> Indexed "reference genome" filename prefix (minus trailing .X.bt2)

-U <r> Files with unpaired reads.

-S <sam> File for SAM output (default: stdout)
you can use 2> to redirect stdout to a file (bowtie writes the summary log files to stdout)



``` bash

#Directories and files

BASEDIR=/media/ubuntu/Herring_aDNA/hybridization_capture/ancient_samples
SAMPLELIST=$BASEDIR/sample_lists/sample_list.txt # Path to a text file with list of prefixes of the fastq files, separated by newline (so, the file name with no extension). 
FASTQDIR=$BASEDIR'/'trimmed_fastq # Path to folder containing fastq files.
GENOMEDIR=/media/ubuntu/Herring_aDNA/atlantic_herring_genome # Path to folder with genome.
GENOME=GCA_900700415 #genome prefix
OUTPUTDIR=$BASEDIR'/'sam # path to folder with sam files (output)


# Command

# Loop over each sample fastq file and align it to genome, then output a sam file

for SAMPLEFILE in `cat $SAMPLELIST`
do
  bowtie2 -q -x $GENOMEDIR'/'$GENOME -U $FASTQDIR'/'$SAMPLEFILE.fastq -S $OUTPUTDIR'/'$SAMPLEFILE.sam
done

```


# Convert sam files to bam format, filter, remove PCR duplicates, and index bam files
Next, I filtered the bam files. I removed any sequences that were shorter than 30 nucleotides long, removed sequences that had a mapping quality below 30, converted the files into bam format, and then sorted and indexed the bam files. All this was done using using *samtools*

- samtools view: prints all alignments in the specified input alignment file (in SAM, BAM,  or  CRAM format) to standard output. Can use the samtools *view* command to convert a sam file to bam file

  -S: input file is in sam format

  -b: output file should be in bam format

  -q <integer>: discards reads whose mapping quality is below this number

  -m <integer>: only outputs alignments with the number of bases greater than or equal to the integer specified

  -h: Include the header in the output.

- samtools sort : Sort alignments by leftmost coordinates

- samtools rmdup : Remove potential PCR duplicates: if multiple read pairs have identical external coordinates, only retain the pair with highest mapping quality.

  -s : Remove duplicates for single-end reads. By default, the command works for paired-end reads only.

- samtools index : Index a coordinate-sorted BAM or CRAM file for fast random access. Index the bam files to quickly extract alignments overlapping particular genomic regions.This index is needed when region arguments are used to limit samtools view and similar commands to particular regions of interest.


``` bash
# Directories and files
BASEDIR=/media/ubuntu/Herring_aDNA/hybridization_capture/ancient_samples
SAMPLELIST=$BASEDIR/sample_lists/sample_list.txt # Path to a text file with list of prefixes of the fastq files, separated by newline (so, the file name with no extension). 
SAMDIR=$BASEDIR'/'sam # path to folder with sam files (input)
BAMDIR=$BASEDIR'/'bam # path to folder with bam files (output)

# Command
for SAMPLEFILE in `cat $SAMPLELIST`
do
  samtools view -S -b -h -q 30 -m 30 $SAMDIR'/'$SAMPLEFILE.sam | \
  # using a unix pipe (input file is taken from previous step and designated by '-')
  samtools sort - | \
  samtools rmdup -s - $BAMDIR'/'${SAMPLEFILE}'_sorted_rd.bam'
  samtools index $BAMDIR'/'${SAMPLEFILE}'_sorted_rd.bam'
  samtools view $BAMDIR'/'${SAMPLEFILE}'_sorted_rd.bam' | wc -l # count the number of alignments
done

```




    
