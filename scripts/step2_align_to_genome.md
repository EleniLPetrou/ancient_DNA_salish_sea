# Download and index Atlantic herring genome

I used the chromosome-level assembly of the Atlantic herring genome (Pettersson et al. 2019) as a reference genome for this project. Here is the NCBI link to the genome : https://www.ncbi.nlm.nih.gov/genome/genomes/15477
Assembly ID : GCA_900700415.1

After I downloaded the genome from NCBI, I indexed it using *bowtie2*. Indexing compresses the size of the file and makes queries fast.

## Index Atlantic herring genome

``` bash

BASEDIR=/media/ubuntu/Herring_aDNA/atlantic_herring_genome_chromosomes # Path to directory with genome.
REF=GCA_900700415.1_Ch_v2.0.2_genomic.fna # Genome fasta file name
BASENAME=GCA_900700415 # Write bt2 data to files with this dir/basename


bowtie2-build $BASEDIR'/'$REF $BASEDIR'/'$BASENAME

```


# Align ancient samples to the genome

## Explanation of terms:

bowtie2 -q -x <bt2-idx> -U <r> -S <sam>
  
  -q query input files are in fastq format

-x <bt2-idx> Indexed "reference genome" filename prefix (minus trailing .X.bt2)

-U <r> Files with unpaired reads.

-S <sam> File for SAM output (default: stdout)

``` bash


```
