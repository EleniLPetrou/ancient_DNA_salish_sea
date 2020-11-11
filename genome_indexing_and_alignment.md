# Download and index Atlantic herring genome
 NCBI link to genome (Pettersson et al. 2019): https://www.ncbi.nlm.nih.gov/genome/genomes/15477
 Assembly ID : GCA_900700415.1

After I downloaded the genome, I indexed it using bowtie2. This compresses the size of the file and makes queries fast.

# Index Genome

``` bash

BASEDIR=~/media/ubuntu/Herring_aDNA/atlantic_herring_genome_chromosomes # Path to the base directory for the project.
REF=GCA_900700415.1_Ch_v2.0.2_genomic.fna # Genome
BASENAME=GCA_900700415 # Write bt2 data to files with this dir/basename


bowtie2-build $BASEDIR'/'$REF $BASENAME

```

# Trim Illumina adapters from raw sequencing data

The raw sequence data looks like this:

 * The barcode (P5 index) is on the 5' end 
 * The p7 adaptor is on the 3' end


## We will do a two-step trimming process using cutadapt

1. Remove p7 adaptor on on 3' end, using the command -a.

-a : The sequence of the adapter is given with the -a option; this trims the full adaptor sequence anywhere on the read , trims partial adapter sequence at 3' end, and trims full adapter sequences at 3' end

-q : Remove low-quality sequences below a certain phred score

-O : regulates how much overlap you allow with your adapter.

-m : sets the minimum allowed sequence length to keep after trimming

-o : specifies name of the output file

Note that the input file has to be in fasta or fastq format
In this case, the adapter sequence is the same as *B03.P7.part1.F* in Meyer and Kircher 2010: doi:10.1101/pdb.prot5448




``` bash
BASEDIR=~/media/ubuntu/Herring_aDNA/hybridization_capture/raw_data
ADAPTER=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
PHREDSCORE=20
OVERLAP=1
MINLENGTH=20
INPUTFASTQ=Undetermined_S0_L003_R1_001.fastq
OUTPUT=Undetermined_S0_L003_R1_001_cut.fastq

cutadapt -a $ADAPTER -q $PHREDSCORE -O $OVERLAP -m $MINLENGTH -o $OUTPUT $INPUTFASTQ

```



# Align ancient samples to the genome

## Explanation of terms:

bowtie2 -q -x <bt2-idx> -U <r> -S <sam>

-q query input files are in fastq format

-x <bt2-idx> Indexed "reference genome" filename prefix (minus trailing .X.bt2)

-U <r> Files with unpaired reads.

-S <sam> File for SAM output (default: stdout)

``` bash

BASEDIR=~/media/ubuntu/Herring_aDNA/Herring_mybaits_aDNA/trimmed_fastq # Path to the base directory

```

'''


