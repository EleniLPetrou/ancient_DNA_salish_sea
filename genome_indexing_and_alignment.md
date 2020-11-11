# Download and index Atlantic herring genome
 NCBI link to genome (Pettersson et al. 2019): https://www.ncbi.nlm.nih.gov/genome/genomes/15477
 Assembly ID : GCA_900700415.1

After I downloaded the genome, I indexed it using *bowtie2*. This compresses the size of the file and makes queries fast.

# Index Atlantic herring genome

``` bash

BASEDIR=~/media/ubuntu/Herring_aDNA/atlantic_herring_genome_chromosomes # Path to the base directory for the project.
REF=GCA_900700415.1_Ch_v2.0.2_genomic.fna # Genome
BASENAME=GCA_900700415 # Write bt2 data to files with this dir/basename


bowtie2-build $BASEDIR'/'$REF $BASENAME

```

# Trim Illumina adapters from raw ancient sequencing data

The raw sequence data looks like this:

 * The barcode (P5 index) is on the 5' end 
 * The p7 adaptor is on the 3' end


## We will do a two-step trimming process using cutadapt

### 1. Remove p7 adaptor on on 3' end, using the command -a.

-a : The sequence of the adapter is given with the -a option; this trims the full adaptor sequence anywhere on the read , trims partial adapter sequence at 3' end, and trims full adapter sequences at 3' end

-q : Remove low-quality sequences below a certain phred score

-O : regulates how much overlap you allow with your adapter.

-m : sets the minimum allowed sequence length to keep after trimming

-o : specifies name of the output file

Note that the input file has to be in fasta or fastq format
In this case, the adapter sequence is the same as *B03.P7.part1.F* in Meyer and Kircher 2010: doi:10.1101/pdb.prot5448

http://cshprotocols.cshlp.org/content/2010/6/pdb.prot5448.abstract



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
### 2. Demultiplex the samples

The barcodes are saved in a fasta file called *barcodes.fasta*. Its format is like this:

``` bash
>2B_01
^ATGACTGC
>2B_06
^CCGGCGAC

```

To remove the barcode on the 3' end, use the -g argument in cutadapt. Here is an explanation of the arguments used:

 -g : The sequence of the barcode is given with the -g option; trims the full barcode sequence anywhere on the read; trims partial barcode sequence at 5' end; and trims full barcode sequences at 5' end.If the barcode sequence starts with the character '^',  the barcode is 'anchored'. An anchored barcode must  appear in its entirety at the 5' end of the read (it  is a prefix of the read)

-- no-indels: specifies that we do not allow indels in the barcodes

-e : Maximum allowed error rate in barcode (default = 0.1)

-m : sets the minimum allowed sequence length to keep after trimming


``` bash

BASEDIR=~/media/ubuntu/Herring_aDNA/hybridization_capture/raw_data
MINLENGTH=20
ERROR=0.125
INPUTFILE=Undetermined_S0_L003_R1_001_cut.fastq

cutadapt -g file:barcodes.fasta -e $ERROR -m $MINLENGTH --no-indels -o "{name}_cut_trim.fastq" $INPUTFILE --discard-untrimmed

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


