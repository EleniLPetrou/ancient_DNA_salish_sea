
# Trim Illumina adapters from raw ancient sequencing data

Once I had raw sequencing data from the hybridization capture of aDNA, I trimmed Illumina adapters.

The raw sequence data looks like this:

 * The barcode (P5 index) is on the 5' end 
 * The p7 adaptor is on the 3' end


## Two-step trimming process using *cutadapt*

### 1. To remove p7 adaptor on on 3' end, use the command -a.
*Note that the input file has to be in fasta or fastq format

Explanation of *cutadapt* arguments used:

-a : The sequence of the adapter is given with the -a option; this trims the full adaptor sequence anywhere on the read , trims partial adapter sequence at 3' end, and trims full adapter sequences at 3' end

-q : Remove low-quality sequences below a certain phred score

-O : regulates how much overlap you allow with your adapter.

-m : sets the minimum allowed sequence length to keep after trimming

-o : specifies name of the output file



In the aDNA samples, the Illumina adapter sequence is the same as *B03.P7.part1.F* sequence reported in Meyer and Kircher 2010: doi:10.1101/pdb.prot5448 - http://cshprotocols.cshlp.org/content/2010/6/pdb.prot5448.abstract



``` bash
BASEDIR=/media/ubuntu/Herring_aDNA/hybridization_capture/ancient_samples/raw_data #directory with raw data
ADAPTER=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC #adapter used in this study
PHREDSCORE=20 #minimum phred score
OVERLAP=1 #overlap
MINLENGTH=20 #minimum sequence length
INPUT=Undetermined_S0_L003_R1_001.fastq.gz #input file
OUTPUT=Undetermined_S0_L003_R1_001_cut.fastq.gz #output file

cutadapt -a $ADAPTER -q $PHREDSCORE -O $OVERLAP -m $MINLENGTH -o $BASEDIR'/'$OUTPUT $BASEDIR'/'$INPUT

```
### 2. Demultiplex the samples

Next, ancient samples were demultiplexed using *cutadapt*. 

I created a directory to store the demultiplexed fastq files

``` bash
BASEDIR=/media/ubuntu/Herring_aDNA/hybridization_capture/ancient_samples # directory for ancient samples

mkdir $BASEDIR'/'trimmed_fastq # make a directory to store demultiplexed fastq files

```

Next, I saved the individual indeces (barcodes) in a fasta file called *barcodes.fasta*. Each sample name is followed by a unique index, like this:

``` bash
>2B_01
^ATGACTGC
>2B_06
^CCGGCGAC

```

To remove the barcode on the 3' end, use the -g argument in cutadapt. Here is an explanation of the arguments used:

 -g : The sequence of the barcode is given with the -g option; trims the full barcode sequence anywhere on the read; trims partial barcode sequence at 5' end; and trims full barcode sequences at 5' end.If the barcode sequence starts with the character '^',  the barcode is 'anchored'. An anchored barcode must  appear in its entirety at the 5' end of the read (it is a prefix of the read).

-- no-indels: specifies that we do not allow indels in the barcodes

-e : Maximum allowed error rate in barcode (default = 0.1)

-m : sets the minimum allowed sequence length to keep after trimming


``` bash

BASEDIR=/media/ubuntu/Herring_aDNA/hybridization_capture/ancient_samples # base directory
BARCODES=raw_data/barcodes.fasta # relative path and file name for barcodes
MINLENGTH=20
ERROR=0.125 # I used a slightly larger error rate than the default
INPUTFILE=raw_data/Undetermined_S0_L003_R1_001_cut.fastq.gz # relative path and file name for sequencing data (adapters removed)
OUTDIR=trimmed_fastq #name of output directory


cutadapt -g file:$BASEDIR'/'$BARCODES -e $ERROR -m $MINLENGTH --no-indels --discard-untrimmed -o $BASEDIR'/'$OUTDIR'/'"{name}_cut_trim.fastq" $BASEDIR'/'$INPUTFILE 

```

# Quality control of sequencing data using *FastQC*

It is always a great idea to look at the demultiplexed fastq files and verify that:
1. data quality looks good
2. there are no adapters still hanging around the sequencing data

I did this using the program *FastQC*

``` bash
BASEDIR=/media/ubuntu/Herring_aDNA/hybridization_capture/ancient_samples #path of base directory
INPUTDIR=$BASEDIR'/'trimmed_fastq #path to fastq files

# Make a directory to hold the fastqc results
mkdir $BASEDIR'/'fastqc 

# Specify path to the folder where you want the fastqc results to be stored
OUTPUTDIR=$BASEDIR'/'fastqc 

for FILE in $INPUTDIR'/'*.fastq # for any file ending in .fastq in this directory
do
fastqc -f fastq -o ${OUTPUTDIR} ${FILE}
done

```

