## Step 1:  Prepare vcf file for plink

Plink is very grumpy about vcf format (it expects human genome data), so I had to tidy up my vcf file before using plink


``` bash
# Specify working directories and file names

BASEDIR=/media/ubuntu/Herring_aDNA/hybridization_capture/merged_analyses #base directory
VCFDIR=$BASEDIR'/'variants_filtered #vcf directory
INFILE=0002.filt.HWE.recode.vcf #name of input vcf 
BASENAME=0002.filt.HWE #name of input vcf (without file extension)

# Get into the working directory where the vcf files are stored
cd $VCFDIR

# If needed, rename the samples because plink takes sample names that are in the format "pop_sample"

# save a text file with the names of samples in your vcf
bcftools query --list-samples $INFILE > $BASENAME'_samples'.txt

# use sed to find and replace the trailing garbage in each sample name, as needed:
# sed -i 's/old-text/new-text/g' input.txt. The s is the substitute command of sed for find and replace.
# It tells sed to find all occurrences of ‘old-text’ and replace with ‘new-text’ in a file named input.txt

sed -i 's/_sorted_rd_realign//g' $BASENAME'_samples'.txt 


# Rename the samples using bcftools reheader

bcftools reheader --samples $BASENAME'_samples'.txt $INFILE > $BASENAME'.tidy.recode'.vcf


# Because the SNP ID field is blank in my vcf (there is no rs iD for my snps), I have to create a unique SNP ID for each of my loci.
# I will do this using bcftools. The command -set-id assign ID on the fly. The format is the same as in the query command (see below). By default all existing IDs are replaced.


bcftools annotate --set-id '%CHROM\_%POS' $BASENAME'.tidy.recode'.vcf --output $BASENAME'.tidy.snpid.recode'.vcf

```
## Step 2:  Prepare vcf file for plink


``` bash 

# Specify the directory names and file names
BASEDIR=/media/ubuntu/Herring_aDNA/hybridization_capture/merged_analyses #base directory
VCFDIR=$BASEDIR'/'variants_filtered #vcf directory
OUTDIR=$BASEDIR'/'plink
INFILE=0002.filt.HWE.tidy.snpid.recode.vcf #name of input vcf (in plink format)
OUTFILE=0002.filt.HWE.tidy.snpid.recode

############
conda activate plink_environment

# Estimate identity by descent for each sample
plink --genome  --allow-extra-chr --vcf $VCFDIR'/'$INFILE --out $OUTDIR'/'$OUTFILE

```

### How to interpret the output of file file.genome

Identical twins or duplicate samples with have a Z2 score of 1. The Z2 score expresses the probability of having two alleles identical by descent. 

