# Check for duplicate samples in the ancient data set, using identity by descent

Each ancient sample was a tiny herring bone. To minimize the probability of sampling multiple bones from each individual fish, we tried to use distinct prootic bones whenever possible. This was not always possible and to achieve decent sample sizes from each archaeological layer we extracted DNA from vertebra as well. In this analysis, I estimate pairwise identity by descent (IBD) in the ancient samples using  *plink*. 

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
## Step 2:  Use *plink* to estimate IBD and identify ancient samples that might be replicate individuals


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

Z0 is the probability that at a given locus 0 alleles are identical by descent. If samples are unrelated,  Z0 will be close to 1.

PI_HAT is a measure of overall alleles that are identical by descent. If samples are unrelated, PI_HAT will be close to 0.

Z0, Z1, and Z2 segregate out the probabilities of having IBD of 0, 1, or 2 over the loci, which gives  a way of discriminating between relationship types. Ideal parent-offspring pair has (Z0, Z1, Z2) = (0, 1, 0), i.e. all loci have one allele identical by descent; ideal full sibling = (1/4, 1/2, 1/4), i.e. 25% of loci have 0 alleles IBD, 50% have 1 allele IBD, 25% have 2 alleles IBD. Identical twins or duplicate samples with have an IBD Z2 score that is equal to 1. 

