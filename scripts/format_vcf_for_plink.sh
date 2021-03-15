# Prepare non-model organism vcf for plink

# Specify working directories and file names

BASEDIR=/media/ubuntu/Herring_aDNA/hybridization_capture/merged_analyses #base directory
VCFDIR=$BASEDIR'/'variants_filtered #vcf directory
INFILE=0003.filt.HWE.recode.vcf #name of input vcf 
BASENAME=0003.filt.HWE #name of input vcf (without file extension)

##################################################################
# Get into the working directory where the vcf files are stored
cd $VCFDIR

# 1- If needed, rename the samples because plink takes sample names that are in the format "pop_sample"

# save a text file with the names of samples in your vcf
bcftools query --list-samples $INFILE > $BASENAME'_samples'.txt

# use sed to find and replace the trailing garbage in each sample name, as needed:
# sed -i 's/old-text/new-text/g' input.txt. The s is the substitute command of sed for find and replace.
# It tells sed to find all occurrences of ‘old-text’ and replace with ‘new-text’ in a file named input.txt

#sed -i 's/_sorted_rd_realign//g' $BASENAME'_samples'.txt # ancient samples only
sed -i 's/_sorted.bam//g' $BASENAME'_samples'.txt # modern samples only

# 3-  rename the samples using bcftools reheader

bcftools reheader --samples $BASENAME'_samples'.txt $INFILE > $BASENAME'.tidy.recode'.vcf


# Because the SNP ID field is blank in my vcf (there is no rs iD for my snps), I have to create a unique SNP ID for each of my loci.
# I will do this using bcftools. The command -set-id assign ID on the fly. The format is the same as in the query command (see below). By default all existing IDs are replaced.


bcftools annotate --set-id '%CHROM\_%POS' $BASENAME'.tidy.recode'.vcf --output $BASENAME'.tidy.snpid.recode'.vcf
