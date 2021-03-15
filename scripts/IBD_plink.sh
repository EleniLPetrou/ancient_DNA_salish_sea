# Specify the directory names and file names
BASEDIR=/media/ubuntu/Herring_aDNA/hybridization_capture/merged_analyses #base directory
VCFDIR=$BASEDIR'/'variants_filtered #vcf directory
OUTDIR=$BASEDIR'/'plink
INFILE=0003.filt.HWE.tidy.snpid.recode.vcf #name of input vcf (in plink format)
OUTFILE=0003.filt.HWE.tidy.snpid.recode


##########################################################
conda activate plink_environment

# Estimate identity by descent for each sample
plink --genome  --allow-extra-chr --vcf $VCFDIR'/'$INFILE --out $OUTDIR'/'$OUTFILE

# How to interpret the output of file file.genome
# Identical twins or duplicate samples with have a Z2 score of 1. The Z2 score expresses the probability of having two alleles identical by descent
