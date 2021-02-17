# Linkage Disequilibrium
As a first stab at this, I estimated LD separately for modern and ancient samples.

## Part 1: In merged_analyses directory ->
- Make folders to hold LD output:LD_modern, LD_ancient

## Part 2: 
 - For each class of samples (ancient vs. modern) split the vcf by chromosome

### Modern samples
``` bash
BASEDIR=/media/ubuntu/Herring_aDNA/hybridization_capture/merged_analyses
VCFDIR=$BASEDIR'/'variants_filtered #directory containing vcf files
OUTDIR=$BASEDIR'/'LD_modern
VCF=0003.filt.HWE.recode.vcf #(with extension)
BASENAME=0003.filt.HWE #(vcf file without extension)

# List of chromosomes in the Atlantic herring genome

chromosomes="
LR535874.1
LR535858.1
LR535869.1
LR535872.1
LR535862.1
LR535871.1
LR535882.1
LR535867.1
LR535863.1
LR535866.1
LR535861.1
LR535873.1
LR535868.1
LR535859.1
LR535865.1
LR535870.1
LR535878.1
LR535864.1
LR535881.1
LR535860.1
LR535877.1
LR535876.1
LR535880.1
LR535857.1
LR535879.1
LR535875.1"

for chromosome in $chromosomes
do 
    vcftools --vcf $VCFDIR'/'$VCF \
    --chr ${chromosome} \
    --recode --recode-INFO-all \
    --out  $OUTDIR'/'$BASENAME.${chromosome};
done

```

Hypothesis:

- The modern samples should display long-range LD on chrom 8, 12, and 15 (as seen with RAdseq paper)
- the ancient samples...I guess if LD is still concentrated on those chromosomes, then that means that the putative inversions are quite old. I wonder what we will find, this is so exciting!!
