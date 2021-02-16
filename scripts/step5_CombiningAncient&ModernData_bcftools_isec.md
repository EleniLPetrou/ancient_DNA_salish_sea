# Combining the ancient and modern genotypes into a data set

The ancient and modern genotype data were combined using *bcftools isec* command, and we retained SNPs that were present in both data sets. 
Subsequently, genotypes with less than a read depth of three sequences were recoded as missing, and all individuals with more than 20% missing genotypes were removed from the data. In the modern samples, loci were tested for deviations from Hardy-Weinberg Equilibium (HWE) using the exact test based on 10,000 Monte Carlo permutations of alleles as implemented in the R package in the R package *HWxtest* (Engels 2019). We used the false discovery rate to correct for multiple testing of HWE; a locus was identified as being out of HWE when the q-value was less than 0.05 and it was removed from all subsequent analyses.

## Merging genotype data with bcftools isec command
 - About:   isec creates intersections, unions and complements of VCF files. 
 - Usage:   bcftools isec [options] <A.vcf.gz> <B.vcf.gz> [...]
 ### Examples:
 - Create intersection and complements of two sets, saving the output in dir/*
 - NB: The vcf files have to be gzipped and indexed before you can use isec command.
 
    ```bcftools isec A.vcf.gz B.vcf.gz -p dir```
  - Arguments used:
  --collapse: Controls how to treat records with duplicate positions and defines compatible records across multiple input files.
              Here by "compatible" we mean records which should be considered as identical by the tools. 
              For example, when performing line intersections, the desire may be to consider as identical all sites with matching positions (bcftools isec --collapse all)
              
   --output-type <b|u|z|v>   b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]
    
  ### Make a new directory and folder that will hold the output of merged analyses: 
  
  /media/ubuntu/hybridization_capture/merged_analyses/variants
   
   ``` bash
   # Directory names
   BASEDIR=/media/ubuntu/Herring_aDNA/hybridization_capture #base directory
   MODERNDIR=$BASEDIR'/'modern_samples/variants
   ANCIENTDIR=$BASEDIR'/'ancient_samples/variants
   MERGEDIR=$BASEDIR'/'merged_analyses/variants
   
   #file names
   ANCIENTFILE=ancient_call_results_filt.vcf
   MODERNFILE=modern_call_results_filt.vcf
   ##########
  

   # Copy and paste the filtered vcf files (from steps 4a and 4b) to the merged_analyses/variants folder, to continue working with them:
   
   cp $ANCIENTDIR'/'$ANCIENTFILE $MERGEDIR
   cp $MODERNDIR'/'$MODERNFILE $MERGEDIR
   
   # Use bgzip command to zip the files:
   
   cd $MERGEDIR
   
   bgzip $MODERNFILE
   bgzip $ANCIENTFILE
   
   # use bcftools index command to index the zipped files
   bcftools index $MODERNFILE'.gz'
   bcftools index $ANCIENTFILE'.gz'
  
  # Use bcftools isec command to merge the two data sets:
  
 bcftools isec $ANCIENTFILE'.gz' $MODERNFILE'.gz' \
 --collapse all \
 --output-type v \
 -p $MERGEDIR
  
  ```
  
  ## Explanation of output files from bcftools isec:
  - 0000.vcf: for records private to $ANCIENTFILE
  - 0001.vcf: for records private to $MODERNFILE
  - 0002.vcf: for records from $ANCIENTFILE shared by both $ANCIENTFILE and $MODERNFILE 
  - 0003.vcf: for records from $MODERNFILE shared by both $ANCIENTFILE and $MODERNFILE

## Filtering ancient data for sequencing depth
After looking at some of the summary statistics plots for the ancient data (0002.vcf), I realized that many of the SNPs identified in the ancient data (N= 7974) were represented by a single read, and could be the result of sequencing error. For this reason, if a SNP had a read depth <3 in the ancient samples, I decided to set that SNP's genotypes to "missing data".

``` bash
BASEDIR=/media/ubuntu/Herring_aDNA/hybridization_capture #base directory
MERGEDIR=$BASEDIR'/'merged_analyses/variants
INFILE=0002.vcf #name of input vcf 
FILTFILE=0002.minDP3.recode.vcf #name of output vcf (with extension)
FILTPREFIX=0002.minDP3 #name of output vcf (without extension)

####

cd $MERGEDIR
vcftools --vcf $INFILE \
--minDP 3 \
--recode --recode-INFO-all \
--out $FILTPREFIX

# calculate some summary statistics using vcftools (these will be used for filtering the individuals and genotypes later on)
vcftools --vcf $FILTFILE --depth --out $FILTPREFIX
vcftools --vcf $FILTFILE  --site-mean-depth --out $FILTPREFIX
vcftools --vcf $FILTFILE  --site-quality --out $FILTPREFIX
vcftools --vcf $FILTFILE  --missing-indv --out $FILTPREFIX
vcftools --vcf $FILTFILE  --missing-site --out $FILTPREFIX
vcftools --vcf $FILTFILE  --get-INFO MQ --get-INFO AD --get-INFO GQ --out $FILTPREFIX
```

## Filtering ancient data: 
 - remove sites with more than 20% missing data and individuals with more than 30% missing data

### Create a "badlist" of ancient SNPs with more than 20% missing data.
``` bash 
BASEDIR=/media/ubuntu/Herring_aDNA/hybridization_capture #base directory
MERGEDIR=$BASEDIR'/'merged_analyses/variants
####

cd $MERGEDIR
head out.lmiss
mawk '$6 > 0.20' 0002.minDP3.lmiss | cut -f1,2 | mawk '!/CHR/' > ancient_bad_loci.txt
head ancient_bad_loci.txt
```
### Create a "badlist" of ancient individuals with more than 30% missing data.
``` bash
BASEDIR=/media/ubuntu/Herring_aDNA/hybridization_capture #base directory
MERGEDIR=$BASEDIR'/'merged_analyses/variants
####

cd $MERGEDIR
head out.imiss
mawk '$5 > 0.30' 0002.minDP3.imiss | cut -f1 | mawk '!/IN/'> ancient_bad_indiv.txt
cat ancient_bad_indiv.txt
```
### Now that we have badlists of ancient sites and individuals, we can remove those sites and samples from our vcf file using vcftools:

``` bash
vcftools --vcf 0002.minDP3.recode.vcf \
--remove ancient_bad_indiv.txt \
--exclude-positions ancient_bad_loci.txt \
--recode --recode-INFO-all \
--out 0002.minDP3.filt

```
After filtering, kept 43 out of 47 Individuals and kept 6601 out of a possible 7974 SNPs in the ancient samples

## Filtering modern data: 
 - If we run the vcftools summary stats for the modern data (0003.vcf from above), we will see that there are no SNPs with more than 20% missing data

``` bash
BASEDIR=/media/ubuntu/Herring_aDNA/hybridization_capture #base directory
MERGEDIR=$BASEDIR'/'merged_analyses/variants
INFILE=0003.vcf #name of input vcf 
PREFIX=0003 #name of output vcf (without extension)

####



# calculate some summary statistics using vcftools (these will be used for filtering the individuals and genotypes later on)
vcftools --vcf $INFILE --depth --out $PREFIX
vcftools --vcf $INFILE  --site-mean-depth --out $PREFIX
vcftools --vcf $INFILE  --site-quality --out $PREFIX
vcftools --vcf $INFILE  --missing-indv --out $PREFIX
vcftools --vcf $INFILE  --missing-site --out $PREFIX
vcftools --vcf $INFILE  --get-INFO MQ --get-INFO AD --get-INFO GQ --out $PREFIX
```

 - However, there are individuals with more than 30% missing data.
 
 ``` bash

# Create a "badlist" of individuals with more than 30% missing data.

mawk '$5 > 0.30' 0003.imiss | cut -f1 | mawk '!/IN/'> modern_bad_indiv.txt
cat modern_bad_indiv.txt

vcftools --vcf 0003.vcf \
--remove modern_bad_indiv.txt \
--exclude-positions ancient_bad_loci.txt \
--recode --recode-INFO-all \
--out 0003.filt
 ```
 
After filtering the modern genotypes, I kept 381 out of 382 Individuals and 6601 out of a possible 7974 Sites









  
