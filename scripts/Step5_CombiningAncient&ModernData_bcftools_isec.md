# Combining the ancient and modern genotypes into a data set

The ancient and modern genotype data were combined using *bcftools isec* command, and we retained SNPs that were present in both data sets. 
Subsequently, genotypes with a read depth below two were recoded as missing, and all individuals with more than 20% missing genotypes were removed from the data. 
In the modern samples, loci were tested for deviations from Hardy-Weinberg Equilibium (HWE) using the exact test based on 10,000 Monte Carlo permutations of alleles as implemented in the R package 
in the R package *HWxtest* (Engels 2019). We used the false discovery rate to correct for multiple testing of HWE anda locus was identified as being out of 
HWE when the q-value was less than 0.05 and it was removed from all subsequent analyses.

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
    
  ### Make a new directory and folder that will hold the output of merged analyses: /media/ubuntu/hybridization_capture/merged_analyses/variants
   
   ``` bash
   # Directory names
   BASEDIR=/media/ubuntu/hybridization_capture #base directory
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
  - 0002.vcf: for records from $ANCIENTFILE shared by both $ANCIENTFILE AND $MODERNFILE
  - 0003.vcf: for records from $MODERNFILE shared by both $ANCIENTFILE AND $MODERNFILE


  
