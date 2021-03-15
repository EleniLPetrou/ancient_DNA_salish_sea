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

sed -i 's/_sorted_rd_realign.bam//g' $BASENAME'_samples'.txt 


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

### How to interpret the output of output file.genome

Z0 is the probability that at a given locus 0 alleles are identical by descent. If samples are unrelated,  Z0 will be close to 1.

PI_HAT is a measure of overall alleles that are identical by descent. If samples are unrelated, PI_HAT will be close to 0.

Z0, Z1, and Z2 segregate out the probabilities of having IBD of 0, 1, or 2 over the loci, which gives  a way of discriminating between relationship types. Ideal parent-offspring pair has (Z0, Z1, Z2) = (0, 1, 0), i.e. all loci have one allele identical by descent; ideal full sibling = (1/4, 1/2, 1/4), i.e. 25% of loci have 0 alleles IBD, 50% have 1 allele IBD, 25% have 2 alleles IBD. Identical twins or duplicate samples with have an IBD Z2 score that is equal to 1. 

## Step 3: Parse the output of IBD analysis in R

The code below can be used to identify samples are are duplicate pairs and plot the distribution of Z2 scores. 

``` r
# The purpose of this script is to make plots of identity by descent, estimated by plink --genome

# Load libraries
library(tidyverse)
library(cowplot)

# Specify the directory containing the data tables:
DATADIR <- "G:/hybridization_capture/merged_analyses/plink"

# Setwd
setwd(DATADIR)
list.files()

# specify input and output file names:
ancient_file <- "0002.filt.HWE.tidy.snpid.recode.genome"
modern_file <- "0003.filt.HWE.tidy.snpid.recode.genome" 

outfile1 <- "identity_by_descent.pdf"
outfile2 <- "identity_by_descent.txt"

# read in the data
ancient_df <- read.delim(ancient_file, sep = "")
modern_df <- read.delim(modern_file, sep = "")  

# Are there any samples whose Z2 score is ~1? These are either identical twins or duplicates


ancient_df %>%
  filter(Z2 > 0.50)

modern_df %>%
  filter(Z2 > 0.50)  

dup_df <- bind_rows(filter(ancient_df, Z2 > 0.5),filter(modern_df, Z2 > 0.5))%>%
  unite(sample1, FID1:IID1) %>%
  unite(sample2, FID2:IID2)

# Plot the Z2 score for each collection


(plot_anc <- ggplot(ancient_df, aes(x = Z2)) +
    geom_histogram(binwidth = 0.01, fill = "#fc8d59") +
    theme_bw() +
    xlab("Z2 score") +
    ggtitle("Ancient herring"))

(plot_mod <- ggplot(modern_df, aes(x = Z2)) +
    geom_histogram(binwidth = 0.01, fill = "#91bfdb") +
    theme_bw()+
    xlab("Z2 score")+
    ggtitle("Modern herring"))


# Merge the plots

(multi_plot <- plot_grid(plot_anc, plot_mod, labels = c('A', 'B'), label_size = 12))


# save plot to pdf file

ggsave(outfile1,
       plot = multi_plot)

# save a tab-delimited text file with the names of duplicate pairs

write.table(dup_df, 
            file = outfile2, 
            append = FALSE, 
            quote = FALSE, 
            sep = "\t",
            row.names = FALSE)


```

Results:

I identified 3 pairs of samples that are likely duplicate individuals in the anceint data, and 1 set of samples in the modern data (surprise! lol). 
All other samples had a Z2 value that was ~0.

| sample1        | sample 2     | Z2           |
| :------------- | :----------: | -----------: |
|  2B_08         | 2B_13        | 0.9095       |
|  2B_10         | 2B_12        | 0.9637       |
|  2B_14         | 2B_19        | 0.9677       |
| SmBy15_005     | SmBy15_006   | 0.9865     |

I will remove the following individuals from the vcf because they have larger amounts of missing genotypes: 
- 2B_08
- 2B_10
- 2B_19
- SmBy15_006

``` bash
BASEDIR=/media/ubuntu/Herring_aDNA/hybridization_capture/merged_analyses #base directory
VCFDIR=$BASEDIR'/'variants_filtered #vcf directory

########
cd $VCFDIR

vcftools --vcf 0003.filt.HWE.tidy.snpid.recode.vcf \
--remove-indv SmBy15_006 \
--recode --recode-INFO-all \
--out 0003.filt.HWE.tidy.snpid.nodup


vcftools --vcf 0002.filt.HWE.tidy.snpid.recode.vcf \
--remove-indv 2B_08 \
--remove-indv 2B_10 \
--remove-indv 2B_19 \
--recode --recode-INFO-all \
--out 0002.filt.HWE.tidy.snpid.nodup

```

## Step 4: Estimate the genotyping error rate between replicate samples 

I estimated the genotyping error rate between the pairs of samples that were replicate individuals in the ancient and modern data sets.
I did this by writing a simple R script that takes a vcf file as input and compares genotypes across specific user-defined pairs of samples. 
The genotyping error rate was defined as the number of genotype mismatches divided by the total number of genotyped loci.

``` r
# The purpose of this script is to take a vcf file and evaluate the genotyping error rate between pairs of individuals
# that are known duplicate samples. I define the genotyping error rate as the number of mismatched genotypes divided by 
# the number of loci genotyped (i.e. excluding missing loci with missing data)


###################
library(vcfR)
library(genetics)
library(ggplot2)
library(tidyr)
library(dplyr)

#To run this code, put all of your vcf files in a single directory

#setwd

setwd("G:/hybridization_capture/merged_analyses/variants_filtered")
list.files()

# set file names

input_fileName <-  "0002.filt.HWE.tidy.snpid.recode.vcf" 


# Specify which samples are the duplicate pairs
pair1 <- c("2B_08", "2B_13")
pair2 <- c("2B_10", "2B_12")
pair3 <- c("2B_14", "2B_19")

################################################################################
#read in vcf files in directory using vcfR, and start data processing

vcf_data <- read.vcfR(input_fileName )

vcf_data #take a peek
head(vcf_data)

#save metadata as dataframe - you will need this later for plotting
vcf_df <- as.data.frame(vcf_data@fix)
head(vcf_df) #check
class(vcf_df) #should be dataframe


#use vcfR to extract the genotypes from the vcf file --> make matrix
gt <- extract.gt(vcf_data, return.alleles = TRUE)
gt[1:4, 1:4] #take a peek
class(gt) #should be matrix

gt_df <- as.data.frame(gt)
head(gt_df)

# pair-wise comparison
list_pairs <- list(pair1, pair2, pair3)

for (pair in list_pairs) {
  pair_df <- gt_df %>%
  select(pair) 
  
  mismatch_df <- pair_df %>%
  mutate(geno_mismatch = if_else(pair_df[1] == pair_df[2], 1, 0)) %>%
  filter(geno_mismatch == 0)
  
  missing_df <- mismatch_df %>%
  filter(mismatch_df[1] == "." | mismatch_df[2] == ".")
  
  
  genotyping_error_rate = (nrow(mismatch_df) - nrow(missing_df))/(nrow(gt_df) - nrow(missing_df))
  print(pair)
  print(genotyping_error_rate)
  }

```

| sample1        | sample 2     | genotyping error rate           |
| :------------- | :----------: | -----------: |
|  2B_08         | 2B_13        | 0.0239       |
|  2B_10         | 2B_12        | 0.0100       |
|  2B_14         | 2B_19        | 0.0085       |
| SmBy15_005     | SmBy15_006   | 0.0036     |
