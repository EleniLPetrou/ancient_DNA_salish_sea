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



## Summary of results

# pair1 genotyping_error_rate = 0.0239 (all het <-> hom discrepancies)
# pair2 genotyping_error_rate = 0.0100 (1 hom <-> hom discrepancies, the rest all het <-> hom discrepancies)
# pair3 genotyping_error_rate = 0.0085 (all het <-> hom discrepancies)
