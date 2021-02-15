# Filter the ancient and modern vcf files for loci that are in HWE

To do this, use the following R script:

``` r
# Eleni Petrou. 20191206

# Purpose of script: 

# This script takes as input data in a vcf file.
# It subsequently allows the user to calculate deviations of loci from HWE, 
# and visualize these outputs and save them to a text document.
# Finally, based on the data visualizations, the user
# can  decide how to filter the data.
######################################################################################

# Load required packages
library(vcfR)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(adegenet)
library(HWxtest)
library(stats)
######################################################################################
#setwd
setwd("/media/ubuntu/Herring_aDNA/hybridization_capture/merged_analyses/variants_filtered")
list.files()

# set output_path

mypath_out <- ("/media/ubuntu/Herring_aDNA/hybridization_capture/merged_analyses/HWE_FIS/") 

# specify input file names:

file_name <- "0003.filt.recode.vcf"
base_name <- "0003.filt"

# specify output file name:
output_file <- paste0("results_snps_outHWE_FDR_",base_name, ".txt" )

######################################################################################
# Read in the data (with vcfR) and save it as a df 
my_vcf <-  read.vcfR(file_name)

# Transform the data into a genind object, that hierfstat can use
my_genind <- vcfR2genind(my_vcf)

######################################################################################
# Do some data processing, in preparation for the analyses
# Save a vector of unique locus names
my_loci<- unique(my_genind$loc.fac)
length(my_loci) #check the number of loci

# Create a vector of individual names
name_vec <- indNames(my_genind)
head(name_vec)
tail(name_vec)

name_vec2 <-gsub("_sorted.bam","", name_vec) #remove the trailing garbage from the name vector
head(name_vec2)

# Create a vector of population names
pop_vec <- sapply(strsplit(name_vec2, "_"), `[`, 1)

# Assign those population names to your genind object
pop(my_genind) <- pop_vec
my_genind@pop #check

######################################################################################
# Do the analyses - HWE

hw_df <- hwdf(hwx.test(my_genind))
head(hw_df)

# The first column shows the population, such as “P1,” and locus, such as 
#“BmC5R700_Y9101” of each sample. The second column is the P value from the LLR criterion by 
# default. N is the number of diploid individuals in the sample and k is the number of alleles.

######################################################################################
# Process the output of the analyses

#1. Set row name as a column_to_rownames()
#2. Separate the column in tidyr using separate function, so there is a population, CHR, and POS column
#3. Name the CHR correctly

hw_df2 <- rownames_to_column(hw_df, var = "temp_col")
head(hw_df2)

hw_df3 <- hw_df2 %>%
  separate(temp_col, c("pop", "chr", "chr2", "pos"))

hw_df4 <- hw_df3 %>%
  unite("CHR", chr, chr2, sep = ".")

head(hw_df4)


# Plot the distribution of p-values per population
plot1 <- ggplot(data = hw_df4, aes(`P-val(LLR)`, fill = pop)) + # specify dataframe to use in plot and set x + y parameters
  geom_histogram(binwidth = 0.01) +
  ggtitle(" Distribution of p-values of exact Hardy-Weinberg test in all populations") +
  facet_wrap(~pop)+
  theme(panel.grid=element_blank(),
        axis.text.x = element_text(angle=90, hjust=1))+
  xlab("p-value of exact Hardy-Weinberg test")

plot1

ggsave(paste0(base_name,"_HWE_pvalue_plot.pdf"), plot = plot1, path = mypath_out)

# count if there are any populations that have many loci out of HWE
nloci <- as.numeric(nrow(myFis))

pop_sign <- hw_df4 %>%
  filter(`P-val(LLR)` < 0.05) %>%
  group_by(pop) %>%
  tally() %>%
  mutate(proportion_sig = n/6601)


# visualize if there are any loci that are out of HWE in many populations
loci_sign <- hw_df4 %>%
  filter(`P-val(LLR)` < 0.05) %>%
  group_by(CHR, pos) %>%
  tally()

head(loci_sign)
nrow(loci_sign) #1965 loci out of HWE in 1 or more pops


# Try and adjust p-values by implementing correction for the FDR using the
# p.adjust function from the stats package

hw_df5<- hw_df4 %>%
  group_by(pop) %>%
  mutate(q_value = p.adjust(`P-val(LLR)`, method = "BH"))

# plot the distribution of q-values 

plot2 <- ggplot(data = hw_df5, aes(q_value, fill = pop)) + # specify dataframe to use in plot and set x + y parameters
  geom_histogram(binwidth = 0.01) +
  ggtitle(" Distribution of q-values of exact Hardy-Weinberg test in all populations") +
  facet_wrap(~pop)+
  theme(panel.grid=element_blank(),
        axis.text.x = element_text(angle=90, hjust=1))+
  xlab("q-value of exact Hardy-Weinberg test")

plot2

ggsave(paste0(base_name,"_HWE_qvalue_plot.pdf"), plot = plot2, path = mypath_out)

#  are any loci that are out of HWE in many populations, using q-values?

snps_outHWE <- hw_df5 %>%
  filter(q_value < 0.05) %>%
  group_by(CHR, pos) %>%
  tally()
  
head(snps_outHWE)
nrow(snps_outHWE) #52 SNPs out of HWE


final_df <- snps_outHWE %>%
  select(CHR, pos)
  
# Write out the results of this analysis
write.table(final_df , file = paste0(mypath_out,output_file) , quote = FALSE, row.names = FALSE )






```


Now that I have identified a list of loci that are out of HWE in the modern samples,
I can exclude those loci from both the ancient and modern vcf files.
First, I will copy the list into the variants_filtered directory

``` bash
cd /media/ubuntu/Herring_aDNA/hybridization_capture/merged_analyses/variants_filtered
vcftools --vcf 0003.filt.recode.vcf \
--exclude-positions results_snps_outHWE_FDR_0003.filt.txt \
--recode --recode-INFO-all \
--out 0003.filt.HWE
```
After filtering the modern samples, kept 381 out of 381 Individuals and 6549 out of a possible 6601 Sites


Now, do the same for the ancient samples:

``` bash

vcftools --vcf 0002.filt.recode.vcf \
--exclude-positions results_snps_outHWE_FDR_0003.filt.txt \
--recode --recode-INFO-all \
--out 0002.filt.HWE

```
After filtering the ancient samples, kept 43 out of 43 Individuals and 6549 out of a possible 6601 Sites
