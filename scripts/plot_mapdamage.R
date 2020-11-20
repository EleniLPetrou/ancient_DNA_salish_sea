# Load libraries
library(tidyverse)
library(cowplot)

##########################################################################
# Specify some directories

base.dir <- "G:/hybridization_capture/ancient_samples" #base directory
mapdamage.dir <- "mapdamage_bam" # folder containing mapdamage results
sample.dir <- "sample_lists" # folder containing list of samples

out.file <- "mapdamage_plot.pdf" #name of plot to save at the end

#Specify the file names containing aggregated mapDamage results

in.file1 <- "5pCtoT_freq_allsamples.txt"
in.file2 <- "3pGtoA_freq_allsamples.txt"
in.file3 <- "sample_list.txt"

##########################################################################
# read in data

fivep.df<- read.delim(paste0(base.dir, "/", mapdamage.dir, "/", in.file1), header = FALSE, col.names = c("pos", "misincorporation"))
threep.df<- read.delim(paste0(base.dir, "/", mapdamage.dir, "/", in.file2), header = FALSE, col.names = c("pos", "misincorporation"))
sample.df <- read.delim(paste0(base.dir, "/", sample.dir, "/", in.file3), header = FALSE, col.names = c("sample"))

# format data for plotting

n.bases <- 25 # mapDamage looks at the first (or last) 25 bases from the 5' or 3' end of sequence

sample <- sample.df$sample %>% #vector of sample names
  rep(n.bases) %>% #repeat this vector by the number of bases in mapDamage
  sort() %>% #sort alphabetically (A to Z)
  as.vector() #save output as vector


fivep.df <- cbind(fivep.df, sample)

fivep.df <- fivep.df %>%
  separate(sample, into = c("layer", "i"), remove = FALSE) #get arch layer data for plotting

threep.df <- cbind(threep.df, sample)

threep.df <- threep.df %>%
  separate(sample, into = c("layer", "i"), remove = FALSE) %>% #get arch layer data for plotting
  mutate(pos2 = -pos) #add an extra column with neg sign in pos, for plotting

# plot

plot1 <- ggplot()+
  geom_line(data = fivep.df, aes(x= pos, y = misincorporation, group = sample, color = layer), alpha = 0.7)+
  ylim(0,0.2)+
  ylab('Frequency of C to T substitution')+
  xlab("distance from 5' end (bp)")+
  theme_bw()+
  scale_color_manual(name = "Archaeological layer", 
                     values = c("#E64B35B2", "#91D1C2B2"), 
                     labels = c("post-1860 CE", "1100-1336 CE"))


plot2 <- ggplot()+
  geom_line(data = threep.df, aes(x= pos2, y = misincorporation, group = sample, color = layer), alpha = 0.7)+
  ylim(0,0.2)+
  ylab('Frequency of G to A substitution')+
  xlab("distance from 3' end (bp)")+
  theme_bw()+
  scale_color_manual(name = "Archaeological layer", 
                     values = c("#E64B35B2", "#91D1C2B2"), 
                     labels = c("post-1860 CE", "1100-1336 CE"))

plot_grid(plot1, plot2, align = "h")



# arrange the  plots in a single row
prow <- plot_grid(
  plot1 + theme(legend.position="none"),
  plot2 + theme(legend.position="none"),
  align = 'vh',
  labels = c("A", "B"),
  hjust = -1,
  nrow = 1
)

prow


# extract the legend from one of the plots
legend <- get_legend(
  # create some space to the left of the legend
  plot1 + theme(legend.box.margin = margin(0, 0, 0, 12))
)

# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
final.plot <- plot_grid(prow, legend, rel_widths = c(3, .75))
final.plot

# save the output

ggsave(out.file, plot = final.plot, path = paste0(base.dir, "/", mapdamage.dir))
