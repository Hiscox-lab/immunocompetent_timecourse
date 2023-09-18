## Plotting sgmrna abundance using output from LeTRS tool (Dong et al., 2022)

library(data.table)
library(cowplot)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyverse)
library(reshape2)

setwd("/home/hannahg/projects/dstl_project/data/nimagen/Raw_fastq/LeTRS/kj_all_25/") # set working directory to where LeTRS output known junction txt files are located

# function to read in known junction files per sample, all located in one directory
letrs_readin_knownjct <- function(filepath ="./", file_pattern, file_list) {
  
  temp <- list.files(filepath, pattern=file_pattern) # creates list with filenames as in the directory
  
  file_list <- list() # creates empty list to add to in the for loop
  
  for (i in 1:length(temp)) { # read in all data files as DFs, split up bracket values to get peak norm counts all and with at least 1 primer
    
    file_list[[i]] <- read.delim(temp[i]) 
    file_list[[i]] <- file_list[[i]][-c(11:17), ] # removes rows with writing in
    file_list[[i]]$Name <- assign(paste( temp[i]), temp[i])
    setDT(file_list[[i]])[, paste0("peak_normalized_count", 1:2) := tstrsplit(peak_normalized_count, "(", type.convert = TRUE, fixed = TRUE)]
    setDT(file_list[[i]])[, paste0("peak_normalized_count2", 1:4) := tstrsplit(peak_normalized_count2, ",", type.convert = TRUE, fixed = TRUE)]
    setDT(file_list[[i]])[, paste0("peak_normalized_count24", 1:2) := tstrsplit(peak_normalized_count24, ")", type.convert = TRUE, fixed = TRUE)]
    names(file_list[[i]])[names(file_list[[i]]) == "peak_normalized_count241"] <- "peak_normcount_atleast1primer"
    names(file_list[[i]])[names(file_list[[i]]) == "peak_normalized_count1"] <- "peak_normcount_all"
    file_list[[i]]$subgenome <- factor(file_list[[i]]$subgenome, levels=c("S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10"))
    file_list[[i]]$proportion1 <- file_list[[i]]$peak_normcount_all / sum(file_list[[i]]$peak_normcount_all) # getting proportion using sum of column to generate %
  }
  
  myfilenames <- gsub(".tab","", temp) # give names to each file, removing unnecessary strings
  names(file_list) <- myfilenames
  return(file_list)
}

aa_list_letrs <- letrs_readin_knownjct(filepath = "./", file_pattern = ".tab", file_list = aa_list)

# for participants with more than one timepoint cbind on S2/S3 onto sample1 so can plot together
# knownjunction_patient_sample == kj_pXsX
# read in known junction files and rename columns so can plot different timepoints on one graph

kj_p1s1 <- aa_list_letrs[["known_junction36"]]
kj_p1s2 <- aa_list_letrs[["known_junction37"]]
names(kj_p1s2)[names(kj_p1s2) == "peak_normcount_all"] <- "peak_normcount_all_S2"
names(kj_p1s2)[names(kj_p1s2) == "proportion1"] <- "proportion2"
kj_p1 <- cbind(kj_p1s1, kj_p1s2[,c("peak_normcount_all_S2", "proportion2"), drop=FALSE])
kj_p1$Patient <- "Participant 01"

kj_p2s1 <- aa_list_letrs[["known_junction39"]]
kj_p2s2 <- aa_list_letrs[["known_junction40"]]
names(kj_p2s2)[names(kj_p2s2) == "peak_normcount_all"] <- "peak_normcount_all_S2"
names(kj_p2s2)[names(kj_p2s2) == "proportion1"] <- "proportion2"
kj_p2 <- cbind(kj_p2s1, kj_p2s2[,c("peak_normcount_all_S2", "proportion2"), drop=FALSE])
kj_p2$Patient <- "Participant 02"

kj_p3s1 <- aa_list_letrs[["known_junction43"]]
kj_p3s2 <- aa_list_letrs[["known_junction44"]]
names(kj_p3s2)[names(kj_p3s2) == "peak_normcount_all"] <- "peak_normcount_all_S2"
names(kj_p3s2)[names(kj_p3s2) == "proportion1"] <- "proportion2"
kj_p3 <- cbind(kj_p3s1, kj_p3s2[,c("peak_normcount_all_S2", "proportion2"), drop=FALSE])
kj_p3$Patient <- "Participant 03"

kj_p4s1 <- aa_list_letrs[["known_junction47"]]
kj_p4s1$Patient <- "Participant 04"
kj_p5s1 <- aa_list_letrs[["known_junction49"]]
kj_p5s2 <- aa_list_letrs[["known_junction50"]]
names(kj_p5s2)[names(kj_p5s2) == "peak_normcount_all"] <- "peak_normcount_all_S2"
names(kj_p5s2)[names(kj_p5s2) == "proportion1"] <- "proportion2"
kj_p5 <- cbind(kj_p1s1, kj_p1s2[,c("peak_normcount_all_S2", "proportion2"), drop=FALSE])
kj_p5$Patient <- "Participant 05"

kj_p6s1 <- aa_list_letrs[["known_junction51"]]
kj_p6s2 <- aa_list_letrs[["known_junction52"]]
names(kj_p6s2)[names(kj_p6s2) == "peak_normcount_all"] <- "peak_normcount_all_S2"
names(kj_p6s2)[names(kj_p6s2) == "proportion1"] <- "proportion2"
kj_p6 <- cbind(kj_p6s1, kj_p6s2[,c("peak_normcount_all_S2", "proportion2"), drop=FALSE])
kj_p6$Patient <- "Participant 06"

kj_p7s1 <- aa_list_letrs[["known_junction53"]]
kj_p7s2 <- aa_list_letrs[["known_junction54"]]
names(kj_p7s2)[names(kj_p7s2) == "peak_normcount_all"] <- "peak_normcount_all_S2"
names(kj_p7s2)[names(kj_p7s2) == "proportion1"] <- "proportion2"
kj_p7 <- cbind(kj_p7s1, kj_p7s2[,c("peak_normcount_all_S2", "proportion2"), drop=FALSE])
kj_p7$Patient <- "Participant 07"

kj_p9s1 <- aa_list_letrs[["known_junction59"]]
kj_p9s2 <- aa_list_letrs[["known_junction60"]]
names(kj_p9s2)[names(kj_p9s2) == "peak_normcount_all"] <- "peak_normcount_all_S2"
names(kj_p9s2)[names(kj_p9s2) == "proportion1"] <- "proportion2"
kj_p9s3 <- aa_list_letrs[["known_junction61"]]
names(kj_p9s3)[names(kj_p9s3) == "peak_normcount_all"] <- "peak_normcount_all_S3"
names(kj_p9s3)[names(kj_p9s3) == "proportion1"] <- "proportion3"
kj_p9 <- cbind(kj_p9s1, kj_p9s2[,c("peak_normcount_all_S2", "proportion2"), drop=FALSE], kj_p9s3[,c("peak_normcount_all_S3", "proportion3"), drop= FALSE])
kj_p9$Patient <- "Participant 09"

kj_p10s1 <- aa_list_letrs[["known_junction62"]]
kj_p10s1$Patient <- "Participant 10"
kj_p11s3 <- aa_list_letrs[["known_junction65"]]
kj_p11s3$Patient <- "Participant 11"
kj_p12s1 <- aa_list_letrs[["known_junction67"]]
kj_p12s2 <- aa_list_letrs[["known_junction68"]]
names(kj_p12s2)[names(kj_p12s2) == "peak_normcount_all"] <- "peak_normcount_all_S2"
names(kj_p12s2)[names(kj_p12s2) == "proportion1"] <- "proportion2"
kj_p12 <- cbind(kj_p12s1, kj_p12s2[,c("peak_normcount_all_S2", "proportion2"), drop=FALSE])
kj_p12$Patient <- "Participant 12"

kj_p13s1 <- aa_list_letrs[["known_junction70"]]
kj_p13s1$Patient <- "Participant 13"
kj_p14s1 <- aa_list_letrs[["known_junction74"]]
kj_p14s1$Patient <- "Participant 14"
kj_p15s2 <- aa_list_letrs[["known_junction79"]]
kj_p15s2$Patient <- "Participant 15"
kj_p16s1 <- aa_list_letrs[["known_junction82"]]
kj_p16s2 <- aa_list_letrs[["known_junction83"]]
names(kj_p16s2)[names(kj_p16s2) == "peak_normcount_all"] <- "peak_normcount_all_S2"
names(kj_p16s2)[names(kj_p16s2) == "proportion1"] <- "proportion2"
kj_p16 <- cbind(kj_p16s1, kj_p16s2[,c("peak_normcount_all_S2", "proportion2"), drop=FALSE])
kj_p16$Patient <- "Participant 16"

rbind_kj <- rbind(kj_p1, kj_p2, kj_p3, kj_p4s1, kj_p5, kj_p6, kj_p7, kj_p9, kj_p10s1, kj_p11s3, kj_p12, kj_p13s1, kj_p14s1, kj_p15s2, kj_p16, fill=TRUE) # rbind all participants dataframes to plot
# need to fill in NA as some patients have 1/2/3 samples etc

melt_rbindkj <- reshape2::melt(rbind_kj, id.vars = c('subgenome', 'Patient'), measure.vars= c('proportion1', 'proportion2', 'proportion3'))
melt_rbindkj_pnc <- reshape2::melt(rbind_kj, id.vars = c('subgenome', 'Patient'), measure.vars= c('peak_normcount_all'))

# plot as bar graph in ggplot
letrs_plot <- ggplot() + geom_bar(data = melt_rbindkj, aes(x = subgenome, y = value, fill = variable), position = "dodge", stat = "identity") +
  labs(x= "Subgenome", y="Proportion normalised count") + scale_fill_discrete(name= "Sample", labels = c("1", "2", "3")) + 
  theme(legend.position = "bottom", 
        axis.title.x = element_text(face="bold", size=14), 
        legend.title = element_text(size=14), 
        axis.title.y = element_text(face="bold", size=14),
        axis.text.x = element_text(angle=90, vjust=0.7)) +
  facet_wrap(~ Patient) 

ggsave(filename = "letrs_nim_propnormcount_participants.tiff", plot= letrs_plot, device='tiff', dpi= 300, width=8, height=8)



