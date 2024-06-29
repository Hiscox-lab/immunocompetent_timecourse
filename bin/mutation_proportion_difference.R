### Proportion difference plots with different colours representing proteins using parse.txt files outputted from Syn_NonSyn_parse_aa_V3.pl post DiversiTools analysis

library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyverse)
library(reshape2)
library('cowplot')

# read in parse.txt file for each participant and add column to differentiate when merging
nonsyn_1 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3636_TAGTCGTCAA-ATGGCTCGGT_L002.final_AA_parse.txt", sep="\t")
nonsyn_1$AAposition <- 1:nrow(nonsyn_1) # column to create amino acid position across the whole genome not just the protein
nonsyn_2 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3637_TGGACTAATT-CAGCGGATAA_L002.final_AA_parse.txt", sep="\t")
nonsyn_2$AAposition <- 1:nrow(nonsyn_2)

nonsyn_1$Proportion1 <-  nonsyn_1$CntNonSyn / nonsyn_1$AAcoverage # work out proportions unfiltered for both timepoints
nonsyn_2$Proportion2 <-  nonsyn_2$CntNonSyn / nonsyn_2$AAcoverage
nonsyn_1$Proportion1<-as.numeric(nonsyn_1$Proportion1) # make proportions numeric for calculations and plotting
nonsyn_2$Proportion2<-as.numeric(nonsyn_2$Proportion2)
names(nonsyn_2)[names(nonsyn_2) == "AAcoverage"] <- "AAcoverage2" # change name for sample2 so can copy to dataframe of sample 1
names(nonsyn_2)[names(nonsyn_2) == "CntNonSyn"] <- "CntNonSyn2" # change name for sample2 so can copy to dataframe of sample 1

# first work out proportion, cbind proportion and coverage
nonsyn_3 <- cbind(nonsyn_1, nonsyn_2[, c("Proportion2", "AAcoverage2", "CntNonSyn2")]) # cbind proportion, cov and cnt to filter
nonsyn_3 <- nonsyn_3 %>% filter(AAcoverage > 20) %>% filter(AAcoverage2 > 20) # filter coverage at 20
nonsyn_3 <- nonsyn_3 %>% filter(CntNonSyn > 3) %>% filter(CntNonSyn2 > 3) # filter non syn counts 3
nonsyn_3$Difference <-  nonsyn_3$Proportion2 - nonsyn_3$Proportion1 # work out difference and make numeric
nonsyn_3$Difference <- as.numeric(nonsyn_3$Difference)
nonsyn_3 <- nonsyn_3 %>% filter(Difference != 0)
nonsyn_3$Name <- 'Participant 01'
participant1df <- nonsyn_3

# Participant 02 repeat the same as above
nonsyn_1 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3639_ATCTTATGAT-GCCTGGCAGT_L002.final_AA_parse.txt", sep="\t")
nonsyn_1$AAposition <- 1:nrow(nonsyn_1)
nonsyn_2 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3640_ACCAACGCTG-TAACTTGGAG_L002.final_AA_parse.txt", sep="\t")
nonsyn_2$AAposition <- 1:nrow(nonsyn_2)

nonsyn_1$Proportion1 <-  nonsyn_1$CntNonSyn / nonsyn_1$AAcoverage # work out proportions unfiltered for both timepoints
nonsyn_2$Proportion2 <-  nonsyn_2$CntNonSyn / nonsyn_2$AAcoverage
nonsyn_1$Proportion1<-as.numeric(nonsyn_1$Proportion1) # make proportions numeric
nonsyn_2$Proportion2<-as.numeric(nonsyn_2$Proportion2)
names(nonsyn_2)[names(nonsyn_2) == "AAcoverage"] <- "AAcoverage2" # change name for sample2 so can copy to df of sample 1
names(nonsyn_2)[names(nonsyn_2) == "CntNonSyn"] <- "CntNonSyn2" # change name for sample2 so can copy to df of sample 1

nonsyn_3 <- cbind(nonsyn_1, nonsyn_2[, c("Proportion2", "AAcoverage2", "CntNonSyn2")]) # cbind proportion, cov and cnt to filter
nonsyn_3 <- nonsyn_3 %>% filter(AAcoverage > 20) %>% filter(AAcoverage2 > 20) # filter cov
nonsyn_3 <- nonsyn_3 %>% filter(CntNonSyn > 3) %>% filter(CntNonSyn2 > 3) # filter non syn counts 3
nonsyn_3$Difference <-  nonsyn_3$Proportion2 - nonsyn_3$Proportion1 # work out difference and make numeric
nonsyn_3$Difference <- as.numeric(nonsyn_3$Difference)
nonsyn_3 <- nonsyn_3 %>% filter(Difference != 0)
nonsyn_3$Name <- 'Participant 02'
participant2df <- nonsyn_3

# Participant 3
nonsyn_1 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3643_TACGGCGCGA-ACGTATGCGC_L002.final_AA_parse.txt", sep="\t")
nonsyn_1$AAposition <- 1:nrow(nonsyn_1)
nonsyn_2 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3644_AGCTATCAAC-TGATTGCTTC_L002.final_AA_parse.txt", sep="\t")
nonsyn_2$AAposition <- 1:nrow(nonsyn_2)

nonsyn_1$Proportion1 <-  nonsyn_1$CntNonSyn / nonsyn_1$AAcoverage # work out proportions unfiltered for both timepoints
nonsyn_2$Proportion2 <-  nonsyn_2$CntNonSyn / nonsyn_2$AAcoverage
nonsyn_1$Proportion1<-as.numeric(nonsyn_1$Proportion1) # make proportions numeric
nonsyn_2$Proportion2<-as.numeric(nonsyn_2$Proportion2)
names(nonsyn_2)[names(nonsyn_2) == "AAcoverage"] <- "AAcoverage2" # change name for sample2 so can copy to df of sample 1
names(nonsyn_2)[names(nonsyn_2) == "CntNonSyn"] <- "CntNonSyn2" # change name for sample2 so can copy to df of sample 1

nonsyn_3 <- cbind(nonsyn_1, nonsyn_2[, c("Proportion2", "AAcoverage2", "CntNonSyn2")]) # cbind proportion, cov and cnt to filter
nonsyn_3 <- nonsyn_3 %>% filter(AAcoverage > 20) %>% filter(AAcoverage2 > 20) # filter cov
nonsyn_3 <- nonsyn_3 %>% filter(CntNonSyn > 3) %>% filter(CntNonSyn2 > 3) # filter non syn counts 3
nonsyn_3$Difference <-  nonsyn_3$Proportion2 - nonsyn_3$Proportion1 # work out difference and make numeric
nonsyn_3$Difference <- as.numeric(nonsyn_3$Difference)
nonsyn_3 <- nonsyn_3 %>% filter(Difference != 0)
nonsyn_3$Name <- 'Participant 03'
participant3df <- nonsyn_3

# Participant 5
nonsyn_1 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3649_CTTAGAGGCA-AGCGTCTGGT_L002.final_AA_parse.txt", sep="\t")
nonsyn_1$AAposition <- 1:nrow(nonsyn_1)
nonsyn_2 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3650_ATAGAGCATT-GAGGACCGAT_L002.final_AA_parse.txt", sep="\t")
nonsyn_2$AAposition <- 1:nrow(nonsyn_2)

nonsyn_1$Proportion1 <-  nonsyn_1$CntNonSyn / nonsyn_1$AAcoverage # work out proportions unfiltered for both timepoints
nonsyn_2$Proportion2 <-  nonsyn_2$CntNonSyn / nonsyn_2$AAcoverage
nonsyn_1$Proportion1<-as.numeric(nonsyn_1$Proportion1) # make proportions numeric
nonsyn_2$Proportion2<-as.numeric(nonsyn_2$Proportion2)
names(nonsyn_2)[names(nonsyn_2) == "AAcoverage"] <- "AAcoverage2" # change name for sample2 so can copy to df of sample 1
names(nonsyn_2)[names(nonsyn_2) == "CntNonSyn"] <- "CntNonSyn2" # change name for sample2 so can copy to df of sample 1

nonsyn_3 <- cbind(nonsyn_1, nonsyn_2[, c("Proportion2", "AAcoverage2", "CntNonSyn2")]) # cbind proportion, cov and cnt to filter
nonsyn_3 <- nonsyn_3 %>% filter(AAcoverage > 20) %>% filter(AAcoverage2 > 20) # filter cov
nonsyn_3 <- nonsyn_3 %>% filter(CntNonSyn > 3) %>% filter(CntNonSyn2 > 3) # filter non syn counts 3
nonsyn_3$Difference <-  nonsyn_3$Proportion2 - nonsyn_3$Proportion1 # work out difference and make numeric
nonsyn_3$Difference <- as.numeric(nonsyn_3$Difference)
nonsyn_3 <- nonsyn_3 %>% filter(Difference != 0)
nonsyn_3$Name <- 'Participant 05'
participant5df <- nonsyn_3

# Participant 6
nonsyn_1 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3651_GGTAGCGCAT-TACTGGATAA_L002.final_AA_parse.txt", sep="\t")
nonsyn_1$AAposition <- 1:nrow(nonsyn_1)
nonsyn_2 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3652_TAATGATATT-CGAGTCGAAG_L002.final_AA_parse.txt", sep="\t")
nonsyn_2$AAposition <- 1:nrow(nonsyn_2)

nonsyn_1$Proportion1 <-  nonsyn_1$CntNonSyn / nonsyn_1$AAcoverage # work out proportions unfiltered for both timepoints
nonsyn_2$Proportion2 <-  nonsyn_2$CntNonSyn / nonsyn_2$AAcoverage
nonsyn_1$Proportion1<-as.numeric(nonsyn_1$Proportion1) # make proportions numeric
nonsyn_2$Proportion2<-as.numeric(nonsyn_2$Proportion2)
names(nonsyn_2)[names(nonsyn_2) == "AAcoverage"] <- "AAcoverage2" # change name for sample2 so can copy to df of sample 1
names(nonsyn_2)[names(nonsyn_2) == "CntNonSyn"] <- "CntNonSyn2" # change name for sample2 so can copy to df of sample 1

nonsyn_3 <- cbind(nonsyn_1, nonsyn_2[, c("Proportion2", "AAcoverage2", "CntNonSyn2")]) # cbind proportion, cov and cnt to filter
nonsyn_3 <- nonsyn_3 %>% filter(AAcoverage > 20) %>% filter(AAcoverage2 > 20) # filter cov
nonsyn_3 <- nonsyn_3 %>% filter(CntNonSyn > 3) %>% filter(CntNonSyn2 > 3) # filter non syn counts 3
nonsyn_3$Difference <-  nonsyn_3$Proportion2 - nonsyn_3$Proportion1 # work out difference and make numeric
nonsyn_3$Difference <- as.numeric(nonsyn_3$Difference)
nonsyn_3 <- nonsyn_3 %>% filter(Difference != 0)
nonsyn_3$Name <- 'Participant 06'
participant6df <- nonsyn_3

# Participant 7
nonsyn_1 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3653_GGCTAAGGCG-TTCTTCTATT_L002.final_AA_parse.txt", sep="\t")
nonsyn_1$AAposition <- 1:nrow(nonsyn_1)
nonsyn_2 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3654_TTAGTACCTT-GACCAATTCT_L002.final_AA_parse.txt", sep="\t")
nonsyn_2$AAposition <- 1:nrow(nonsyn_2)

nonsyn_1$Proportion1 <-  nonsyn_1$CntNonSyn / nonsyn_1$AAcoverage # work out proportions unfiltered for both timepoints
nonsyn_2$Proportion2 <-  nonsyn_2$CntNonSyn / nonsyn_2$AAcoverage
nonsyn_1$Proportion1<-as.numeric(nonsyn_1$Proportion1) # make proportions numeric
nonsyn_2$Proportion2<-as.numeric(nonsyn_2$Proportion2)
names(nonsyn_2)[names(nonsyn_2) == "AAcoverage"] <- "AAcoverage2" # change name for sample2 so can copy to df of sample 1
names(nonsyn_2)[names(nonsyn_2) == "CntNonSyn"] <- "CntNonSyn2" # change name for sample2 so can copy to df of sample 1

nonsyn_3 <- cbind(nonsyn_1, nonsyn_2[, c("Proportion2", "AAcoverage2", "CntNonSyn2")]) # cbind proportion, cov and cnt to filter
nonsyn_3 <- nonsyn_3 %>% filter(AAcoverage > 20) %>% filter(AAcoverage2 > 20) # filter cov
nonsyn_3 <- nonsyn_3 %>% filter(CntNonSyn > 3) %>% filter(CntNonSyn2 > 3) # filter non syn counts 3
nonsyn_3$Difference <-  nonsyn_3$Proportion2 - nonsyn_3$Proportion1 # work out difference and make numeric
nonsyn_3$Difference <- as.numeric(nonsyn_3$Difference)
nonsyn_3 <- nonsyn_3 %>% filter(Difference != 0)
nonsyn_3$Name <- 'Participant 07'
participant7df <- nonsyn_3

# Participant 9
nonsyn_1 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3659_GTAGCTACCG-TCAGAACGAC_L002.final_AA_parse.txt", sep="\t")
nonsyn_1$AAposition <- 1:nrow(nonsyn_1)
nonsyn_2 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3660_ATTAATAGAA-GGCAGAGGAG_L002.final_AA_parse.txt", sep="\t")
nonsyn_2$AAposition <- 1:nrow(nonsyn_2)

nonsyn_1$Proportion1 <-  nonsyn_1$CntNonSyn / nonsyn_1$AAcoverage # work out proportions unfiltered for both timepoints
nonsyn_2$Proportion2 <-  nonsyn_2$CntNonSyn / nonsyn_2$AAcoverage
nonsyn_1$Proportion1<-as.numeric(nonsyn_1$Proportion1) # make proportions numeric
nonsyn_2$Proportion2<-as.numeric(nonsyn_2$Proportion2)
names(nonsyn_2)[names(nonsyn_2) == "AAcoverage"] <- "AAcoverage2" # change name for sample2 so can copy to df of sample 1
names(nonsyn_2)[names(nonsyn_2) == "CntNonSyn"] <- "CntNonSyn2" # change name for sample2 so can copy to df of sample 1

nonsyn_3 <- cbind(nonsyn_1, nonsyn_2[, c("Proportion2", "AAcoverage2", "CntNonSyn2")]) # cbind proportion, cov and cnt to filter
nonsyn_3 <- nonsyn_3 %>% filter(AAcoverage > 20) %>% filter(AAcoverage2 > 20) # filter cov
nonsyn_3 <- nonsyn_3 %>% filter(CntNonSyn > 3) %>% filter(CntNonSyn2 > 3) # filter non syn counts 3
nonsyn_3$Difference <-  nonsyn_3$Proportion1 - nonsyn_3$Proportion2 # work out difference and make numeric
nonsyn_3$Difference <- as.numeric(nonsyn_3$Difference)
nonsyn_3 <- nonsyn_3 %>% filter(Difference != 0)
nonsyn_3$Name <- 'Participant 09 S1, S2'
participant9s21df <- nonsyn_3

# Participant 9 s3-s2; participant with three timepoints
nonsyn_4 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3661_ATATGATGAA-ATGAAGGAGG_L002.final_AA_parse.txt", sep="\t")
nonsyn_4$AAposition <- 1:nrow(nonsyn_4)
nonsyn_4$Proportion4 <-  nonsyn_4$CntNonSyn / nonsyn_4$AAcoverage
nonsyn_4$Proportion4<-as.numeric(nonsyn_4$Proportion4)
names(nonsyn_4)[names(nonsyn_4) == "AAcoverage"] <- "AAcoverage3"
names(nonsyn_4)[names(nonsyn_4) == "CntNonSyn"] <- "CntNonSyn3"

nonsyn_5 <- cbind(nonsyn_2, nonsyn_4[, c("Proportion4", "AAcoverage3", "CntNonSyn3")]) # cbind proportion, cov and cnt to filter
nonsyn_5 <- nonsyn_5 %>% filter(AAcoverage2 > 20) %>% filter(AAcoverage3 > 20)
nonsyn_3 <- nonsyn_3 %>% filter(CntNonSyn > 3) %>% filter(CntNonSyn2 > 3) # filter non syn counts 3
nonsyn_5$Difference <-  nonsyn_5$Proportion4 - nonsyn_5$Proportion2 # work out difference and make numeric
nonsyn_5$Difference <- as.numeric(nonsyn_5$Difference)
names(nonsyn_5)[names(nonsyn_5) == "AAcoverage2"] <- "AAcoverage" #renaming to like others so can rbind
names(nonsyn_5)[names(nonsyn_5) == "CntNonSyn2"] <- "CntNonSyn"
names(nonsyn_5)[names(nonsyn_5) == "CntNonSyn3"] <- "CntNonSyn2"
names(nonsyn_5)[names(nonsyn_5) == "AAcoverage3"] <- "AAcoverage2"
names(nonsyn_5)[names(nonsyn_5) == "Proportion2"] <- "Proportion1"
names(nonsyn_5)[names(nonsyn_5) == "Proportion4"] <- "Proportion2"
nonsyn_3 <- nonsyn_3 %>% filter(Difference != 0)
nonsyn_5$Name <- 'Participant 09 S2, S3'
participant9s32df <- nonsyn_5

# Participant 11
nonsyn_1 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3664_ACCATAAGCG-TGAATAACTA_L002.final_AA_parse.txt", sep="\t")
nonsyn_1$AAposition <- 1:nrow(nonsyn_1)
nonsyn_2 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3665_CCGCTGCATT-GAGTCCTGCA_L002.final_AA_parse.txt", sep="\t")
nonsyn_2$AAposition <- 1:nrow(nonsyn_2)

nonsyn_1$Proportion1 <-  nonsyn_1$CntNonSyn / nonsyn_1$AAcoverage # work out proportions unfiltered for both timepoints
nonsyn_2$Proportion2 <-  nonsyn_2$CntNonSyn / nonsyn_2$AAcoverage
nonsyn_1$Proportion1<-as.numeric(nonsyn_1$Proportion1) # make proportions numeric
nonsyn_2$Proportion2<-as.numeric(nonsyn_2$Proportion2)
names(nonsyn_2)[names(nonsyn_2) == "AAcoverage"] <- "AAcoverage2" # change name for sample2 so can copy to df of sample 1
names(nonsyn_2)[names(nonsyn_2) == "CntNonSyn"] <- "CntNonSyn2" # change name for sample2 so can copy to df of sample 1

nonsyn_3 <- cbind(nonsyn_1, nonsyn_2[, c("Proportion2", "AAcoverage2", "CntNonSyn2")]) # cbind proportion, cov and cnt to filter
nonsyn_3 <- nonsyn_3 %>% filter(AAcoverage > 20) %>% filter(AAcoverage2 > 20) # filter cov
nonsyn_3 <- nonsyn_3 %>% filter(CntNonSyn > 3) %>% filter(CntNonSyn2 > 3) # filter non syn counts 3
nonsyn_3$Difference <-  nonsyn_3$Proportion2 - nonsyn_3$Proportion1 # work out difference and make numeric
nonsyn_3$Difference <- as.numeric(nonsyn_3$Difference)
nonsyn_3 <- nonsyn_3 %>% filter(Difference != 0)
nonsyn_3$Name <- 'Participant 11'
participant11df <- nonsyn_3

# Participant 12
nonsyn_1 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3667_GGAACGATTC-CGGCGAGTCC_L002.final_AA_parse.txt", sep="\t")
nonsyn_1$AAposition <- 1:nrow(nonsyn_1)
nonsyn_2 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3668_AAGCATCTCG-CTCATGATTG_L002.final_AA_parse.txt", sep="\t")
nonsyn_2$AAposition <- 1:nrow(nonsyn_2)

nonsyn_1$Proportion1 <-  nonsyn_1$CntNonSyn / nonsyn_1$AAcoverage # work out proportions unfiltered for both timepoints
nonsyn_2$Proportion2 <-  nonsyn_2$CntNonSyn / nonsyn_2$AAcoverage
nonsyn_1$Proportion1<-as.numeric(nonsyn_1$Proportion1) # make proportions numeric
nonsyn_2$Proportion2<-as.numeric(nonsyn_2$Proportion2)
names(nonsyn_2)[names(nonsyn_2) == "AAcoverage"] <- "AAcoverage2" # change name for sample2 so can copy to df of sample 1
names(nonsyn_2)[names(nonsyn_2) == "CntNonSyn"] <- "CntNonSyn2" # change name for sample2 so can copy to df of sample 1

nonsyn_3 <- cbind(nonsyn_1, nonsyn_2[, c("Proportion2", "AAcoverage2", "CntNonSyn2")]) # cbind proportion, cov and cnt to filter
nonsyn_3 <- nonsyn_3 %>% filter(AAcoverage > 20) %>% filter(AAcoverage2 > 20) # filter cov
nonsyn_3 <- nonsyn_3 %>% filter(CntNonSyn > 3) %>% filter(CntNonSyn2 > 3) # filter non syn counts 3
nonsyn_3$Difference <-  nonsyn_3$Proportion2 - nonsyn_3$Proportion1 # work out difference and make numeric
nonsyn_3$Difference <- as.numeric(nonsyn_3$Difference)
nonsyn_3 <- nonsyn_3 %>% filter(Difference != 0)
nonsyn_3$Name <- 'Participant 12'
participant12df <- nonsyn_3

# Participant 15
nonsyn_1 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3678_ATACTATATT-GGTATCATCT_L002.final_AA_parse.txt", sep="\t")
nonsyn_1$AAposition <- 1:nrow(nonsyn_1)
nonsyn_2 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3679_TAGTTATGCG-CTGGCTTAGT_L002.final_AA_parse.txt", sep="\t")
nonsyn_2$AAposition <- 1:nrow(nonsyn_2)

nonsyn_1$Proportion1 <-  nonsyn_1$CntNonSyn / nonsyn_1$AAcoverage # work out proportions unfiltered for both timepoints
nonsyn_2$Proportion2 <-  nonsyn_2$CntNonSyn / nonsyn_2$AAcoverage
nonsyn_1$Proportion1<-as.numeric(nonsyn_1$Proportion1) # make proportions numeric
nonsyn_2$Proportion2<-as.numeric(nonsyn_2$Proportion2)
names(nonsyn_2)[names(nonsyn_2) == "AAcoverage"] <- "AAcoverage2" # change name for sample2 so can copy to df of sample 1
names(nonsyn_2)[names(nonsyn_2) == "CntNonSyn"] <- "CntNonSyn2" # change name for sample2 so can copy to df of sample 1

nonsyn_3 <- cbind(nonsyn_1, nonsyn_2[, c("Proportion2", "AAcoverage2", "CntNonSyn2")]) # cbind proportion, cov and cnt to filter
nonsyn_3 <- nonsyn_3 %>% filter(AAcoverage > 20) %>% filter(AAcoverage2 > 20) # filter cov
nonsyn_3 <- nonsyn_3 %>% filter(CntNonSyn > 3) %>% filter(CntNonSyn2 > 3) # filter non syn counts 3
nonsyn_3$Difference <-  nonsyn_3$Proportion2 - nonsyn_3$Proportion1 # work out difference and make numeric
nonsyn_3$Difference <- as.numeric(nonsyn_3$Difference)
nonsyn_3 <- nonsyn_3 %>% filter(Difference != 0)
nonsyn_3$Name <- 'Participant 15'
participant15df <- nonsyn_3

#  Participant 16
nonsyn_1 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3682_GTCCGTAAGG-CGACCTATAC_L002.final_AA_parse.txt", sep="\t")
nonsyn_1$AAposition <- 1:nrow(nonsyn_1)
nonsyn_2 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3683_CCGGTCGACG-GCTGCATAGC_L002.final_AA_parse.txt", sep="\t")
nonsyn_2$AAposition <- 1:nrow(nonsyn_2)

nonsyn_1$Proportion1 <-  nonsyn_1$CntNonSyn / nonsyn_1$AAcoverage # work out proportions unfiltered for both timepoints
nonsyn_2$Proportion2 <-  nonsyn_2$CntNonSyn / nonsyn_2$AAcoverage
nonsyn_1$Proportion1<-as.numeric(nonsyn_1$Proportion1) # make proportions numeric
nonsyn_2$Proportion2<-as.numeric(nonsyn_2$Proportion2)
names(nonsyn_2)[names(nonsyn_2) == "AAcoverage"] <- "AAcoverage2" # change name for sample2 so can copy to df of sample 1
names(nonsyn_2)[names(nonsyn_2) == "CntNonSyn"] <- "CntNonSyn2" # change name for sample2 so can copy to df of sample 1

nonsyn_3 <- cbind(nonsyn_1, nonsyn_2[, c("Proportion2", "AAcoverage2", "CntNonSyn2")]) # cbind proportion, cov and cnt to filter
nonsyn_3 <- nonsyn_3 %>% filter(AAcoverage > 20) %>% filter(AAcoverage2 > 20) # filter cov
nonsyn_3 <- nonsyn_3 %>% filter(CntNonSyn > 3) %>% filter(CntNonSyn2 > 3) # filter non syn counts 3
nonsyn_3$Difference <-  nonsyn_3$Proportion2 - nonsyn_3$Proportion1 # work out difference and make numeric
nonsyn_3$Difference <- as.numeric(nonsyn_3$Difference)
nonsyn_3 <- nonsyn_3 %>% filter(Difference != 0)
nonsyn_3$Name <- 'Participant 16'
participant16df <- nonsyn_3

#facet wrap and plot

propdiff_rbind_15N <- rbind(participant1df, participant2df, participant3df, participant5df, participant6df, participant7df, participant9s21df, participant9s32df, participant12df, participant16df) # rbind all participants

# change protein names for plot
propdiff_rbind_15N$Protein[propdiff_rbind_15N$Protein == "nsp12_1"] <- "NSP12"
propdiff_rbind_15N$Protein[propdiff_rbind_15N$Protein == "nsp12_2"] <- "NSP12"
propdiff_rbind_15N$Protein[propdiff_rbind_15N$Protein == "nsp1"] <- "NSP1"
propdiff_rbind_15N$Protein[propdiff_rbind_15N$Protein == "nsp2"] <- "NSP2"
propdiff_rbind_15N$Protein[propdiff_rbind_15N$Protein == "nsp3"] <- "NSP3"
propdiff_rbind_15N$Protein[propdiff_rbind_15N$Protein == "nsp4"] <- "NSP4"
propdiff_rbind_15N$Protein[propdiff_rbind_15N$Protein == "nsp5"] <- "NSP5"
propdiff_rbind_15N$Protein[propdiff_rbind_15N$Protein == "nsp6"] <- "NSP6"
propdiff_rbind_15N$Protein[propdiff_rbind_15N$Protein == "nsp7"] <- "NSP7"
propdiff_rbind_15N$Protein[propdiff_rbind_15N$Protein == "nsp8"] <- "NSP8"
propdiff_rbind_15N$Protein[propdiff_rbind_15N$Protein == "nsp9"] <- "NSP9"
propdiff_rbind_15N$Protein[propdiff_rbind_15N$Protein == "nsp10"] <- "NSP10"
propdiff_rbind_15N$Protein[propdiff_rbind_15N$Protein == "nsp13"] <- "NSP13"
propdiff_rbind_15N$Protein[propdiff_rbind_15N$Protein == "nsp14"] <- "NSP14"
propdiff_rbind_15N$Protein[propdiff_rbind_15N$Protein == "nsp15"] <- "NSP15"
propdiff_rbind_15N$Protein[propdiff_rbind_15N$Protein == "nsp16"] <- "NSP16"

propdiff_rbind_15N$Protein <- factor(propdiff_rbind_15N$Protein, levels=c("NSP1", "NSP2", "NSP3", "NSP4", "NSP5", "NSP6", "NSP7", "NSP8", "NSP9", "NSP10", "NSP12", "NSP13", "NSP14", "NSP15", "NSP16", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10"))

# generate colour scheme
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# melt data and plot
melt_dfrbind <- reshape2::melt(propdiff_rbind_15N, id.vars = c('AAposition', 'Name', 'Protein'), measure.vars= 'Difference') # melt data for ggplotting
dfplot <- ggplot(data= melt_dfrbind, aes(x=AAposition, y=value, colour=Protein)) + 
  geom_point(size=1.25) + ylab("Difference in proportion") + xlab("Amino acid position") + theme_bw() +
  scale_colour_manual(values= col_vector) +
  ylim(-1,1) +
  theme(legend.position = "bottom", axis.title.x = element_text(face ="bold", size=14), 
        legend.title = element_text(face = "bold", size=14), 
        axis.title.y = element_text(face="bold", size = 14))

yyy <- dfplot + facet_wrap(vars(Name))

ggsave(filename= "propdiff_colours_20X_participants.tiff", plot = yyy, device = 'tiff', dpi= 300, width=13, height=8)

