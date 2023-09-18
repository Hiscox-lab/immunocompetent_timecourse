### Proportion difference plots with different colours representing proteins 

library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyverse)
library(reshape2)
library('cowplot')

# could R bind all sample 1 and then cbind on rbind sample 2?

nonsyn_1 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3636_TAGTCGTCAA-ATGGCTCGGT_L002.final_AA_parse.txt", sep="\t")
nonsyn_1$AAposition <- 1:nrow(nonsyn_1)
nonsyn_2 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3637_TGGACTAATT-CAGCGGATAA_L002.final_AA_parse.txt", sep="\t")
nonsyn_2$AAposition <- 1:nrow(nonsyn_2)

nonsyn_1$Proportion1 <-  nonsyn_1$CntNonSyn / nonsyn_1$AAcoverage # work out proportions unfiltered for both timepoints
nonsyn_2$Proportion2 <-  nonsyn_2$CntNonSyn / nonsyn_2$AAcoverage
nonsyn_1$Proportion1<-as.numeric(nonsyn_1$Proportion1) # make proportions numeric
nonsyn_2$Proportion2<-as.numeric(nonsyn_2$Proportion2)
names(nonsyn_2)[names(nonsyn_2) == "AAcoverage"] <- "AAcoverage2" # change name for sample2 so can copy to df of sample 1
names(nonsyn_2)[names(nonsyn_2) == "CntNonSyn"] <- "CntNonSyn2" # change name for sample2 so can copy to df of sample 1

# first work out proportion, cbind proportion and coverage
nonsyn_3 <- cbind(nonsyn_1, nonsyn_2[, c("Proportion2", "AAcoverage2", "CntNonSyn2")]) # cbind proportion, cov and cnt to filter
nonsyn_3 <- nonsyn_3 %>% filter(AAcoverage > 20) %>% filter(AAcoverage2 > 20) # filter cov
#nonsyn_3 <- nonsyn_3 %>% filter(CntNonSyn > 4) %>% filter(CntNonSyn2 > 4) # filter non syn counts at arbitrary value
nonsyn_3$Difference <-  nonsyn_3$Proportion2 - nonsyn_3$Proportion1 # work out difference and make numeric
nonsyn_3$Difference <- as.numeric(nonsyn_3$Difference)
nonsyn_3 <- nonsyn_3 %>% filter(Difference != 0)
nonsyn_3$Name <- 'Participant 01'
Patient1df <- nonsyn_3

#patient 02
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
#nonsyn_3 <- nonsyn_3 %>% filter(CntNonSyn > 4) %>% filter(CntNonSyn2 > 4) # filter non syn counts at arbitrary value
nonsyn_3$Difference <-  nonsyn_3$Proportion2 - nonsyn_3$Proportion1 # work out difference and make numeric
nonsyn_3$Difference <- as.numeric(nonsyn_3$Difference)
nonsyn_3 <- nonsyn_3 %>% filter(Difference != 0)
nonsyn_3$Name <- 'Participant 02'
Patient2df <- nonsyn_3

#patient 3
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
#nonsyn_3 <- nonsyn_3 %>% filter(CntNonSyn > 4) %>% filter(CntNonSyn2 > 4) # filter non syn counts at arbitrary value
nonsyn_3$Difference <-  nonsyn_3$Proportion2 - nonsyn_3$Proportion1 # work out difference and make numeric
nonsyn_3$Difference <- as.numeric(nonsyn_3$Difference)
nonsyn_3 <- nonsyn_3 %>% filter(Difference != 0)
nonsyn_3$Name <- 'Participant 03'
Patient3df <- nonsyn_3

#patient 5
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
#nonsyn_3 <- nonsyn_3 %>% filter(CntNonSyn > 2) %>% filter(CntNonSyn2 > 2) # filter non syn counts at arbitrary value
nonsyn_3$Difference <-  nonsyn_3$Proportion2 - nonsyn_3$Proportion1 # work out difference and make numeric
nonsyn_3$Difference <- as.numeric(nonsyn_3$Difference)
nonsyn_3 <- nonsyn_3 %>% filter(Difference != 0)
nonsyn_3$Name <- 'Participant 05'
Patient5df <- nonsyn_3

#patient 6
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
#nonsyn_3 <- nonsyn_3 %>% filter(CntNonSyn > 4) %>% filter(CntNonSyn2 > 4) # filter non syn counts at arbitrary value
nonsyn_3$Difference <-  nonsyn_3$Proportion2 - nonsyn_3$Proportion1 # work out difference and make numeric
nonsyn_3$Difference <- as.numeric(nonsyn_3$Difference)
nonsyn_3 <- nonsyn_3 %>% filter(Difference != 0)
nonsyn_3$Name <- 'Participant 06'
Patient6df <- nonsyn_3

#patient 7
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
#nonsyn_3 <- nonsyn_3 %>% filter(CntNonSyn > 4) %>% filter(CntNonSyn2 > 4) # filter non syn counts at arbitrary value
nonsyn_3$Difference <-  nonsyn_3$Proportion2 - nonsyn_3$Proportion1 # work out difference and make numeric
nonsyn_3$Difference <- as.numeric(nonsyn_3$Difference)
nonsyn_3 <- nonsyn_3 %>% filter(Difference != 0)
nonsyn_3$Name <- 'Participant 07'
Patient7df <- nonsyn_3

#patient 9
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
#nonsyn_3 <- nonsyn_3 %>% filter(CntNonSyn > 4) %>% filter(CntNonSyn2 > 4) # filter non syn counts at arbitrary value
nonsyn_3$Difference <-  nonsyn_3$Proportion1 - nonsyn_3$Proportion2 # work out difference and make numeric
nonsyn_3$Difference <- as.numeric(nonsyn_3$Difference)
nonsyn_3 <- nonsyn_3 %>% filter(Difference != 0)
nonsyn_3$Name <- 'Participant 09 S1, S2'
Patient9s21df <- nonsyn_3

#patient 9 s3-s2
nonsyn_4 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3661_ATATGATGAA-ATGAAGGAGG_L002.final_AA_parse.txt", sep="\t")
nonsyn_4$AAposition <- 1:nrow(nonsyn_4)
nonsyn_4$Proportion4 <-  nonsyn_4$CntNonSyn / nonsyn_4$AAcoverage
nonsyn_4$Proportion4<-as.numeric(nonsyn_4$Proportion4)
names(nonsyn_4)[names(nonsyn_4) == "AAcoverage"] <- "AAcoverage3"
names(nonsyn_4)[names(nonsyn_4) == "CntNonSyn"] <- "CntNonSyn3"

nonsyn_5 <- cbind(nonsyn_2, nonsyn_4[, c("Proportion4", "AAcoverage3", "CntNonSyn3")]) # cbind proportion, cov and cnt to filter
nonsyn_5 <- nonsyn_5 %>% filter(AAcoverage2 > 20) %>% filter(AAcoverage3 > 20)
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
Patient9s32df <- nonsyn_5

#patient 11
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
#nonsyn_3 <- nonsyn_3 %>% filter(CntNonSyn > 4) %>% filter(CntNonSyn2 > 4) # filter non syn counts at arbitrary value
nonsyn_3$Difference <-  nonsyn_3$Proportion2 - nonsyn_3$Proportion1 # work out difference and make numeric
nonsyn_3$Difference <- as.numeric(nonsyn_3$Difference)
nonsyn_3 <- nonsyn_3 %>% filter(Difference != 0)
nonsyn_3$Name <- 'Participant 11'
Patient11df <- nonsyn_3

# patient 12
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
#nonsyn_3 <- nonsyn_3 %>% filter(CntNonSyn > 2) %>% filter(CntNonSyn2 > 2) # filter non syn counts at arbitrary value
nonsyn_3$Difference <-  nonsyn_3$Proportion2 - nonsyn_3$Proportion1 # work out difference and make numeric
nonsyn_3$Difference <- as.numeric(nonsyn_3$Difference)
nonsyn_3 <- nonsyn_3 %>% filter(Difference != 0)
nonsyn_3$Name <- 'Participant 12'
Patient12df <- nonsyn_3

#patient 15
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
#nonsyn_3 <- nonsyn_3 %>% filter(CntNonSyn > 4) %>% filter(CntNonSyn2 > 4) # filter non syn counts at arbitrary value
nonsyn_3$Difference <-  nonsyn_3$Proportion2 - nonsyn_3$Proportion1 # work out difference and make numeric
nonsyn_3$Difference <- as.numeric(nonsyn_3$Difference)
nonsyn_3 <- nonsyn_3 %>% filter(Difference != 0)
nonsyn_3$Name <- 'Participant 15'
Patient15df <- nonsyn_3

# patient 16
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
#nonsyn_3 <- nonsyn_3 %>% filter(CntNonSyn > 4) %>% filter(CntNonSyn2 > 4) # filter non syn counts at arbitrary value
nonsyn_3$Difference <-  nonsyn_3$Proportion2 - nonsyn_3$Proportion1 # work out difference and make numeric
nonsyn_3$Difference <- as.numeric(nonsyn_3$Difference)
nonsyn_3 <- nonsyn_3 %>% filter(Difference != 0)
nonsyn_3$Name <- 'Participant 16'
Patient16df <- nonsyn_3

#facet wrap and plot

#propdiff_rbind <- rbind(Patient1df, Patient2df, Patient3df, Patient5df, Patient6df, Patient7df, Patient9df, Patient11df, Patient12df, Patient15df, Patient16df)
propdiff_rbind_15N <- rbind(Patient1df, Patient2df, Patient3df, Patient5df, Patient6df, Patient7df, Patient9s21df, Patient9s32df, Patient12df, Patient16df)

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

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

melt_dfrbind <- reshape2::melt(propdiff_rbind_15N, id.vars = c('AAposition', 'Name', 'Protein'), measure.vars= 'Difference') 
dfplot <- ggplot(data= melt_dfrbind, aes(x=AAposition, y=value, colour=Protein)) + 
  geom_point(size=1.25) + ylab("Difference in proportion") + xlab("Amino acid position") + theme_bw() +
  scale_colour_manual(values= col_vector) +
  ylim(-1,1) +
  theme(legend.position = "bottom", axis.title.x = element_text(face ="bold", size=14), 
        legend.title = element_text(face = "bold", size=14), 
        axis.title.y = element_text(face="bold", size = 14))

yyy <- dfplot + facet_wrap(vars(Name))
# make a decision of whether to filter cntnonsyn 

ggsave(filename= "propdiff_colours_20X_participants.tiff", plot = yyy, device = 'tiff', dpi= 300, width=13, height=8)

