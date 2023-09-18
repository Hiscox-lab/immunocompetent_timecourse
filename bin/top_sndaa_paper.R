# top and snd aa for nimagen data
# read in AA.txt separately and add in AAposition then can rbind and do all analysis at once

library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyverse)
library(reshape2) 

tsaa_36 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3636_TAGTCGTCAA-ATGGCTCGGT_L002.final_AA.txt", sep="\t")
tsaa_36$AAposition <- 1:nrow(tsaa_36)
tsaa_36$Name <- 'P 01'
tsaa_36$Sample <- 'S1'
tsaa_37 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3637_TGGACTAATT-CAGCGGATAA_L002.final_AA.txt", sep="\t")
tsaa_37$AAposition <- 1:nrow(tsaa_37)
tsaa_37$Name <- 'P 01'
tsaa_37$Sample <- 'S2'
tsaa_39 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3639_ATCTTATGAT-GCCTGGCAGT_L002.final_AA.txt", sep="\t")
tsaa_39$AAposition <- 1:nrow(tsaa_39)
tsaa_39$Name <- 'P 02'
tsaa_39$Sample <- 'S1'
tsaa_40 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3640_ACCAACGCTG-TAACTTGGAG_L002.final_AA.txt", sep="\t")
tsaa_40$AAposition <- 1:nrow(tsaa_40)
tsaa_40$Name <- 'P 02'
tsaa_40$Sample <- 'S2'
tsaa_43 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3643_TACGGCGCGA-ACGTATGCGC_L002.final_AA.txt", sep="\t")
tsaa_43$AAposition <- 1:nrow(tsaa_43)
tsaa_43$Name <- 'P 03'
tsaa_43$Sample <- 'S1'
tsaa_44 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3644_AGCTATCAAC-TGATTGCTTC_L002.final_AA.txt", sep="\t")
tsaa_44$AAposition <- 1:nrow(tsaa_44)
tsaa_44$Name <- 'P 03'
tsaa_44$Sample <- 'S2'
tsaa_47 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3647_AGAGGTCATT-GGCATATATA_L002.final_AA.txt", sep="\t")
tsaa_47$AAposition <- 1:nrow(tsaa_47)
tsaa_47$Name <- 'P 04'
tsaa_47$Sample <- 'S1'
tsaa_49 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3649_CTTAGAGGCA-AGCGTCTGGT_L002.final_AA.txt", sep="\t")
tsaa_49$AAposition <- 1:nrow(tsaa_49)
tsaa_49$Name <- 'P 05'
tsaa_49$Sample <- 'S1'
tsaa_50 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3650_ATAGAGCATT-GAGGACCGAT_L002.final_AA.txt", sep="\t")
tsaa_50$AAposition <- 1:nrow(tsaa_50)
tsaa_50$Name <- 'P 05'
tsaa_50$Sample <- 'S2'
tsaa_51 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3651_GGTAGCGCAT-TACTGGATAA_L002.final_AA.txt", sep="\t")
tsaa_51$AAposition <- 1:nrow(tsaa_51)
tsaa_51$Name <- 'P 06'
tsaa_51$Sample <- 'S1'
tsaa_52 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3652_TAATGATATT-CGAGTCGAAG_L002.final_AA.txt", sep="\t")
tsaa_52$AAposition <- 1:nrow(tsaa_52)
tsaa_52$Name <- 'P 06'
tsaa_52$Sample <- 'S2'
tsaa_53 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3653_GGCTAAGGCG-TTCTTCTATT_L002.final_AA.txt", sep="\t")
tsaa_53$AAposition <- 1:nrow(tsaa_53)
tsaa_53$Name <- 'P 07'
tsaa_53$Sample <- 'S1'
tsaa_54 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3654_TTAGTACCTT-GACCAATTCT_L002.final_AA.txt", sep="\t")
tsaa_54$AAposition <- 1:nrow(tsaa_54)
tsaa_54$Name <- 'P 07'
tsaa_54$Sample <- 'S2'
tsaa_59 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3659_GTAGCTACCG-TCAGAACGAC_L002.final_AA.txt", sep="\t")
tsaa_59$AAposition <- 1:nrow(tsaa_59)
tsaa_59$Name <- 'P 09'
tsaa_59$Sample <- 'S1'
tsaa_60 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3660_ATTAATAGAA-GGCAGAGGAG_L002.final_AA.txt", sep="\t")
tsaa_60$AAposition <- 1:nrow(tsaa_60)
tsaa_60$Name <- 'P 09'
tsaa_60$Sample <- 'S2'
tsaa_61 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3661_ATATGATGAA-ATGAAGGAGG_L002.final_AA.txt", sep="\t")
tsaa_61$AAposition <- 1:nrow(tsaa_61)
tsaa_61$Name <- 'P 09'
tsaa_61$Sample <- 'S3'
tsaa_62 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3662_AATGGACGTA-GAAGGTTCGG_L002.final_AA.txt", sep="\t")
tsaa_62$AAposition <- 1:nrow(tsaa_62)
tsaa_62$Name <- 'P 10'
tsaa_62$Sample <- 'S1'
tsaa_64 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3664_ACCATAAGCG-TGAATAACTA_L002.final_AA.txt", sep="\t")
tsaa_64$AAposition <- 1:nrow(tsaa_64)
tsaa_64$Name <- 'P 11'
tsaa_64$Sample <- 'S2'
tsaa_65 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3665_CCGCTGCATT-GAGTCCTGCA_L002.final_AA.txt", sep="\t")
tsaa_65$AAposition <- 1:nrow(tsaa_65)
tsaa_65$Name <- 'P 11'
tsaa_65$Sample <- 'S3'
tsaa_66 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3666_TTGGTTGGTA-GCCATAGTAG_L002.final_AA.txt", sep="\t")
tsaa_66$AAposition <- 1:nrow(tsaa_66)
tsaa_66$Name <- 'P 11'
tsaa_66$Sample <- 'S4'
tsaa_67 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3667_GGAACGATTC-CGGCGAGTCC_L002.final_AA.txt", sep="\t")
tsaa_67$AAposition <- 1:nrow(tsaa_67)
tsaa_67$Name <- 'P 12'
tsaa_67$Sample <- 'S1'
tsaa_68 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3668_AAGCATCTCG-CTCATGATTG_L002.final_AA.txt", sep="\t")
tsaa_68$AAposition <- 1:nrow(tsaa_68)
tsaa_68$Name <- 'P 12'
tsaa_68$Sample <- 'S2'
tsaa_70 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3670_TATAGATGGC-GCGAGAGAAG_L002.final_AA.txt", sep="\t")
tsaa_70$AAposition <- 1:nrow(tsaa_70)
tsaa_70$Name <- 'P 13'
tsaa_70$Sample <- 'S1'
tsaa_74 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3674_TCCGTATTCC-TCTAATTATT_L002.final_AA.txt", sep="\t")
tsaa_74$AAposition <- 1:nrow(tsaa_74)
tsaa_74$Name <- 'P 14'
tsaa_74$Sample <- 'S1'
tsaa_78 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3678_ATACTATATT-GGTATCATCT_L002.final_AA.txt", sep="\t")
tsaa_78$AAposition <- 1:nrow(tsaa_78)
tsaa_78$Name <- 'P 15'
tsaa_78$Sample <- 'S1'
tsaa_79 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3679_TAGTTATGCG-CTGGCTTAGT_L002.final_AA.txt", sep="\t")
tsaa_79$AAposition <- 1:nrow(tsaa_79)
tsaa_79$Name <- 'P 15'
tsaa_79$Sample <- 'S2'
tsaa_80 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3680_AACTTACGCT-TAACCTGCAA_L002.final_AA.txt", sep="\t")
tsaa_80$AAposition <- 1:nrow(tsaa_80)
tsaa_80$Name <- 'P 15'
tsaa_80$Sample <- 'S3'
tsaa_82 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3682_GTCCGTAAGG-CGACCTATAC_L002.final_AA.txt", sep="\t")
tsaa_82$AAposition <- 1:nrow(tsaa_82)
tsaa_82$Name <- 'P 16'
tsaa_82$Sample <- 'S1'
tsaa_83 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3683_CCGGTCGACG-GCTGCATAGC_L002.final_AA.txt", sep="\t")
tsaa_83$AAposition <- 1:nrow(tsaa_83)
tsaa_83$Name <- 'P 16'
tsaa_83$Sample <- 'S2'

rbind_nim_tsaa <- rbind(tsaa_36, tsaa_37, tsaa_39, tsaa_40, tsaa_43, tsaa_44, tsaa_47, tsaa_49, tsaa_50, tsaa_51, tsaa_52, tsaa_53, tsaa_54, tsaa_59, tsaa_60, tsaa_61, tsaa_62, tsaa_64, tsaa_65, tsaa_66, tsaa_67, tsaa_68, tsaa_70, tsaa_74, tsaa_78, tsaa_79, tsaa_80, tsaa_82, tsaa_83)
rbind_nim_tsaa_15N <- rbind(tsaa_36, tsaa_37, tsaa_39, tsaa_40, tsaa_43, tsaa_44, tsaa_47, tsaa_49, tsaa_50, tsaa_51, tsaa_52, tsaa_53, tsaa_54, tsaa_59, tsaa_60, tsaa_61, tsaa_62, tsaa_65, tsaa_67, tsaa_68, tsaa_70, tsaa_74, tsaa_79, tsaa_82, tsaa_83)

#for all 15N
rbind_nim_tsaa_15N$AAcoverage<-as.numeric(rbind_nim_tsaa_15N$AAcoverage) # make all the values you might want to plot as numeric
rbind_nim_tsaa_15N$TopAAcnt <- as.numeric(rbind_nim_tsaa_15N$TopAAcnt)
rbind_nim_tsaa_15N$SndAAcnt <- as.numeric(rbind_nim_tsaa_15N$SndAAcnt)
rbind_nim_tsaa_15N$TrdAAcnt <- as.numeric(rbind_nim_tsaa_15N$TrdAAcnt)
rbind_nim_tsaa_15N$ProportionTopAA <-  rbind_nim_tsaa_15N$TopAAcnt / rbind_nim_tsaa_15N$AAcoverage # work out proportions using coverage
rbind_nim_tsaa_15N$Proportion2ndAA <-  rbind_nim_tsaa_15N$SndAAcnt / rbind_nim_tsaa_15N$AAcoverage
rbind_nim_tsaa_15N$Proportion3rdAA <-  rbind_nim_tsaa_15N$TrdAAcnt / rbind_nim_tsaa_15N$AAcoverage
rbind_nim_tsaa_15N$ProportionTopAA <- as.numeric(rbind_nim_tsaa_15N$ProportionTopAA) # make proportions numeric
rbind_nim_tsaa_15N$Proportion2ndAA <- as.numeric(rbind_nim_tsaa_15N$Proportion2ndAA)
rbind_nim_tsaa_15N$Proportion3rdAA <- as.numeric(rbind_nim_tsaa_15N$Proportion3rdAA)
rbind_nim_tsaa_15N$AAposition <- as.numeric(rbind_nim_tsaa_15N$AAposition) # make positions numeric
rbind_nim_tsaa_15N$AAPosition <- as.numeric(rbind_nim_tsaa_15N$AAPosition) 
rbind_nim_tsaa_15N <- rbind_nim_tsaa_15N %>% filter(Sample != "Sample") %>% filter(AAcoverage > 20) # filter to get variation from refAA and for coverage of 20
cols <- c("RefAA", "AAPosition", "TopAA") # make vector of column names to make mutant column
rbind_nim_tsaa_15N$mutant <- do.call(paste, c(rbind_nim_tsaa_15N[cols])) # add mutant column to dataframe
cols2AA <- c("RefAA", "AAPosition", "SndAA") # repeat for SndAA
rbind_nim_tsaa_15N$MutantSndAA <- do.call(paste, c(rbind_nim_tsaa_15N[cols2AA]))
cols3AA <- c("RefAA", "AAPosition", "TrdAA") # repeat for trdAA
rbind_nim_tsaa_15N$MutantTrdAA <- do.call(paste, c(rbind_nim_tsaa_15N[cols3AA]))

#write.csv(rbind_nim_tsaa_15N,"./proportions_aa_nim_15N.csv", row.names = FALSE) # print csv of results of prop AAs

melt_rbindall <- reshape2::melt(rbind_nim_tsaa_15N, id.vars = c('Name', 'Sample', 'AAposition'), measure.vars=  c('ProportionTopAA', 'Proportion2ndAA'))
melt_rbindall$variable <- as.factor(melt_rbindall$variable) # change to allow plotting of graph

tsaa_all <- ggplot(data= melt_rbindall, aes(x=AAposition, y=value)) + geom_point(aes(colour=factor(variable)), size=1) + 
  xlim(0,9755) + ylim(0,1) + ylab("Proportion") + xlab("Amino acid position") +
  theme(text = element_text(size = 11)) + theme(legend.position = "bottom") + scale_colour_discrete(name= "Mutant", labels = c("Top amino acid count", "Second amino acid count")) 
 
x <- tsaa_all + facet_grid(Sample ~ Name) +theme_bw()# looks better for single participant
xx <- tsaa_all + facet_grid(Name ~ Sample) # better for visualising many participants
ggsave(filename= "topsndaa_nim_participants_filt20.tiff", plot = xx, device = 'tiff', dpi= 300, width=7, height=10)





