# heatmap to compare syn and non syn mutations across proteins

library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyverse) 
library(reshape2)
library(RColorBrewer)
library("pheatmap")


### SYN NONSYN HEATMAP ###
# use cat parse.txt file 

parse_all_nim <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/nimnan25_parse.txt", sep="\t")
parse15N_nim <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/nim15%_parse.txt", sep="\t")

parse15N_nim$CntNonSyn<-as.numeric(parse15N_nim$CntNonSyn)
parse15N_nim$CntSyn<-as.numeric(parse15N_nim$CntSyn)
parse15N_nim$AAcoverage<-as.numeric(parse15N_nim$AAcoverage)
parse15N_nim <- parse15N_nim %>% filter(AAcoverage > 20)
parse15N_nim$prop_syn <- parse15N_nim$CntSyn / parse15N_nim$AAcoverage
parse15N_nim$prop_nonsyn <- parse15N_nim$CntNonSyn / parse15N_nim$AAcoverage
parse15N_nim$Protein[parse15N_nim$Protein == "nsp12_1"] <- "NSP12"
parse15N_nim$Protein[parse15N_nim$Protein == "nsp12_2"] <- "NSP12"
parse15N_nim$Protein[parse15N_nim$Protein == "nsp1"] <- "NSP1"
parse15N_nim$Protein[parse15N_nim$Protein == "nsp2"] <- "NSP2"
parse15N_nim$Protein[parse15N_nim$Protein == "nsp3"] <- "NSP3"
parse15N_nim$Protein[parse15N_nim$Protein == "nsp4"] <- "NSP4"
parse15N_nim$Protein[parse15N_nim$Protein == "nsp5"] <- "NSP5"
parse15N_nim$Protein[parse15N_nim$Protein == "nsp6"] <- "NSP6"
parse15N_nim$Protein[parse15N_nim$Protein == "nsp7"] <- "NSP7"
parse15N_nim$Protein[parse15N_nim$Protein == "nsp8"] <- "NSP8"
parse15N_nim$Protein[parse15N_nim$Protein == "nsp9"] <- "NSP9"
parse15N_nim$Protein[parse15N_nim$Protein == "nsp10"] <- "NSP10"
parse15N_nim$Protein[parse15N_nim$Protein == "nsp13"] <- "NSP13"
parse15N_nim$Protein[parse15N_nim$Protein == "nsp14"] <- "NSP14"
parse15N_nim$Protein[parse15N_nim$Protein == "nsp15"] <- "NSP15"
parse15N_nim$Protein[parse15N_nim$Protein == "nsp16"] <- "NSP16"

# aim to calculate mean of each per protein 

Mean_syn <- as.data.frame(parse15N_nim %>% group_by(Protein, Sample) %>% summarize(Mean_syn=mean(prop_syn, na.rm=T))) 
Mean_nonsyn <- as.data.frame(parse15N_nim %>% group_by(Protein, Sample) %>% summarize(Mean_nonsyn=mean(prop_nonsyn, na.rm=T))) 

Mean_nonsyn$Sample[Mean_nonsyn$Sample == "3636_TAGTCGTCAA-ATGGCTCGGT_L002.final"] <- "01 S1"
Mean_nonsyn$Sample[Mean_nonsyn$Sample == "3637_TGGACTAATT-CAGCGGATAA_L002.final"] <- "01 S2"
Mean_nonsyn$Sample[Mean_nonsyn$Sample == "3639_ATCTTATGAT-GCCTGGCAGT_L002.final"] <- "02 S1"
Mean_nonsyn$Sample[Mean_nonsyn$Sample == "3640_ACCAACGCTG-TAACTTGGAG_L002.final"] <- "02 S2"
Mean_nonsyn$Sample[Mean_nonsyn$Sample == "3643_TACGGCGCGA-ACGTATGCGC_L002.final"] <- "03 S1"
Mean_nonsyn$Sample[Mean_nonsyn$Sample == "3644_AGCTATCAAC-TGATTGCTTC_L002.final"] <- "03 S2"
Mean_nonsyn$Sample[Mean_nonsyn$Sample == "3647_AGAGGTCATT-GGCATATATA_L002.final"] <- "04 S1"
Mean_nonsyn$Sample[Mean_nonsyn$Sample == "3649_CTTAGAGGCA-AGCGTCTGGT_L002.final"] <- "05 S1"
Mean_nonsyn$Sample[Mean_nonsyn$Sample == "3650_ATAGAGCATT-GAGGACCGAT_L002.final"] <- "05 S2"
Mean_nonsyn$Sample[Mean_nonsyn$Sample == "3651_GGTAGCGCAT-TACTGGATAA_L002.final"] <- "06 S1"
Mean_nonsyn$Sample[Mean_nonsyn$Sample == "3652_TAATGATATT-CGAGTCGAAG_L002.final"] <- "06 S2"
Mean_nonsyn$Sample[Mean_nonsyn$Sample == "3653_GGCTAAGGCG-TTCTTCTATT_L002.final"] <- "07 S1"
Mean_nonsyn$Sample[Mean_nonsyn$Sample == "3654_TTAGTACCTT-GACCAATTCT_L002.final"] <- "07 S2"
Mean_nonsyn$Sample[Mean_nonsyn$Sample == "3659_GTAGCTACCG-TCAGAACGAC_L002.final"] <- "09 S1"
Mean_nonsyn$Sample[Mean_nonsyn$Sample == "3660_ATTAATAGAA-GGCAGAGGAG_L002.final"] <- "09 S2"
Mean_nonsyn$Sample[Mean_nonsyn$Sample == "3661_ATATGATGAA-ATGAAGGAGG_L002.final"] <- "09 S3"
Mean_nonsyn$Sample[Mean_nonsyn$Sample == "3662_AATGGACGTA-GAAGGTTCGG_L002.final"] <- "10 S1"
#Mean_nonsyn$Sample[Mean_nonsyn$Sample == "3664_ACCATAAGCG-TGAATAACTA_L002.final"] <- "11 S2"
Mean_nonsyn$Sample[Mean_nonsyn$Sample == "3665_CCGCTGCATT-GAGTCCTGCA_L002.final"] <- "11 S3"
#Mean_nonsyn$Sample[Mean_nonsyn$Sample == "3666_TTGGTTGGTA-GCCATAGTAG_L002.final"] <- "11 S4"
Mean_nonsyn$Sample[Mean_nonsyn$Sample == "3667_GGAACGATTC-CGGCGAGTCC_L002.final"] <- "12 S1"
Mean_nonsyn$Sample[Mean_nonsyn$Sample == "3668_AAGCATCTCG-CTCATGATTG_L002.final"] <- "12 S2"
Mean_nonsyn$Sample[Mean_nonsyn$Sample == "3670_TATAGATGGC-GCGAGAGAAG_L002.final"] <- "13 S1"
Mean_nonsyn$Sample[Mean_nonsyn$Sample == "3674_TCCGTATTCC-TCTAATTATT_L002.final"] <- "14 S1"
#Mean_nonsyn$Sample[Mean_nonsyn$Sample == "3678_ATACTATATT-GGTATCATCT_L002.final"] <- "15 S1"
Mean_nonsyn$Sample[Mean_nonsyn$Sample == "3679_TAGTTATGCG-CTGGCTTAGT_L002.final"] <- "15 S2"
#Mean_nonsyn$Sample[Mean_nonsyn$Sample == "3680_AACTTACGCT-TAACCTGCAA_L002.final"] <- "15 S3"
Mean_nonsyn$Sample[Mean_nonsyn$Sample == "3682_GTCCGTAAGG-CGACCTATAC_L002.final"] <- "16 S1"
Mean_nonsyn$Sample[Mean_nonsyn$Sample == "3683_CCGGTCGACG-GCTGCATAGC_L002.final"] <- "16 S2"

Mean_syn$Sample[Mean_syn$Sample == "3636_TAGTCGTCAA-ATGGCTCGGT_L002.final"] <- "01 S1"
Mean_syn$Sample[Mean_syn$Sample == "3637_TGGACTAATT-CAGCGGATAA_L002.final"] <- "01 S2"
Mean_syn$Sample[Mean_syn$Sample == "3639_ATCTTATGAT-GCCTGGCAGT_L002.final"] <- "02 S1"
Mean_syn$Sample[Mean_syn$Sample == "3640_ACCAACGCTG-TAACTTGGAG_L002.final"] <- "02 S2"
Mean_syn$Sample[Mean_syn$Sample == "3643_TACGGCGCGA-ACGTATGCGC_L002.final"] <- "03 S1"
Mean_syn$Sample[Mean_syn$Sample == "3644_AGCTATCAAC-TGATTGCTTC_L002.final"] <- "03 S2"
Mean_syn$Sample[Mean_syn$Sample == "3647_AGAGGTCATT-GGCATATATA_L002.final"] <- "04 S1"
Mean_syn$Sample[Mean_syn$Sample == "3649_CTTAGAGGCA-AGCGTCTGGT_L002.final"] <- "05 S1"
Mean_syn$Sample[Mean_syn$Sample == "3650_ATAGAGCATT-GAGGACCGAT_L002.final"] <- "05 S2"
Mean_syn$Sample[Mean_syn$Sample == "3651_GGTAGCGCAT-TACTGGATAA_L002.final"] <- "06 S1"
Mean_syn$Sample[Mean_syn$Sample == "3652_TAATGATATT-CGAGTCGAAG_L002.final"] <- "06 S2"
Mean_syn$Sample[Mean_syn$Sample == "3653_GGCTAAGGCG-TTCTTCTATT_L002.final"] <- "07 S1"
Mean_syn$Sample[Mean_syn$Sample == "3654_TTAGTACCTT-GACCAATTCT_L002.final"] <- "07 S2"
Mean_syn$Sample[Mean_syn$Sample == "3659_GTAGCTACCG-TCAGAACGAC_L002.final"] <- "09 S1"
Mean_syn$Sample[Mean_syn$Sample == "3660_ATTAATAGAA-GGCAGAGGAG_L002.final"] <- "09 S2"
Mean_syn$Sample[Mean_syn$Sample == "3661_ATATGATGAA-ATGAAGGAGG_L002.final"] <- "09 S3"
Mean_syn$Sample[Mean_syn$Sample == "3662_AATGGACGTA-GAAGGTTCGG_L002.final"] <- "10 S1"
#Mean_syn$Sample[Mean_syn$Sample == "3664_ACCATAAGCG-TGAATAACTA_L002.final"] <- "Patient 11 S2"
Mean_syn$Sample[Mean_syn$Sample == "3665_CCGCTGCATT-GAGTCCTGCA_L002.final"] <- "11 S3"
#Mean_syn$Sample[Mean_syn$Sample == "3666_TTGGTTGGTA-GCCATAGTAG_L002.final"] <- "Patient 11 S4"
Mean_syn$Sample[Mean_syn$Sample == "3667_GGAACGATTC-CGGCGAGTCC_L002.final"] <- "12 S1"
Mean_syn$Sample[Mean_syn$Sample == "3668_AAGCATCTCG-CTCATGATTG_L002.final"] <- "12 S2"
Mean_syn$Sample[Mean_syn$Sample == "3670_TATAGATGGC-GCGAGAGAAG_L002.final"] <- "13 S1"
Mean_syn$Sample[Mean_syn$Sample == "3674_TCCGTATTCC-TCTAATTATT_L002.final"] <- "14 S1"
#Mean_syn$Sample[Mean_syn$Sample == "3678_ATACTATATT-GGTATCATCT_L002.final"] <- "Patient 15 S1"
Mean_syn$Sample[Mean_syn$Sample == "3679_TAGTTATGCG-CTGGCTTAGT_L002.final"] <- "15 S2"
#Mean_syn$Sample[Mean_syn$Sample == "3680_AACTTACGCT-TAACCTGCAA_L002.final"] <- "Patient 15 S3"
Mean_syn$Sample[Mean_syn$Sample == "3682_GTCCGTAAGG-CGACCTATAC_L002.final"] <- "16 S1"
Mean_syn$Sample[Mean_syn$Sample == "3683_CCGGTCGACG-GCTGCATAGC_L002.final"] <- "16 S2"

Mean_syn$Protein <- as.factor(Mean_syn$Protein)
Mean_syn$Protein <- factor(Mean_syn$Protein, levels= c("5' UTR", "NSP1", "NSP2", "NSP3", "NSP4", "NSP5", "NSP6", "NSP7", "NSP8", "NSP9", "NSP10", "NSP12", "NSP13", "NSP14", "NSP15", "NSP16", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10", "3' UTR"))
Mean_nonsyn$Protein <- as.factor(Mean_nonsyn$Protein)
Mean_nonsyn$Protein <- factor(Mean_nonsyn$Protein, levels= c("5' UTR", "NSP1", "NSP2", "NSP3", "NSP4", "NSP5", "NSP6", "NSP7", "NSP8", "NSP9", "NSP10", "NSP12", "NSP13", "NSP14", "NSP15", "NSP16", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10", "3' UTR"))


## dN/dS ratio using tstv ratio script dNonsyn / dSyn
# use parse15N_nim

metadata_dstl <- read.csv("/home/hannahg/projects/dstl_project/data/nimagen/nimagen_metadata_bam.csv") # metadata sheet with sample info read in
parse15N_nim <- merge(metadata_dstl, parse15N_nim, by = "Sample")
parse15N_nim$Protein <- factor(parse15N_nim$Protein, levels = c("5'UTR", "NSP1", "NSP2", "NSP3", "NSP4", "NSP5", "NSP6", "NSP7", "NSP8", "NSP9", "NSP10", "NSP12", "NSP13", "NSP14", "NSP15", "NSP16", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10", "3'UTR"))

protlengths <- parse15N_nim %>% group_by(Sample, Protein) %>% count(AAPosition) %>% tally
protlengths<- rename(protlengths, ProtLength=n)
parse15N_nim<- merge(parse15N_nim, protlengths, by=c("Sample","Protein"))

parse15N_nim$CntNonSyn<-as.numeric(parse15N_nim$CntNonSyn)
parse15N_nim$CntSyn<-as.numeric(parse15N_nim$CntSyn)
parse15N_nim$AAcoverage<-as.numeric(parse15N_nim$AAcoverage)
parse15N_nim <- parse15N_nim %>% filter(AAcoverage > 20)
parse15N_nim <- parse15N_nim %>% filter(Protein != "NA")
parse15N_nim$prop_nonsyn <- parse15N_nim$CntNonSyn / parse15N_nim$AAcoverage
parse15N_nim$prop_syn <- parse15N_nim$CntSyn / parse15N_nim$AAcoverage

## Iahs mean tstv per protein 
mean_nonsyn_perprot <- parse15N_nim %>% filter(AAcoverage >20) %>% group_by(Protein, Sample, ProtLength) %>% 
  summarize(Meannonsyn=mean(prop_nonsyn, na.rm=T ), stdev=sd(prop_nonsyn, na.rm=T), across())
mean_nonsyn_perprot$Meannonsyn_pnorm <- mean_nonsyn_perprot$Meannonsyn / mean_nonsyn_perprot$ProtLength
#mean_ts_perprot$Meants_pnorm_scale <- scale(mean_ts_perprot$Meants_pnorm)
#mean_ts_perprot$MeantsNorm <- scale(mean_ts_perprot$Meants)

mean_syn_perprot <- parse15N_nim %>% filter(AAcoverage >20) %>% group_by(Protein, Sample, ProtLength) %>% 
  summarize(Meansyn=mean(prop_syn, na.rm=T ), stdev=sd(prop_syn, na.rm=T), across())
mean_syn_perprot$Meansyn_pnorm <- mean_syn_perprot$Meansyn / mean_syn_perprot$ProtLength
#mean_tv_perprot$Meantv_pnorm_scale <- scale(mean_tv_perprot$Meantv_pnorm)
#mean_tv_perprot$MeantvNorm <- scale(mean_tv_perprot$Meantv)

meansns_perprot <- merge(mean_nonsyn_perprot, mean_syn_perprot, by = c("Sample", "Protein", "ProtLength", "AAPosition", "Timepoint", "Participant"))

meansns_perprot$snsratio_norm <- meansns_perprot$MeannonsynNorm / meansns_perprot$MeansynNorm
meansns_perprot$snsratio_norm <- meansns_perprot$Meannonsyn_pnorm / meansns_perprot$Meansyn_pnorm # for unscaled
names(meansns_perprot)[names(meansns_perprot) == "snsratio_norm"] <- "Ratio"

## use this one ?? ##
dnds<- ggplot(meansns_perprot, aes(Protein, Participant, fill= Ratio)) + 
  geom_tile(colour="black") + scale_fill_gradient(low= "white", high= "blue") + ylab("Participant") + xlab("Protein") +
  ggtitle("dN/dS Ratio") + theme_bw() + theme(axis.text.x = element_text(angle=90, vjust= 0.7)) + facet_wrap(~Timepoint) 

ggsave(filename= "dndsratio_nim_hm_20X.tiff", plot= dnds, device='tiff', dpi=300, width=11.5, height=7)



### INDELS HEATMAP ###
# getting position on entropy file for nsps by reading in individually and then rbinding them all at the end

nsps <- read.csv("/home/hannahg/projects/dstl_project/data/below_15%_N/protein_position.csv")

ent_36 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3636_TAGTCGTCAA-ATGGCTCGGT_L002.final_entropy.txt", sep="\t")
ent_36_prot <- cbind(ent_36, nsps["Protein"])
ent_37 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3637_TGGACTAATT-CAGCGGATAA_L002.final_entropy.txt", sep="\t")
ent_37_prot <- cbind(ent_37, nsps["Protein"])
ent_39 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3639_ATCTTATGAT-GCCTGGCAGT_L002.final_entropy.txt", sep="\t")
ent_39_prot <- cbind(ent_39, nsps["Protein"])
ent_40 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3640_ACCAACGCTG-TAACTTGGAG_L002.final_entropy.txt", sep="\t")
ent_40_prot <- cbind(ent_40, nsps["Protein"])
ent_43 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3643_TACGGCGCGA-ACGTATGCGC_L002.final_entropy.txt", sep="\t")
ent_43_prot <- cbind(ent_43, nsps["Protein"])
ent_44 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3644_AGCTATCAAC-TGATTGCTTC_L002.final_entropy.txt", sep="\t")
ent_44_prot <- cbind(ent_44, nsps["Protein"])
ent_47 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3647_AGAGGTCATT-GGCATATATA_L002.final_entropy.txt", sep="\t")
ent_47_prot <- cbind(ent_47, nsps["Protein"])
ent_49 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3649_CTTAGAGGCA-AGCGTCTGGT_L002.final_entropy.txt", sep="\t")
ent_49_prot <- cbind(ent_49, nsps["Protein"])
ent_50 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3650_ATAGAGCATT-GAGGACCGAT_L002.final_entropy.txt", sep="\t")
ent_50_prot <- cbind(ent_50, nsps["Protein"])
ent_51 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3651_GGTAGCGCAT-TACTGGATAA_L002.final_entropy.txt", sep="\t")
ent_51_prot <- cbind(ent_51, nsps["Protein"])
ent_52 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3652_TAATGATATT-CGAGTCGAAG_L002.final_entropy.txt", sep="\t")
ent_52_prot <- cbind(ent_52, nsps["Protein"])
ent_53 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3653_GGCTAAGGCG-TTCTTCTATT_L002.final_entropy.txt", sep="\t")
ent_53_prot <- cbind(ent_53, nsps["Protein"])
ent_54 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3654_TTAGTACCTT-GACCAATTCT_L002.final_entropy.txt", sep="\t")
ent_54_prot <- cbind(ent_54, nsps["Protein"])
ent_59 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3659_GTAGCTACCG-TCAGAACGAC_L002.final_entropy.txt", sep="\t")
ent_59_prot <- cbind(ent_59, nsps["Protein"])
ent_60 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3660_ATTAATAGAA-GGCAGAGGAG_L002.final_entropy.txt", sep="\t")
ent_60_prot <- cbind(ent_60, nsps["Protein"])
ent_61 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3661_ATATGATGAA-ATGAAGGAGG_L002.final_entropy.txt", sep="\t")
ent_61_prot <- cbind(ent_61, nsps["Protein"])
ent_62 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3662_AATGGACGTA-GAAGGTTCGG_L002.final_entropy.txt", sep="\t")
ent_62_prot <- cbind(ent_62, nsps["Protein"])
ent_64 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3664_ACCATAAGCG-TGAATAACTA_L002.final_entropy.txt", sep="\t")
ent_64_prot <- cbind(ent_64, nsps["Protein"])
ent_65 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3665_CCGCTGCATT-GAGTCCTGCA_L002.final_entropy.txt", sep="\t")
ent_65_prot <- cbind(ent_65, nsps["Protein"])
ent_66 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3666_TTGGTTGGTA-GCCATAGTAG_L002.final_entropy.txt", sep="\t")
ent_66_prot <- cbind(ent_66, nsps["Protein"])
ent_67 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3667_GGAACGATTC-CGGCGAGTCC_L002.final_entropy.txt", sep="\t")
ent_67_prot <- cbind(ent_67, nsps["Protein"])
ent_68 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3668_AAGCATCTCG-CTCATGATTG_L002.final_entropy.txt", sep="\t")
ent_68_prot <- cbind(ent_68, nsps["Protein"])
ent_70 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3670_TATAGATGGC-GCGAGAGAAG_L002.final_entropy.txt", sep="\t")
ent_70_prot <- cbind(ent_70, nsps["Protein"])
ent_74 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3674_TCCGTATTCC-TCTAATTATT_L002.final_entropy.txt", sep="\t")
ent_74_prot <- cbind(ent_74, nsps["Protein"])
ent_78 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3678_ATACTATATT-GGTATCATCT_L002.final_entropy.txt", sep="\t")
ent_78_prot <- cbind(ent_78, nsps["Protein"])
ent_79 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3679_TAGTTATGCG-CTGGCTTAGT_L002.final_entropy.txt", sep="\t")
ent_79_prot <- cbind(ent_79, nsps["Protein"])
ent_80 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3680_AACTTACGCT-TAACCTGCAA_L002.final_entropy.txt", sep="\t")
ent_80_prot <- cbind(ent_80, nsps["Protein"])
ent_82 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3682_GTCCGTAAGG-CGACCTATAC_L002.final_entropy.txt", sep="\t")
ent_82_prot <- cbind(ent_82, nsps["Protein"])
ent_83 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3683_CCGGTCGACG-GCTGCATAGC_L002.final_entropy.txt", sep="\t")
ent_83_prot <- cbind(ent_83, nsps["Protein"])

all_ent <- rbind(ent_36_prot, ent_37_prot, ent_39_prot, ent_40_prot, ent_43_prot, ent_44_prot, ent_47_prot, ent_49_prot, ent_50_prot, ent_51_prot, ent_52_prot, ent_53_prot, ent_54_prot, ent_59_prot, ent_60_prot, ent_61_prot, ent_62_prot, ent_64_prot, ent_65_prot, ent_66_prot, ent_67_prot, ent_68_prot, ent_70_prot, ent_74_prot, ent_78_prot, ent_79_prot, ent_80_prot, ent_82_prot, ent_83_prot)
nim15N_ent <- rbind(ent_36_prot, ent_37_prot, ent_39_prot, ent_40_prot, ent_43_prot, ent_44_prot, ent_47_prot, ent_49_prot, ent_50_prot, ent_51_prot, ent_52_prot, ent_53_prot, ent_54_prot, ent_59_prot, ent_60_prot, ent_61_prot, ent_62_prot, ent_65_prot, ent_67_prot, ent_68_prot, ent_70_prot, ent_74_prot, ent_79_prot, ent_82_prot, ent_83_prot)

nim15N_ent$Protein[nim15N_ent$Protein == "nsp12"] <- "NSP12"
nim15N_ent$Protein[nim15N_ent$Protein == "nsp1"] <- "NSP1"
nim15N_ent$Protein[nim15N_ent$Protein == "nsp2"] <- "NSP2"
nim15N_ent$Protein[nim15N_ent$Protein == "nsp3"] <- "NSP3"
nim15N_ent$Protein[nim15N_ent$Protein == "nsp4"] <- "NSP4"
nim15N_ent$Protein[nim15N_ent$Protein == "nsp5"] <- "NSP5"
nim15N_ent$Protein[nim15N_ent$Protein == "nsp6"] <- "NSP6"
nim15N_ent$Protein[nim15N_ent$Protein == "nsp7"] <- "NSP7"
nim15N_ent$Protein[nim15N_ent$Protein == "nsp8"] <- "NSP8"
nim15N_ent$Protein[nim15N_ent$Protein == "nsp9"] <- "NSP9"
nim15N_ent$Protein[nim15N_ent$Protein == "nsp10"] <- "NSP10"
nim15N_ent$Protein[nim15N_ent$Protein == "nsp13"] <- "NSP13"
nim15N_ent$Protein[nim15N_ent$Protein == "nsp14"] <- "NSP14"
nim15N_ent$Protein[nim15N_ent$Protein == "nsp15"] <- "NSP15"
nim15N_ent$Protein[nim15N_ent$Protein == "nsp16"] <- "NSP16"


### TSTV ### 
# use all_ent or nim15N_ent from above 

### for nim15N_ent creating mean per protein of tstv 25/10/22 ###

metadata_dstl <- read.csv("/home/hannahg/projects/dstl_project/data/nimagen/nimagen_metadata.csv") # metadata sheet with sample info read in
nim15N_ent <- merge(metadata_dstl, nim15N_ent, by = "Sample")
nim15N_ent$Protein <- factor(nim15N_ent$Protein, levels = c("5'UTR", "NSP1", "NSP2", "NSP3", "NSP4", "NSP5", "NSP6", "NSP7", "NSP8", "NSP9", "NSP10", "NSP12", "NSP13", "NSP14", "NSP15", "NSP16", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10", "3'UTR"))

protlengths <- nim15N_ent %>% group_by(Sample, Protein) %>% count(Position) %>% tally
protlengths<- rename(protlengths, ProtLength=n)
nim15N_ent<- merge(nim15N_ent, protlengths, by=c("Sample","Protein"))

nim15N_ent$CntTs<-as.numeric(nim15N_ent$CntTs)
nim15N_ent$CntTv<-as.numeric(nim15N_ent$CntTv)
nim15N_ent$Coverage<-as.numeric(nim15N_ent$Coverage)
nim15N_ent <- nim15N_ent %>% filter(Coverage > 20)
nim15N_ent <- nim15N_ent %>% filter(Protein != "NA")
nim15N_ent$prop_ts <- nim15N_ent$CntTs / nim15N_ent$Coverage
nim15N_ent$prop_tv <- nim15N_ent$CntTv / nim15N_ent$Coverage

## Iahs mean tstv per protein 
mean_ts_perprot <- nim15N_ent %>% filter(Coverage >20) %>% group_by(Protein, Sample, ProtLength) %>% 
  summarize(Meants=mean(prop_ts, na.rm=T ), stdev=sd(prop_ts, na.rm=T), across())
mean_ts_perprot$Meants_pnorm <- mean_ts_perprot$Meants / mean_ts_perprot$ProtLength
#mean_ts_perprot$Meants_pnorm_scale <- scale(mean_ts_perprot$Meants_pnorm)
#mean_ts_perprot$MeantsNorm <- scale(mean_ts_perprot$Meants)

mean_tv_perprot <- nim15N_ent %>% filter(Coverage >20) %>% group_by(Protein, Sample, ProtLength) %>% 
  summarize(Meantv=mean(prop_tv, na.rm=T ), stdev=sd(prop_tv, na.rm=T), across())
mean_tv_perprot$Meantv_pnorm <- mean_tv_perprot$Meantv / mean_tv_perprot$ProtLength
#mean_tv_perprot$Meantv_pnorm_scale <- scale(mean_tv_perprot$Meantv_pnorm)
#mean_tv_perprot$MeantvNorm <- scale(mean_tv_perprot$Meantv)

meantstv_perprot <- merge(mean_ts_perprot, mean_tv_perprot, by = c("Sample", "Protein", "ProtLength", "Position", "Timepoint", "Patient", "Ct_value", "DPI", "Participant"))

meantstv_perprot$tstvratio_norm <- meantstv_perprot$MeantsNorm / meantstv_perprot$MeantvNorm
meantstv_perprot$tstvratio_norm <- meantstv_perprot$Meants_pnorm / meantstv_perprot$Meantv_pnorm # for unscaled
names(meantstv_perprot)[names(meantstv_perprot) == "tstvratio_norm"] <- "Ratio"

                                      ## use this one ?? ##
tstvvv <- ggplot(meantstv_perprot, aes(Protein, Participant, fill= Ratio)) + 
  geom_tile(colour="black") + scale_fill_gradient(low= "white", high= "blue") + ylab("Participant") + xlab("Protein") +
  ggtitle("Ts/Tv Ratio") + theme_bw() + theme(axis.text.x = element_text(angle=90, vjust= 0.7)) + facet_wrap(~Timepoint) 


ggsave(filename= "tstvratio_nim_hm_20X.tiff", plot= tstvvv, device='tiff', dpi=300, width=11.5, height=7)


