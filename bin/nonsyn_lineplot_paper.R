# Script to plot linegraph for synonymous and non-synonymous mutations across the genome using parse.txt output from Syn_NonSyn_parse_aa_V3.pl post DiversiTools analysis
# Add column of nrow amino acid position and then rbind and rename column for end plot

nsl_36 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3636_TAGTCGTCAA-ATGGCTCGGT_L002.final_AA_parse.txt", sep="\t")
nsl_36$AAposition <- 1:nrow(nsl_36)
nsl_36$Name <- 'P 01'
nsl_36$Sample <- 'S1'
nsl_37 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3637_TGGACTAATT-CAGCGGATAA_L002.final_AA_parse.txt", sep="\t")
nsl_37$AAposition <- 1:nrow(nsl_37)
nsl_37$Name <- 'P 01'
nsl_37$Sample <- 'S2'
nsl_39 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3639_ATCTTATGAT-GCCTGGCAGT_L002.final_AA_parse.txt", sep="\t")
nsl_39$AAposition <- 1:nrow(nsl_39)
nsl_39$Name <- 'P 02'
nsl_39$Sample <- 'S1'
nsl_40 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3640_ACCAACGCTG-TAACTTGGAG_L002.final_AA_parse.txt", sep="\t")
nsl_40$AAposition <- 1:nrow(nsl_40)
nsl_40$Name <- 'P 02'
nsl_40$Sample <- 'S2'
nsl_43 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3643_TACGGCGCGA-ACGTATGCGC_L002.final_AA_parse.txt", sep="\t")
nsl_43$AAposition <- 1:nrow(nsl_43)
nsl_43$Name <- 'P 03'
nsl_43$Sample <- 'S1'
nsl_44 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3644_AGCTATCAAC-TGATTGCTTC_L002.final_AA_parse.txt", sep="\t")
nsl_44$AAposition <- 1:nrow(nsl_44)
nsl_44$Name <- 'P 03'
nsl_44$Sample <- 'S2'
nsl_47 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3647_AGAGGTCATT-GGCATATATA_L002.final_AA_parse.txt", sep="\t")
nsl_47$AAposition <- 1:nrow(nsl_47)
nsl_47$Name <- 'P 04'
nsl_47$Sample <- 'S1'
nsl_49 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3649_CTTAGAGGCA-AGCGTCTGGT_L002.final_AA_parse.txt", sep="\t")
nsl_49$AAposition <- 1:nrow(nsl_49)
nsl_49$Name <- 'P 05'
nsl_49$Sample <- 'S1'
nsl_50 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3650_ATAGAGCATT-GAGGACCGAT_L002.final_AA_parse.txt", sep="\t")
nsl_50$AAposition <- 1:nrow(nsl_50)
nsl_50$Name <- 'P 05'
nsl_50$Sample <- 'S2'
nsl_51 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3651_GGTAGCGCAT-TACTGGATAA_L002.final_AA_parse.txt", sep="\t")
nsl_51$AAposition <- 1:nrow(nsl_51)
nsl_51$Name <- 'P 06'
nsl_51$Sample <- 'S1'
nsl_52 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3652_TAATGATATT-CGAGTCGAAG_L002.final_AA_parse.txt", sep="\t")
nsl_52$AAposition <- 1:nrow(nsl_52)
nsl_52$Name <- 'P 06'
nsl_52$Sample <- 'S2'
nsl_53 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3653_GGCTAAGGCG-TTCTTCTATT_L002.final_AA_parse.txt", sep="\t")
nsl_53$AAposition <- 1:nrow(nsl_53)
nsl_53$Name <- 'P 07'
nsl_53$Sample <- 'S1'
nsl_54 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3654_TTAGTACCTT-GACCAATTCT_L002.final_AA_parse.txt", sep="\t")
nsl_54$AAposition <- 1:nrow(nsl_54)
nsl_54$Name <- 'P 07'
nsl_54$Sample <- 'S2'
nsl_59 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3659_GTAGCTACCG-TCAGAACGAC_L002.final_AA_parse.txt", sep="\t")
nsl_59$AAposition <- 1:nrow(nsl_59)
nsl_59$Name <- 'P 09'
nsl_59$Sample <- 'S1'
nsl_60 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3660_ATTAATAGAA-GGCAGAGGAG_L002.final_AA_parse.txt", sep="\t")
nsl_60$AAposition <- 1:nrow(nsl_60)
nsl_60$Name <- 'P 09'
nsl_60$Sample <- 'S2'
nsl_61 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3661_ATATGATGAA-ATGAAGGAGG_L002.final_AA_parse.txt", sep="\t")
nsl_61$AAposition <- 1:nrow(nsl_61)
nsl_61$Name <- 'P 09'
nsl_61$Sample <- 'S3'
nsl_62 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3662_AATGGACGTA-GAAGGTTCGG_L002.final_AA_parse.txt", sep="\t")
nsl_62$AAposition <- 1:nrow(nsl_62)
nsl_62$Name <- 'P 10'
nsl_62$Sample <- 'S1'
nsl_64 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3664_ACCATAAGCG-TGAATAACTA_L002.final_AA_parse.txt", sep="\t")
nsl_64$AAposition <- 1:nrow(nsl_64)
nsl_64$Name <- 'P 11'
nsl_64$Sample <- 'S2'
nsl_65 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3665_CCGCTGCATT-GAGTCCTGCA_L002.final_AA_parse.txt", sep="\t")
nsl_65$AAposition <- 1:nrow(nsl_65)
nsl_65$Name <- 'P 11'
nsl_65$Sample <- 'S3'
nsl_66 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3666_TTGGTTGGTA-GCCATAGTAG_L002.final_AA_parse.txt", sep="\t")
nsl_66$AAposition <- 1:nrow(nsl_66)
nsl_66$Name <- 'P 11'
nsl_66$Sample <- 'S4'
nsl_67 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3667_GGAACGATTC-CGGCGAGTCC_L002.final_AA_parse.txt", sep="\t")
nsl_67$AAposition <- 1:nrow(nsl_67)
nsl_67$Name <- 'P 12'
nsl_67$Sample <- 'S1'
nsl_68 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3668_AAGCATCTCG-CTCATGATTG_L002.final_AA_parse.txt", sep="\t")
nsl_68$AAposition <- 1:nrow(nsl_68)
nsl_68$Name <- 'P 12'
nsl_68$Sample <- 'S2'
nsl_70 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3670_TATAGATGGC-GCGAGAGAAG_L002.final_AA_parse.txt", sep="\t")
nsl_70$AAposition <- 1:nrow(nsl_70)
nsl_70$Name <- 'P 13'
nsl_70$Sample <- 'S1'
nsl_74 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3674_TCCGTATTCC-TCTAATTATT_L002.final_AA_parse.txt", sep="\t")
nsl_74$AAposition <- 1:nrow(nsl_74)
nsl_74$Name <- 'P 14'
nsl_74$Sample <- 'S1'
nsl_78 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3678_ATACTATATT-GGTATCATCT_L002.final_AA_parse.txt", sep="\t")
nsl_78$AAposition <- 1:nrow(nsl_78)
nsl_78$Name <- 'P 15'
nsl_78$Sample <- 'S1'
nsl_79 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3679_TAGTTATGCG-CTGGCTTAGT_L002.final_AA_parse.txt", sep="\t")
nsl_79$AAposition <- 1:nrow(nsl_79)
nsl_79$Name <- 'P 15'
nsl_79$Sample <- 'S2'
nsl_80 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3680_AACTTACGCT-TAACCTGCAA_L002.final_AA_parse.txt", sep="\t")
nsl_80$AAposition <- 1:nrow(nsl_80)
nsl_80$Name <- 'P 15'
nsl_80$Sample <- 'S3'
nsl_82 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3682_GTCCGTAAGG-CGACCTATAC_L002.final_AA_parse.txt", sep="\t")
nsl_82$AAposition <- 1:nrow(nsl_82)
nsl_82$Name <- 'P 16'
nsl_82$Sample <- 'S1'
nsl_83 <- read.delim("/home/hannahg/projects/dstl_project/data/nimagen/DiversiTools/3683_CCGGTCGACG-GCTGCATAGC_L002.final_AA_parse.txt", sep="\t")
nsl_83$AAposition <- 1:nrow(nsl_83)
nsl_83$Name <- 'P 16'
nsl_83$Sample <- 'S2'

rbind_nim_15N <- rbind(nsl_36, nsl_37, nsl_39, nsl_40, nsl_43, nsl_44, nsl_47, nsl_49, nsl_50, nsl_51, nsl_52, nsl_53, nsl_54, nsl_59, nsl_60, nsl_61, nsl_62, nsl_65, nsl_67, nsl_68, nsl_70, nsl_74, nsl_79, nsl_82, nsl_83) # rbind all dataframes for plotting

rbind_nim_15N$AAcoverage<-as.numeric(rbind_nim_15N$AAcoverage) # change variables to numeric for calculations
rbind_nim_15N$CntNonSyn<-as.numeric(rbind_nim_15N$CntNonSyn) 
rbind_nim_15N$AAposition<-as.numeric(rbind_nim_15N$AAposition) 
rbind_nim_15N$Proportion <-  rbind_nim_15N$CntNonSyn / rbind_nim_15N$AAcoverage # new column with proportion of cntnonsyn
rbind_nim_15N <- rbind_nim_15N %>% filter(AAcoverage > 20) # filter coverage at 20
rbind_nim_15N$Proportion<-as.numeric(rbind_nim_15N$Proportion) # make variable numeric to allow plotting
rbind_nim_15N$Propsyn <- rbind_nim_15N$CntSyn / rbind_nim_15N$AAcoverage # calculate proportions
rbind_nim_15N$allprop <- (rbind_nim_15N$CntNonSyn + rbind_nim_15N$CntSyn) / rbind_nim_15N$AAcoverage

names(rbind_nim_15N)[names(rbind_nim_15N) == "Proportion"] <- "Non-synonymous nucleotide mutations" # change names for plot
names(rbind_nim_15N)[names(rbind_nim_15N) == "Propsyn"] <- "Synonymous nucleotide mutations"

cbp1 <- c("#E69F00", "#56B4E9", "#009E73", "#7CFC00") # make vector colour schemes for plot
cbp2 <- c("#E69F00", "#56B4E9", "#009E73")
cbp3 <- c("#009E73", "#E69F00")

#line plot with different colours for synonymous and non-synonymous changes using above vectors
melt_lineplotrbind <- reshape2::melt(rbind_nim_15N, id.vars = c('AAposition', 'Sample', 'Name', 'Non-synonymous nucleotide mutations', 'Synonymous nucleotide mutations'), measure.vars= c('Non-synonymous nucleotide mutations', 'Synonymous nucleotide mutations')) # melt data for ggplot
lpp <- ggplot(data= melt_lineplotrbind, aes(x=AAposition, y=value, colour=variable)) + geom_line() +scale_colour_manual(values=cbp3) +
  ylab("Proportion") + xlab("Amino acid position") + ylim(0,1) + 
  theme(axis.title.x = element_text(face ="bold", size=13), 
        legend.title = element_blank(), 
        legend.position = "bottom",
        axis.title.y = element_text(face="bold", size = 13)) 

gglpp<- lpp + facet_grid(Name ~ Sample) # facet plots to see all data on one plot
ggsave(filename= "synnonsyncolour_lineplot_nim_participants.tiff", plot = gglpp, device='tiff', dpi= 300, width=8, height=12)



