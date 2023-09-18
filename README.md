# SARS-CoV-2 population dynamics in immunocompetent individuals in a closed transmission chain shows genomic diversity over the course of infection
\
Hannah Goldswain<sup>1</sup>, Rebekah Penrice-Randal<sup>1</sup>, I’ah Donovan-Banfield<sup>1</sup>, Xiaofeng Dong<sup>1</sup>, Nadine Randle<sup>1</sup>, Yan Ryan<sup>1</sup>, Alex Rzeszutek<sup>2</sup>, Jack Pilgrim<sup>2</sup>, Emma Keyser<sup>3</sup>, Simon A. Weller<sup>3</sup>, Emma J. Hutley<sup>4</sup>, Catherine Hartley<sup>1</sup>, Tessa Prince<sup>1</sup>, Alistair C. Darby<sup>1</sup>, Niall Aye Maung<sup>5</sup>, Henry Nwume<sup>3</sup>, Julian A. Hiscox<sup>1,6</sup>, and Stevan R. Emmett<sup>3</sup>.

<sup>1</sup>Institute for Infection, Veterinary and Ecological Sciences, University of Liverpool, UK. \
<sup>2</sup>Centre for Genomic Research, University of Liverpool, UK.\
<sup>3</sup>Defence Science Technology Laboratory, Porton Down, UK.\
<sup>4</sup>Centre for Defence Pathology, OCT Centre, Royal Centre for Defence Medicine, Birmingham, B15 2WB.\
<sup>5</sup>British Army, Hunter House, St Omer Barracks, Aldershot, Hampshire, GU11 2BG\
<sup>6</sup>A*STAR Infectious Diseases Laboratories (A*STAR ID Labs), Agency for Science, Technology and Research (A*STAR), Singapore. \


## Sequencing of longitudinal nasopharyngeal swab samples from immunocompetent participants in an isolated geographical location

1. 50 participant samples sequenced via NimaGen-Illumina sequencing (Coolen et al., 2021).
2. EasySeq pipeline used to analyse data (https://github.com/JordyCoolen/easyseq_covid19).
3. Sequences filtered at 85% coverage across the genome.
4. Further downstream analysis using DiversiTools (https://josephhughes.github.io/DiversiTools/tutorial.html), Pangolin (https://github.com/cov-lineages/pangolin), Snipit (https://github.com/aineniamh/snipit) and Nextclade (https://clades.nextstrain.org/).
5. Plotting of outputs using custom R scripts in this repository.

## Scripts included in repository
```Syn_NonSyn_parse_aa_V3.pl # developed by XD/RPR/IDB/HG 
heatmap_updated_paper.R # developed by HG/IDB
letrs_paper.R # developed by HG
nonsyn_lineplot_paper.R # developed by HG/IDB
prop_diff_colours_paper.R # developed by HG/IDB
top_sndaa_paper.R # developed by HG/IDB 
```
Contributors: XD- Xiaofeng Dong, RPR- Rebekah Penrice-Randal, IDB- I'ah Donovan Banfield, HG- Hannah Goldswain

## Raw data availability
Raw fastq files have been deposited in the National Center for Biotechnology Information (NCBI) Short Read Archive (SRA) under the bioproject: . Accessible here: 

## Contact us
Hannah Goldswain at hannah.goldswain@liverpool.ac.uk, Julian Hiscox at julianh@liverpool.ac.uk.

