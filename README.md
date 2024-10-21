# Phylogenomics and diversification analysis of the NRC helper clade from 123 Solanaceae proteomes.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13362063.svg)](https://doi.org/10.5281/zenodo.13362063)
[![DOI](https://img.shields.io/badge/bioRxiv-doi.org/10.1101/2023.12.17.572070-BE2634.svg)](https://doi.org/10.1101/2023.12.17.572070)
[![DOI](https://img.shields.io/badge/plos_biology-doi.org/10.1371/journal.pbio.3002868-1334C6.svg)](https://doi.org/10.1371/journal.pbio.3002868)

# Supporting scripts and material for "Activation of plant immunity through conversion of a helper NLR homodimer into a resistosome"
Muniyandi Selvaraj, AmirAli Toghani, Hsuan Pai, Yu Sugihara, Jiorgos Kourelis, Enoch Lok Him Yuen, Tarhan Ibrahim, He Zhao, Rongrong Xie, Abbas Maqbool, Juan Carlos De la Concepcion, Mark J Banfield, Lida Derevnina, Benjamin Petre, David M Lawson, Tolga O Bozkurt, Chih-Hang Wu, Sophien Kamoun, Mauricio P Contreras


Resources:
Software                            | Source
------------------------------------| ------------------------------------
*MAFFT v7.520*                      | (https://github.com/GSLBiotech/mafft)
*FastTree v2.1.11*                  | (http://www.microbesonline.org/fasttree/)
*Dendroscope v3.8.8*                | (https://software-ab.cs.uni-tuebingen.de/download/dendroscope3/welcome.html)
*R v4.3.1*                          | (https://cran.r-project.org/)
*ClipKIT v2.0.1*                    | (https://github.com/JLSteenwyk/ClipKIT)
*NLRtracker*                        | (https://github.com/slt666666/NLRtracker)

R packages:
```R
install.packages("tidyverse")
install.packages("ggseqlogo")
install.packages("entropy")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biostrings")
BiocManager::install("ggtree")
```
