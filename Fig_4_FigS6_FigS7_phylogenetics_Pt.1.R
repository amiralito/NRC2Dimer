library(tidyverse)
library(Biostrings)
library(readxl)
library(ggtree)

# Set working directory to where the NLRtracker output of the 124 Solanaceae genomes are (from DOI: 10.5281/zenodo.10354350)
setwd("/path/to/NLRtracker")

# Make a function to add the file name to the proteins
read_delim_name <- function(flnm) {
  read_delim(flnm) %>%
    mutate(filename = flnm)
}


# Import the metadata
metadata <- list.files(pattern = "*_Domains.tsv", # metadata
                       full.names = T, recursive = TRUE) %>%
  map_df(~read_delim_name(.)) # bind all .tsv files together in a single dataframe
names(metadata)[11] <- "file_name" # correct the column name


# Modify the filename to RefSeq ID
metadata$file_name <- gsub("*_Domains.tsv", "", metadata$file_name)
metadata$file_name <- gsub(".*/", "", metadata$file_name)

metadata <- filter(metadata, !(metadata$file_name == "Niben261_genome.annotation.proteins")) # remove the additional benthamiana genome


# import genome data
genome_data <- read_excel("/path/to/DataS2.xlsx")
genome_meta <- genome_data[,c(1,2,3,19,23)]







# Subset based on "Status"
NLR_meta <- filter(metadata, metadata$Status == "NLR") # NLRs
NLR_meta <- filter(NLR_meta, NLR_meta$type == "CHAIN") # full length NLRs

write_csv(NLR_meta, "/path/to/NLR_meta.csv")

NLR_raw_seq <- AAStringSet(NLR_meta$sequence)
NLR_raw_seq@ranges@NAMES <- NLR_meta$seqname

writeXStringSet(NLR_raw_seq, "~/path/to/NLR_raw.fasta")

# filter the unwanted architectures
NLR_meta_filtered <- NLR_meta[NLR_meta$Simple %in% c("CNL","CN"),]

# merge with NLR meta

NLR_meta_filtered <- NLR_meta_filtered %>% left_join(genome_meta, by = "file_name")

# subset NBARC
NBARC_meta_filterd <- metadata[metadata$seqname %in% NLR_meta_filtered$seqname,]
NBARC_meta_filterd <- NBARC_meta_filterd[NBARC_meta_filterd$description == "NBARC",]


# NLR sequences dataframe
NLR_meta_filtered_seq <- NLR_meta_filtered[,c(1,10,11)]
NLR_seq <- AAStringSet(NLR_meta_filtered_seq$sequence)
NLR_seq@ranges@NAMES <- NLR_meta_filtered_seq$seqname

NBARC_meta_filtered_seq <- NBARC_meta_filterd[,c(1,10,11)]
NBARC_seq <- AAStringSet(NBARC_meta_filtered_seq$sequence)
NBARC_seq@ranges@NAMES <- NBARC_meta_filtered_seq$seqname

# import reference NRC sequences

## refplantNLR
RefPlantNLR_NBARC <- readAAStringSet("/path/to/RefPlantNLR_NBARC.fasta")

RefPlantNLR <- readAAStringSet("/path/to/RefPlantNLR.fasta")


## reference NRCs

Ca_NRC_NBARC <- readAAStringSet("/path/to/Ca_NRC_NBARC.fasta")
Sl_NRC_NBARC <- readAAStringSet("/path/to/Sl_NRC_NBARC.fasta")
St_NRC_NBARC <- readAAStringSet("/path/to/St_NRC_NBARC.fasta")
Nb_NRC_NBARC <- readAAStringSet("/path/to/Nb_NRC_NBARC.fasta")
NRC_ref_NBARC <- readAAStringSet("/path/to/NRC_ref_NBARC.fasta")


Ca_NRC <- readAAStringSet("/path/to/Ca_NRC.fasta")
Sl_NRC <- readAAStringSet("/path/to/Sl_NRC.fasta")
St_NRC <- readAAStringSet("/path/to/St_NRC.fasta")
Nb_NRC <- readAAStringSet("/path/to/Nb_NRC.fasta")
NRC_ref <- readAAStringSet("/path/to/NRC_ref.fasta")


NBARC_ref <- c(NBARC_seq, RefPlantNLR_NBARC, Ca_NRC_NBARC, Sl_NRC_NBARC) # combine them with main sequences

### St, Nb, and NRC_ref will be used for later stages

# filter out the truncated NBARC domains
NBARC_ref_filtered <- NBARC_ref[NBARC_ref@ranges@width > 300,]
NBARC_ref_filtered <- NBARC_ref_filtered[NBARC_ref_filtered@ranges@width < 400,]

NBARC_seq_filtered <- NBARC_seq[NBARC_seq@ranges@NAMES %in% NBARC_ref_filtered@ranges@NAMES]

NBARC_seq_filtered_meta <- NLR_meta_filtered[NLR_meta_filtered$seqname %in% NBARC_seq_filtered@ranges@NAMES,]
write_csv(NBARC_seq_filtered_meta, "/path/to/NBARC_seq_filtered_meta.csv")




# import and extract NRC sequences
NRC_tree <- read.tree("/path/to/NRC.tree")
NRC_NBARC <- NBARC_ref_filtered[NBARC_ref_filtered@ranges@NAMES %in% NRC_tree$tip.label]

NRC <- NLR_seq[NLR_seq@ranges@NAMES %in% NRC_tree$tip.label]

NRC_meta <- NLR_meta_filtered[NLR_meta_filtered$seqname %in% NRC@ranges@NAMES,]
write_csv(NRC_meta, "/path/to/NRC_meta.csv")

# exctract the NRCH
NRCH_tree <- read.tree("/path/to/NRCH.tree")

NRCH_meta <- NLR_meta_filtered[NLR_meta_filtered$seqname %in% NRCH_tree$tip.label,]
write_csv(NRCH_meta,"/path/to/NRCH_meta.csv")

NRCH_NBARC <- NBARC_seq_filtered[NBARC_seq_filtered@ranges@NAMES %in% NRCH_tree$tip.label]
writeXStringSet(NRCH_NBARC, "/path/to/NRCH_NBARC.fasta")
NRCH <- NLR_seq[NLR_seq@ranges@NAMES %in% NRCH_tree$tip.label]
writeXStringSet(NRCH, "/path/to/NRCH.fasta")


# import the deduplicated NRC sequences from CD-HIT output

NRCH_clu <- readAAStringSet("/path/to/NRCH_clu.fasta")
NRCH_clu_NBARC <- NBARC_seq_filtered[NBARC_seq_filtered@ranges@NAMES %in% NRCH_clu@ranges@NAMES]


# only keep sequences within 750-950 AA range

NRCH_clu_filtered <- NRCH_clu[NRCH_clu@ranges@width < 950]
NRCH_clu_filtered <- NRCH_clu_filtered[NRCH_clu_filtered@ranges@width > 750]
writeXStringSet(NRCH_clu_filtered, "/path/to/NRCH_clu_filtered.fasta")

NRCH_clu_NBARC_filtered <- NRCH_clu_NBARC[NRCH_clu_NBARC@ranges@NAMES %in% NRCH_clu_filtered@ranges@NAMES]
writeXStringSet(NRCH_clu_NBARC_filtered, "/path/to/NRCH_clu_NBARC_filtered.fasta")


# merge it with reference sequences
NRCH_clu_filtered_ref <- c(NRCH_clu_filtered, Ca_NRC, St_NRC, Sl_NRC, Nb_NRC, NRC_ref)
NRCH_clu_NBARC_filtered_ref <- c(NRCH_clu_NBARC_filtered, 
                                 Ca_NRC_NBARC, 
                                 St_NRC_NBARC, 
                                 Sl_NRC_NBARC, 
                                 Nb_NRC_NBARC, 
                                 NRC_ref_NBARC)

writeXStringSet(NRCH_clu_filtered_ref, "/path/to/NRCH_clu_filtered_ref.fasta")
writeXStringSet(NRCH_clu_NBARC_filtered_ref, "/path/to/NRCH_clu_NBARC_filtered_ref.fasta")

