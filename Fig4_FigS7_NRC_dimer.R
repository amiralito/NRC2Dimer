library(tidyverse)
library(Biostrings)
library(ggseqlogo)
library(ggtree)


# import the alignment file (local alignment > this is the file we're gonna use for the paper)
NRCH_clu_filtered_ref_trimmed <- readAAStringSet("/path/to/NRCH_clu_filtered_ref_local_trimmed.afa")


# extract the three surfaces
NRCH_S1a <- subseq(NRCH_clu_filtered_ref_trimmed, start = 234, end = 237)
NRCH_S1b <- subseq(NRCH_clu_filtered_ref_trimmed, start = 256, end = 262)
NRCH_S1c <- subseq(NRCH_clu_filtered_ref_trimmed, start = 288, end = 292)
NRCH_S1d <- subseq(NRCH_clu_filtered_ref_trimmed, start = 522, end = 529)
NRCH_S1e <- subseq(NRCH_clu_filtered_ref_trimmed, start = 550, end = 553)
NRCH_S1f <- subseq(NRCH_clu_filtered_ref_trimmed, start = 578, end = 582)

NRCH_S2a <- subseq(NRCH_clu_filtered_ref_trimmed, start = 544, end = 549)
NRCH_S2b <- subseq(NRCH_clu_filtered_ref_trimmed, start = 566, end = 573)

# extract the conserved motifs
NRCH_MADA <- subseq(NRCH_clu_filtered_ref_trimmed, start = 12, end = 32)
NRCH_Ploop <- subseq(NRCH_clu_filtered_ref_trimmed, start = 199, end = 209)
NRCH_MHD <- subseq(NRCH_clu_filtered_ref_trimmed, start = 495, end = 497)

# convert to dataframe
NRCH_S1a_df <- as.data.frame(NRCH_S1a)
NRCH_S1b_df <- as.data.frame(NRCH_S1b)
NRCH_S1c_df <- as.data.frame(NRCH_S1c)
NRCH_S1d_df <- as.data.frame(NRCH_S1d)
NRCH_S1e_df <- as.data.frame(NRCH_S1e)
NRCH_S1f_df <- as.data.frame(NRCH_S1f)

NRCH_S2a_df <- as.data.frame(NRCH_S2a)
NRCH_S2b_df <- as.data.frame(NRCH_S2b)


NRCH_MADA_df <- as.data.frame(NRCH_MADA)
NRCH_Ploop_df <- as.data.frame(NRCH_Ploop)
NRCH_MHD_df <- as.data.frame(NRCH_MHD)


# import the subclades based on tree
NRC0_tree <- read.tree("/path/to/NRC0.tree")
NRC1_tree <- read.tree("/path/to/NRC1.tree")
NRC2_tree <- read.tree("/path/to/NRC2.tree")
NRC3_tree <- read.tree("/path/to/NRC3.tree")
NRC4a_tree <- read.tree("/path/to/NRC4a.tree")
NRC4b_tree <- read.tree("/path/to/NRC4b.tree")
NRC4c_tree <- read.tree("/path/to/NRC4c.tree")
NRC4t_tree <- read.tree("/path/to/NRC4t.tree")
NRC4other_tree <- read.tree("/path/to/NRC4other.tree")
NRC5_tree <- read.tree("/path/to/NRC5.tree")
NRC6_tree <- read.tree("/path/to/NRC6.tree")
NRC7_tree <- read.tree("/path/to/NRC7.tree")
NRC8_tree <- read.tree("/path/to/NRC8.tree")
NRC9_tree <- read.tree("/path/to/NRC9.tree")
NRCX_tree <- read.tree("/path/to/NRCX.tree")



## extract the motifs for each subclade

# NRC0
NRC0_S1a <- NRCH_S1a_df[rownames(NRCH_S1a_df) %in% NRC0_tree$tip.label,] %>% as.data.frame()
NRC0_S1b <- NRCH_S1b_df[rownames(NRCH_S1b_df) %in% NRC0_tree$tip.label,] %>% as.data.frame()
NRC0_S1c <- NRCH_S1c_df[rownames(NRCH_S1c_df) %in% NRC0_tree$tip.label,] %>% as.data.frame()
NRC0_S1d <- NRCH_S1d_df[rownames(NRCH_S1d_df) %in% NRC0_tree$tip.label,] %>% as.data.frame()
NRC0_S1e <- NRCH_S1e_df[rownames(NRCH_S1e_df) %in% NRC0_tree$tip.label,] %>% as.data.frame()
NRC0_S1f <- NRCH_S1f_df[rownames(NRCH_S1f_df) %in% NRC0_tree$tip.label,] %>% as.data.frame()

NRC0_S2a <- NRCH_S2a_df[rownames(NRCH_S2a_df) %in% NRC0_tree$tip.label,] %>% as.data.frame()
NRC0_S2b <- NRCH_S2b_df[rownames(NRCH_S2b_df) %in% NRC0_tree$tip.label,] %>% as.data.frame()

NRC0_MADA <- NRCH_MADA_df[rownames(NRCH_MADA_df) %in% NRC0_tree$tip.label,] %>% as.data.frame()
NRC0_Ploop <- NRCH_Ploop_df[rownames(NRCH_Ploop_df) %in% NRC0_tree$tip.label,] %>% as.data.frame()
NRC0_MHD <- NRCH_MHD_df[rownames(NRCH_MHD_df) %in% NRC0_tree$tip.label,] %>% as.data.frame()


# NRC1
NRC1_S1a <- NRCH_S1a_df[rownames(NRCH_S1a_df) %in% NRC1_tree$tip.label,] %>% as.data.frame()
NRC1_S1b <- NRCH_S1b_df[rownames(NRCH_S1b_df) %in% NRC1_tree$tip.label,] %>% as.data.frame()
NRC1_S1c <- NRCH_S1c_df[rownames(NRCH_S1c_df) %in% NRC1_tree$tip.label,] %>% as.data.frame()
NRC1_S1d <- NRCH_S1d_df[rownames(NRCH_S1d_df) %in% NRC1_tree$tip.label,] %>% as.data.frame()
NRC1_S1e <- NRCH_S1e_df[rownames(NRCH_S1e_df) %in% NRC1_tree$tip.label,] %>% as.data.frame()
NRC1_S1f <- NRCH_S1f_df[rownames(NRCH_S1f_df) %in% NRC1_tree$tip.label,] %>% as.data.frame()

NRC1_S2a <- NRCH_S2a_df[rownames(NRCH_S2a_df) %in% NRC1_tree$tip.label,] %>% as.data.frame()
NRC1_S2b <- NRCH_S2b_df[rownames(NRCH_S2b_df) %in% NRC1_tree$tip.label,] %>% as.data.frame()

NRC1_MADA <- NRCH_MADA_df[rownames(NRCH_MADA_df) %in% NRC1_tree$tip.label,] %>% as.data.frame()
NRC1_Ploop <- NRCH_Ploop_df[rownames(NRCH_Ploop_df) %in% NRC1_tree$tip.label,] %>% as.data.frame()
NRC1_MHD <- NRCH_MHD_df[rownames(NRCH_MHD_df) %in% NRC1_tree$tip.label,] %>% as.data.frame()


# NRC2
NRC2_S1a <- NRCH_S1a_df[rownames(NRCH_S1a_df) %in% NRC2_tree$tip.label,] %>% as.data.frame()
NRC2_S1b <- NRCH_S1b_df[rownames(NRCH_S1b_df) %in% NRC2_tree$tip.label,] %>% as.data.frame()
NRC2_S1c <- NRCH_S1c_df[rownames(NRCH_S1c_df) %in% NRC2_tree$tip.label,] %>% as.data.frame()
NRC2_S1d <- NRCH_S1d_df[rownames(NRCH_S1d_df) %in% NRC2_tree$tip.label,] %>% as.data.frame()
NRC2_S1e <- NRCH_S1e_df[rownames(NRCH_S1e_df) %in% NRC2_tree$tip.label,] %>% as.data.frame()
NRC2_S1f <- NRCH_S1f_df[rownames(NRCH_S1f_df) %in% NRC2_tree$tip.label,] %>% as.data.frame()

NRC2_S2a <- NRCH_S2a_df[rownames(NRCH_S2a_df) %in% NRC2_tree$tip.label,] %>% as.data.frame()
NRC2_S2b <- NRCH_S2b_df[rownames(NRCH_S2b_df) %in% NRC2_tree$tip.label,] %>% as.data.frame()

NRC2_MADA <- NRCH_MADA_df[rownames(NRCH_MADA_df) %in% NRC2_tree$tip.label,] %>% as.data.frame()
NRC2_Ploop <- NRCH_Ploop_df[rownames(NRCH_Ploop_df) %in% NRC2_tree$tip.label,] %>% as.data.frame()
NRC2_MHD <- NRCH_MHD_df[rownames(NRCH_MHD_df) %in% NRC2_tree$tip.label,] %>% as.data.frame()


# NRC3
NRC3_S1a <- NRCH_S1a_df[rownames(NRCH_S1a_df) %in% NRC3_tree$tip.label,] %>% as.data.frame()
NRC3_S1b <- NRCH_S1b_df[rownames(NRCH_S1b_df) %in% NRC3_tree$tip.label,] %>% as.data.frame()
NRC3_S1c <- NRCH_S1c_df[rownames(NRCH_S1c_df) %in% NRC3_tree$tip.label,] %>% as.data.frame()
NRC3_S1d <- NRCH_S1d_df[rownames(NRCH_S1d_df) %in% NRC3_tree$tip.label,] %>% as.data.frame()
NRC3_S1e <- NRCH_S1e_df[rownames(NRCH_S1e_df) %in% NRC3_tree$tip.label,] %>% as.data.frame()
NRC3_S1f <- NRCH_S1f_df[rownames(NRCH_S1f_df) %in% NRC3_tree$tip.label,] %>% as.data.frame()

NRC3_S2a <- NRCH_S2a_df[rownames(NRCH_S2a_df) %in% NRC3_tree$tip.label,] %>% as.data.frame()
NRC3_S2b <- NRCH_S2b_df[rownames(NRCH_S2b_df) %in% NRC3_tree$tip.label,] %>% as.data.frame()

NRC3_MADA <- NRCH_MADA_df[rownames(NRCH_MADA_df) %in% NRC3_tree$tip.label,] %>% as.data.frame()
NRC3_Ploop <- NRCH_Ploop_df[rownames(NRCH_Ploop_df) %in% NRC3_tree$tip.label,] %>% as.data.frame()
NRC3_MHD <- NRCH_MHD_df[rownames(NRCH_MHD_df) %in% NRC3_tree$tip.label,] %>% as.data.frame()


# NRC4a
NRC4a_S1a <- NRCH_S1a_df[rownames(NRCH_S1a_df) %in% NRC4a_tree$tip.label,] %>% as.data.frame()
NRC4a_S1b <- NRCH_S1b_df[rownames(NRCH_S1b_df) %in% NRC4a_tree$tip.label,] %>% as.data.frame()
NRC4a_S1c <- NRCH_S1c_df[rownames(NRCH_S1c_df) %in% NRC4a_tree$tip.label,] %>% as.data.frame()
NRC4a_S1d <- NRCH_S1d_df[rownames(NRCH_S1d_df) %in% NRC4a_tree$tip.label,] %>% as.data.frame()
NRC4a_S1e <- NRCH_S1e_df[rownames(NRCH_S1e_df) %in% NRC4a_tree$tip.label,] %>% as.data.frame()
NRC4a_S1f <- NRCH_S1f_df[rownames(NRCH_S1f_df) %in% NRC4a_tree$tip.label,] %>% as.data.frame()

NRC4a_S2a <- NRCH_S2a_df[rownames(NRCH_S2a_df) %in% NRC4a_tree$tip.label,] %>% as.data.frame()
NRC4a_S2b <- NRCH_S2b_df[rownames(NRCH_S2b_df) %in% NRC4a_tree$tip.label,] %>% as.data.frame()

NRC4a_MADA <- NRCH_MADA_df[rownames(NRCH_MADA_df) %in% NRC4a_tree$tip.label,] %>% as.data.frame()
NRC4a_Ploop <- NRCH_Ploop_df[rownames(NRCH_Ploop_df) %in% NRC4a_tree$tip.label,] %>% as.data.frame()
NRC4a_MHD <- NRCH_MHD_df[rownames(NRCH_MHD_df) %in% NRC4a_tree$tip.label,] %>% as.data.frame()


# NRC4b
NRC4b_S1a <- NRCH_S1a_df[rownames(NRCH_S1a_df) %in% NRC4b_tree$tip.label,] %>% as.data.frame()
NRC4b_S1b <- NRCH_S1b_df[rownames(NRCH_S1b_df) %in% NRC4b_tree$tip.label,] %>% as.data.frame()
NRC4b_S1c <- NRCH_S1c_df[rownames(NRCH_S1c_df) %in% NRC4b_tree$tip.label,] %>% as.data.frame()
NRC4b_S1d <- NRCH_S1d_df[rownames(NRCH_S1d_df) %in% NRC4b_tree$tip.label,] %>% as.data.frame()
NRC4b_S1e <- NRCH_S1e_df[rownames(NRCH_S1e_df) %in% NRC4b_tree$tip.label,] %>% as.data.frame()
NRC4b_S1f <- NRCH_S1f_df[rownames(NRCH_S1f_df) %in% NRC4b_tree$tip.label,] %>% as.data.frame()

NRC4b_S2a <- NRCH_S2a_df[rownames(NRCH_S2a_df) %in% NRC4b_tree$tip.label,] %>% as.data.frame()
NRC4b_S2b <- NRCH_S2b_df[rownames(NRCH_S2b_df) %in% NRC4b_tree$tip.label,] %>% as.data.frame()

NRC4b_MADA <- NRCH_MADA_df[rownames(NRCH_MADA_df) %in% NRC4b_tree$tip.label,] %>% as.data.frame()
NRC4b_Ploop <- NRCH_Ploop_df[rownames(NRCH_Ploop_df) %in% NRC4b_tree$tip.label,] %>% as.data.frame()
NRC4b_MHD <- NRCH_MHD_df[rownames(NRCH_MHD_df) %in% NRC4b_tree$tip.label,] %>% as.data.frame()


# NRC4c
NRC4c_S1a <- NRCH_S1a_df[rownames(NRCH_S1a_df) %in% NRC4c_tree$tip.label,] %>% as.data.frame()
NRC4c_S1b <- NRCH_S1b_df[rownames(NRCH_S1b_df) %in% NRC4c_tree$tip.label,] %>% as.data.frame()
NRC4c_S1c <- NRCH_S1c_df[rownames(NRCH_S1c_df) %in% NRC4c_tree$tip.label,] %>% as.data.frame()
NRC4c_S1d <- NRCH_S1d_df[rownames(NRCH_S1d_df) %in% NRC4c_tree$tip.label,] %>% as.data.frame()
NRC4c_S1e <- NRCH_S1e_df[rownames(NRCH_S1e_df) %in% NRC4c_tree$tip.label,] %>% as.data.frame()
NRC4c_S1f <- NRCH_S1f_df[rownames(NRCH_S1f_df) %in% NRC4c_tree$tip.label,] %>% as.data.frame()

NRC4c_S2a <- NRCH_S2a_df[rownames(NRCH_S2a_df) %in% NRC4c_tree$tip.label,] %>% as.data.frame()
NRC4c_S2b <- NRCH_S2b_df[rownames(NRCH_S2b_df) %in% NRC4c_tree$tip.label,] %>% as.data.frame()

NRC4c_MADA <- NRCH_MADA_df[rownames(NRCH_MADA_df) %in% NRC4c_tree$tip.label,] %>% as.data.frame()
NRC4c_Ploop <- NRCH_Ploop_df[rownames(NRCH_Ploop_df) %in% NRC4c_tree$tip.label,] %>% as.data.frame()
NRC4c_MHD <- NRCH_MHD_df[rownames(NRCH_MHD_df) %in% NRC4c_tree$tip.label,] %>% as.data.frame()


# NRC4t
NRC4t_S1a <- NRCH_S1a_df[rownames(NRCH_S1a_df) %in% NRC4t_tree$tip.label,] %>% as.data.frame()
NRC4t_S1b <- NRCH_S1b_df[rownames(NRCH_S1b_df) %in% NRC4t_tree$tip.label,] %>% as.data.frame()
NRC4t_S1c <- NRCH_S1c_df[rownames(NRCH_S1c_df) %in% NRC4t_tree$tip.label,] %>% as.data.frame()
NRC4t_S1d <- NRCH_S1d_df[rownames(NRCH_S1d_df) %in% NRC4t_tree$tip.label,] %>% as.data.frame()
NRC4t_S1e <- NRCH_S1e_df[rownames(NRCH_S1e_df) %in% NRC4t_tree$tip.label,] %>% as.data.frame()
NRC4t_S1f <- NRCH_S1f_df[rownames(NRCH_S1f_df) %in% NRC4t_tree$tip.label,] %>% as.data.frame()

NRC4t_S2a <- NRCH_S2a_df[rownames(NRCH_S2a_df) %in% NRC4t_tree$tip.label,] %>% as.data.frame()
NRC4t_S2b <- NRCH_S2b_df[rownames(NRCH_S2b_df) %in% NRC4t_tree$tip.label,] %>% as.data.frame()

NRC4t_MADA <- NRCH_MADA_df[rownames(NRCH_MADA_df) %in% NRC4t_tree$tip.label,] %>% as.data.frame()
NRC4t_Ploop <- NRCH_Ploop_df[rownames(NRCH_Ploop_df) %in% NRC4t_tree$tip.label,] %>% as.data.frame()
NRC4t_MHD <- NRCH_MHD_df[rownames(NRCH_MHD_df) %in% NRC4t_tree$tip.label,] %>% as.data.frame()


# NRC4other
NRC4other_S1a <- NRCH_S1a_df[rownames(NRCH_S1a_df) %in% NRC4other_tree$tip.label,] %>% as.data.frame()
NRC4other_S1b <- NRCH_S1b_df[rownames(NRCH_S1b_df) %in% NRC4other_tree$tip.label,] %>% as.data.frame()
NRC4other_S1c <- NRCH_S1c_df[rownames(NRCH_S1c_df) %in% NRC4other_tree$tip.label,] %>% as.data.frame()
NRC4other_S1d <- NRCH_S1d_df[rownames(NRCH_S1d_df) %in% NRC4other_tree$tip.label,] %>% as.data.frame()
NRC4other_S1e <- NRCH_S1e_df[rownames(NRCH_S1e_df) %in% NRC4other_tree$tip.label,] %>% as.data.frame()
NRC4other_S1f <- NRCH_S1f_df[rownames(NRCH_S1f_df) %in% NRC4other_tree$tip.label,] %>% as.data.frame()

NRC4other_S2a <- NRCH_S2a_df[rownames(NRCH_S2a_df) %in% NRC4other_tree$tip.label,] %>% as.data.frame()
NRC4other_S2b <- NRCH_S2b_df[rownames(NRCH_S2b_df) %in% NRC4other_tree$tip.label,] %>% as.data.frame()

NRC4other_MADA <- NRCH_MADA_df[rownames(NRCH_MADA_df) %in% NRC4other_tree$tip.label,] %>% as.data.frame()
NRC4other_Ploop <- NRCH_Ploop_df[rownames(NRCH_Ploop_df) %in% NRC4other_tree$tip.label,] %>% as.data.frame()
NRC4other_MHD <- NRCH_MHD_df[rownames(NRCH_MHD_df) %in% NRC4other_tree$tip.label,] %>% as.data.frame()


# NRC5
NRC5_S1a <- NRCH_S1a_df[rownames(NRCH_S1a_df) %in% NRC5_tree$tip.label,] %>% as.data.frame()
NRC5_S1b <- NRCH_S1b_df[rownames(NRCH_S1b_df) %in% NRC5_tree$tip.label,] %>% as.data.frame()
NRC5_S1c <- NRCH_S1c_df[rownames(NRCH_S1c_df) %in% NRC5_tree$tip.label,] %>% as.data.frame()
NRC5_S1d <- NRCH_S1d_df[rownames(NRCH_S1d_df) %in% NRC5_tree$tip.label,] %>% as.data.frame()
NRC5_S1e <- NRCH_S1e_df[rownames(NRCH_S1e_df) %in% NRC5_tree$tip.label,] %>% as.data.frame()
NRC5_S1f <- NRCH_S1f_df[rownames(NRCH_S1f_df) %in% NRC5_tree$tip.label,] %>% as.data.frame()

NRC5_S2a <- NRCH_S2a_df[rownames(NRCH_S2a_df) %in% NRC5_tree$tip.label,] %>% as.data.frame()
NRC5_S2b <- NRCH_S2b_df[rownames(NRCH_S2b_df) %in% NRC5_tree$tip.label,] %>% as.data.frame()

NRC5_MADA <- NRCH_MADA_df[rownames(NRCH_MADA_df) %in% NRC5_tree$tip.label,] %>% as.data.frame()
NRC5_Ploop <- NRCH_Ploop_df[rownames(NRCH_Ploop_df) %in% NRC5_tree$tip.label,] %>% as.data.frame()
NRC5_MHD <- NRCH_MHD_df[rownames(NRCH_MHD_df) %in% NRC5_tree$tip.label,] %>% as.data.frame()


# NRC6
NRC6_S1a <- NRCH_S1a_df[rownames(NRCH_S1a_df) %in% NRC6_tree$tip.label,] %>% as.data.frame()
NRC6_S1b <- NRCH_S1b_df[rownames(NRCH_S1b_df) %in% NRC6_tree$tip.label,] %>% as.data.frame()
NRC6_S1c <- NRCH_S1c_df[rownames(NRCH_S1c_df) %in% NRC6_tree$tip.label,] %>% as.data.frame()
NRC6_S1d <- NRCH_S1d_df[rownames(NRCH_S1d_df) %in% NRC6_tree$tip.label,] %>% as.data.frame()
NRC6_S1e <- NRCH_S1e_df[rownames(NRCH_S1e_df) %in% NRC6_tree$tip.label,] %>% as.data.frame()
NRC6_S1f <- NRCH_S1f_df[rownames(NRCH_S1f_df) %in% NRC6_tree$tip.label,] %>% as.data.frame()

NRC6_S2a <- NRCH_S2a_df[rownames(NRCH_S2a_df) %in% NRC6_tree$tip.label,] %>% as.data.frame()
NRC6_S2b <- NRCH_S2b_df[rownames(NRCH_S2b_df) %in% NRC6_tree$tip.label,] %>% as.data.frame()

NRC6_MADA <- NRCH_MADA_df[rownames(NRCH_MADA_df) %in% NRC6_tree$tip.label,] %>% as.data.frame()
NRC6_Ploop <- NRCH_Ploop_df[rownames(NRCH_Ploop_df) %in% NRC6_tree$tip.label,] %>% as.data.frame()
NRC6_MHD <- NRCH_MHD_df[rownames(NRCH_MHD_df) %in% NRC6_tree$tip.label,] %>% as.data.frame()


# NRC7
NRC7_S1a <- NRCH_S1a_df[rownames(NRCH_S1a_df) %in% NRC7_tree$tip.label,] %>% as.data.frame()
NRC7_S1b <- NRCH_S1b_df[rownames(NRCH_S1b_df) %in% NRC7_tree$tip.label,] %>% as.data.frame()
NRC7_S1c <- NRCH_S1c_df[rownames(NRCH_S1c_df) %in% NRC7_tree$tip.label,] %>% as.data.frame()
NRC7_S1d <- NRCH_S1d_df[rownames(NRCH_S1d_df) %in% NRC7_tree$tip.label,] %>% as.data.frame()
NRC7_S1e <- NRCH_S1e_df[rownames(NRCH_S1e_df) %in% NRC7_tree$tip.label,] %>% as.data.frame()
NRC7_S1f <- NRCH_S1f_df[rownames(NRCH_S1f_df) %in% NRC7_tree$tip.label,] %>% as.data.frame()

NRC7_S2a <- NRCH_S2a_df[rownames(NRCH_S2a_df) %in% NRC7_tree$tip.label,] %>% as.data.frame()
NRC7_S2b <- NRCH_S2b_df[rownames(NRCH_S2b_df) %in% NRC7_tree$tip.label,] %>% as.data.frame()

NRC7_MADA <- NRCH_MADA_df[rownames(NRCH_MADA_df) %in% NRC7_tree$tip.label,] %>% as.data.frame()
NRC7_Ploop <- NRCH_Ploop_df[rownames(NRCH_Ploop_df) %in% NRC7_tree$tip.label,] %>% as.data.frame()
NRC7_MHD <- NRCH_MHD_df[rownames(NRCH_MHD_df) %in% NRC7_tree$tip.label,] %>% as.data.frame()


# NRC8
NRC8_S1a <- NRCH_S1a_df[rownames(NRCH_S1a_df) %in% NRC8_tree$tip.label,] %>% as.data.frame()
NRC8_S1b <- NRCH_S1b_df[rownames(NRCH_S1b_df) %in% NRC8_tree$tip.label,] %>% as.data.frame()
NRC8_S1c <- NRCH_S1c_df[rownames(NRCH_S1c_df) %in% NRC8_tree$tip.label,] %>% as.data.frame()
NRC8_S1d <- NRCH_S1d_df[rownames(NRCH_S1d_df) %in% NRC8_tree$tip.label,] %>% as.data.frame()
NRC8_S1e <- NRCH_S1e_df[rownames(NRCH_S1e_df) %in% NRC8_tree$tip.label,] %>% as.data.frame()
NRC8_S1f <- NRCH_S1f_df[rownames(NRCH_S1f_df) %in% NRC8_tree$tip.label,] %>% as.data.frame()

NRC8_S2a <- NRCH_S2a_df[rownames(NRCH_S2a_df) %in% NRC8_tree$tip.label,] %>% as.data.frame()
NRC8_S2b <- NRCH_S2b_df[rownames(NRCH_S2b_df) %in% NRC8_tree$tip.label,] %>% as.data.frame()

NRC8_MADA <- NRCH_MADA_df[rownames(NRCH_MADA_df) %in% NRC8_tree$tip.label,] %>% as.data.frame()
NRC8_Ploop <- NRCH_Ploop_df[rownames(NRCH_Ploop_df) %in% NRC8_tree$tip.label,] %>% as.data.frame()
NRC8_MHD <- NRCH_MHD_df[rownames(NRCH_MHD_df) %in% NRC8_tree$tip.label,] %>% as.data.frame()


# NRC9
NRC9_S1a <- NRCH_S1a_df[rownames(NRCH_S1a_df) %in% NRC9_tree$tip.label,] %>% as.data.frame()
NRC9_S1b <- NRCH_S1b_df[rownames(NRCH_S1b_df) %in% NRC9_tree$tip.label,] %>% as.data.frame()
NRC9_S1c <- NRCH_S1c_df[rownames(NRCH_S1c_df) %in% NRC9_tree$tip.label,] %>% as.data.frame()
NRC9_S1d <- NRCH_S1d_df[rownames(NRCH_S1d_df) %in% NRC9_tree$tip.label,] %>% as.data.frame()
NRC9_S1e <- NRCH_S1e_df[rownames(NRCH_S1e_df) %in% NRC9_tree$tip.label,] %>% as.data.frame()
NRC9_S1f <- NRCH_S1f_df[rownames(NRCH_S1f_df) %in% NRC9_tree$tip.label,] %>% as.data.frame()

NRC9_S2a <- NRCH_S2a_df[rownames(NRCH_S2a_df) %in% NRC9_tree$tip.label,] %>% as.data.frame()
NRC9_S2b <- NRCH_S2b_df[rownames(NRCH_S2b_df) %in% NRC9_tree$tip.label,] %>% as.data.frame()

NRC9_MADA <- NRCH_MADA_df[rownames(NRCH_MADA_df) %in% NRC9_tree$tip.label,] %>% as.data.frame()
NRC9_Ploop <- NRCH_Ploop_df[rownames(NRCH_Ploop_df) %in% NRC9_tree$tip.label,] %>% as.data.frame()
NRC9_MHD <- NRCH_MHD_df[rownames(NRCH_MHD_df) %in% NRC9_tree$tip.label,] %>% as.data.frame()


# NRCX
NRCX_S1a <- NRCH_S1a_df[rownames(NRCH_S1a_df) %in% NRCX_tree$tip.label,] %>% as.data.frame()
NRCX_S1b <- NRCH_S1b_df[rownames(NRCH_S1b_df) %in% NRCX_tree$tip.label,] %>% as.data.frame()
NRCX_S1c <- NRCH_S1c_df[rownames(NRCH_S1c_df) %in% NRCX_tree$tip.label,] %>% as.data.frame()
NRCX_S1d <- NRCH_S1d_df[rownames(NRCH_S1d_df) %in% NRCX_tree$tip.label,] %>% as.data.frame()
NRCX_S1e <- NRCH_S1e_df[rownames(NRCH_S1e_df) %in% NRCX_tree$tip.label,] %>% as.data.frame()
NRCX_S1f <- NRCH_S1f_df[rownames(NRCH_S1f_df) %in% NRCX_tree$tip.label,] %>% as.data.frame()

NRCX_S2a <- NRCH_S2a_df[rownames(NRCH_S2a_df) %in% NRCX_tree$tip.label,] %>% as.data.frame()
NRCX_S2b <- NRCH_S2b_df[rownames(NRCH_S2b_df) %in% NRCX_tree$tip.label,] %>% as.data.frame()

NRCX_MADA <- NRCH_MADA_df[rownames(NRCH_MADA_df) %in% NRCX_tree$tip.label,] %>% as.data.frame()
NRCX_Ploop <- NRCH_Ploop_df[rownames(NRCH_Ploop_df) %in% NRCX_tree$tip.label,] %>% as.data.frame()
NRCX_MHD <- NRCH_MHD_df[rownames(NRCH_MHD_df) %in% NRCX_tree$tip.label,] %>% as.data.frame()





ggseqlogo(NRCH_MADA_df)



ggplot() +
  geom_logo(NRC6_MADA, col_scheme = color_scheme) +
  theme_light() +
  theme_custom



theme_custom <- theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  axis.title = element_blank(),
  legend.position = "none",
  panel.background = element_rect(fill = "transparent", colour = NA), # No fill, no border
  plot.background = element_rect(fill = "transparent", colour = NA), # No fill, no border
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  plot.margin = unit(c(0,0,0,0), "lines"),
  plot.title = element_blank(),
  plot.subtitle = element_blank(),
  plot.caption = element_blank(),
  panel.border = element_blank()
)



color_scheme <- make_col_scheme(chars = c("K","R","H",
                                          "D","E",
                                          "A","V","L","I","W","F","M","Y",
                                          "S","T","N","Q",
                                          "P","C","G"),
                                groups = c("1","1","1",
                                           "2","2",
                                           "3","3","3","3","3","3","3","3",
                                           "4","4","4","4",
                                           "5","5","5"),
                                cols = c("#FF1F1F","#FF1F1F","#FF1F1F",
                                         "#CC00CC","#CC00CC",
                                         "#3333FF","#3333FF","#3333FF","#3333FF","#3333FF","#3333FF","#3333FF","#3333FF",
                                         "#1FFF1F","#1FFF1F","#1FFF1F","#1FFF1F",
                                         "#FFAE00","#FFAE00","#FFAE00"))





## export the motif sequence logo
setwd("~/Desktop/NRC2_dimer/analysis/v3/figures/motifs/v3/")

# loop through each dataframe within the list
for (name in names(motifs)) {
  df <- motifs[[name]]
  
  # create the plot
  plot <- ggplot() +
    geom_logo(df, col_scheme = color_scheme) +
    theme_light() +
    theme_custom
  
  ggsave(filename = paste0(name, ".png"), plot = plot, width = 10, height = 2)
}







# each motif separately

S1a_motifs <- list(
  NRCH_S1a = NRCH_S1a_df,
  NRC0_S1a = 	NRC0_S1a,
  NRC1_S1a = 	NRC1_S1a,
  NRC2_S1a = 	NRC2_S1a,
  NRC3_S1a = 	NRC3_S1a,
  NRC4a_S1a = 	NRC4a_S1a,
  NRC4b_S1a = 	NRC4b_S1a,
  NRC4c_S1a = 	NRC4c_S1a,
  NRC4t_S1a = 	NRC4t_S1a,
  NRC4other_S1a = 	NRC4other_S1a,
  NRC5_S1a = 	NRC5_S1a,
  NRC6_S1a = 	NRC6_S1a,
  NRC7_S1a = 	NRC7_S1a,
  NRC8_S1a = 	NRC8_S1a,
  NRC9_S1a = 	NRC9_S1a,
  NRCX_S1a = 	NRCX_S1a
)



S1b_motifs <- list(
  NRCH_S1b = NRCH_S1b_df,
  NRC0_S1b = 	NRC0_S1b,
  NRC1_S1b = 	NRC1_S1b,
  NRC2_S1b = 	NRC2_S1b,
  NRC3_S1b = 	NRC3_S1b,
  NRC4a_S1b = 	NRC4a_S1b,
  NRC4b_S1b = 	NRC4b_S1b,
  NRC4c_S1b = 	NRC4c_S1b,
  NRC4t_S1b = 	NRC4t_S1b,
  NRC4other_S1b = 	NRC4other_S1b,
  NRC5_S1b = 	NRC5_S1b,
  NRC6_S1b = 	NRC6_S1b,
  NRC7_S1b = 	NRC7_S1b,
  NRC8_S1b = 	NRC8_S1b,
  NRC9_S1b = 	NRC9_S1b,
  NRCX_S1b = 	NRCX_S1b
)



S1c_motifs <- list(
  NRCH_S1c = NRCH_S1c_df,
  NRC0_S1c = 	NRC0_S1c,
  NRC1_S1c = 	NRC1_S1c,
  NRC2_S1c = 	NRC2_S1c,
  NRC3_S1c = 	NRC3_S1c,
  NRC4a_S1c = 	NRC4a_S1c,
  NRC4b_S1c = 	NRC4b_S1c,
  NRC4c_S1c = 	NRC4c_S1c,
  NRC4t_S1c = 	NRC4t_S1c,
  NRC4other_S1c = 	NRC4other_S1c,
  NRC5_S1c = 	NRC5_S1c,
  NRC6_S1c = 	NRC6_S1c,
  NRC7_S1c = 	NRC7_S1c,
  NRC8_S1c = 	NRC8_S1c,
  NRC9_S1c = 	NRC9_S1c,
  NRCX_S1c = 	NRCX_S1c
)



S1d_motifs <- list(
  NRCH_S1d = NRCH_S1d_df,
  NRC0_S1d = 	NRC0_S1d,
  NRC1_S1d = 	NRC1_S1d,
  NRC2_S1d = 	NRC2_S1d,
  NRC3_S1d = 	NRC3_S1d,
  NRC4a_S1d = 	NRC4a_S1d,
  NRC4b_S1d = 	NRC4b_S1d,
  NRC4c_S1d = 	NRC4c_S1d,
  NRC4t_S1d = 	NRC4t_S1d,
  NRC4other_S1d = 	NRC4other_S1d,
  NRC5_S1d = 	NRC5_S1d,
  NRC6_S1d = 	NRC6_S1d,
  NRC7_S1d = 	NRC7_S1d,
  NRC8_S1d = 	NRC8_S1d,
  NRC9_S1d = 	NRC9_S1d,
  NRCX_S1d = 	NRCX_S1d
)



S1e_motifs <- list(
  NRCH_S1e = NRCH_S1e_df,
  NRC0_S1e = 	NRC0_S1e,
  NRC1_S1e = 	NRC1_S1e,
  NRC2_S1e = 	NRC2_S1e,
  NRC3_S1e = 	NRC3_S1e,
  NRC4a_S1e = 	NRC4a_S1e,
  NRC4b_S1e = 	NRC4b_S1e,
  NRC4c_S1e = 	NRC4c_S1e,
  NRC4t_S1e = 	NRC4t_S1e,
  NRC4other_S1e = 	NRC4other_S1e,
  NRC5_S1e = 	NRC5_S1e,
  NRC6_S1e = 	NRC6_S1e,
  NRC7_S1e = 	NRC7_S1e,
  NRC8_S1e = 	NRC8_S1e,
  NRC9_S1e = 	NRC9_S1e,
  NRCX_S1e = 	NRCX_S1e
)



S1f_motifs <- list(
  NRCH_S1f = NRCH_S1f_df,
  NRC0_S1f = 	NRC0_S1f,
  NRC1_S1f = 	NRC1_S1f,
  NRC2_S1f = 	NRC2_S1f,
  NRC3_S1f = 	NRC3_S1f,
  NRC4a_S1f = 	NRC4a_S1f,
  NRC4b_S1f = 	NRC4b_S1f,
  NRC4c_S1f = 	NRC4c_S1f,
  NRC4t_S1f = 	NRC4t_S1f,
  NRC4other_S1f = 	NRC4other_S1f,
  NRC5_S1f = 	NRC5_S1f,
  NRC6_S1f = 	NRC6_S1f,
  NRC7_S1f = 	NRC7_S1f,
  NRC8_S1f = 	NRC8_S1f,
  NRC9_S1f = 	NRC9_S1f,
  NRCX_S1f = 	NRCX_S1f
)



S2a_motifs <- list(
  NRCH_S2a = NRCH_S2a_df,
  NRC0_S2a = 	NRC0_S2a,
  NRC1_S2a = 	NRC1_S2a,
  NRC2_S2a = 	NRC2_S2a,
  NRC3_S2a = 	NRC3_S2a,
  NRC4a_S2a = 	NRC4a_S2a,
  NRC4b_S2a = 	NRC4b_S2a,
  NRC4c_S2a = 	NRC4c_S2a,
  NRC4t_S2a = 	NRC4t_S2a,
  NRC4other_S2a = 	NRC4other_S2a,
  NRC5_S2a = 	NRC5_S2a,
  NRC6_S2a = 	NRC6_S2a,
  NRC7_S2a = 	NRC7_S2a,
  NRC8_S2a = 	NRC8_S2a,
  NRC9_S2a = 	NRC9_S2a,
  NRCX_S2a = 	NRCX_S2a
)



S2b_motifs <- list(
  NRCH_S2b = NRCH_S2b_df,
  NRC0_S2b = 	NRC0_S2b,
  NRC1_S2b = 	NRC1_S2b,
  NRC2_S2b = 	NRC2_S2b,
  NRC3_S2b = 	NRC3_S2b,
  NRC4a_S2b = 	NRC4a_S2b,
  NRC4b_S2b = 	NRC4b_S2b,
  NRC4c_S2b = 	NRC4c_S2b,
  NRC4t_S2b = 	NRC4t_S2b,
  NRC4other_S2b = 	NRC4other_S2b,
  NRC5_S2b = 	NRC5_S2b,
  NRC6_S2b = 	NRC6_S2b,
  NRC7_S2b = 	NRC7_S2b,
  NRC8_S2b = 	NRC8_S2b,
  NRC9_S2b = 	NRC9_S2b,
  NRCX_S2b = 	NRCX_S2b
)




MADA_motifs <- list(
  NRCH_MADA = NRCH_MADA_df,
  NRC0_MADA = 	NRC0_MADA,
  NRC1_MADA = 	NRC1_MADA,
  NRC2_MADA = 	NRC2_MADA,
  NRC3_MADA = 	NRC3_MADA,
  NRC4a_MADA = 	NRC4a_MADA,
  NRC4b_MADA = 	NRC4b_MADA,
  NRC4c_MADA = 	NRC4c_MADA,
  NRC4t_MADA = 	NRC4t_MADA,
  NRC4other_MADA = 	NRC4other_MADA,
  NRC5_MADA = 	NRC5_MADA,
  NRC6_MADA = 	NRC6_MADA,
  NRC7_MADA = 	NRC7_MADA,
  NRC8_MADA = 	NRC8_MADA,
  NRC9_MADA = 	NRC9_MADA,
  NRCX_MADA = 	NRCX_MADA
)



Ploop_motifs <- list(
  NRCH_Ploop = NRCH_Ploop_df,
  NRC0_Ploop = 	NRC0_Ploop,
  NRC1_Ploop = 	NRC1_Ploop,
  NRC2_Ploop = 	NRC2_Ploop,
  NRC3_Ploop = 	NRC3_Ploop,
  NRC4a_Ploop = 	NRC4a_Ploop,
  NRC4b_Ploop = 	NRC4b_Ploop,
  NRC4c_Ploop = 	NRC4c_Ploop,
  NRC4t_Ploop = 	NRC4t_Ploop,
  NRC4other_Ploop = 	NRC4other_Ploop,
  NRC5_Ploop = 	NRC5_Ploop,
  NRC6_Ploop = 	NRC6_Ploop,
  NRC7_Ploop = 	NRC7_Ploop,
  NRC8_Ploop = 	NRC8_Ploop,
  NRC9_Ploop = 	NRC9_Ploop,
  NRCX_Ploop = 	NRCX_Ploop
)



MHD_motifs <- list(
  NRCH_MHD = NRCH_MHD_df,
  NRC0_MHD = 	NRC0_MHD,
  NRC1_MHD = 	NRC1_MHD,
  NRC2_MHD = 	NRC2_MHD,
  NRC3_MHD = 	NRC3_MHD,
  NRC4a_MHD = 	NRC4a_MHD,
  NRC4b_MHD = 	NRC4b_MHD,
  NRC4c_MHD = 	NRC4c_MHD,
  NRC4t_MHD = 	NRC4t_MHD,
  NRC4other_MHD = 	NRC4other_MHD,
  NRC5_MHD = 	NRC5_MHD,
  NRC6_MHD = 	NRC6_MHD,
  NRC7_MHD = 	NRC7_MHD,
  NRC8_MHD = 	NRC8_MHD,
  NRC9_MHD = 	NRC9_MHD,
  NRCX_MHD = 	NRCX_MHD
)






# loop through Stretch 1a
for (name in names(S1a_motifs)) {
  df <- S1a_motifs[[name]]
  
  # create the plot
  plot <- ggplot() +
    geom_logo(df, col_scheme = color_scheme) +
    theme_light() +
    theme_custom
  
  ggsave(filename = paste0(name, ".png"), plot = plot, width = 4, height = 2)
}


# loop through Stretch 1b
for (name in names(S1b_motifs)) {
  df <- S1b_motifs[[name]]
  
  # create the plot
  plot <- ggplot() +
    geom_logo(df, col_scheme = color_scheme) +
    theme_light() +
    theme_custom
  
  ggsave(filename = paste0(name, ".png"), plot = plot, width = 7, height = 2)
}


# loop through Stretch 1c
for (name in names(S1c_motifs)) {
  df <- S1c_motifs[[name]]
  
  # create the plot
  plot <- ggplot() +
    geom_logo(df, col_scheme = color_scheme) +
    theme_light() +
    theme_custom
  
  ggsave(filename = paste0(name, ".png"), plot = plot, width = 5, height = 2)
}


# loop through Stretch 1d
for (name in names(S1d_motifs)) {
  df <- S1d_motifs[[name]]
  
  # create the plot
  plot <- ggplot() +
    geom_logo(df, col_scheme = color_scheme) +
    theme_light() +
    theme_custom
  
  ggsave(filename = paste0(name, ".png"), plot = plot, width = 8, height = 2)
}



# loop through Stretch 1e
for (name in names(S1e_motifs)) {
  df <- S1e_motifs[[name]]
  
  # create the plot
  plot <- ggplot() +
    geom_logo(df, col_scheme = color_scheme) +
    theme_light() +
    theme_custom
  
  ggsave(filename = paste0(name, ".png"), plot = plot, width = 4, height = 2)
}



# loop through Stretch 1f
for (name in names(S1f_motifs)) {
  df <- S1f_motifs[[name]]
  
  # create the plot
  plot <- ggplot() +
    geom_logo(df, col_scheme = color_scheme) +
    theme_light() +
    theme_custom
  
  ggsave(filename = paste0(name, ".png"), plot = plot, width = 5, height = 2)
}



# loop through Stretch 2a
for (name in names(S2a_motifs)) {
  df <- S2a_motifs[[name]]
  
  # create the plot
  plot <- ggplot() +
    geom_logo(df, col_scheme = color_scheme) +
    theme_light() +
    theme_custom
  
  ggsave(filename = paste0(name, ".png"), plot = plot, width = 6, height = 2)
}



# loop through Stretch 2b
for (name in names(S2b_motifs)) {
  df <- S2b_motifs[[name]]
  
  # create the plot
  plot <- ggplot() +
    geom_logo(df, col_scheme = color_scheme) +
    theme_light() +
    theme_custom
  
  ggsave(filename = paste0(name, ".png"), plot = plot, width = 8, height = 2)
}



# loop through MADA
for (name in names(MADA_motifs)) {
  df <- MADA_motifs[[name]]
  
  # create the plot
  plot <- ggplot() +
    geom_logo(df, col_scheme = color_scheme) +
    theme_light() +
    theme_custom
  
  ggsave(filename = paste0(name, ".png"), plot = plot, width = 21, height = 2)
}


# loop through Ploop
for (name in names(Ploop_motifs)) {
  df <- Ploop_motifs[[name]]
  
  # create the plot
  plot <- ggplot() +
    geom_logo(df, col_scheme = color_scheme) +
    theme_light() +
    theme_custom
  
  ggsave(filename = paste0(name, ".png"), plot = plot, width = 11, height = 2)
}



# loop through MHD
for (name in names(MHD_motifs)) {
  df <- MHD_motifs[[name]]
  
  # create the plot
  plot <- ggplot() +
    geom_logo(df, col_scheme = color_scheme) +
    theme_light() +
    theme_custom
  
  ggsave(filename = paste0(name, ".png"), plot = plot, width = 3, height = 2)
}




### let's calculate the genus distribution in each clade

# import the NLR metadata file
NLR_meta <- read_excel("/path/to/DataS6.xlsx")

# import the genome metadata file
Genome_data_v3 <- read_excel("/path/to/DataS2.xlsx")

# merge both
NLR_metadata <- NLR_meta %>% left_join(Genome_data_v3, by = "file_name")


# extract the metadata for each clade
NRC0_meta <- NLR_metadata[NLR_metadata$seqname %in% NRC0_tree$tip.label,]
NRC1_meta <- NLR_metadata[NLR_metadata$seqname %in% NRC1_tree$tip.label,]
NRC2_meta <- NLR_metadata[NLR_metadata$seqname %in% NRC2_tree$tip.label,]
NRC3_meta <- NLR_metadata[NLR_metadata$seqname %in% NRC3_tree$tip.label,]
NRCX_meta <- NLR_metadata[NLR_metadata$seqname %in% NRCX_tree$tip.label,]
NRC4a_meta <- NLR_metadata[NLR_metadata$seqname %in% NRC4a_tree$tip.label,]
NRC4b_meta <- NLR_metadata[NLR_metadata$seqname %in% NRC4b_tree$tip.label,]
NRC4c_meta <- NLR_metadata[NLR_metadata$seqname %in% NRC4c_tree$tip.label,]
NRC4t_meta <- NLR_metadata[NLR_metadata$seqname %in% NRC4t_tree$tip.label,]
NRC4other_meta <- NLR_metadata[NLR_metadata$seqname %in% NRC4other_tree$tip.label,]
NRC5_meta <- NLR_metadata[NLR_metadata$seqname %in% NRC5_tree$tip.label,]
NRC6_meta <- NLR_metadata[NLR_metadata$seqname %in% NRC6_tree$tip.label,]
NRC7_meta <- NLR_metadata[NLR_metadata$seqname %in% NRC7_tree$tip.label,]
NRC8_meta <- NLR_metadata[NLR_metadata$seqname %in% NRC8_tree$tip.label,]
NRC9_meta <- NLR_metadata[NLR_metadata$seqname %in% NRC9_tree$tip.label,]






