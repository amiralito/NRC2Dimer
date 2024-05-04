# the outputs are going to be exported here
setwd("/path/to/your/figure/directory")

theme_custom2 <- theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  legend.position = "none",
  panel.background = element_blank(),
  plot.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  plot.title = element_blank(),
  plot.subtitle = element_blank(),
  plot.caption = element_blank(),
  axis.line = element_line(colour = "black", size = 0.25)
)

# extract the alignment for each clade
# calculate and visualize the entropy

# NRCH
NRCH_full_pos <- convert_to_dataframe(NRCH_clu_filtered_ref_trimmed)
NRCH_full_ent <- calculate_entropy_df(NRCH_full_pos)

# annotate the file based on motifs and stretches
NRCH_full_ent <- NRCH_full_ent %>%
  mutate(color_segment = case_when(
    Position >= 234 & Position <= 237 |
      Position >= 256 & Position <= 262 |
      Position >= 288 & Position <= 292 |
      Position >= 522 & Position <= 529 |
      Position >= 550 & Position <= 553 | 
      Position >= 578 & Position <= 582 ~ "stretch 1",
    Position >= 544 & Position <= 549 |
      Position >= 566 & Position <= 573 ~ "stretch 2",
    Position >= 1 & Position <= 32 | 
      Position >= 199 & Position <= 209 | 
      Position >= 495 & Position <= 497 ~ "motifs",
    TRUE ~ "normal"
  ))

NRCH_full_ent_p <- ggplot(NRCH_full_ent, aes(x = Position, y = Entropy)) +
  geom_segment(aes(x = Position, xend = Position, y = 0, yend = Entropy, color = color_segment)) +
  scale_color_manual(values = c("stretch 1" = "#45B64B", "stretch 2" = "#EAC435", "motifs" = "#AFAFAF", "normal" = "black")) +
  ylim(c(0,4.32)) +
  xlim(c(0,923)) +
  geom_hline(yintercept = 1.5, linetype = "dotted", color = "#6C6C6C", size = 0.5, alpha = 1) +
  theme_custom2

ggsave(filename = "NRCH_full.png", NRCH_full_ent_p, width = 10, height = 1.5, dpi = "retina")




# NRC0
NRC0_align <- NRCH_clu_filtered_ref_trimmed[NRCH_clu_filtered_ref_trimmed@ranges@NAMES %in% NRC0_tree$tip.label]

NRC0_full_pos <- convert_to_dataframe(NRC0_align)
NRC0_full_ent <- calculate_entropy_df(NRC0_full_pos)

# annotate the file based on motifs and stretches
NRC0_full_ent <- NRC0_full_ent %>%
  mutate(color_segment = case_when(
    Position >= 234 & Position <= 237 |
      Position >= 256 & Position <= 262 |
      Position >= 288 & Position <= 292 |
      Position >= 522 & Position <= 529 |
      Position >= 550 & Position <= 553 | 
      Position >= 578 & Position <= 582 ~ "stretch 1",
    Position >= 544 & Position <= 549 |
      Position >= 566 & Position <= 573 ~ "stretch 2",
    Position >= 1 & Position <= 32 | 
      Position >= 199 & Position <= 209 | 
      Position >= 495 & Position <= 497 ~ "motifs",
    TRUE ~ "normal"
  ))

NRC0_full_ent_p <- ggplot(NRC0_full_ent, aes(x = Position, y = Entropy)) +
  geom_segment(aes(x = Position, xend = Position, y = 0, yend = Entropy, color = color_segment)) +
  scale_color_manual(values = c("stretch 1" = "#45B64B", "stretch 2" = "#EAC435", "motifs" = "#AFAFAF", "normal" = "black")) +
  ylim(c(0,4.32)) +
  xlim(c(0,923)) +
  geom_hline(yintercept = 1.5, linetype = "dotted", color = "#6C6C6C", size = 0.5, alpha = 1) +
  theme_custom2

ggsave(filename = "NRC0_full.png", NRC0_full_ent_p, width = 10, height = 1.5, dpi = "retina")




# NRC1
NRC1_align <- NRCH_clu_filtered_ref_trimmed[NRCH_clu_filtered_ref_trimmed@ranges@NAMES %in% NRC1_tree$tip.label]

NRC1_full_pos <- convert_to_dataframe(NRC1_align)
NRC1_full_ent <- calculate_entropy_df(NRC1_full_pos)

# annotate the file based on motifs and stretches
NRC1_full_ent <- NRC1_full_ent %>%
  mutate(color_segment = case_when(
    Position >= 234 & Position <= 237 |
      Position >= 256 & Position <= 262 |
      Position >= 288 & Position <= 292 |
      Position >= 522 & Position <= 529 |
      Position >= 550 & Position <= 553 | 
      Position >= 578 & Position <= 582 ~ "stretch 1",
    Position >= 544 & Position <= 549 |
      Position >= 566 & Position <= 573 ~ "stretch 2",
    Position >= 1 & Position <= 32 | 
      Position >= 199 & Position <= 209 | 
      Position >= 495 & Position <= 497 ~ "motifs",
    TRUE ~ "normal"
  ))

NRC1_full_ent_p <- ggplot(NRC1_full_ent, aes(x = Position, y = Entropy)) +
  geom_segment(aes(x = Position, xend = Position, y = 0, yend = Entropy, color = color_segment)) +
  scale_color_manual(values = c("stretch 1" = "#45B64B", "stretch 2" = "#EAC435", "motifs" = "#AFAFAF", "normal" = "black")) +
  ylim(c(0,4.32)) +
  xlim(c(0,923)) +
  geom_hline(yintercept = 1.5, linetype = "dotted", color = "#6C6C6C", size = 0.5, alpha = 1) +
  theme_custom2

ggsave(filename = "NRC1_full.png", NRC1_full_ent_p, width = 10, height = 1.5, dpi = "retina")




# NRC2
NRC2_align <- NRCH_clu_filtered_ref_trimmed[NRCH_clu_filtered_ref_trimmed@ranges@NAMES %in% NRC2_tree$tip.label]

NRC2_full_pos <- convert_to_dataframe(NRC2_align)
NRC2_full_ent <- calculate_entropy_df(NRC2_full_pos)

# annotate the file based on motifs and stretches
NRC2_full_ent <- NRC2_full_ent %>%
  mutate(color_segment = case_when(
    Position >= 234 & Position <= 237 |
      Position >= 256 & Position <= 262 |
      Position >= 288 & Position <= 292 |
      Position >= 522 & Position <= 529 |
      Position >= 550 & Position <= 553 | 
      Position >= 578 & Position <= 582 ~ "stretch 1",
    Position >= 544 & Position <= 549 |
      Position >= 566 & Position <= 573 ~ "stretch 2",
    Position >= 1 & Position <= 32 | 
      Position >= 199 & Position <= 209 | 
      Position >= 495 & Position <= 497 ~ "motifs",
    TRUE ~ "normal"
  ))

NRC2_full_ent_p <- ggplot(NRC2_full_ent, aes(x = Position, y = Entropy)) +
  geom_segment(aes(x = Position, xend = Position, y = 0, yend = Entropy, color = color_segment)) +
  scale_color_manual(values = c("stretch 1" = "#45B64B", "stretch 2" = "#EAC435", "motifs" = "#AFAFAF", "normal" = "black")) +
  ylim(c(0,4.32)) +
  xlim(c(0,923)) +
  geom_hline(yintercept = 1.5, linetype = "dotted", color = "#6C6C6C", size = 0.5, alpha = 1) +
  theme_custom2

ggsave(filename = "NRC2_full.png", NRC2_full_ent_p, width = 10, height = 1.5, dpi = "retina")




# NRC3
NRC3_align <- NRCH_clu_filtered_ref_trimmed[NRCH_clu_filtered_ref_trimmed@ranges@NAMES %in% NRC3_tree$tip.label]

NRC3_full_pos <- convert_to_dataframe(NRC3_align)
NRC3_full_ent <- calculate_entropy_df(NRC3_full_pos)


# annotate the file based on motifs and stretches
NRC3_full_ent <- NRC3_full_ent %>%
  mutate(color_segment = case_when(
    Position >= 234 & Position <= 237 |
      Position >= 256 & Position <= 262 |
      Position >= 288 & Position <= 292 |
      Position >= 522 & Position <= 529 |
      Position >= 550 & Position <= 553 | 
      Position >= 578 & Position <= 582 ~ "stretch 1",
    Position >= 544 & Position <= 549 |
      Position >= 566 & Position <= 573 ~ "stretch 2",
    Position >= 1 & Position <= 32 | 
      Position >= 199 & Position <= 209 | 
      Position >= 495 & Position <= 497 ~ "motifs",
    TRUE ~ "normal"
  ))

NRC3_full_ent_p <- ggplot(NRC3_full_ent, aes(x = Position, y = Entropy)) +
  geom_segment(aes(x = Position, xend = Position, y = 0, yend = Entropy, color = color_segment)) +
  scale_color_manual(values = c("stretch 1" = "#45B64B", "stretch 2" = "#EAC435", "motifs" = "#AFAFAF", "normal" = "black")) +
  ylim(c(0,4.32)) +
  xlim(c(0,923)) +
  geom_hline(yintercept = 1.5, linetype = "dotted", color = "#6C6C6C", size = 0.5, alpha = 1) +
  theme_custom2

ggsave(filename = "NRC3_full.png", NRC3_full_ent_p, width = 10, height = 1.5, dpi = "retina")




# NRCX
NRCX_align <- NRCH_clu_filtered_ref_trimmed[NRCH_clu_filtered_ref_trimmed@ranges@NAMES %in% NRCX_tree$tip.label]

NRCX_full_pos <- convert_to_dataframe(NRCX_align)
NRCX_full_ent <- calculate_entropy_df(NRCX_full_pos)

# annotate the file based on motifs and stretches
NRCX_full_ent <- NRCX_full_ent %>%
  mutate(color_segment = case_when(
    Position >= 234 & Position <= 237 |
      Position >= 256 & Position <= 262 |
      Position >= 288 & Position <= 292 |
      Position >= 522 & Position <= 529 |
      Position >= 550 & Position <= 553 | 
      Position >= 578 & Position <= 582 ~ "stretch 1",
    Position >= 544 & Position <= 549 |
      Position >= 566 & Position <= 573 ~ "stretch 2",
    Position >= 1 & Position <= 32 | 
      Position >= 199 & Position <= 209 | 
      Position >= 495 & Position <= 497 ~ "motifs",
    TRUE ~ "normal"
  ))

NRCX_full_ent_p <- ggplot(NRCX_full_ent, aes(x = Position, y = Entropy)) +
  geom_segment(aes(x = Position, xend = Position, y = 0, yend = Entropy, color = color_segment)) +
  scale_color_manual(values = c("stretch 1" = "#45B64B", "stretch 2" = "#EAC435", "motifs" = "#AFAFAF", "normal" = "black")) +
  ylim(c(0,4.32)) +
  xlim(c(0,923)) +
  geom_hline(yintercept = 1.5, linetype = "dotted", color = "#6C6C6C", size = 0.5, alpha = 1) +
  theme_custom2

ggsave(filename = "NRCX_full.png", 
       NRCX_full_ent_p, width = 10, height = 1.5, dpi = "retina")





# NRC9
NRC9_align <- NRCH_clu_filtered_ref_trimmed[NRCH_clu_filtered_ref_trimmed@ranges@NAMES %in% NRC9_tree$tip.label]

NRC9_full_pos <- convert_to_dataframe(NRC9_align)
NRC9_full_ent <- calculate_entropy_df(NRC9_full_pos)

# annotate the file based on motifs and stretches
NRC9_full_ent <- NRC9_full_ent %>%
  mutate(color_segment = case_when(
    Position >= 234 & Position <= 237 |
      Position >= 256 & Position <= 262 |
      Position >= 288 & Position <= 292 |
      Position >= 522 & Position <= 529 |
      Position >= 550 & Position <= 553 | 
      Position >= 578 & Position <= 582 ~ "stretch 1",
    Position >= 544 & Position <= 549 |
      Position >= 566 & Position <= 573 ~ "stretch 2",
    Position >= 1 & Position <= 32 | 
      Position >= 199 & Position <= 209 | 
      Position >= 495 & Position <= 497 ~ "motifs",
    TRUE ~ "normal"
  ))

NRC9_full_ent_p <- ggplot(NRC9_full_ent, aes(x = Position, y = Entropy)) +
  geom_segment(aes(x = Position, xend = Position, y = 0, yend = Entropy, color = color_segment)) +
  scale_color_manual(values = c("stretch 1" = "#45B64B", "stretch 2" = "#EAC435", "motifs" = "#AFAFAF", "normal" = "black")) +
  ylim(c(0,4.32)) +
  xlim(c(0,923)) +
  geom_hline(yintercept = 1.5, linetype = "dotted", color = "#6C6C6C", size = 0.5, alpha = 1) +
  theme_custom2

ggsave(filename = "NRC9_full.png", 
       NRC9_full_ent_p, width = 10, height = 1.5, dpi = "retina")





# NRC8
NRC8_align <- NRCH_clu_filtered_ref_trimmed[NRCH_clu_filtered_ref_trimmed@ranges@NAMES %in% NRC8_tree$tip.label]

NRC8_full_pos <- convert_to_dataframe(NRC8_align)
NRC8_full_ent <- calculate_entropy_df(NRC8_full_pos)

# annotate the file based on motifs and stretches
NRC8_full_ent <- NRC8_full_ent %>%
  mutate(color_segment = case_when(
    Position >= 234 & Position <= 237 |
      Position >= 256 & Position <= 262 |
      Position >= 288 & Position <= 292 |
      Position >= 522 & Position <= 529 |
      Position >= 550 & Position <= 553 | 
      Position >= 578 & Position <= 582 ~ "stretch 1",
    Position >= 544 & Position <= 549 |
      Position >= 566 & Position <= 573 ~ "stretch 2",
    Position >= 1 & Position <= 32 | 
      Position >= 199 & Position <= 209 | 
      Position >= 495 & Position <= 497 ~ "motifs",
    TRUE ~ "normal"
  ))

NRC8_full_ent_p <- ggplot(NRC8_full_ent, aes(x = Position, y = Entropy)) +
  geom_segment(aes(x = Position, xend = Position, y = 0, yend = Entropy, color = color_segment)) +
  scale_color_manual(values = c("stretch 1" = "#45B64B", "stretch 2" = "#EAC435", "motifs" = "#AFAFAF", "normal" = "black")) +
  ylim(c(0,4.32)) +
  xlim(c(0,923)) +
  geom_hline(yintercept = 1.5, linetype = "dotted", color = "#6C6C6C", size = 0.5, alpha = 1) +
  theme_custom2

ggsave(filename = "NRC8_full.png", 
       NRC8_full_ent_p, width = 10, height = 1.5, dpi = "retina")





# NRC6
NRC6_align <- NRCH_clu_filtered_ref_trimmed[NRCH_clu_filtered_ref_trimmed@ranges@NAMES %in% NRC6_tree$tip.label]

NRC6_full_pos <- convert_to_dataframe(NRC6_align)
NRC6_full_ent <- calculate_entropy_df(NRC6_full_pos)

# annotate the file based on motifs and stretches
NRC6_full_ent <- NRC6_full_ent %>%
  mutate(color_segment = case_when(
    Position >= 234 & Position <= 237 |
      Position >= 256 & Position <= 262 |
      Position >= 288 & Position <= 292 |
      Position >= 522 & Position <= 529 |
      Position >= 550 & Position <= 553 | 
      Position >= 578 & Position <= 582 ~ "stretch 1",
    Position >= 544 & Position <= 549 |
      Position >= 566 & Position <= 573 ~ "stretch 2",
    Position >= 1 & Position <= 32 | 
      Position >= 199 & Position <= 209 | 
      Position >= 495 & Position <= 497 ~ "motifs",
    TRUE ~ "normal"
  ))

NRC6_full_ent_p <- ggplot(NRC6_full_ent, aes(x = Position, y = Entropy)) +
  geom_segment(aes(x = Position, xend = Position, y = 0, yend = Entropy, color = color_segment)) +
  scale_color_manual(values = c("stretch 1" = "#45B64B", "stretch 2" = "#EAC435", "motifs" = "#AFAFAF", "normal" = "black")) +
  ylim(c(0,4.32)) +
  xlim(c(0,923)) +
  geom_hline(yintercept = 1.5, linetype = "dotted", color = "#6C6C6C", size = 0.5, alpha = 1) +
  theme_custom2

ggsave(filename = "NRC6_full.png", 
       NRC6_full_ent_p, width = 10, height = 1.5, dpi = "retina")






# NRC7
NRC7_align <- NRCH_clu_filtered_ref_trimmed[NRCH_clu_filtered_ref_trimmed@ranges@NAMES %in% NRC7_tree$tip.label]

NRC7_full_pos <- convert_to_dataframe(NRC7_align)
NRC7_full_ent <- calculate_entropy_df(NRC7_full_pos)

# annotate the file based on motifs and stretches
NRC7_full_ent <- NRC7_full_ent %>%
  mutate(color_segment = case_when(
    Position >= 234 & Position <= 237 |
      Position >= 256 & Position <= 262 |
      Position >= 288 & Position <= 292 |
      Position >= 522 & Position <= 529 |
      Position >= 550 & Position <= 553 | 
      Position >= 578 & Position <= 582 ~ "stretch 1",
    Position >= 544 & Position <= 549 |
      Position >= 566 & Position <= 573 ~ "stretch 2",
    Position >= 1 & Position <= 32 | 
      Position >= 199 & Position <= 209 | 
      Position >= 495 & Position <= 497 ~ "motifs",
    TRUE ~ "normal"
  ))

NRC7_full_ent_p <- ggplot(NRC7_full_ent, aes(x = Position, y = Entropy)) +
  geom_segment(aes(x = Position, xend = Position, y = 0, yend = Entropy, color = color_segment)) +
  scale_color_manual(values = c("stretch 1" = "#45B64B", "stretch 2" = "#EAC435", "motifs" = "#AFAFAF", "normal" = "black")) +
  ylim(c(0,4.32)) +
  xlim(c(0,923)) +
  geom_hline(yintercept = 1.5, linetype = "dotted", color = "#6C6C6C", size = 0.5, alpha = 1) +
  theme_custom2

ggsave(filename = "NRC7_full.png", 
       NRC7_full_ent_p, width = 10, height = 1.5, dpi = "retina")





# NRC5
NRC5_align <- NRCH_clu_filtered_ref_trimmed[NRCH_clu_filtered_ref_trimmed@ranges@NAMES %in% NRC5_tree$tip.label]

NRC5_full_pos <- convert_to_dataframe(NRC5_align)
NRC5_full_ent <- calculate_entropy_df(NRC5_full_pos)

# annotate the file based on motifs and stretches
NRC5_full_ent <- NRC5_full_ent %>%
  mutate(color_segment = case_when(
    Position >= 234 & Position <= 237 |
      Position >= 256 & Position <= 262 |
      Position >= 288 & Position <= 292 |
      Position >= 522 & Position <= 529 |
      Position >= 550 & Position <= 553 | 
      Position >= 578 & Position <= 582 ~ "stretch 1",
    Position >= 544 & Position <= 549 |
      Position >= 566 & Position <= 573 ~ "stretch 2",
    Position >= 1 & Position <= 32 | 
      Position >= 199 & Position <= 209 | 
      Position >= 495 & Position <= 497 ~ "motifs",
    TRUE ~ "normal"
  ))

NRC5_full_ent_p <- ggplot(NRC5_full_ent, aes(x = Position, y = Entropy)) +
  geom_segment(aes(x = Position, xend = Position, y = 0, yend = Entropy, color = color_segment)) +
  scale_color_manual(values = c("stretch 1" = "#45B64B", "stretch 2" = "#EAC435", "motifs" = "#AFAFAF", "normal" = "black")) +
  ylim(c(0,4.32)) +
  xlim(c(0,923)) +
  geom_hline(yintercept = 1.5, linetype = "dotted", color = "#6C6C6C", size = 0.5, alpha = 1) +
  theme_custom2

ggsave(filename = "NRC5_full.png", 
       NRC5_full_ent_p, width = 10, height = 1.5, dpi = "retina")





# NRC4c
NRC4c_align <- NRCH_clu_filtered_ref_trimmed[NRCH_clu_filtered_ref_trimmed@ranges@NAMES %in% NRC4c_tree$tip.label]

NRC4c_full_pos <- convert_to_dataframe(NRC4c_align)
NRC4c_full_ent <- calculate_entropy_df(NRC4c_full_pos)

# annotate the file based on motifs and stretches
NRC4c_full_ent <- NRC4c_full_ent %>%
  mutate(color_segment = case_when(
    Position >= 234 & Position <= 237 |
      Position >= 256 & Position <= 262 |
      Position >= 288 & Position <= 292 |
      Position >= 522 & Position <= 529 |
      Position >= 550 & Position <= 553 | 
      Position >= 578 & Position <= 582 ~ "stretch 1",
    Position >= 544 & Position <= 549 |
      Position >= 566 & Position <= 573 ~ "stretch 2",
    Position >= 1 & Position <= 32 | 
      Position >= 199 & Position <= 209 | 
      Position >= 495 & Position <= 497 ~ "motifs",
    TRUE ~ "normal"
  ))

NRC4c_full_ent_p <- ggplot(NRC4c_full_ent, aes(x = Position, y = Entropy)) +
  geom_segment(aes(x = Position, xend = Position, y = 0, yend = Entropy, color = color_segment)) +
  scale_color_manual(values = c("stretch 1" = "#45B64B", "stretch 2" = "#EAC435", "motifs" = "#AFAFAF", "normal" = "black")) +
  ylim(c(0,4.32)) +
  xlim(c(0,923)) +
  geom_hline(yintercept = 1.5, linetype = "dotted", color = "#6C6C6C", size = 0.5, alpha = 1) +
  theme_custom2

ggsave(filename = "NRC4c_full.png", 
       NRC4c_full_ent_p, width = 10, height = 1.5, dpi = "retina")





# NRC4b
NRC4b_align <- NRCH_clu_filtered_ref_trimmed[NRCH_clu_filtered_ref_trimmed@ranges@NAMES %in% NRC4b_tree$tip.label]

NRC4b_full_pos <- convert_to_dataframe(NRC4b_align)
NRC4b_full_ent <- calculate_entropy_df(NRC4b_full_pos)

# annotate the file based on motifs and stretches
NRC4b_full_ent <- NRC4b_full_ent %>%
  mutate(color_segment = case_when(
    Position >= 234 & Position <= 237 |
      Position >= 256 & Position <= 262 |
      Position >= 288 & Position <= 292 |
      Position >= 522 & Position <= 529 |
      Position >= 550 & Position <= 553 | 
      Position >= 578 & Position <= 582 ~ "stretch 1",
    Position >= 544 & Position <= 549 |
      Position >= 566 & Position <= 573 ~ "stretch 2",
    Position >= 1 & Position <= 32 | 
      Position >= 199 & Position <= 209 | 
      Position >= 495 & Position <= 497 ~ "motifs",
    TRUE ~ "normal"
  ))

NRC4b_full_ent_p <- ggplot(NRC4b_full_ent, aes(x = Position, y = Entropy)) +
  geom_segment(aes(x = Position, xend = Position, y = 0, yend = Entropy, color = color_segment)) +
  scale_color_manual(values = c("stretch 1" = "#45B64B", "stretch 2" = "#EAC435", "motifs" = "#AFAFAF", "normal" = "black")) +
  ylim(c(0,4.32)) +
  xlim(c(0,923)) +
  geom_hline(yintercept = 1.5, linetype = "dotted", color = "#6C6C6C", size = 0.5, alpha = 1) +
  theme_custom2

ggsave(filename = "NRC4b_full.png", 
       NRC4b_full_ent_p, width = 10, height = 1.5, dpi = "retina")





# NRC4other
NRC4other_align <- NRCH_clu_filtered_ref_trimmed[NRCH_clu_filtered_ref_trimmed@ranges@NAMES %in% NRC4other_tree$tip.label]

NRC4other_full_pos <- convert_to_dataframe(NRC4other_align)
NRC4other_full_ent <- calculate_entropy_df(NRC4other_full_pos)

# annotate the file based on motifs and stretches
NRC4other_full_ent <- NRC4other_full_ent %>%
  mutate(color_segment = case_when(
    Position >= 234 & Position <= 237 |
      Position >= 256 & Position <= 262 |
      Position >= 288 & Position <= 292 |
      Position >= 522 & Position <= 529 |
      Position >= 550 & Position <= 553 | 
      Position >= 578 & Position <= 582 ~ "stretch 1",
    Position >= 544 & Position <= 549 |
      Position >= 566 & Position <= 573 ~ "stretch 2",
    Position >= 1 & Position <= 32 | 
      Position >= 199 & Position <= 209 | 
      Position >= 495 & Position <= 497 ~ "motifs",
    TRUE ~ "normal"
  ))

NRC4other_full_ent_p <- ggplot(NRC4other_full_ent, aes(x = Position, y = Entropy)) +
  geom_segment(aes(x = Position, xend = Position, y = 0, yend = Entropy, color = color_segment)) +
  scale_color_manual(values = c("stretch 1" = "#45B64B", "stretch 2" = "#EAC435", "motifs" = "#AFAFAF", "normal" = "black")) +
  ylim(c(0,4.32)) +
  xlim(c(0,923)) +
  geom_hline(yintercept = 1.5, linetype = "dotted", color = "#6C6C6C", size = 0.5, alpha = 1) +
  theme_custom2

ggsave(filename = "NRC4other_full.png", 
       NRC4other_full_ent_p, width = 10, height = 1.5, dpi = "retina")





# NRC4t
NRC4t_align <- NRCH_clu_filtered_ref_trimmed[NRCH_clu_filtered_ref_trimmed@ranges@NAMES %in% NRC4t_tree$tip.label]

NRC4t_full_pos <- convert_to_dataframe(NRC4t_align)
NRC4t_full_ent <- calculate_entropy_df(NRC4t_full_pos)

# annotate the file based on motifs and stretches
NRC4t_full_ent <- NRC4t_full_ent %>%
  mutate(color_segment = case_when(
    Position >= 234 & Position <= 237 |
      Position >= 256 & Position <= 262 |
      Position >= 288 & Position <= 292 |
      Position >= 522 & Position <= 529 |
      Position >= 550 & Position <= 553 | 
      Position >= 578 & Position <= 582 ~ "stretch 1",
    Position >= 544 & Position <= 549 |
      Position >= 566 & Position <= 573 ~ "stretch 2",
    Position >= 1 & Position <= 32 | 
      Position >= 199 & Position <= 209 | 
      Position >= 495 & Position <= 497 ~ "motifs",
    TRUE ~ "normal"
  ))

NRC4t_full_ent_p <- ggplot(NRC4t_full_ent, aes(x = Position, y = Entropy)) +
  geom_segment(aes(x = Position, xend = Position, y = 0, yend = Entropy, color = color_segment)) +
  scale_color_manual(values = c("stretch 1" = "#45B64B", "stretch 2" = "#EAC435", "motifs" = "#AFAFAF", "normal" = "black")) +
  ylim(c(0,4.32)) +
  xlim(c(0,923)) +
  geom_hline(yintercept = 1.5, linetype = "dotted", color = "#6C6C6C", size = 0.5, alpha = 1) +
  theme_custom2

ggsave(filename = "NRC4t_full.png", 
       NRC4t_full_ent_p, width = 10, height = 1.5, dpi = "retina")





# NRC4a
NRC4a_align <- NRCH_clu_filtered_ref_trimmed[NRCH_clu_filtered_ref_trimmed@ranges@NAMES %in% NRC4a_tree$tip.label]

NRC4a_full_pos <- convert_to_dataframe(NRC4a_align)
NRC4a_full_ent <- calculate_entropy_df(NRC4a_full_pos)

# annotate the file based on motifs and stretches
NRC4a_full_ent <- NRC4a_full_ent %>%
  mutate(color_segment = case_when(
    Position >= 234 & Position <= 237 |
      Position >= 256 & Position <= 262 |
      Position >= 288 & Position <= 292 |
      Position >= 522 & Position <= 529 |
      Position >= 550 & Position <= 553 | 
      Position >= 578 & Position <= 582 ~ "stretch 1",
    Position >= 544 & Position <= 549 |
      Position >= 566 & Position <= 573 ~ "stretch 2",
    Position >= 1 & Position <= 32 | 
      Position >= 199 & Position <= 209 | 
      Position >= 495 & Position <= 497 ~ "motifs",
    TRUE ~ "normal"
  ))

NRC4a_full_ent_p <- ggplot(NRC4a_full_ent, aes(x = Position, y = Entropy)) +
  geom_segment(aes(x = Position, xend = Position, y = 0, yend = Entropy, color = color_segment)) +
  scale_color_manual(values = c("stretch 1" = "#45B64B", "stretch 2" = "#EAC435", "motifs" = "#AFAFAF", "normal" = "black")) +
  ylim(c(0,4.32)) +
  xlim(c(0,923)) +
  geom_hline(yintercept = 1.5, linetype = "dotted", color = "#6C6C6C", size = 0.5, alpha = 1) +
  theme_custom2

ggsave(filename = "NRC4a_full.png", 
       NRC4a_full_ent_p, width = 10, height = 1.5, dpi = "retina")





# make a list of the full length

full_length_ent <- list(NRCH = NRCH_full_ent,
                        NRC0 = NRC0_full_ent,
                        NRCX = NRCX_full_ent,
                        NRC3 = NRC3_full_ent,
                        NRC2 = NRC2_full_ent,
                        NRC1 = NRC1_full_ent,
                        NRC9 = NRC9_full_ent,
                        NRC8 = NRC8_full_ent,
                        NRC6 = NRC6_full_ent,
                        NRC7 = NRC7_full_ent,
                        NRC5 = NRC5_full_ent,
                        NRC4c = NRC4c_full_ent,
                        NRC4b = NRC4b_full_ent,
                        NRC4other = NRC4other_full_ent,
                        NRC4t = NRC4t_full_ent,
                        NRC4a = NRC4a_full_ent)

full_length_df <- merge_dataframes(full_length_ent)
write_csv(full_length_df,"/path/to/full_length_entropy.csv", col_names = TRUE) # later this was integrated as part of DataS4




