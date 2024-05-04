library(tidyverse)
library(entropy)


myColors <- c("#0d47a1","#ffffff","#D10000")


# Function to convert a Biostrings/dataframe object to a dataframe and break the sequence into characters
convert_to_dataframe <- function(alignment) {
  # Extract sequences as strings
  sequences <- sapply(alignment, as.character)
  # Split each sequence into individual characters and convert to a matrix
  alignment_matrix <- do.call(rbind, strsplit(sequences, ""))
  # Convert the matrix to a dataframe
  alignment_df <- as.data.frame(alignment_matrix, stringsAsFactors = FALSE)
  # Naming columns as positions
  colnames(alignment_df) <- paste(seq_len(ncol(alignment_df)))
  
  return(alignment_df)
}

# Wrapper function to process a list of Biostrings/dataframe objects
convert_list_to_dataframes <- function(alignment_list) {
  # Apply the convert_to_dataframe function to each element in the list
  lapply(alignment_list, convert_to_dataframe)
}



# Function to calculate entropy for multiple columns and return as a dataframe
calculate_entropy_df <- function(dataframe) {
  # Function to calculate entropy for a single column
  calculate_column_entropy <- function(column) {
    # Remove any NA values (if there are any)
    column <- na.omit(column)
    # Calculate frequencies
    freq_table <- table(column)
    # Convert to proportions
    proportions <- freq_table / sum(freq_table)
    # Calculate and return entropy
    entropy.empirical(proportions, unit = "log2")
  }
  
  # Ensure that the dataframe is processed column-wise
  entropy_values <- sapply(dataframe, calculate_column_entropy, simplify = FALSE)
  
  # Flatten the list if necessary
  if (is.list(entropy_values)) {
    entropy_values <- unlist(entropy_values)
  }
  
  # Create a dataframe with positions and entropy values
  entropy_df <- data.frame(Position = seq_along(entropy_values), Entropy = entropy_values)
  
  return(entropy_df)
}




# Wrapper function to process a list of dataframes
calculate_entropy_for_list <- function(list_of_dataframes) {
  # Apply the calculate_entropy_df function to each dataframe in the list
  lapply(list_of_dataframes, calculate_entropy_df)
}


# the outputs are going to be exported here
setwd("/path/to/your/figure/directory")



## calculate entropy and export entropy plot
# set the breakpoints for plot coloring
breakpoints <- c(0, 1.5, 4.32)

# calculate entropy for S1a
S1a_pos <- convert_list_to_dataframes(S1a_motifs)
S1a_entropy <- calculate_entropy_for_list(S1a_pos)


for (name in names(S1a_entropy)) {
  df <- S1a_entropy[[name]]
  
  # create the plot
  plot <- ggplot(df, aes(x = Position, y = 1, fill = Entropy)) +
    geom_tile(color = NA) +
    scale_fill_gradientn(colors = myColors, 
                         values = scales::rescale(breakpoints, to = c(0, 1)), 
                         limits = c(0, 4.32)) +
    theme_minimal() +
    theme_custom
  
  ggsave(filename = paste0(name, ".png"), plot = plot, width = 4, height = 0.5, dpi = "retina")
}




# calculate entropy for S1b
S1b_pos <- convert_list_to_dataframes(S1b_motifs)
S1b_entropy <- calculate_entropy_for_list(S1b_pos)


for (name in names(S1b_entropy)) {
  df <- S1b_entropy[[name]]
  
  # create the plot
  plot <- ggplot(df, aes(x = Position, y = 1, fill = Entropy)) +
    geom_tile(color = NA) +
    scale_fill_gradientn(colors = myColors, 
                         values = scales::rescale(breakpoints, to = c(0, 1)), 
                         limits = c(0, 4.32)) +
    theme_minimal() +
    theme_custom
  
  ggsave(filename = paste0(name, ".png"), plot = plot, width = 7, height = 0.5, dpi = "retina")
}




# calculate entropy for S1c
S1c_pos <- convert_list_to_dataframes(S1c_motifs)
S1c_entropy <- calculate_entropy_for_list(S1c_pos)


for (name in names(S1c_entropy)) {
  df <- S1c_entropy[[name]]
  
  # create the plot
  plot <- ggplot(df, aes(x = Position, y = 1, fill = Entropy)) +
    geom_tile(color = NA) +
    scale_fill_gradientn(colors = myColors, 
                         values = scales::rescale(breakpoints, to = c(0, 1)), 
                         limits = c(0, 4.32)) +
    theme_minimal() +
    theme_custom
  
  ggsave(filename = paste0(name, ".png"), plot = plot, width = 5, height = 0.5, dpi = "retina")
}




# calculate entropy for S1d
S1d_pos <- convert_list_to_dataframes(S1d_motifs)
S1d_entropy <- calculate_entropy_for_list(S1d_pos)


for (name in names(S1d_entropy)) {
  df <- S1d_entropy[[name]]
  
  # create the plot
  plot <- ggplot(df, aes(x = Position, y = 1, fill = Entropy)) +
    geom_tile(color = NA) +
    scale_fill_gradientn(colors = myColors, 
                         values = scales::rescale(breakpoints, to = c(0, 1)), 
                         limits = c(0, 4.32)) +
    theme_minimal() +
    theme_custom
  
  ggsave(filename = paste0(name, ".png"), plot = plot, width = 8, height = 0.5, dpi = "retina")
}




# calculate entropy for S1e
S1e_pos <- convert_list_to_dataframes(S1e_motifs)
S1e_entropy <- calculate_entropy_for_list(S1e_pos)


for (name in names(S1e_entropy)) {
  df <- S1e_entropy[[name]]
  
  # create the plot
  plot <- ggplot(df, aes(x = Position, y = 1, fill = Entropy)) +
    geom_tile(color = NA) +
    scale_fill_gradientn(colors = myColors, 
                         values = scales::rescale(breakpoints, to = c(0, 1)), 
                         limits = c(0, 4.32)) +
    theme_minimal() +
    theme_custom
  
  ggsave(filename = paste0(name, ".png"), plot = plot, width = 4, height = 0.5, dpi = "retina")
}



# calculate entropy for S1f
S1f_pos <- convert_list_to_dataframes(S1f_motifs)
S1f_entropy <- calculate_entropy_for_list(S1f_pos)


for (name in names(S1f_entropy)) {
  df <- S1f_entropy[[name]]
  
  # create the plot
  plot <- ggplot(df, aes(x = Position, y = 1, fill = Entropy)) +
    geom_tile(color = NA) +
    scale_fill_gradientn(colors = myColors, 
                         values = scales::rescale(breakpoints, to = c(0, 1)), 
                         limits = c(0, 4.32)) +
    theme_minimal() +
    theme_custom
  
  ggsave(filename = paste0(name, ".png"), plot = plot, width = 5, height = 0.5, dpi = "retina")
}



# calculate entropy for S2a
S2a_pos <- convert_list_to_dataframes(S2a_motifs)
S2a_entropy <- calculate_entropy_for_list(S2a_pos)


for (name in names(S2a_entropy)) {
  df <- S2a_entropy[[name]]
  
  # create the plot
  plot <- ggplot(df, aes(x = Position, y = 1, fill = Entropy)) +
    geom_tile(color = NA) +
    scale_fill_gradientn(colors = myColors, 
                         values = scales::rescale(breakpoints, to = c(0, 1)), 
                         limits = c(0, 4.32)) +
    theme_minimal() +
    theme_custom
  
  ggsave(filename = paste0(name, ".png"), plot = plot, width = 6, height = 0.5, dpi = "retina")
}



# calculate entropy for S2b
S2b_pos <- convert_list_to_dataframes(S2b_motifs)
S2b_entropy <- calculate_entropy_for_list(S2b_pos)


for (name in names(S2b_entropy)) {
  df <- S2b_entropy[[name]]
  
  # create the plot
  plot <- ggplot(df, aes(x = Position, y = 1, fill = Entropy)) +
    geom_tile(color = NA) +
    scale_fill_gradientn(colors = myColors, 
                         values = scales::rescale(breakpoints, to = c(0, 1)), 
                         limits = c(0, 4.32)) +
    theme_minimal() +
    theme_custom
  
  ggsave(filename = paste0(name, ".png"), plot = plot, width = 8, height = 0.5, dpi = "retina")
}







# calculate entropy for MADA
MADA_pos <- convert_list_to_dataframes(MADA_motifs)
MADA_entropy <- calculate_entropy_for_list(MADA_pos)


for (name in names(MADA_entropy)) {
  df <- MADA_entropy[[name]]
  
  # create the plot
  plot <- ggplot(df, aes(x = Position, y = 1, fill = Entropy)) +
    geom_tile(color = NA) +
    scale_fill_gradientn(colors = myColors, limits = c(0,4.32)) +
    theme_minimal() +
    theme_custom
  
  ggsave(filename = paste0(name, ".png"), plot = plot, width = 21, height = 0.5, dpi = "retina")
}




# calculate entropy for MHD
MHD_pos <- convert_list_to_dataframes(MHD_motifs)
MHD_entropy <- calculate_entropy_for_list(MHD_pos)


for (name in names(MHD_entropy)) {
  df <- MHD_entropy[[name]]
  
  # create the plot
  plot <- ggplot(df, aes(x = Position, y = 1, fill = Entropy)) +
    geom_tile(color = NA) +
    scale_fill_gradientn(colors = myColors, limits = c(0,4.32)) +
    theme_minimal() +
    theme_custom
  
  ggsave(filename = paste0(name, ".png"), plot = plot, width = 3, height = 0.5, dpi = "retina")
}




# calculate entropy for Ploop
Ploop_pos <- convert_list_to_dataframes(Ploop_motifs)
Ploop_entropy <- calculate_entropy_for_list(Ploop_pos)


for (name in names(Ploop_entropy)) {
  df <- Ploop_entropy[[name]]
  
  # create the plot
  plot <- ggplot(df, aes(x = Position, y = 1, fill = Entropy)) +
    geom_tile(color = NA) +
    scale_fill_gradientn(colors = myColors, limits = c(0,4.32)) +
    theme_minimal() +
    theme_custom
  
  ggsave(filename = paste0(name, ".png"), plot = plot, width = 11, height = 0.5, dpi = "retina")
}








# export the quantitative data as tables

# function to merge the dataframes in a list
merge_dataframes <- function(df_list) {
  # Initialize an empty dataframe to store the merged result
  merged_df <- data.frame(Position = df_list[[1]]$Position)
  
  # Loop through each dataframe and merge
  for (name in names(df_list)) {
    # Create a temporary dataframe with renamed 'Entropy' column
    temp_df <- data.frame(Position = df_list[[name]]$Position)
    temp_df[[paste0(name)]] <- df_list[[name]]$Entropy
    
    # Merge with the main dataframe
    merged_df <- merge(merged_df, temp_df, by = "Position", all = TRUE)
  }
  
  return(merged_df)
}




# export the entropy scores of the motifs and stretches
S1a_df <- merge_dataframes(S1a_entropy)
write_csv(S1a_df,"/path/to/S1a_entropy.csv", col_names = TRUE)

S1b_df <- merge_dataframes(S1b_entropy)
write_csv(S1b_df,"/path/to/S1b_entropy.csv", col_names = TRUE)

S1c_df <- merge_dataframes(S1c_entropy)
write_csv(S1c_df,"~/path/to/S1c_entropy.csv", col_names = TRUE)

S1d_df <- merge_dataframes(S1d_entropy)
write_csv(S1d_df,"~/path/to/S1d_entropy.csv", col_names = TRUE)

S1e_df <- merge_dataframes(S1e_entropy)
write_csv(S1e_df,"/path/to/S1e_entropy.csv", col_names = TRUE)

S1f_df <- merge_dataframes(S1f_entropy)
write_csv(S1f_df,"/path/to/S1f_entropy.csv", col_names = TRUE)

S2a_df <- merge_dataframes(S2a_entropy)
write_csv(S2a_df,"/path/to/S2a_entropy.csv", col_names = TRUE)

S2b_df <- merge_dataframes(S2b_entropy)
write_csv(S2b_df,"/path/to/S2b_entropy.csv", col_names = TRUE)





MADA_df <- merge_dataframes(MADA_entropy) %>% as.data.frame()
write_csv(MADA_df,"/path/to/MADA_entropy.csv", col_names = TRUE)

MHD_df <- merge_dataframes(MHD_entropy) %>% as.data.frame()
write_csv(MHD_df,"/path/to/MHD_entropy.csv", col_names = TRUE)

Ploop_df <- merge_dataframes(Ploop_entropy) %>% as.data.frame()
write_csv(Ploop_df,"/path/to/Ploop_entropy.csv", col_names = TRUE)


