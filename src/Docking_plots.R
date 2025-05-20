## Analyzing IPSAE Docking Data on AF3 predicted structures
## JRR 3/14/2025
# This is mostly just compiling the data from the output txt files
#############################################################################
## Load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
## Set working directory
setwd('~/Desktop/AlphaFold3')
###############################################################################
## Functions
load_data <- function(filename, model){
  data <- read.table(paste0(filename, '/', filename, '_model_', model, '_15_15.txt'), 
                     header = T)
  return(data)
}

pull_scores <- function(scores_df, input_data, model, chain_key, pair_list){
  for (pair in pair_list){
    protein_a <- pair[1]
    a_name <- chain_key$Protein[which(chain_key$Chain == protein_a)]
    protein_b <- pair[2]
    b_name <- chain_key$Protein[which(chain_key$Chain == protein_b)]
    reduced_data <- subset(input_data, Chn1 == protein_a & Chn2 == protein_b 
                           & Type == 'max')
    dq <- reduced_data[[11]] 
    dq2 <- reduced_data[[12]]
    lis <- reduced_data[[13]]
    new_data <- data.frame(
      Chn1 = a_name,
      Chn2 = b_name,
      Model = model,
      DockQ = dq,
      DockQ2 = dq2,
      LIS = lis
    )
    scores_df <- rbind(scores_df, new_data)
  }
  return(scores_df)
}

generate_plots <- function(scores_df, path){
  colnames(scores_df) <- c("Chn1", 'Chn2', 'Model', 'DockQ', 'DockQ2', 'LIS')
  scores_df <- scores_df %>%
    mutate(
      norm_Chn1 = pmin(Chn1, Chn2),
      norm_Chn2 = pmax(Chn1, Chn2)
    )
  scores_means <- scores_df %>%
    group_by(norm_Chn1, norm_Chn2) %>%
    summarize(
      mean_DockQ = mean(DockQ), 
      sd_DockQ = sd(DockQ),
      mean_DockQ2 = mean(DockQ2),
      sd_DockQ2 = sd(DockQ2),
      mean_LIS = mean(LIS), 
      sd_LIS = sd(LIS), 
      .groups = "drop"
    ) %>%
    mutate(norm_Chn1 = norm_Chn1, norm_Chn2 = norm_Chn2)
  
  plot_data<- scores_means %>%
    pivot_longer(
      cols = starts_with("mean_"),
      names_to = "score_type",
      values_to = "mean_value"
    ) %>%
    mutate(
      sd_value = case_when(
        score_type == "mean_DockQ" ~ sd_DockQ, 
        score_type == "mean_DockQ2" ~ sd_DockQ2,
        score_type == "mean_LIS" ~ sd_LIS
      ), 
      score_color = case_when(
        score_type == "mean_DockQ" ~ "blue3",
        score_type == "mean_DockQ2" ~ "tomato",
        score_type == "mean_LIS" ~ "gold"
      )
    )
  point_data <- scores_df %>%
    pivot_longer(
      cols = c(DockQ, DockQ2, LIS), 
      names_to = "score_type", 
      values_to = "score_value"
    )
  
  p <- ggplot(plot_data, aes(x = paste(norm_Chn1, norm_Chn2, sep = "-"),
                             y = mean_value, fill = score_type))+
    geom_col(position = position_dodge(), width = 0.7, alpha = 0.9, color = 'black') +
    geom_point(data = point_data, aes(x = paste(norm_Chn1, norm_Chn2, sep = "-"), y = score_value), 
               color = "black", position = position_dodge(width = 0.7), alpha = 0.7, size = .5)+
    geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value),
                  position = position_dodge(width = 0.7), width = 0.2) +
    labs(x = 'Chain Pairings', y = "Scores", fill = "Score Type")+
    scale_y_continuous(expand = c(0, 0), limits = c(0, 0.7))+
    geom_hline(yintercept = 0.23, linetype = 'dotted')+
    scale_fill_manual(values = c("mean_DockQ" = 'blue3', "mean_DockQ2" = "tomato", 
                                 "mean_LIS" = "gold"))+
    ggtitle(path)+
    theme_bw()+
    theme(axis.line = element_line(color='black'),
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.ticks.length = unit(10, 'pt'),
          axis.minor.ticks.length = rel(0.5),
          axis.title.y = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          plot.title = element_text(size = 14),
          text = element_text(size = 10, family = 'Helvetica'),
          legend.position = 'right',
          axis.text.x = element_text(angle = 45, hjust = 1))
  print(p)
  ggsave(filename = paste0(path, 'pairings_plot.pdf'), plot = p,
         device = 'pdf', height = 5, width = 7)
}

run_analysis <- function(path, n_model, chain_key){
  chains <- c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K')
  pairs_list1 <- list(c('A', 'B'), c('C', 'D'), c('C', 'H'), c('D', 'E'),
                     c('E', 'I'), c('I', 'F'), c('I', 'J'), c('F', 'G'), c('F', 'J'), 
                     c('G', 'H'), c('G', 'J'), c('H', 'K'), c('K', 'J'), c('A', 'I'),
                     c('A', 'F'), c('B', 'G'), c('B', 'H'), c('B', 'F'), c('C', 'K'),
                     c('E', 'J'), c('F', 'K'), c('G', 'K')
  )
  pairs_list2 <- list(c('K', 'A'), c('B', 'C'), c('B', 'G'), c('C', 'D'),
                      c('D', 'H'), c('H', 'E'), c('H', 'I'), c('E', 'F'), c('E', 'I'), 
                      c('F', 'G'), c('F', 'I'), c('G', 'J'), c('J', 'I'), c('K', 'H'),
                      c('K', 'E'), c('A', 'F'), c('A', 'G'), c('A', 'E'), c('B', 'J'),
                      c('D', 'I'), c('E', 'J'), c('F', 'J')
  )
  if (chain_key == 'chain_key_1'){
    chain_key_df <- read.table('Chain_key_1.txt', header = T)
    normalized_pairs <-lapply(pairs_list1, function(x) sort(x))
    unique_pairs <- unique(sapply(normalized_pairs, paste, collapse = '-'))
    fpl <- strsplit(unique_pairs, '-')
  }
  else if (chain_key == 'chain_key_2'){
    chain_key_df <- read.table('Chain_key_2.txt', header = T)
    normalized_pairs <-lapply(pairs_list2, function(x) sort(x))
    unique_pairs <- unique(sapply(normalized_pairs, paste, collapse = '-'))
    fpl <- strsplit(unique_pairs, '-')
  }
  
  model_list <- seq(0, n_model - 1)
  scores_df <- data.frame(
    Chn1 = character(),
    Chn2 = character(),
    Model = integer(), 
    DockQ = numeric(), 
    DockQ2 = numeric(), 
    LIS = numeric(), 
    stringsAsFactors = F
  )
  for (i in model_list){
    df <- load_data(filename = path, i)
    scores_df <- pull_scores(scores_df, df, i, chain_key_df, fpl)
  }
  generate_plots(scores_df, path)
}
get_data <- function(path, n_model, chain_key){
  chains <- c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K')
  pairs_list1 <- list(c('A', 'B'), c('C', 'D'), c('C', 'H'), c('D', 'E'),
                      c('E', 'I'), c('I', 'F'), c('I', 'J'), c('F', 'G'), c('F', 'J'), 
                      c('G', 'H'), c('G', 'J'), c('H', 'K'), c('K', 'J'), c('A', 'I'),
                      c('A', 'F'), c('B', 'G'), c('B', 'H'), c('B', 'F'), c('C', 'K'),
                      c('E', 'J'), c('F', 'K'), c('G', 'K')
  )
  pairs_list2 <- list(c('K', 'A'), c('B', 'C'), c('B', 'G'), c('C', 'D'),
                      c('D', 'H'), c('H', 'E'), c('H', 'I'), c('E', 'F'), c('E', 'I'), 
                      c('F', 'G'), c('F', 'I'), c('G', 'J'), c('J', 'I'), c('K', 'H'),
                      c('K', 'E'), c('A', 'F'), c('A', 'G'), c('A', 'E'), c('B', 'J'),
                      c('D', 'I'), c('E', 'J'), c('F', 'J')
  )
  if (chain_key == 'chain_key_1'){
    chain_key_df <- read.table('Chain_key_1.txt', header = T)
    normalized_pairs <-lapply(pairs_list1, function(x) sort(x))
    unique_pairs <- unique(sapply(normalized_pairs, paste, collapse = '-'))
    fpl <- strsplit(unique_pairs, '-')
  }
  else if (chain_key == 'chain_key_2'){
    chain_key_df <- read.table('Chain_key_2.txt', header = T)
    normalized_pairs <-lapply(pairs_list2, function(x) sort(x))
    unique_pairs <- unique(sapply(normalized_pairs, paste, collapse = '-'))
    fpl <- strsplit(unique_pairs, '-')
  }
  
  model_list <- seq(0, n_model - 1)
  scores_df <- data.frame(
    Chn1 = character(),
    Chn2 = character(),
    Model = integer(), 
    DockQ = numeric(), 
    DockQ2 = numeric(), 
    LIS = numeric(), 
    stringsAsFactors = F
  )
  for (i in model_list){
    df <- load_data(filename = path, i)
    scores_df <- pull_scores(scores_df, df, i, chain_key_df, fpl)
  }
  return(scores_df)
}

#############################################################################
## Running functions that extracts data into a tidy format
## Input dataframe scores are the results tables from running the ipsae.py script
run_analysis('fold_proteasome_alpha_1_4', 5, 'chain_key_1')
run_analysis('fold_sacch_proteasome', 5, 'chain_key_1')
run_analysis('fold_proteasome_alpha_scpoc4', 5, 'chain_key_2')
run_analysis('fold_proteasome_alpha_cd421', 5, 'chain_key_2')
run_analysis('fold_proteasome_alpha_cd422', 5, 'chain_key_2')
run_analysis('fold_proteasome_alpha_cd423', 5, 'chain_key_2')
run_analysis('fold_proteasome_alpha_cd424', 5, 'chain_key_2')
run_analysis('fold_proteasome_alpha_cd425', 5, 'chain_key_2')

############################################################################
## Pulls data from tidy format for chains that would be expected to interact based on the 3D structure
cau_wt <- get_data('fold_proteasome_alpha_1_4', 5, 'chain_key_1')
sc_wt <- get_data('fold_sacch_proteasome', 5, 'chain_key_1')
cau_scpoc4 <- get_data('fold_proteasome_alpha_scpoc4', 5, 'chain_key_2')
cau_cd421 <- get_data('fold_proteasome_alpha_cd421', 5, 'chain_key_2')
cau_cd422 <- get_data('fold_proteasome_alpha_cd422', 5, 'chain_key_2')
cau_cd423 <- get_data('fold_proteasome_alpha_cd423', 5, 'chain_key_2')
cau_cd424 <- get_data('fold_proteasome_alpha_cd424', 5, 'chain_key_2')
cau_cd425 <- get_data('fold_proteasome_alpha_cd425', 5, 'chain_key_2')

# filter for IRC25-POC4 interactions
protein_a <- 'Poc4'
protein_b <- 'Irc25'
protein_c <- 'chap_des'
cau_wt_filt <- subset(cau_wt, Chn1 == protein_a & Chn2 == protein_b)
cau_wt_filt$structure <- 'Cau_wt'

sc_wt_filt <- subset(sc_wt, Chn1 == protein_a & Chn2 == protein_b)
sc_wt_filt$structure <- 'Sc_wt'


cau_scpoc4_filt <- subset(cau_scpoc4, Chn2 == protein_c & Chn1 == protein_b)
cau_scpoc4_filt$structure <- 'Cau_ScPoc4'

cau_cd421_filt <- subset(cau_cd421, Chn2 == protein_c & Chn1 == protein_b)
cau_cd421_filt$structure <- 'Cau_cd421'

cau_cd422_filt <- subset(cau_cd422, Chn2 == protein_c & Chn1 == protein_b)
cau_cd422_filt$structure <- 'Cau_cd422'

cau_cd423_filt <- subset(cau_cd423, Chn2 == protein_c & Chn1 == protein_b)
cau_cd423_filt$structure <- 'Cau_cd423'

cau_cd424_filt <- subset(cau_cd424, Chn2 == protein_c & Chn1 == protein_b)
cau_cd424_filt$structure <- 'Cau_cd424'

cau_cd425_filt <- subset(cau_cd425, Chn2 == protein_c & Chn1 == protein_b)
cau_cd425_filt$structure <- 'Cau_cd425'

poc4_irc25_df <- rbind(cau_wt_filt, sc_wt_filt, cau_scpoc4_filt, 
                       cau_cd421_filt, cau_cd422_filt, cau_cd423_filt, 
                       cau_cd424_filt, cau_cd425_filt)
poc4_irc25_clean <- poc4_irc25_df %>%
  select(-Chn1, -Chn2)

poc4_irc25_means <- poc4_irc25_clean %>%
  group_by(structure) %>%
  summarize(
    mean_DockQ = mean(DockQ), 
    sd_DockQ = sd(DockQ),
    mean_DockQ2 = mean(DockQ2),
    sd_DockQ2 = sd(DockQ2),
    mean_LIS = mean(LIS), 
    sd_LIS = sd(LIS), 
    .groups = "drop"
  )
poc4_irc25_means$structure <- fct_relevel(poc4_irc25_means$structure, 
                                          "Cau_wt", "Sc_wt", 'Cau_ScPoc4', 
                                          "Cau_cd421", "Cau_cd422", "Cau_cd423", 
                                          "Cau_cd424", 'Cau_cd425')
poc4_irc25_means$activity <- c("y", "y", "n", "n", "n", "n", "y", "y")
poc4_irc25_clean$structure <- fct_relevel(poc4_irc25_clean$structure, 
                                          "Cau_wt", "Sc_wt", 'Cau_ScPoc4', 
                                          "Cau_cd421", "Cau_cd422", "Cau_cd423", 
                                          "Cau_cd424", 'Cau_cd425')

poc4_irc25_clean <- merge(poc4_irc25_clean, poc4_irc25_means[, c("structure", "activity")], by = "structure")

ggplot(poc4_irc25_means, aes(x = structure, y = mean_DockQ2, fill = activity))+
  geom_col(width = 0.5, color = 'black')+
  geom_point(data = poc4_irc25_clean, aes(x = structure, y = DockQ2), color = 'black',
             alpha = 1, position = position_jitter(width = 0.15), size = 2)+
  geom_errorbar(aes(ymin = mean_DockQ2 - sd_DockQ2, ymax = mean_DockQ2 + sd_DockQ2), 
                width = 0.2)+
  scale_fill_manual(values = c("y" = 'seagreen2', "n" = 'gray48'))+
  ylab('DockQ2 Score')+
  xlab('Model')+
  scale_y_continuous(expand = c(0,0), limits = c(0, 0.6),
                     minor_breaks = seq(0, 0.6, 0.02))+
  geom_hline(yintercept = 0.23, linetype = 'dotted')+
  ggtitle('Poc4-Irc25 DockQ2 Scores')+
  guides(y = guide_axis(minor.ticks = T))+
  theme_bw()+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks.length = unit(10, 'pt'),
        axis.minor.ticks.length = rel(0.5),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        plot.title = element_text(size = 14),
        text = element_text(size = 10, family = 'Helvetica'),
        legend.position = 'right',
  )
## Poc4-Pup2 interaction
# filter for Pup2-POC4 interactions
protein_a <- 'Poc4'
protein_b <- 'Pup2'
protein_c <- 'chap_des'
cau_wt_filt <- subset(cau_wt, Chn1 == protein_a & Chn2 == protein_b)
cau_wt_filt$structure <- 'Cau_wt'

sc_wt_filt <- subset(sc_wt, Chn1 == protein_a & Chn2 == protein_b)
sc_wt_filt$structure <- 'Sc_wt'

cau_scpoc4_filt <- subset(cau_scpoc4, Chn2 == protein_c & Chn1 == protein_b)
cau_scpoc4_filt$structure <- 'Cau_ScPoc4'

cau_cd421_filt <- subset(cau_cd421, Chn2 == protein_c & Chn1 == protein_b)
cau_cd421_filt$structure <- 'Cau_cd421'

cau_cd422_filt <- subset(cau_cd422, Chn2 == protein_c & Chn1 == protein_b)
cau_cd422_filt$structure <- 'Cau_cd422'

cau_cd423_filt <- subset(cau_cd423, Chn2 == protein_c & Chn1 == protein_b)
cau_cd423_filt$structure <- 'Cau_cd423'

cau_cd424_filt <- subset(cau_cd424, Chn2 == protein_c & Chn1 == protein_b)
cau_cd424_filt$structure <- 'Cau_cd424'

cau_cd425_filt <- subset(cau_cd425, Chn2 == protein_c & Chn1 == protein_b)
cau_cd425_filt$structure <- 'Cau_cd425'

poc4_Pup2_df <- rbind(cau_wt_filt, sc_wt_filt, cau_scpoc4_filt, 
                       cau_cd421_filt, cau_cd422_filt, cau_cd423_filt, 
                       cau_cd424_filt, cau_cd425_filt)
poc4_Pup2_clean <- poc4_Pup2_df %>%
  select(-Chn1, -Chn2)

poc4_Pup2_means <- poc4_Pup2_clean %>%
  group_by(structure) %>%
  summarize(
    mean_DockQ = mean(DockQ), 
    sd_DockQ = sd(DockQ),
    mean_DockQ2 = mean(DockQ2),
    sd_DockQ2 = sd(DockQ2),
    mean_LIS = mean(LIS), 
    sd_LIS = sd(LIS), 
    .groups = "drop"
  )
poc4_Pup2_means$structure <- fct_relevel(poc4_Pup2_means$structure, 
                                          "Cau_wt", "Sc_wt", 'Cau_ScPoc4', 
                                          "Cau_cd421", "Cau_cd422", "Cau_cd423", 
                                          "Cau_cd424", 'Cau_cd425')
poc4_Pup2_means$activity <- c("y", "y", "n", "n", "n", "n", "y", "y")
poc4_Pup2_clean$structure <- fct_relevel(poc4_Pup2_clean$structure, 
                                          "Cau_wt", "Sc_wt", 'Cau_ScPoc4', 
                                          "Cau_cd421", "Cau_cd422", "Cau_cd423", 
                                          "Cau_cd424", 'Cau_cd425')

poc4_Pup2_clean <- merge(poc4_Pup2_clean, poc4_Pup2_means[, c("structure", "activity")], by = "structure")

ggplot(poc4_Pup2_means, aes(x = structure, y = mean_DockQ2, fill = activity))+
  geom_col(width = 0.5, color = 'black')+
  geom_point(data = poc4_Pup2_clean, aes(x = structure, y = DockQ2), color = 'black',
             alpha = 1, position = position_jitter(width = 0.15), size = 2)+
  geom_errorbar(aes(ymin = mean_DockQ2 - sd_DockQ2, ymax = mean_DockQ2 + sd_DockQ2), 
                width = 0.2)+
  scale_fill_manual(values = c("y" = 'seagreen2', "n" = 'gray48'))+
  ylab('DockQ2 Score')+
  xlab('Model')+
  scale_y_continuous(expand = c(0,0), limits = c(0, 0.6),
                     minor_breaks = seq(0, 0.6, 0.02))+
  geom_hline(yintercept = 0.23, linetype = 'dotted')+
  ggtitle('Poc4-Pup2 DockQ2 Scores')+
  guides(y = guide_axis(minor.ticks = T))+
  theme_bw()+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks.length = unit(10, 'pt'),
        axis.minor.ticks.length = rel(0.5),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        plot.title = element_text(size = 14),
        text = element_text(size = 10, family = 'Helvetica'),
        legend.position = 'right',
  )