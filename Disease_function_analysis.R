#Disease function analysis

library(tidyverse)
library(purrr)
library(pheatmap)

dis_myeloid_c0 <- read.csv('data/Disease_c0.csv')
dis_myeloid_c1 <- read.csv('data/Disease_c1.csv')
dis_myeloid_c2 <- read.csv('data/Disease_c2.csv')
dis_myeloid_c3 <- read.csv('data/Disease_c3.csv')
dis_myeloid_c4 <- read.csv('data/Disease_c4.csv')
dis_myeloid_c5 <- read.csv('data/Disease_c5.csv')
dis_myeloid_c6 <- read.csv('data/Disease_c6.csv')
dis_myeloid_c7 <- read.csv('data/Disease_c7.csv')
dis_myeloid_c8 <- read.csv('data/Disease_c8.csv')
dis_myeloid_c10 <- read.csv('data/Disease_c10.csv')


#Function to obtain siginifcant data without na for z-score
clean_data <- function(df) {
  df_final <- df |>
    filter(!is.na(Activation.z.score)) |>
    filter(p.value < 0.05) |>
    mutate_at('Activation.z.score', as.numeric) |>
    select(-Predicted.Activation.State)
  colnames(df_final)[1] <- 'Categories'
  return(df_final)
}

#Function to get avg activation of each category
activation_avg <- function(df, colname) {
  df_final <- clean_data(df) |>
    group_by(Categories) |>
    summarise(za = mean(Activation.z.score))
  colnames(df_final)[2] <- colname
  rownames(df_final) <- df_final$Categories
  return(df_final)
}


c0 <- activation_avg(dis_myeloid_c0, 'c0')
c1 <- activation_avg(dis_myeloid_c1, 'c1')
c2 <- activation_avg(dis_myeloid_c2, 'c2')
c3 <- activation_avg(dis_myeloid_c3, 'c3')
c4 <- activation_avg(dis_myeloid_c4, 'c4')
c5 <- activation_avg(dis_myeloid_c5, 'c5')
c6 <- activation_avg(dis_myeloid_c6, 'c6')
c7 <- activation_avg(dis_myeloid_c7, 'c7')
c8 <- activation_avg(dis_myeloid_c8, 'c8')
c10 <- activation_avg(dis_myeloid_c10, 'c10')


Exp_table <- c0 |>
  inner_join(c1, by = 'Categories') |>
  inner_join(c2, by = 'Categories') |>
  inner_join(c3, by = 'Categories') |>
  inner_join(c4, by = 'Categories') |>
  inner_join(c5, by = 'Categories') |>
  inner_join(c6, by = 'Categories') |>
  inner_join(c7, by = 'Categories') |>
  inner_join(c8, by = 'Categories') |>
  inner_join(c10, by = 'Categories')


#Change row name of Exp_table
Categories <- Exp_table$Categories
Exp_table <- Exp_table|>
  select(-Categories)
rownames(Exp_table) <- Categories



#Function to generate heatmap with sig label
Generate_heatmap <- function(df){
  heat <- t(scale(t(as.matrix(df))))
  pheatmap(heat,cluster_rows = F,cluster_cols = T,border_color = NA, cellwidth = 30, cellheight = 16)
}

Generate_heatmap(Exp_table)

