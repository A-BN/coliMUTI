# Loading stuff ----
library(tidyverse)

all_strain <- 
  sapply(X = list.files(path = './gd/', pattern = '*.gd', full.names = TRUE), 
       FUN = function(x) read.csv(x, sep = '\t', stringsAsFactors = FALSE))

all_strain <-
  do.call(bind_rows, all_strain)

# Unique ID ----
all_strain |>
  rowwise() |> 
  mutate(mut_id = paste(mutation_category, position_start,position_end,ref_seq, new_seq, sep = '|')) |> 
  ungroup() |> 
  group_by() |> 
  identity() -> all_strain

### Mutations that are too close (on the same read) are a bit worrisome
all_strain$n_neighbor <- NA
read_size <- 150

for (i in 1:nrow(all_strain)) {
  curr_pos <- all_strain$position[i]
  curr_id <- all_strain$mut_id[i]
  curr_strain <- all_strain$title[i]
  curr_range <- c(curr_pos - read_size,curr_pos + read_size)
  is_do <- all_strain$position >= curr_range[1]
  is_up <- all_strain$position <= curr_range[2]
  curr_n_mut <- length(unique(all_strain$mut_id[is_up & is_do])) - 1 
  all_strain$n_neighbor[i] <- (sum(curr_n_mut))
}

# Filtering out neighbors ----
all_strain |> 
  ungroup() |> 
  filter(n_neighbor < 1) |> 
  group_by(mut_id) |> 
  mutate(n_strains = length(unique(title))) |> 
  identity() -> all_strain_filt_1

# Presence Absence ----
## Populating a presence absence matrix of mutations by strain

pre_ab <- 
  matrix(data = 0, nrow = length(unique(all_strain_filt_1$mut_id)), 
         ncol = length(unique(all_strain_filt_1$title)))
colnames(pre_ab) <- unique(all_strain_filt_1$title)
rownames(pre_ab) <- unique(all_strain_filt_1$mut_id)

for (i in 1:ncol(pre_ab)) {
  curr_strain <- colnames(pre_ab)[i]
  curr_muts <- all_strain_filt_1$mut_id[all_strain$title == curr_strain]
  pre_ab[match(curr_muts, rownames(pre_ab)), i] <- 1  
}
# Filtering  ---- 
## mutations recovered by mapping the strain on itself: CONTROL and PAS
not_mut <- 
  rownames(pre_ab)[rowSums(pre_ab) == ncol(pre_ab)]

# Control is a the PAS133 strain after cultivation on solid medium
mut_spurious_ctl <- 
  rownames(pre_ab)[pre_ab[,"CONTROL"] == 1]
# PAS is obtained after the variant calling using the very reads used to construct the reference.   
mut_spurious_pas <-
  rownames(pre_ab)[pre_ab[,"PAS"] == 1]
mut_spurious_all <- 
  unique(c(not_mut, mut_spurious_ctl, mut_spurious_pas))

pre_ab <- 
  pre_ab[! rownames(pre_ab) %in% mut_spurious_all, ! 
           colnames(pre_ab) %in% c('CONTROL','PAS')]

all_strain_filt_1 |> 
  filter(mut_id %in% rownames(pre_ab)) |> 
  identity() -> all_strain_filt 

# Everything just disappears here...
# No mutations, no cry.

# # Distance ----
# dist(t(pre_ab), method = 'man')
# heatmap(as.matrix(dist(t(pre_ab), method = 'man')), symm = TRUE)
# plot(hclust(d = dist(t(pre_ab), method = 'man'), method = 'com'))
# 
# # Gene level convergence ----
# all_genes <- all_strain_filt$locus_tag
# all_genes <- 
#   str_replace_all(string = all_genes, pattern = "\\|", replacement = "/")
# all_genes <-
#   unlist(str_split(string = all_genes, pattern = '/'))
# table(all_genes) 
# # it seems that there is no convergence besides identical mutations in different
# # strain. I believe them to be completely artefactual.
# 
# write_tsv(x = all_strain_filt, file = 'all_strain_filt', col_names = TRUE)
