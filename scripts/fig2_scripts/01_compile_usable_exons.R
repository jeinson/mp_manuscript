# This script compiles all exons that pass filtering across all tissues, 
# and exports it in a format that can be used for haplotype calling. 

# The first in a series of analysis scripts

setwd("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance")
source("~/myPackages.R")

tissues <- read_lines("sQTL_v8_anno/completed_tissues.txt")

psi_fp <- function(tiss){
  sprintf("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/sQTL_v8_anno/bed/%s_psi_v8_protocol.bed.gz", tiss)
}

tiss_list <- list()
for(tiss in tissues){
  message(tiss)
  tiss_list[[tiss]] <- read_tsv(psi_fp(tiss))[,1:4]
}

all_exons <- bind_rows(tiss_list)
all_exons <- distinct(all_exons)

all_exons <- arrange(all_exons, nchar(`#Chr`), `#Chr`, start)

# Add the exon and gene info
exon_id_map <- readRDS("../../data/gtex_stuff/gtex_v8_exon_id_map.rds")
id_exon_map <- names(exon_id_map)
names(id_exon_map) <- exon_id_map

all_exons$exon_id <- id_exon_map[all_exons$ID]
all_exons$gene_id <- all_exons$exon_id %>% str_remove("_.*") %>% remove_trailing_digit()

write_tsv(all_exons, "gtex_lof_in_exons/all_exons_passing_filtering.bed")
