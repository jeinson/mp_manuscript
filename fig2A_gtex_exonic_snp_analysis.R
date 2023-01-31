# This analyzes exon inclusion across all tissues available, to try to do
# something as similar to the analysis in Castel. et. al as possible. 
# 
# combines everything together into one file, and leaves the individual
# tissue results in the supplement. 

rm(list = ls())

# Set wd to wherever mp_manuscript is located
library(tidyverse)
library(magrittr)
source("tompen_utility_functions_manuscript.R")
library(ggplot2)

rv_bytiss <- readRDS("data/fig2_data/rare_variants_exon_psi_zscore_bytissue.rds")
tissues <- names(rv_bytiss)

rv_bytiss <- 
  rv_bytiss %>%
  map(~ .x %>% filter(!is.na(psi) & !is.na(cadd_phred)))

#### Filtering and Data Processing ####
# Combine all rare variants into one file
rvs <- bind_rows(rv_bytiss, .id = "Tissue")

# Flag variants that are within 10bp of the start or end of an exon. 
rvs <- 
  rvs %>%
  mutate(exon_start = exon_coord %>% str_split("_") %>% map_chr(2) %>% as.numeric()) %>%
  mutate(exon_end   = exon_coord %>% str_split("_") %>% map_chr(3) %>% as.numeric()) %>%
  mutate(snp_flag = pos - exon_start < 10 | exon_end - pos < 10)

# Filter variants with a gnomad_af < .005
rvs <- rvs %>% filter(gnomad_af < .005)

# Clean up the annotations
snp_cons_order <- read_tsv("data/fig2_data/consdetail_rank_copy.tsv")$consdetail
rvs$consequence %<>% 
  str_split(",") %>%
  map_chr(~ names(which.min(sapply(.x, match, snp_cons_order))))

# Add pretty annotations
consequence_pretty_map <- read_tsv("data/fig2_data/consdetail_rank_pretty_names.csv") %>%
  select(consdetail, pretty_name) %>%
  deframe()
rvs %<>%
  mutate(pretty_name = consequence_pretty_map[consequence])

# Further classify the variants

rvs %<>%
  mutate(variant_class = ifelse(consequence == "synonymous", "synonymous", 
                                ifelse(consequence == "missense", 
                                       ifelse(cadd_phred > 15, 
                                              "pathogenic_missense", 
                                              "benign_missense"),
                                       "other")))

rvs %<>% mutate(deleterious_status = ifelse(cadd_phred > 15, "deleterious", "non_deleterious"))

#### Run statistical tests for difference ----
# Look at overall changes in inclusion
rvs_cross <-
  rvs %>%
  group_by(indv, var) %>%
  mutate(exon_psi_Zscore = mean(exon_psi_Zscore)) %>%
  slice(1) %>%
  ungroup %>%
  filter(!snp_flag)

# Split by deleterious vs. non-deleterious
table(rvs_cross$deleterious_status)
with(rvs_cross, tapply(exon_psi_Zscore, deleterious_status, mean))

# Test if there is a difference
wilcox.test(exon_psi_Zscore ~ deleterious_status, data = rvs_cross)

#### Summary plots of all exons and PSI ----
# What do these distributions look like in general?
psi_density_plt <- 
  rvs %>% select(exon_id, psi, deleterious_status) %>% distinct %>%
  ggplot(aes(psi, fill = deleterious_status)) + 
  geom_histogram(bins = 20, col = "black") + 
  scale_fill_manual(values = c("red", "blue")) +
  theme_classic()  +
  theme(legend.position = "none") 

psi_Zscore_plt <- 
  rvs %>% select(exon_id, exon_psi_Zscore) %>% distinct %>%
  ggplot(aes(exon_psi_Zscore)) +
  geom_histogram(fill = "cornflowerblue", col = "black") + 
  xlim(-5,5) +
  theme_classic() 

psi_Zscore_rank_plt <- 
  ggplot(rvs, aes(rank(exon_psi_Zscore), fill = deleterious_status)) +
  geom_density(alpha = .5) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_fill_manual(values = c("red", "blue")) +
  xlab("PSI Z-Score Rank") +
  theme_classic() +
  theme(legend.position = "none") 

# Save plots 
save_plot("fig2S1_psi_density_plt.svg", width = 4, height = 3)
psi_density_plt
dev.off()

save_plot("fig2S2_psi_Zscore_plt.svg", width = 4, height = 3)
psi_Zscore_plt
dev.off()

save_plot("fig2S3_psi_density_rank_plt.svg", width = 4, height = 3)
psi_Zscore_rank_plt
dev.off()

# psi-Z over tissue type ----
# Calculate bootstraped confidence intervals
# https://stackoverflow.com/questions/35953394/calculating-length-of-95-ci-using-dplyr
library(boot)
library(tidyr)
df <- 
  rvs %>%
  filter(!snp_flag) %>%
  group_by(Tissue, deleterious_status) %>%
  #filter(consequence != "stop_gained") %>%
  nest() %>% 
  mutate(boot_res = map(data,
                        ~ boot(data = .$exon_psi_Zscore,
                               statistic = function(x, i) mean(x[i]),
                               R = 1000)),
         boot_res_ci = map(boot_res, boot.ci, type = "perc"),
         mean = map(boot_res_ci, ~ .$t0),
         lower_ci = map(boot_res_ci, ~ .$percent[[4]]),
         upper_ci = map(boot_res_ci, ~ .$percent[[5]]),
         n =  map(data, nrow)) %>% 
  select(-data, -boot_res, -boot_res_ci) %>% 
  unnest(cols = c(n, mean, lower_ci, upper_ci)) %>% 
  ungroup(deleterious_status) %>%
  arrange(deleterious_status) %>%
  mutate(Tissue = paste0(str_replace(Tissue, "_", " "), "\n", 
                         paste0(prettyNum(n, big.mark = ","), collapse = " | "))) %>%
  ungroup

plt <- 
  ggplot(df, aes(Tissue, mean, fill = deleterious_status)) +
  geom_bar(stat = "identity", position = position_dodge(width = .9)) +
  geom_pointrange(aes(ymin = lower_ci, ymax = upper_ci), 
                  size = .5,
                  position = position_dodge(width = .8)) + 
  scale_fill_manual(values = c("red", "blue")) + 
  geom_hline(yintercept = 0, lty = 2) +
  ylab("Mean PSI Zscore (+/- 95% Bootstrap CI)") + 
  gtex_v8_figure_theme() + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) 

save_plot("fig2A_psiZ_by_tissue.svg", width = 8, height = 4)
plt
dev.off()

### Plot Across Variant Annotations ####
# How does this difference look across variant annotations?
usable_cons <- names(sort(table(rvs$consequence),decreasing = T)[1:13])

non_splice_cons_df <- 
  rvs %>%
  filter(!snp_flag & consequence != "splice") %>%
  group_by(indv, var) %>%
  filter(abs(exon_psi_Zscore) == max(abs(exon_psi_Zscore))) %>%
  ungroup %>%
  group_by(pretty_name, deleterious_status) %>%
  nest() %>% 
  mutate(boot_res = map(data,
                        ~ boot(data = .$exon_psi_Zscore,
                               statistic = function(x, i) mean(x[i]),
                               R = 1000)),
         boot_res_ci = map(boot_res, boot.ci, type = "perc"),
         mean = map(boot_res_ci, ~ .$t0),
         lower_ci = map(boot_res_ci, ~ .$percent[[4]]),
         upper_ci = map(boot_res_ci, ~ .$percent[[5]]),
         n =  map(data, nrow)) %>% 
  select(-data, -boot_res, -boot_res_ci) %>% 
  unnest(cols = c(n, mean, lower_ci, upper_ci)) %>% 
  ungroup()

non_splice_cons_df %<>%
  group_by(pretty_name) %>%
  arrange(deleterious_status) %>%
  mutate(pretty_name = paste0(pretty_name, "\n", paste(prettyNum(n,big.mark = ","), collapse = " | ")))

cons_order <- non_splice_cons_df %>% 
  group_by(pretty_name) %>% 
  summarize(total = sum(n)) %>% 
  arrange(desc(total)) %$% 
  pretty_name

non_splice_cons_df %<>%
  filter(pretty_name %in% cons_order[1:13]) %>%
  mutate(pretty_name = factor(pretty_name, level = cons_order[1:13]))

non_splice_variant_plot <- 
  ggplot(non_splice_cons_df, aes(pretty_name, mean, fill = deleterious_status)) +
  geom_bar(stat = "identity", position = position_dodge(width = .9)) +
  geom_pointrange(aes(ymin = lower_ci, ymax = upper_ci), 
                  size = .5,
                  position = position_dodge(width = .9)) + 
  scale_fill_manual(values = c("red", "blue")) +  
  geom_hline(yintercept = 0, lty = 2) +
  ylab("Mean PSI Zscore (+/- 95% Bootstrap CI)") + 
  xlab("Rare variant consequence\nN Deleterious | N Non-Deleterious") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1), 
        legend.position = 'none') 

non_splice_variant_plot

# Now plot just splice-variants
splice_cons_df <- 
  rvs %>%
  filter(snp_flag | consequence == "splice") %>%
  group_by(indv, var) %>%
  filter(abs(exon_psi_Zscore) == max(abs(exon_psi_Zscore))) %>%
  ungroup %>%
  group_by(deleterious_status) %>%
  nest() %>% 
  mutate(boot_res = map(data,
                        ~ boot(data = .$exon_psi_Zscore,
                               statistic = function(x, i) mean(x[i]),
                               R = 1000)),
         boot_res_ci = map(boot_res, boot.ci, type = "perc"),
         mean = map(boot_res_ci, ~ .$t0),
         lower_ci = map(boot_res_ci, ~ .$percent[[4]]),
         upper_ci = map(boot_res_ci, ~ .$percent[[5]]),
         n =  map(data, nrow)) %>% 
  select(-data, -boot_res, -boot_res_ci) %>% 
  unnest(cols = c(n, mean, lower_ci, upper_ci)) %>% 
  ungroup()

splice_cons_df$pretty_name <- "Within 10bp of Splice Site"

splice_cons_df %<>%
  group_by(pretty_name) %>%
  arrange(deleterious_status) %>%
  mutate(pretty_name = paste0(pretty_name, "\n", paste(prettyNum(n,big.mark = ","), collapse = " | "))) %>%
  mutate(pretty_name = factor(pretty_name))

# splice_variant_plot <- 
#   ggplot(splice_cons_df, aes(pretty_name, mean, fill = deleterious_status)) +
#   geom_bar(stat = "identity", position = position_dodge(width = .9)) +
#   geom_pointrange(aes(ymin = lower_ci, ymax = upper_ci), 
#                   size = .5,
#                   position = position_dodge(width = .9)) + 
#   scale_fill_manual(values = c("red", "blue")) +  
#   geom_hline(yintercept = 0, lty = 2) +
#   ylab("Mean PSI Zscore (+/- 95% Bootstrap CI)") + 
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) 
# 
# splice_variant_plot


# Combine the two together
# library(cowplot)
# plot_grid(non_splice_variant_plot, splice_variant_plot, align = "h", 
#           rel_widths = c(5, 5))

save_plot("fig2B_mean_psiZ_cross_annotations.svg", width = 6, height = 3)
rbind(non_splice_cons_df, splice_cons_df) %>%
ggplot(aes(pretty_name, mean, fill = deleterious_status)) +
  geom_bar(stat = "identity", position = position_dodge(width = .9)) +
  geom_pointrange(aes(ymin = lower_ci, ymax = upper_ci), 
                  size = .3,
                  position = position_dodge(width = .9)) + 
  scale_fill_manual(values = c("red", "blue")) +  
  geom_hline(yintercept = 0, lty = 2) +
  ylab("Mean PSI Zscore\n(+/- 95% Bootstrap CI)") + 
  xlab("Rare variant consequence\nN Deleterious | N Non-Deleterious") +
  gtex_v8_figure_theme() + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) 
dev.off()

message(
"R version 4.2.1 (2022-06-23)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS Ventura 13.0

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] boot_1.3-28     magrittr_2.0.3  forcats_0.5.2   stringr_1.4.1   dplyr_1.0.10   
 [6] purrr_0.3.5     readr_2.1.3     tidyr_1.2.1     tibble_3.1.8    ggplot2_3.3.6  
[11] tidyverse_1.3.2

loaded via a namespace (and not attached):
 [1] tidyselect_1.2.0    lattice_0.20-45     splines_4.2.1       haven_2.5.1        
 [5] gargle_1.2.1        colorspace_2.0-3    vctrs_0.5.0         generics_0.1.3     
 [9] mgcv_1.8-41         utf8_1.2.2          rlang_1.0.6         pillar_1.8.1       
[13] glue_1.6.2          withr_2.5.0         DBI_1.1.3           bit64_4.0.5        
[17] dbplyr_2.2.1        modelr_0.1.9        readxl_1.4.1        lifecycle_1.0.3    
[21] munsell_0.5.0       gtable_0.3.1        cellranger_1.1.0    rvest_1.0.3        
[25] labeling_0.4.2      tzdb_0.3.0          parallel_4.2.1      fansi_1.0.3        
[29] broom_1.0.1         backports_1.4.1     scales_1.2.1        googlesheets4_1.0.1
[33] vroom_1.6.0         jsonlite_1.8.3      farver_2.1.1        fs_1.5.2           
[37] bit_4.0.4           hms_1.1.2           digest_0.6.30       stringi_1.7.8      
[41] grid_4.2.1          cli_3.4.1           tools_4.2.1         crayon_1.5.2       
[45] pkgconfig_2.0.3     Matrix_1.5-1        ellipsis_0.3.2      xml2_1.3.3         
[49] reprex_2.0.2        googledrive_2.0.0   lubridate_1.8.0     assertthat_0.2.1   
[53] httr_1.4.4          rstudioapi_0.14     R6_2.5.1            nlme_3.1-160       
[57] compiler_4.2.1")


