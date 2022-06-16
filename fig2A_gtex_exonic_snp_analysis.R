# This analyzes exon inclusion across all tissues available, to try to do
# something as similar to the analysis in Castel. et. al as possible. 
# 
# combines everything together into one file, and leaves the individual
# tissue results in the supplement. 

rm(list = ls())
setwd("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance")
source("~/myPackages.R")
source("mp_manuscript/tompen_utility_functions_manuscript.R")
library(ggplot2)

rv_bytiss <- readRDS("gtex_lof_in_exons/data/rare_variants_exon_psi_zscore_bytissue.rds")
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
snp_cons_order <- read_tsv("../../data/consdetail_rank_copy.tsv")$consdetail
rvs$consequence %<>% 
  str_split(",") %>%
  map_chr(~ names(which.min(sapply(.x, match, snp_cons_order))))

# Add pretty annotations
consequence_pretty_map <- read_tsv("../../data/consdetail_rank_pretty_names.csv") %>%
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

#### Median PSI Z-Score ####
# How about the median PSI Z-score?
df <- 
  rvs %>%
  group_by(Tissue, deleterious_status) %>%
  #filter(consequence != "stop_gained") %>%
  nest() %>% 
  mutate(boot_res = map(data,
                        ~ boot(data = .$exon_psi_Zscore,
                               statistic = function(x, i) median(x[i]),
                               R = 1000)),
         boot_res_ci = map(boot_res, boot.ci, type = "perc"),
         mean = map(boot_res_ci, ~ .$t0),
         lower_ci = map(boot_res_ci, ~ .$percent[[4]]),
         upper_ci = map(boot_res_ci, ~ .$percent[[5]]),
         n =  map(data, nrow)) %>% 
  select(-data, -boot_res, -boot_res_ci) %>% 
  unnest(cols = c(n, mean, lower_ci, upper_ci)) %>% 
  ungroup()

plt <- 
  ggplot(df, aes(Tissue, mean, fill = deleterious_status)) +
  geom_bar(stat = "identity", position = position_dodge(width = .9)) +
  geom_pointrange(aes(ymin = lower_ci, ymax = upper_ci), 
                  size = .5,
                  position = position_dodge(width = .8)) + 
  geom_text(aes(y = 0.4, label = prettyNum(unlist(n), big.mark = ",")), 
            angle = 90, hjust = 1,
            position = position_dodge(width = .9)) + 
  scale_fill_manual(values = c("red", "blue")) + 
  geom_hline(yintercept = 0, lty = 2) +
  ylab("Median PSI Zscore (+/- 95% Bootstrap CI)") + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) 

plt


# Plot everything together, for the sake of looking nice. 
ggplot(rvs, aes(deleterious_status, exon_psi_Zscore)) + 
  geom_boxplot()

# Phred vs. psi-Zscore
ggplot(rvs, aes(cadd_phred, exon_psi_Zscore)) + 
  geom_point() + 
  geom_smooth()
