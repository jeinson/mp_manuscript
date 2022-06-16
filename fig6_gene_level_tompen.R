# This script is used to run the poisson-binomial version of the modified 
# penetrance test on all gene-level haplotypes in the dataset. 

# I split this into its own script because the poisson-binomial thing
# has a tendency to crash...
# 
# The idea here is to use the exact same set of variants as in the cross gene
# analysis, and run the modified penetrance test one gene at a time. 

rm(list = ls())
setwd("/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/")
source("~/myPackages.R")
source("mp_manuscript/tompen_utility_functions_manuscript.R")
library(tompen)
library(ggplot2)

# Read in the filtered haplotype file
topmed_hap <- read_rds("topmed/topmed_haplotypes_v8_anno_full_with_details.rds")

topmed_hap %<>%
  filter(csnp_count <= 10 & csnp_count > 1) 

# Haplotypes per gene filter (The final way I'm doing this)
haps_per_gene <- table(topmed_hap$gene)
usable_genes <- true_names(!flag_outliers(haps_per_gene))
topmed_hap <- topmed_hap %>% filter(gene %in% usable_genes)

# Split the data frame of all haplotypes into a list of data frames where
# each one contains haplotype records from one gene only 
# Also filter down to csnps with MAF < 1e-3. (Exon length filtering 
# is already done, so don't do it again)
topmed_all_hap_by_gene <- 
  topmed_hap %>%
  split(.$gene) 

genes_with_enough_haps <- map_dbl(topmed_all_hap_by_gene, nrow) > 10
genes_with_enough_haps <- true_names(genes_with_enough_haps)

# How many haplotypes do we get per gene after filtering?
hist(map_dbl(topmed_all_hap_by_gene[genes_with_enough_haps], nrow))

# Make a new pb test that doesn't suck
poisson_binomial_test_new <- 
  function(haps){
    
    x <- haps
    epsilon <- mean((as.numeric(x$beta) - x$exp_beta))
    
    p1 <- poibin::ppoibin(sum(x$beta), x$exp_beta)
    p2 <- 1 - p1
    p <- 2 * min(p1, p2)
    
    return(c(poison_binomial_p = p,
             epsilon = epsilon,
             n_haplotypes = nrow(x)
    )
    )
  }

# Run the test gene by gene ----
tot = length(genes_with_enough_haps)

#--- Make a progress bar
library(progress)
pb <- progress_bar$new(total = tot)
gene_level_penetrance_tests <- list()

for(i in 1:tot){
  this_gene <- genes_with_enough_haps[i]
  x <- topmed_all_hap_by_gene[[this_gene]]
  x$exp_beta <- x$exp_beta_empir
  
  # First test everything together
  full_result <- 
    tryCatch(
      poisson_binomial_test_new(x), 
      error = function(e) c(NA, NA, nrow(x))
    )
  
  # Second split the haplotypes split into deleterious / non-deleterious groups
  del <- x %>% filter(CADD > 15)
  ndel <- x %>% filter(CADD < 15)
  
  #--- Deleterious First
  ifelse(
    nrow(del) > 10,
    
    deleterious_result <- 
      tryCatch(
        poisson_binomial_test_new(del)
        , 
        error = function(e) c(NA, NA, nrow(del))
      ), 
    # otherwise:
    deleterious_result <- c(NA, NA, nrow(del))
  )
  
  #--- Non-deleterious second
  ifelse(
    nrow(ndel) > 10,
    control_result <- 
      tryCatch(
        poisson_binomial_test_new(ndel), 
        error = function(e) c(NA, NA, nrow(ndel))
      ), 
    # otherwise:
    control_result <- c(NA, NA, nrow(ndel))
  )
  
  # Now combine everything back together
  gene_level_penetrance_tests[[this_gene]] <- 
    rbind(full_result, deleterious_result, control_result)
  
  # Tick the progress bar
  pb$tick()
}


### Analysis ####
# Bind things together and run similar tests to the other script
ptests_df <- 
  gene_level_penetrance_tests %>%
  map(as.data.frame) %>%
  map(rownames_to_column, var = "test_type") %>%
  bind_rows(.id = "gene")

ptest_df_wide <-   
  ptests_df %>% 
  pivot_wider(
    names_from = "test_type", 
    values_from = c("poison_binomial_p", "epsilon", "n_haplotypes")) %>% 
  select(-starts_with("V"))

# Join sQTL information
sqtls <- read_tsv("sQTL_v8_anno/cross_tissue_top_sQTLs/top_sQTLs_MAF05.tsv")
ptest_df_wide %<>% 
  left_join(sqtls, by = c("gene" = "group")) %>%
  left_join(distinct(topmed_hap %>% select(gene, oe_lof)))

ptest_df_wide %>%
  select(starts_with("epsilon")) %>%
  colMeans(na.rm = T)

### Look at sQTL allele frequencies to see high frequency sQTLs cause instability
#All haplotypes
save_plot("fig6_all_variants_volcano.svg", width = 6, height = 4)
ggplot(ptest_df_wide, 
       aes(epsilon_full_result, 
           -log10(poison_binomial_p_full_result), 
           color = log10(n_haplotypes_full_result), 
           shape = sQTL_af < .5
           )) + 
  geom_point() +
  scale_color_gradient(low = "#5D3A9B", high = "#E66100") + 
  geom_hline(yintercept = -log10(.05)) +
  geom_hline(yintercept = -log10(.05/nrow(ptest_df_wide)), lty = 2, lwt = .5) + 
  theme_linedraw()
dev.off()

#### SCRATCH ZONE ####
# Can I just figure out what's going on with these weird genes with pval=0, and 
# all with epsilon > 0???
weird_genes <- 
  ptest_df_wide %>% 
  filter(-log10(poison_binomial_p_full_result) > 10) %$%
  gene
gene_level_penetrance_tests

weird_haps <- 
  topmed_hap %>% 
  filter(gene %in% weird_genes)

table(weird_haps$beta)
table(weird_haps$gene)
filter(topmed_hap, gene == "ENSG00000151632") %>% View

####
# Deleterious haplotypes
save_plot("fig6_volcano_deleterious_csnps.svg", width = 5, height = 4)
ggplot(ptest_df_wide, 
       aes(epsilon_deleterious_result, 
           -log10(poison_binomial_p_deleterious_result), 
           #color = log10(n_haplotypes_deleterious_result), 
           color = -log(oe_lof), 
           shape = sQTL_af < .5
       )) + 
  geom_point() +
  scale_color_gradient(low = "#5D3A9B", high = "#E66100") + 
  geom_hline(yintercept = -log10(.05)) +
  geom_hline(yintercept = -log10(.05/nrow(ptest_df_wide)), lty = 2, lwt = .5) + 
  theme_linedraw()
dev.off()

# Is the epsilon value correlated with constraint?
ggplot(ptest_df_wide, 
       aes(oe_lof, epsilon_deleterious_result, 
           color = delta_psi)
       ) + 
  scale_color_gradient(low = "#5D3A9B", high = "#E66100") + 
  geom_hline(yintercept = 0, col = "red", lty = 2) +
  geom_point() + 
  gtex_v8_figure_theme()

# Are there even more genes with negative epsilon values?
binom.test(sum(ptest_df_wide$epsilon_deleterious_result < 0, na.rm = T), 
           sum(!is.na(ptest_df_wide$epsilon_deleterious_result)))

####

# Non-deleterious haplotypes
ggplot(ptest_df_wide, aes(epsilon_control_result, -log10(poison_binomial_p_control_result))) + 
  geom_point(aes(color = sQTL_af, size = n_haplotypes_control_result), alpha = .5) +
  scale_color_gradient(low = "cornflowerblue", high = "firebrick") + 
  geom_hline(yintercept = -log10(.05)) +
  geom_hline(yintercept = -log10(.05/nrow(ptest_df_wide)), lty = 2, lwt = .5) + 
  theme_linedraw()

### Test genes that are depleted for high-penetrance haplotypes with deleterious
### rare variants. 
genes_w_hp_lof_depletion <- 
  ptest_df_wide %>%
  filter(p.adjust(poison_binomial_p_deleterious_result, method = 'fdr') < .05) %>%
  filter(epsilon_deleterious_result < 0) %$%
  gene

message(paste(genes_w_hp_lof_depletion, "\n"))

sqtls %>% 
  filter(group %in% genes_w_hp_lof_depletion) %>%
  left_join(topmed_hap %>% select(gene, oe_lof) %>% distinct, 
            by = c("group" = "gene")) %$%
  oe_lof %>% hist

sqtls_with_label <- sqtls %>%
  mutate(hp_depleted = group %in% genes_w_hp_lof_depletion) %>%
  left_join(topmed_hap %>% select(gene, oe_lof) %>% distinct, 
            by = c("group" = "gene")) 

ggplot(sqtls_with_label, aes(oe_lof, fill = hp_depleted)) + 
  geom_density(alpha = .5)

# Comparison between deleterious and non-deleterious variants
call_significance_code <- function(del_pval, con_pval){
  if(!is.na(del_pval) & !is.na(con_pval)){
    ifelse(del_pval < .05, 
           ifelse(con_pval < .05, "Del *, Con *", "Del *, Con -"),
           ifelse(con_pval < .05, "Del -, Con *", "Del -, Con -")
    )
  } else {
    NA
  }
}

# Add a column for FDR corrected p-values
ptest_df_wide_fdr <- ptest_df_wide %>%
  mutate(across(.cols = starts_with('poison_binomial_p'), 
         .fns = ~ p.adjust(.x, method = 'fdr'), 
         .names = "{col}_adj"))

ptest_df_wide_fdr$test_results <- 
  with(ptest_df_wide_fdr, 
       mapply(call_significance_code, 
              poison_binomial_p_deleterious_result_adj, 
              poison_binomial_p_control_result_adj)
  )

pval_comparison_plt <- 
  ggplot(ptest_df_wide_fdr, aes(epsilon_deleterious_result, epsilon_control_result)) + 
  theme_gray() + 
  geom_hline(yintercept = 0, color = "red", lty = 2) + 
  geom_vline(xintercept = 0, color = "red", lty = 2) + 
  geom_smooth(method = "lm") +
  geom_point(aes(color = test_results)) +
  #scale_size_manual(values = c(1, 2, 2, 2)) + 
  scale_color_manual(values = c("grey", "blue", "darkgreen", "red")) +
  xlab("HP haplotype depletion in genes with deleterious RVs") + 
  ylab("HP haplotype depletion in genes non-deleterious RVs") +
  theme_classic()

save_plot("fig6B_epsilon_comparisons.svg", width = 5, height = 4)
pval_comparison_plt
dev.off()

# What's the significance of this correlation?
del_mod <- lm(epsilon_control_result ~ epsilon_deleterious_result, 
              data = ptest_df_wide)
summary(del_mod)


save_plot("fig6B_epsilon_comparisons_inset.svg", width = 5, height = 5)
ggplot(ptest_df_wide_fdr %>% filter(!is.na(test_results)), aes(test_results, fill = test_results)) + 
  geom_bar() + 
  scale_fill_manual(values = c("grey", "blue", "darkgreen", "red")) + 
  theme_classic() +
  theme(legend.position = "none")
dev.off()

# Make a contingency table for this
result_ctable <- 
  ptest_df_wide_fdr %>%
  filter(!is.na(poison_binomial_p_control_result) & !is.na(poison_binomial_p_deleterious_result)) %$%
  table(poison_binomial_p_control_result_adj < .05, poison_binomial_p_deleterious_result_adj < .05)

fisher.test(result_ctable)

plot(result_ctable, frame.plot = T, xlab = "Control", ylab = "Deleterious")
