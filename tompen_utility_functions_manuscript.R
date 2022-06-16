# These are the tompen plotting functions used in the manuscript

# This function takes in a haplotype call file and adds a clolumn for Beta,
# which is based on the sQTL or eQTL hypothesis. 
# It also computes Exp[Beta] using the Hardy-Weinberg assumption, and the 
# empirical method, which simply counts heterozygotes. 
read_haps <- function(base_dir, tiss_name){
  haps_fp <- list.files(base_dir, pattern = tiss_name, full.names = T)
  haps_fp <- haps_fp[order(nchar(haps_fp), haps_fp)]
  out = tibble()
  for(i in haps_fp) {
    out <- rbind(
      stfu(read_tsv(i, 
                    col_types = list(
                      indv = col_character(),
                      gene = col_character(),
                      exon_coord = col_character(),
                      top_exon_coord = col_character(),
                      chr = col_character(),
                      csnp = col_character(),
                      effect = col_character(),
                      cadd_effect = col_character(),
                      CADD = col_double(),
                      csnp_af = col_double(),
                      csnp_af_gnomad = col_double(),
                      csnp_type = col_character(),
                      csnp_ref_allele = col_character(),
                      csnp_alt_allele = col_character(),
                      csnp_anc_allele = col_character(),
                      dummy = col_character(),
                      esnp = col_character(),
                      esnp_hi_inc_allele = col_character(),
                      esnp_low_inc_allele = col_character(),
                      esnp_anc_allele = col_character(),
                      slope = col_double(),
                      esnp_af = col_double(),
                      recom_esnp_af = col_double(),
                      n_hom_lower = col_double(),
                      n_het = col_double(),
                      n_hom_higher = col_double(),
                      chisq_p = col_double(),
                      haplotype = col_character()
                    )
      )),
      out
    )
  }
  
  # Fix compound homozygotes
  out %<>%
    group_by(indv, gene) %>%
    filter(CADD == max(CADD)) %>%
    ungroup
  
  # Add beta homozygous and beta flags
  hom <- as.logical(c(0, 0, 1, 1, 0, 0, 1, 1))
  names(hom) <- names(beta_config_sqtl)
  out$hom <- hom[out$haplotype]
  out$beta <- beta_config_sqtl[out$haplotype]
  
  out <- out[!is.na(out$hom),]
  # Get the expectation of beta for the homozygotes if we assume HWE
  out$exp_beta_hardy <- 
    ifelse(out$hom, 
           out$recom_esnp_af^2 / 
             ((1-out$recom_esnp_af)^2+out$recom_esnp_af^2), 
           .5
    )
  
  # And if we just use the empirical counts of higher-inclusion homozygotes
  out$exp_beta_empir <- 
    ifelse(out$hom, 
           (out$n_hom_higher + 1) / # Add a pseudocount to prevent the thing from breaking
             (out$n_hom_higher + out$n_hom_lower), 
           .5
    )
  
  out
}

run_enrichment <- function(tweaked_haplotypes) {
  rbind(
    bootstrap_test(filter(tweaked_haplotypes, CADD > 15), exp_beta = "exp_beta_empir"),
    bootstrap_test(filter(tweaked_haplotypes, CADD < 15), exp_beta = "exp_beta_empir")
  )
}

run_comparison_test <- function(tweaked_haplotypes){
  haps1 <- filter(tweaked_haplotypes, CADD > 15)
  haps2 <- filter(tweaked_haplotypes, CADD < 15)
  bootstrap_comparison_test(haps1, haps2)
}

plot_epsilon_plot <- function(x, comp_pval, pval_pos = .01, n_pos = .02) {
  x %>%
    as_tibble %>%
    mutate(ndel = c("del", "non_del")) %>%
    ggplot(aes(y = epsilon, x = "baseline", col = ndel,ymin = lower, ymax = upper)) + 
    geom_pointrange(position=position_dodge(width=.5)) + 
    geom_hline(yintercept = 0, lty = 3) +
    
    # The P-value label
    geom_text(mapping = aes(y = pval_pos, label = bootstrap_p), 
              position = position_dodge(width = .5)) +
    
    # The N-sample label
    geom_text(mapping = aes(y = n_pos, label = prettyNum(n_haplotypes, big.mark = ",")), 
              position = position_dodge(width = .5)) +
    
    # The comparison test p-value
    geom_text(y = -(mean(n_pos, pval_pos)), 
              x = 1, 
              col = "black", 
              label = comp_pval) +
    ylim(-n_pos, n_pos) + 
    scale_color_manual(values = c("red", "blue")) +
    ylab("ε") +
    theme_classic() + 
    theme(axis.text.y = element_blank(), 
          axis.title.y = element_blank(),
          legend.position = "none") +
    coord_flip() +
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(), 
          axis.line.y = element_blank())
}

plot_epsilon_plot_new <-  function(x, comp_pval, pval_pos = .01, n_pos = .02) {
  df <-  x %>%
    as_tibble %>%
    mutate(ndel = c("del", "non_del"))
  
  main_plt <- 
    ggplot(df, aes(y = epsilon, x = "baseline", col = ndel,ymin = lower, ymax = upper)) + 
    geom_pointrange(position=position_dodge(width=.5)) + 
    geom_hline(yintercept = 0, lty = 3) +
    ylim(-n_pos, n_pos) + 
    scale_color_manual(values = c("red", "blue")) +
    ylab("ε") +
    theme_classic() + 
    theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
          legend.position = "none") +
    coord_flip() +
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(), 
          axis.line.y = element_blank(),
          panel.border = element_rect(size = 1, fill = NA) ,
          plot.margin = unit(c(0, 0, 0, 0), "cm")
          )
  
  # Separate plot for just the p-vlaues
  pval_plt <- ggplot(df, aes(y = epsilon, x = "baseline", col = ndel,ymin = lower, ymax = upper)) + 
    geom_text(mapping = aes(y = .1, label = bootstrap_p), 
              position = position_dodge(width = .5)) +
    geom_text(mapping = aes(y = .2, label = n_haplotypes), 
              position = position_dodge(width = .5)) +
    scale_color_manual(values = c("red", "blue")) +
    ylim(.05, .25)+
    theme_classic() + 
    theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
          legend.position = "none") +
    coord_flip() +
    ggtitle("Comparison\np-value") +
    theme(axis.text.y = element_blank(), 
          axis.line.y = element_blank(),
          axis.title.x = element_blank(), 
          axis.text.x = element_blank(), 
          axis.ticks = element_blank(),
          panel.border = element_rect(size = 1, fill = NA), 
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          plot.title = element_text(hjust = .5))
  
  # Combine everything with plotgrid
  plot_grid(main_plt, NULL, pval_plt, ncol = 3, rel_widths = c(1, 0, .2), align = "h")
}

# A big function that just does everything
run_nonsynonymous_comparison <- function(tweaked_haplotypes, ...){
  
  # run each test individually
  message("running individual tests")
  res1 <- 
    rbind(
      bootstrap_test(tweaked_haplotypes %>% filter(effect != "synonymous_variant"), 
                     exp_beta = "exp_beta_empir"), 
      bootstrap_test(tweaked_haplotypes %>% filter(effect == "synonymous_variant"), 
                     exp_beta = "exp_beta_empir")
  )
  
  # run the comparison test
  message("running comparison tests")
  res2 <- 
    bootstrap_comparison_test(
      tweaked_haplotypes %>% filter(effect != "synonymous_variant"), 
      tweaked_haplotypes %>% filter(effect == "synonymous_variant"), 
      exp_beta = "exp_beta_empir"
    )['bootstrap_p']
  
  # plotting result
  plot_epsilon_plot(res1, res2, pval_pos = .04, n_pos = .05)
}

save_plot <- function(fn, save_fn = svg, ...){
  full_path <- paste0(
    "/gpfs/commons/groups/lappalainen_lab/jeinson/projects/modified_penetrance/mp_manuscript/figures/",
    fn)
  save_fn(full_path, ...)
  
}
