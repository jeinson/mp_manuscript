# Figure 2: Plot some figures related to the Poisson-binomial test

library(tompen)
library(ggplot2)

source("tompen_utility_functions_manuscript.R")

# parameter DF

params <- data.frame(
  alpha = c(1, 2, .3, 50, 1.387),
  beta = c(1, 6, .3, 50, 0.954)
)

params$label <- paste0("⍺ = ", as.character(params$alpha), ", β = ", as.character(params$beta))

base <- ggplot(lwd = 2) + xlim(0, 1)

base_w_fns <- 
  base + 
  geom_function(aes(colour = params[1, 3]), fun = dbeta, args = list(shape1 = params[1,1], shape2 = params[1,2])) +
  geom_function(aes(colour = params[2, 3]), fun = dbeta, args = list(shape1 = params[2,1], shape2 = params[2,2])) +
  geom_function(aes(colour = params[3, 3]), fun = dbeta, args = list(shape1 = params[3,1], shape2 = params[3,2])) +
  geom_function(aes(colour = params[4, 3]), fun = dbeta, args = list(shape1 = params[4,1], shape2 = params[4,2])) +
  geom_function(aes(colour = params[5, 3]), fun = dbeta, args = list(shape1 = params[5,1], shape2 = params[5,2])) +
  scale_color_discrete(name = "Beta Parameters") + 
  gtex_v8_figure_theme()

save_plot("fig4A_binomial_distributions.svg", width = 4, height = 2)
base_w_fns
dev.off()
