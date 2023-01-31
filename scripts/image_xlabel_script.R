# Plot haplotype configuration chart
library(EBImage)
library(ggplot2)
library(magrittr)
library(grid)
library(ggthemes)
library(scales)

# Source the script used to add images to a barplot
#devtools::source_gist("1d1bdb00a7b3910d62bf3eec8a77b4a7")

# Function for adding images to plot, from 
# https://jcarroll.com.au/2016/06/02/images-as-x-axis-labels/
add_images_as_xlabels <- function(gg, pics) {
  
  ## ensure that the input is a ggplot
  if(!inherits(gg, "ggplot")) stop("Requires a valid ggplot to attach images to.")
  
  ## extract the components of the ggplot
  gb   <- ggplot_build(gg)
  #xpos <- gb$layout$panel_params[[1]]$x.major
  xpos <- gb$layout$panel_params[[1]]$x.sec$break_positions()
  yrng <- gb$layout$panel_params[[1]]$y.range
  
  ## ensure that the number of pictures to use for labels 
  ## matches the number of x categories
  if(length(xpos) != length(pics)) stop("Detected a different number of pictures to x categories")
  
  ## create a new grob of the images aligned to the x-axis
  ## at the categorical x positions
  my_g <- do.call("grobTree", Map(rasterGrob, pics, x=xpos, y=0))
  
  ## annotate the original ggplot with the new grob
  gg <- gg + annotation_custom(my_g,
                               xmin = -Inf, 
                               xmax =  Inf,
                               ymax = yrng[1] + 0.25*(yrng[2]-yrng[1])/npoints, 
                               ymin = yrng[1] - 1.25*(yrng[2]-yrng[1])/npoints) +
    xlab(NULL)
  
  ## turn off clipping to allow plotting outside of the plot area
  gg2 <- ggplotGrob(gg)
  gg2$layout$clip[gg2$layout$name=="panel"] <- "off"
  
  ## produce the final, combined grob
  grid.newpage()
  grid.draw(gg2)
  
  return(invisible(NULL))
  
}

plot_haplotype_pileup <- function(x){
  # rearrange the data to make to condense it and add error bars
  x %>%
    group_by(pathogenic, haplotype) %>%
    summarise(count = n()) %>%
    ungroup(haplotype) %>%
    mutate(percent = count / sum(count)) %>%
    ungroup() %>%
    split(.$pathogenic) %>%
    lapply( function(x) {
      x$group <- c(1,1,2,2,3,3,4,4)
      x %>% 
        group_by(group) %>%
        mutate(percent = sum(percent)) %>%
        mutate(count = sum(count)) %>%
        slice(1) %>%
        ungroup -> x
      
      tot <- sum(x$count)
      Sp <- sqrt((x$percent*(1-x$percent)) / (x$count+.1))
      x$upper <- x$percent + Sp
      x$lower <- x$percent - Sp
      return(x)
    }
    ) %>%
    do.call("rbind", .) %>%
    
    # plot it, add the axis labels
    ggplot(aes(haplotype, percent, fill = pathogenic)) + 
    geom_bar(stat = "identity", position = "dodge") + 
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.5,
                  position=position_dodge(.9)) +
    geom_text(aes(y = percent + .02, label = count), position = position_dodge(width = 1)) + 
    scale_fill_manual(values = c(pathogenic="#E23A28", benign="#F1988E", control="#5978B9")) +
    theme_minimal() +
    theme(plot.margin = unit(c(0.5,0.5,5,0.5), "lines"),
          axis.text.x = element_blank()) -> 
    
    plt
  
  add_images_as_xlabels(plt, hap_pics[c(1,3,5,7)])
}
