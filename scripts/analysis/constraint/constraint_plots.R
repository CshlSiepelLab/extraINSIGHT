library(ggplot2)
library(cowplot)
source("extract_result_lib.R")


#############################################
## Top 5% of expressed genes in each tissue
#############################################
result_dir <- "../../../results/grch38/gnomad_v3.0/constraint/gene_annotations"
directory_list <- dir(result_dir, full.names = TRUE)
plot_dir <- file.path(result_dir, "plots")
directory_list <- directory_list[!grepl("log|plots", directory_list)]

out_abs <- extract_conservation_estimates(directory_list, FALSE)
#out_resid <- extract_conservation_estimates(directory_list, TRUE)

ord <- out_abs[name == "strong_selection"][order(absolute, decreasing = TRUE)]$label
dir.create(plot_dir, FALSE, TRUE)

## Absolute selection
pos <- position_dodge(width=0.9)
ga <- ggplot(out_abs, aes(x= factor(label, levels = ord), y = absolute, fill = name))+
    geom_bar(stat = "identity", position = pos)+
    geom_errorbar(aes(ymin = absolute - 1.96*se, ymax = absolute + 1.96*se), width = 0.2, position = pos)+
    theme_cowplot()+
    ylim(-0.2,1)+
    scale_fill_brewer()+
    theme(legend.position=c(0.65,0.85), axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_blank())+
    xlab("Annotation")+
    ylab("Parameter estimate")

ggsave("gene_annotation_contraint.pdf", plot = ga, path = plot_dir, width = 6, height = 6)

#############################################
## Top 5% of expressed genes in each tissue
#############################################
result_dir <- "../../../results/grch38/gnomad_v3.0/constraint/top_5_expressed"
directory_list <- dir(result_dir, full.names = TRUE)
plot_dir <- file.path(result_dir, "plots")

directory_list <- directory_list[!grepl("log", directory_list)]

out_abs <- extract_conservation_estimates(directory_list, FALSE)
out_resid <- extract_conservation_estimates(directory_list, TRUE)

ord <- out_abs[name == "strong_selection"][order(absolute, decreasing = TRUE)]$label
dir.create(plot_dir, FALSE, TRUE)

## Absolute selection
pos <- position_dodge(width=0.9)
ga <- ggplot(out_abs, aes(x= factor(label, levels = ord), y = absolute, fill = name))+
    geom_bar(stat = "identity", position = pos)+
    geom_errorbar(aes(ymin = absolute - 1.96*se, ymax = absolute + 1.96*se), width = 0.2, position = pos)+
    theme_cowplot()+
    ylim(0,1)+
    scale_fill_brewer()+
    theme(legend.position=c(0.65,0.85), axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_blank())+
    xlab("Annotation")+
    ylab("Parameter estimate")

ggsave("top_five_percent_tissue_expressed_contraint.pdf", plot = ga, path = plot_dir, width = 16, height = 7)

## TBD comparisons:
##  - Ultra-conserved non-coding sequence
##  - Pleiotropy based on gene expression (both promoters and cds)
##  - Haplo-insufficient vs. not based on pLI (both promoters and cds)

#############################################
## Top 5% of expressed genes in each tissue
#############################################
