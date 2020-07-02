library(ggplot2)
library(cowplot)
library(stringr)
source("extract_result_lib.R")

plot_dir <- file.path("../../../results/grch38/gnomad_v3.0/constraint", "plots")
dir.create(plot_dir, F, T)
#############################################
## Annotations of genic regions
#############################################
result_dir <- "../../../results/grch38/gnomad_v3.0/constraint/gene_annotations"
directory_list <- dir(result_dir, full.names = TRUE)
directory_list <- directory_list[!grepl("log|plots", directory_list)]

out_abs <- extract_conservation_estimates(directory_list, FALSE)
#out_resid <- extract_conservation_estimates(directory_list, TRUE)
out_abs[,label := gsub("_gencode33","", label)]

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

ggsave("gene_annotation_contraint.pdf", plot = ga, path = plot_dir, width = 6, height = 5)

#############################################
## Top 5% of expressed genes in each tissue
#############################################
result_dir <- "../../../results/grch38/gnomad_v3.0/constraint/top_5_expressed"
directory_list <- dir(result_dir, full.names = TRUE)
#plot_dir <- file.path(result_dir, "plots")

directory_list <- directory_list[!grepl("log|plots", directory_list)]

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

ggsave("top_five_percent_tissue_expressed_contraint.pdf", plot = ga, path = plot_dir, width = 17, height = 6)

#############################################
## Ultra-conserved non-coding sequence
#############################################
result_dir <- "~/Projects/extraINSIGHT/results/grch38/gnomad_v3.0/constraint/ucne"
directory_list <- dir(result_dir, full.names = TRUE)
#plot_dir <- file.path(result_dir, "plots")

directory_list <- directory_list[!grepl("log|plots", directory_list)]

out_abs <- extract_conservation_estimates(directory_list, FALSE)
# out_resid <- extract_conservation_estimates(directory_list, TRUE)

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

ggsave("ucne_contraint.pdf", plot = ga, path = plot_dir, width = 6, height = 6)


#############################################
## pLI based haploinsufficiency
#############################################
result_dir <- "~/Projects/extraINSIGHT/results/grch38/gnomad_v3.0/constraint/pli"
directory_list <- dir(result_dir, full.names = TRUE)
#plot_dir <- file.path(result_dir, "plots")

directory_list <- directory_list[!grepl("log|plots", directory_list)]

out_abs <- extract_conservation_estimates(directory_list, FALSE)
label_matrix <- str_match(out_abs$label, "^([a-z]+)_([a-z]+)_(.+)")
out_abs[, region := label_matrix[,4]]
out_abs[, pLI := label_matrix[,3]]
out_abs[, pLI := factor(pLI, levels = c("high", "low"))]

ord <- out_abs[name == "strong_selection"][order(absolute, decreasing = TRUE)]$label
dir.create(plot_dir, FALSE, TRUE)

## Absolute selection
pos <- position_dodge(width=0.9)
ga <- ggplot(out_abs, aes(x = pLI, y = absolute, fill = name))+
    facet_wrap(.~region, nrow = 1)+
    geom_bar(stat = "identity", position = pos)+
    geom_errorbar(aes(ymin = absolute - 1.96*se, ymax = absolute + 1.96*se), width = 0.2, position = pos)+
    theme_cowplot()+
    scale_y_continuous(breaks = pretty(seq(-0.3,1,0.1)), limits = c(-0.3, 1))+
    scale_fill_brewer()+
    theme(legend.position=c(0.80,0.85), axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_blank())+
    xlab("Annotation")+
    ylab("Parameter estimate")

ggsave("pLI_contraint.pdf", plot = ga, path = plot_dir, width = 10, height = 5)

####################################################
## pleiotropy based on cross tissue gene expression
####################################################
result_dir <- "~/Projects/extraINSIGHT/results/grch38/gnomad_v3.0/constraint/expression_pleiotropy"
directory_list <- dir(result_dir, full.names = TRUE)
#plot_dir <- file.path(result_dir, "plots")

directory_list <- directory_list[!grepl("log|plots", directory_list)]

out_abs <- extract_conservation_estimates(directory_list, FALSE)
label_matrix <- str_match(out_abs$label, "(.+)_([a-z]+)$")
out_abs[, region := label_matrix[,3]]
out_abs[, expression_class:= label_matrix[,2]]
express_lvls <- c("broadly_expressed", "intermediate", "tissue_specific", "unexpressed")
out_abs[, expression_class := factor(expression_class, levels = express_lvls)]

ord <- out_abs[name == "strong_selection"][order(absolute, decreasing = TRUE)]$label
dir.create(plot_dir, FALSE, TRUE)

## Absolute selection
pos <- position_dodge(width=0.9)
ga <- ggplot(out_abs, aes(x = expression_class, y = absolute, fill = name))+
    facet_wrap(.~region, nrow = 1)+
    geom_bar(stat = "identity", position = pos)+
    geom_errorbar(aes(ymin = absolute - 1.96*se, ymax = absolute + 1.96*se), width = 0.2, position = pos)+
    theme_cowplot()+
    scale_y_continuous(breaks = pretty(seq(-0.3,1,0.1)), limits = c(-0.3, 1))+
    scale_fill_brewer()+
    theme(legend.position=c(0.80,0.85), axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_blank())+
    xlab("Annotation")+
    ylab("Parameter estimate")

ggsave("expression_pleiotropy_contraint.pdf", plot = ga, path = plot_dir, width = 10, height = 5)

####################################################
## Top level reactome plots
###################################################
result_dir <- "~/Projects/extraINSIGHT/results/grch38/gnomad_v3.0/constraint/reactome"
directory_list <- dir(result_dir, full.names = TRUE)
#plot_dir <- file.path(result_dir, "plots")

directory_list <- directory_list[!grepl("log|plots", directory_list)]

term_mapping <- fread("~/Projects/extraINSIGHT/data/grch38/annotation_bed/reactome/id_pathway_mapping.txt.gz",
                      key = "reactome_id")

out_abs <- extract_conservation_estimates(directory_list, FALSE)
label_matrix <- str_match(out_abs$label, "(.+)_([a-z]+)$")
out_abs[, region := label_matrix[,3]]
out_abs[, pathway:= term_mapping[toupper(label_matrix[,2])]$pathway]
ord_pathway <- out_abs[name == "strong_selection" & region == "cds"][order(absolute, decreasing = TRUE)]$pathway
out_abs[, pathway:= factor(pathway, levels = ord_pathway)]
ord <- out_abs[name == "strong_selection"][order(absolute, decreasing = TRUE)]$label
dir.create(plot_dir, FALSE, TRUE)

## Absolute selection
pos <- position_dodge(width=0.9)
ga <- ggplot(out_abs[region == "cds"], aes(x = pathway, y = absolute, fill = name))+
    facet_wrap(.~region, nrow = 1)+
    geom_bar(stat = "identity", position = pos)+
    geom_errorbar(aes(ymin = absolute - 1.96*se, ymax = absolute + 1.96*se), width = 0.2, position = pos)+
    theme_cowplot()+
    scale_y_continuous(breaks = pretty(seq(-0.3,1,0.1)), limits = c(-0.3, 1))+
    scale_fill_brewer()+
    theme(legend.position=c(0.80,0.85), axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_blank())+
    xlab("Annotation")+
    ylab("Parameter estimate")

ggsave("reactome_contraint.pdf", plot = ga, path = plot_dir, width = 9, height = 6)

####################################################
## MicroRNAs
###################################################
result_dir <- "~/Projects/extraINSIGHT/results/grch38/gnomad_v3.0/constraint/mirna"
directory_list <- dir(result_dir, full.names = TRUE)
#plot_dir <- file.path(result_dir, "plots")

directory_list <- directory_list[!grepl("log|plots", directory_list)]

out_abs <- extract_conservation_estimates(directory_list, FALSE)
label_matrix <- str_match(out_abs$label, "([a-z]+)_([a-z]+)_seed_hg19")
out_abs[, site := label_matrix[,3]]
out_abs[, family:= label_matrix[,2]]
#ord_pathway <- out_abs[name == "strong_selection" & region == "cds"][order(absolute, decreasing = TRUE)]$
ord <- out_abs[name == "strong_selection"][order(absolute, decreasing = TRUE)]$label
dir.create(plot_dir, FALSE, TRUE)

## Absolute selection
pos <- position_dodge(width=0.9)
ga <- ggplot(out_abs[family != "otherfam"], aes(x = site, y = absolute, fill = name))+
    facet_wrap(.~family, nrow = 1)+
    geom_bar(stat = "identity", position = pos)+
    geom_errorbar(aes(ymin = absolute - 1.96*se, ymax = absolute + 1.96*se), width = 0.2, position = pos)+
    theme_cowplot()+
    scale_y_continuous(breaks = pretty(seq(-0.3,1,0.1)))+
    scale_fill_brewer()+
    theme(legend.position=c(0.80,0.85), axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_blank())+
    xlab("Annotation")+
    ylab("Parameter estimate")

ggsave("microRNA_contraint.pdf", plot = ga, path = plot_dir, width = 9, height = 6)

####################################################
## MicroRNAs in 3'-UTRs ------- WIP
###################################################
result_dir <- "~/Projects/extraINSIGHT/results/grch38/gnomad_v3.0/constraint/three_utr_decomposition/"
directory_list <- dir(result_dir, full.names = TRUE)
#plot_dir <- file.path(result_dir, "plots")

directory_list <- directory_list[!grepl("log|plots", directory_list)]

out_abs <- extract_conservation_estimates(directory_list, FALSE)
label_matrix <- str_match(out_abs$label, "([a-z]+)_([a-z]+)_seed_hg19")
out_abs[, site := label_matrix[,3]]
out_abs[, family:= label_matrix[,2]]
#ord_pathway <- out_abs[name == "strong_selection" & region == "cds"][order(absolute, decreasing = TRUE)]$
ord <- out_abs[name == "strong_selection"][order(absolute, decreasing = TRUE)]$label
dir.create(plot_dir, FALSE, TRUE)

## Absolute selection
pos <- position_dodge(width=0.9)
ga <- ggplot(out_abs[family != "otherfam"], aes(x = site, y = absolute, fill = name))+
    facet_wrap(.~family, nrow = 1)+
    geom_bar(stat = "identity", position = pos)+
    geom_errorbar(aes(ymin = absolute - 1.96*se, ymax = absolute + 1.96*se), width = 0.2, position = pos)+
    theme_cowplot()+
    scale_y_continuous(breaks = pretty(seq(-0.3,1,0.1)))+
    scale_fill_brewer()+
    theme(legend.position=c(0.80,0.85), axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_blank())+
    xlab("Annotation")+
    ylab("Parameter estimate")

ggsave("microRNA_contraint.pdf", plot = ga, path = plot_dir, width = 9, height = 6)


####################################################
## MicroRNAs in 3'-UTRs
###################################################
result_dir <- "~/Projects/extraINSIGHT/results/grch38/gnomad_v3.0/constraint/noncoding_rna"
directory_list <- dir(result_dir, full.names = TRUE)

directory_list <- directory_list[!grepl("log|plots", directory_list)]

out_abs <- extract_conservation_estimates(directory_list, FALSE)
#ord_pathway <- out_abs[name == "strong_selection" & region == "cds"][order(absolute, decreasing = TRUE)]$
ord <- out_abs[name == "strong_selection"][order(absolute, decreasing = TRUE)]$label
dir.create(plot_dir, FALSE, TRUE)

## Absolute selection
pos <- position_dodge(width=0.9)
ga <- ggplot(out_abs, aes(x = label, y = absolute, fill = name))+
    geom_bar(stat = "identity", position = pos)+
    geom_errorbar(aes(ymin = absolute - 1.96*se, ymax = absolute + 1.96*se), width = 0.2, position = pos)+
    theme_cowplot()+
    scale_y_continuous(breaks = pretty(seq(-0.3,1,0.1)))+
    scale_fill_brewer()+
    theme(legend.position=c(0.80,0.85), axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_blank())+
    xlab("Annotation")+
    ylab("Parameter estimate")

ggsave("noncoding_rna_contraint.pdf", plot = ga, path = plot_dir, width = 9, height = 6)
