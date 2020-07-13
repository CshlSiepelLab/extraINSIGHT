library(ggplot2)
library(cowplot)
library(stringr)
source("extract_result_lib.R")
source("ratio_uncertainty_propagation.R")

plot_dir <- file.path("../../../results/grch38/gnomad_v3.0/constraint", "ratio_plots")
dir.create(plot_dir, F, T)

ratio_table <- function(sel_tab, annotation_cols = "label", sel_type_col = "name",
                        estimate_col = "absolute", se_col = "se"){
    ## Copy table to avoid modify by reference
    sel_tab <- copy(sel_tab)
    ## Create joint key
    sel_tab[,uid := do.call(paste0,.SD), .SDcols=annotation_cols]
    setkeyv(sel_tab, c("uid", sel_type_col))

    ## Extract paired lists of estimated values
    rho <- unlist(sel_tab[.(unique(uid), "rho")][, ..estimate_col])
    strong_sel <- unlist(sel_tab[.(unique(uid), "strong_selection")][, ..estimate_col])
    se_rho <- unlist(sel_tab[.(unique(uid), "rho")][, ..se_col])
    se_strong_sel <- unlist(sel_tab[.(unique(uid), "strong_selection")][, ..se_col])

    ## Create a table for the ratios
    ratios <- ratio_estimate(strong_sel, se_strong_sel, rho, se_rho)
    ratio_tab <- sel_tab[.(unique(uid), "rho")]
    ratio_tab[[sel_type_col]] <- "ratio"
    ratio_tab[[estimate_col]] <- ratios[[1]]
    ratio_tab[[se_col]] <- ratios[[2]]

    ## Rbind ratio table into selection table
    fin_tab <- rbind(sel_tab, ratio_tab)
        
    ## Remove added column and precompute some values
    fin_tab[,uid:=NULL]
    fin_tab[[sel_type_col]] <- factor(fin_tab[[sel_type_col]],
                                      levels = c("rho","strong_selection", "ratio"))
    return(fin_tab)
}

####
## Calculate genomewide coding EI/I ratio
####
gw_nonc_ei <- fread("../../../results/genomewide/ExtRaINSIGHT/noncoding/strong_selection_estimate.txt")
gw_strong_sel <- gw_nonc_ei$val[1]
gw_nonc_i <-
    fread("grep -A 10 PARAMETERS ../../../results/genomewide/INSIGHT/noncoding/noncoding.model | grep Rho")
gw_rho <- gw_nonc_i$V2 

gw_ratio <- gw_strong_sel / gw_rho
gw_tab <- data.table(name = factor(c("rho", "strong_selection", "ratio"),
                         levels = c("rho", "strong_selection", "ratio")) ,
                     absolute = c(gw_rho, gw_strong_sel, gw_ratio))

####################################################
## MicroRNAs - target regions
###################################################
result_dir <- "~/Projects/extraINSIGHT/results/grch38/gnomad_v3.0/constraint/mirna"
directory_list <- dir(result_dir, full.names = TRUE)
#plot_dir <- file.path(result_dir, "plots")

directory_list <- directory_list[!grepl("log|plots", directory_list)]

out_abs <- extract_conservation_estimates(directory_list, FALSE)
label_matrix <- str_match(out_abs$label, "([a-z]+)_([a-z]+)_seed_hg19")
out_abs[, site := label_matrix[,3]]
out_abs[, family:= label_matrix[,2]]

r_tab <- ratio_table(sel_tab = out_abs)
r_tab <- r_tab[family != "otherfam"]
ord <- r_tab[name == "ratio"][, mean(absolute), by = "site"]$site

pos <- position_dodge(0.2)
ga <- ggplot(data = r_tab, aes(x=factor(site, levels = ord), y = absolute, color = name))+
    facet_wrap(~family, nrow = 1)+
    geom_point(stat = "identity", position = pos)+ 
    geom_errorbar(aes(ymin = pmax(0, absolute - 1.96 * se), ymax = pmin(1, absolute + 1.96 * se)), width = 0.2,
                  position =pos)+
    geom_hline(data = gw_tab, aes(yintercept = absolute, color = name), linetype = 2)+
    theme_cowplot()+
    ylim(0,1)+
    scale_color_brewer(palette = "Dark2", labels = expression(rho, sel[strong], sel[strong]/rho))+
    theme(legend.position=c(0.85,0.85), axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_blank())+
    xlab("Annotation")+
    ylab("Estimate")

ggsave("microRNA_target_seeds_ratio.pdf", plot = ga, path = plot_dir, width = 9, height = 4)

####################################################
## MicroRNAs - from mirbase and targetscan annotations
###################################################
result_dir <- "~/Projects/extraINSIGHT/results/grch38/gnomad_v3.0/constraint/micro_rnas"
directory_list <- dir(result_dir, full.names = TRUE)
#plot_dir <- file.path(result_dir, "plots")

directory_list <- directory_list[!grepl("log|plots", directory_list)]

out_abs <- extract_conservation_estimates(directory_list, FALSE)

r_tab <- ratio_table(sel_tab = out_abs)
ord <- r_tab[name == "ratio"][, mean(absolute), by = "label"][order(V1)]$label

pos <- position_dodge(0.2)
ga <- ggplot(data = r_tab, aes(x=factor(label, levels = ord), y = absolute, color = name))+
    geom_point(stat = "identity", position = pos)+ 
    geom_errorbar(aes(ymin = pmax(0, absolute - 1.96 * se), ymax = pmin(1, absolute + 1.96 * se)), width = 0.2,
                  position =pos)+
    geom_hline(data = gw_tab, aes(yintercept = absolute, color = name), linetype = 2)+
    theme_cowplot()+
    ylim(0,1)+
    scale_color_brewer(palette = "Dark2", labels = expression(rho, sel[strong], sel[strong]/rho))+
    theme(legend.position=c(0.85,0.85), axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_blank())+
    xlab("Annotation")+
    ylab("Estimate")

ggsave("microRNA_gene_ratio.pdf", plot = ga, path = plot_dir, width = 9, height = 4)

####################################################
## UCNEs
###################################################
result_dir <- "~/Projects/extraINSIGHT/results/grch38/gnomad_v3.0/constraint/ucne"
directory_list <- dir(result_dir, full.names = TRUE)
#plot_dir <- file.path(result_dir, "plots")

directory_list <- directory_list[!grepl("log|plots", directory_list)]

out_abs <- extract_conservation_estimates(directory_list, FALSE)
complete_analysis <- names(table(out_abs$label)[table(out_abs$label) == 2])
out_abs <- out_abs[label %in% complete_analysis]

r_tab <- ratio_table(sel_tab = out_abs)
ord <- r_tab[name == "ratio"][, mean(absolute), by = "label"][order(V1)]$label

pos <- position_dodge(0.2)
ga <- ggplot(data = r_tab, aes(x=factor(label, levels = ord), y = absolute, color = name))+
    geom_point(stat = "identity", position = pos)+ 
    geom_errorbar(aes(ymin = pmax(0, absolute - 1.96 * se), ymax = pmin(1, absolute + 1.96 * se)), width = 0.2,
                  position =pos)+
    geom_hline(data = gw_tab, aes(yintercept = absolute, color = name), linetype = 2)+
    theme_cowplot()+
    scale_color_brewer(palette = "Dark2", labels = expression(rho, sel[strong], sel[strong]/rho))+
    theme(legend.position=c(0.7,0.88), axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_blank())+
    xlab("Annotation")+
    ylab("Estimate")

ggsave("ucne_ratio.pdf", plot = ga, path = plot_dir, width = 5, height = 5)


####################################################
## RBPs
###################################################
result_dir <- "~/Projects/extraINSIGHT/results/grch38/gnomad_v3.0/constraint/rna_binding_proteins"
directory_list <- dir(result_dir, full.names = TRUE)
#plot_dir <- file.path(result_dir, "plots")

directory_list <- directory_list[!grepl("log|plots", directory_list)]

out_abs <- extract_conservation_estimates(directory_list, FALSE)
complete_analysis <- names(table(out_abs$label)[table(out_abs$label) == 2])
out_abs <- out_abs[label %in% complete_analysis]

r_tab <- ratio_table(sel_tab = out_abs)
ord <- r_tab[name == "ratio"][, mean(absolute), by = "label"][order(V1)]$label

pos <- position_dodge(0)
ga <- ggplot(data = r_tab, aes(x=factor(label, levels = ord), y = absolute, color = name))+
    geom_point(stat = "identity", position = pos)+ 
    geom_errorbar(aes(ymin = pmax(0, absolute - 1.96 * se), ymax = pmin(1, absolute + 1.96 * se)), width = 0.2,
                  position =pos)+
    geom_hline(data = gw_tab, aes(yintercept = absolute, color = name), linetype = 2)+
    theme_cowplot()+
    scale_color_brewer(palette = "Dark2", labels = expression(rho, sel[strong], sel[strong]/rho))+
    theme(legend.position=c(0.85,0.85), axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_blank())+
    xlab("Annotation")+
    ylab("Estimate")

ggsave("RNABP_binding_sites_ratio.pdf", plot = ga, path = plot_dir, width = 18, height = 6)
