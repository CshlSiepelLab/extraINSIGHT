library(ggplot2)
library(cowplot)
library(stringr)
source("extract_result_lib.R")
source("ratio_uncertainty_propagation.R")

plot_dir <- file.path("../../../results/grch38/gnomad_v3.0/constraint", "ratio_plots")
dir.create(plot_dir, F, T)

## Allows for 1-3 annotation cols, if more than 1, the rest are faceted on

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
gw_cds_ei <- fread("../../../results/genomewide/ExtRaINSIGHT/coding/strong_selection_estimate.txt")
gw_strong_sel <- gw_cds_ei$val[1]
gw_cds_i <-
    fread("grep -A 10 PARAMETERS ../../../results/genomewide/INSIGHT/coding/coding.model | grep Rho")
gw_rho <- gw_cds_i$V2 

gw_ratio <- gw_strong_sel / gw_rho
gw_tab <- data.table(name = factor(c("rho", "strong_selection", "ratio"),
                         levels = c("rho", "strong_selection", "ratio")) ,
                     absolute = c(gw_rho, gw_strong_sel, gw_ratio))

#############################################
## Top 5% of expressed genes in each tissue
#############################################
result_dir <- "../../../results/grch38/gnomad_v3.0/constraint/top_5_expressed"
directory_list <- dir(result_dir, full.names = TRUE)
#plot_dir <- file.path(result_dir, "plots")

directory_list <- directory_list[!grepl("log|plots", directory_list)]

out_abs <- extract_conservation_estimates(directory_list, FALSE)
out_abs[,label:=gsub("_cds", "", label)]
r_tab <- ratio_table(sel_tab = out_abs)
ord <- r_tab[name == "ratio"][order(absolute, decreasing = TRUE)]$label

ga <- ggplot(data = r_tab, aes(x=factor(label, levels = ord), y = absolute, color = name))+
         geom_point(stat = "identity")+ 
         geom_errorbar(aes(ymin = absolute - 1.96 * se, ymax = absolute + 1.96 * se), width = 0.2)+
         geom_hline(data = gw_tab, aes(yintercept = absolute, color = name), linetype = 2)+
         theme_cowplot()+
         ylim(0,1)+
         scale_color_brewer(palette = "Dark2", labels = expression(rho, sel[strong], sel[strong]/rho))+
         theme(legend.position=c(0.65,0.85), axis.text.x = element_text(angle = 45, hjust = 1),
               legend.title = element_blank())+
         xlab("Annotation")+
         ylab("Estimate")

ggsave("top_five_percent_tissue_expressed_ratio.pdf", plot = ga, path = plot_dir, width = 17, height = 6)


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

r_tab <- ratio_table(sel_tab = out_abs)
ord <- r_tab[name == "ratio"][order(absolute, decreasing = TRUE)]$label

ga <- ggplot(data = r_tab[region == "cds"], aes(x=pLI, y = absolute, color = name))+
         geom_point(stat = "identity")+ 
         geom_errorbar(aes(ymin = absolute - 1.96 * se, ymax = absolute + 1.96 * se), width = 0.2)+
         geom_hline(data = gw_tab, aes(yintercept = absolute, color = name), linetype = 2)+
         theme_cowplot()+
         ylim(0,1)+
         scale_color_brewer(palette = "Dark2", labels = expression(rho, sel[strong], sel[strong]/rho))+
         theme(legend.position=c(0.65,0.85), axis.text.x = element_text(angle = 45, hjust = 1),
               legend.title = element_blank())+
         xlab("Annotation")+
         ylab("Estimate")

ggsave("pLI_ratio.pdf", plot = ga, path = plot_dir, width = 4, height = 3)

#############################################
## Check against subsampled cds
#############################################
result_dir <- "../../../results/grch38/gnomad_v3.0/constraint/gene_annotations"
directory_list <- dir(result_dir, full.names = TRUE)
directory_list <- directory_list[!grepl("log|plots", directory_list)]

out_abs <- extract_conservation_estimates(directory_list, FALSE)
#out_resid <- extract_conservation_estimates(directory_list, TRUE)
out_abs[,label := gsub("_gencode33","", label)]
out_abs <- out_abs[label == "cds"]

r_tab <- ratio_table(sel_tab = out_abs)
ord <- r_tab[name == "ratio"][order(absolute, decreasing = TRUE)]$label

ga <- ggplot(data = r_tab, aes(x=factor(label, levels = ord), y = absolute, color = name))+
         geom_point(stat = "identity")+ 
         geom_errorbar(aes(ymin = absolute - 1.96 * se, ymax = absolute + 1.96 * se), width = 0.2)+
         geom_hline(data = gw_tab, aes(yintercept = absolute, color = name), linetype = 2)+
         theme_cowplot()+
         ylim(0,1)+
         scale_color_brewer(palette = "Dark2", labels = expression(rho, sel[strong], sel[strong]/rho))+
         theme(legend.position=c(0.65,0.85), axis.text.x = element_text(angle = 45, hjust = 1),
               legend.title = element_blank())+
         xlab("Annotation")+
         ylab("Estimate")

ggsave("cds_subsampled_ratio.pdf", plot = ga, path = plot_dir, width = 4, height = 3)


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

r_tab <- ratio_table(sel_tab = out_abs[region == "cds"])
ord <- r_tab[name == "ratio"][order(absolute, decreasing = TRUE)]$pathway

ga <- ggplot(data = r_tab, aes(x=factor(pathway, levels = ord), y = absolute, color = name))+
         geom_point(stat = "identity")+ 
         geom_errorbar(aes(ymin = absolute - 1.96 * se, ymax = absolute + 1.96 * se), width = 0.2)+
         geom_hline(data = gw_tab, aes(yintercept = absolute, color = name), linetype = 2)+
         theme_cowplot()+
         ylim(0,1)+
         scale_color_brewer(palette = "Dark2", labels = expression(rho, sel[strong], sel[strong]/rho))+
         theme(legend.position=c(0.65,0.85), axis.text.x = element_text(angle = 45, hjust = 1),
               legend.title = element_blank())+
         xlab("Annotation")+
         ylab("Estimate")

ggsave("reactome_cds_ratio.pdf", plot = ga, path = plot_dir, width = 9, height = 6)
