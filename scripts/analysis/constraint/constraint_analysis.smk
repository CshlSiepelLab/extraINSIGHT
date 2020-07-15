import os.path
import glob

## Important directories
anno_dir_grch38="../../../data/grch38/annotation_bed/"
anno_dir_grch37="../../../data/grch37/annotation_bed/"
out_dir_grch38="../../../results/grch38/gnomad_v3.0/constraint/"

## Parameters
max_sites=1000000 # max number of sites to subsample for analysis
exclude_chr="X Y" # chromosomes to exclude, space delimited
max_proc=5 # max num. of analysis that will be run by batch script simultaineously
default_ei_genome="hg38"

def adaptive_proc(path, subdir, m):
    bed_files = len(glob.glob(os.path.join(path, subdir, "*.bed.*")))
    return(min(bed_files, m))

## Run all analysis
rule all:
    input:
        os.path.join(out_dir_grch38, "gene_annotations"),
        os.path.join(out_dir_grch38, "top_5_expressed"),
        os.path.join(out_dir_grch38, "expression_pleiotropy"),
        os.path.join(out_dir_grch38, "reactome"),
        os.path.join(out_dir_grch38, "noncoding_rna"),
        os.path.join(out_dir_grch38, "rna_binding_proteins"),
        os.path.join(out_dir_grch38, "pli"),
        os.path.join(out_dir_grch38, "ucne"),
        os.path.join(out_dir_grch38, "mirna"),
        os.path.join(out_dir_grch38, "three_utr_decomposition"),
        os.path.join(out_dir_grch38, "tissue_specific_expression"),
        os.path.join(out_dir_grch38, "tissue_group_exclusivity"),
        os.path.join(out_dir_grch38, "phastcons_elements")

####################
## Bed file on hg38
####################
rule gene_annotations:
    input:
        analysis_bed_dir=os.path.join(anno_dir_grch38, "gene_annotations")
    output:
        analysis_out_dir=directory(os.path.join(out_dir_grch38, "gene_annotations"))
    params:
        ei_genome=default_ei_genome,
        bed_genome="hg38"
    threads:
        adaptive_proc(anno_dir_grch38, "gene_annotations", max_proc)
    shell:
        """
        ./extrainsight_vs_insight_batch.R -b {input.analysis_bed_dir} \
        -o {output.analysis_out_dir} --max-sites {max_sites} -g {params.bed_genome} \
        --extrainsight-genome {params.ei_genome} --max-processes {threads} --exclude-chr {exclude_chr}
        """

rule top_five_percent_expressed_genes_tissue:
    input:
        analysis_bed_dir=os.path.join(anno_dir_grch38, "top_5_expressed")
    output:
        analysis_out_dir=directory(os.path.join(out_dir_grch38, "top_5_expressed"))
    params:
        ei_genome=default_ei_genome,
        bed_genome="hg38"
    threads:
        adaptive_proc(anno_dir_grch38, "top_5_expressed", max_proc)
    shell:
        """
        ./extrainsight_vs_insight_batch.R -b {input.analysis_bed_dir} \
        -o {output.analysis_out_dir} --max-sites {max_sites} -g {params.bed_genome} \
        --extrainsight-genome {params.ei_genome} --max-processes {threads} --exclude-chr {exclude_chr}
        """

rule expression_pleiotropy:
    input:
        analysis_bed_dir=os.path.join(anno_dir_grch38, "expression_pleiotropy")
    output:
        analysis_out_dir=directory(os.path.join(out_dir_grch38, "expression_pleiotropy"))
    params:
        ei_genome=default_ei_genome,
        bed_genome="hg38"
    threads:
        adaptive_proc(anno_dir_grch38, "expression_pleiotropy", max_proc)
    shell:
        """
        ./extrainsight_vs_insight_batch.R -b {input.analysis_bed_dir} \
        -o {output.analysis_out_dir} --max-sites {max_sites} -g {params.bed_genome} \
        --extrainsight-genome {params.ei_genome} --max-processes {threads} --exclude-chr {exclude_chr}
        """

rule reactome_cds:
    input:
        analysis_bed_dir=os.path.join(anno_dir_grch38, "reactome")
    output:
        analysis_out_dir=directory(os.path.join(out_dir_grch38, "reactome"))
    params:
        ei_genome=default_ei_genome,
        bed_genome="hg38"
    threads:
        adaptive_proc(anno_dir_grch38, "reactome", max_proc)
    shell:
        """
        ./extrainsight_vs_insight_batch.R -b {input.analysis_bed_dir} \
        -o {output.analysis_out_dir} --max-sites {max_sites} -g {params.bed_genome} \
        --extrainsight-genome {params.ei_genome} --max-processes {threads} --exclude-chr {exclude_chr}
        """

rule noncoding_rna:
    input:
        analysis_bed_dir=os.path.join(anno_dir_grch38, "noncoding_rna")
    output:
        analysis_out_dir=directory(os.path.join(out_dir_grch38, "noncoding_rna"))
    params:
        ei_genome=default_ei_genome,
        bed_genome="hg38"
    threads:
        adaptive_proc(anno_dir_grch38, "noncoding_rna", max_proc)
    shell:
        """
        ./extrainsight_vs_insight_batch.R -b {input.analysis_bed_dir} \
        -o {output.analysis_out_dir} --max-sites {max_sites} -g {params.bed_genome} \
        --extrainsight-genome {params.ei_genome} --max-processes {threads} --exclude-chr {exclude_chr}
        """

rule rna_binding_proteins:
    input:
        analysis_bed_dir=os.path.join(anno_dir_grch38, "rna_binding_proteins")
    output:
        analysis_out_dir=directory(os.path.join(out_dir_grch38, "rna_binding_proteins"))
    params:
        ei_genome=default_ei_genome,
        bed_genome="hg38"
    threads:
        adaptive_proc(anno_dir_grch38, "rna_binding_proteins", max_proc)
    shell:
        """
        ./extrainsight_vs_insight_batch.R -b {input.analysis_bed_dir} \
        -o {output.analysis_out_dir} --max-sites {max_sites} -g {params.bed_genome} \
        --extrainsight-genome {params.ei_genome} --max-processes {threads} --exclude-chr {exclude_chr}
        """

rule phastcons_100way:
    input:
        analysis_bed_dir=os.path.join(anno_dir_grch38, "phastcons_elements")
    output:
        analysis_out_dir=directory(os.path.join(out_dir_grch38, "phastcons_elements"))
    params:
        ei_genome=default_ei_genome,
        bed_genome="hg38"
    threads:
        adaptive_proc(anno_dir_grch38, "phastcons_elements", max_proc)
    shell:
        """
        ./extrainsight_vs_insight_batch.R -b {input.analysis_bed_dir} \
        -o {output.analysis_out_dir} --max-sites {max_sites} -g {params.bed_genome} \
        --extrainsight-genome {params.ei_genome} --max-processes {threads} --exclude-chr {exclude_chr}
        """

rule micro_rna_genes:
    input:
        analysis_bed_dir=os.path.join(anno_dir_grch38, "micro_rnas")
    output:
        analysis_out_dir=directory(os.path.join(out_dir_grch38, "micro_rnas"))
    params:
        ei_genome=default_ei_genome,
        bed_genome="hg38"
    threads:
        adaptive_proc(anno_dir_grch38, "micro_rnas", max_proc)
    shell:
        """
        ./extrainsight_vs_insight_batch.R -b {input.analysis_bed_dir} \
        -o {output.analysis_out_dir} --max-sites {max_sites} -g {params.bed_genome} \
        --extrainsight-genome {params.ei_genome} --max-processes {threads} --exclude-chr {exclude_chr}
        """
        
####################
## Bed file on hg19
####################
rule pli:
    input:
        analysis_bed_dir=os.path.join(anno_dir_grch37, "pli_annotations")
    output:
        analysis_out_dir=directory(os.path.join(out_dir_grch38, "pli"))
    params:
        ei_genome=default_ei_genome,
        bed_genome="hg19"
    threads:
        adaptive_proc(anno_dir_grch37, "pli_annotations", max_proc)
    shell:
        """
        ./extrainsight_vs_insight_batch.R -b {input.analysis_bed_dir} \
        -o {output.analysis_out_dir} --max-sites {max_sites} -g {params.bed_genome} \
        --extrainsight-genome {params.ei_genome} --max-processes {threads} --exclude-chr {exclude_chr}
        """

rule ucne:
    input:
        analysis_bed_dir=os.path.join(anno_dir_grch37, "ucne")
    output:
        analysis_out_dir=directory(os.path.join(out_dir_grch38, "ucne"))
    params:
        ei_genome=default_ei_genome,
        bed_genome="hg19"
    threads:
        adaptive_proc(anno_dir_grch37, "ucne", max_proc)
    shell:
        """
        ./extrainsight_vs_insight_batch.R -b {input.analysis_bed_dir} \
        -o {output.analysis_out_dir} --max-sites {max_sites} -g {params.bed_genome} \
        --extrainsight-genome {params.ei_genome} --max-processes {threads} --exclude-chr {exclude_chr}
        """

rule mirna:
    input:
        analysis_bed_dir=os.path.join(anno_dir_grch37, "mirna")
    output:
        analysis_out_dir=directory(os.path.join(out_dir_grch38, "mirna"))
    params:
        ei_genome=default_ei_genome,
        bed_genome="hg19"
    threads:
        adaptive_proc(anno_dir_grch37, "mirna", max_proc)
    shell:
        """
        ./extrainsight_vs_insight_batch.R -b {input.analysis_bed_dir} \
        -o {output.analysis_out_dir} --max-sites {max_sites} -g {params.bed_genome} \
        --extrainsight-genome {params.ei_genome} --max-processes {threads} --exclude-chr {exclude_chr}
        """

rule three_utr_mirna_decomposition:
    input:
        analysis_bed_dir=os.path.join(anno_dir_grch37, "three_utr_decomposition")
    output:
        analysis_out_dir=directory(os.path.join(out_dir_grch38, "three_utr_decomposition"))
    params:
        ei_genome=default_ei_genome,
        bed_genome="hg19"
    threads:
        adaptive_proc(anno_dir_grch37, "three_utr_decomposition", max_proc)
    shell:
        """
        ./extrainsight_vs_insight_batch.R -b {input.analysis_bed_dir} \
        -o {output.analysis_out_dir} --max-sites {max_sites} -g {params.bed_genome} \
        --extrainsight-genome {params.ei_genome} --max-processes {threads} --exclude-chr {exclude_chr}
        """

rule tissue_specific_expression:
    input:
        analysis_bed_dir=os.path.join(anno_dir_grch37, "tissue_specific_expression")
    output:
        analysis_out_dir=directory(os.path.join(out_dir_grch38, "tissue_specific_expression"))
    params:
        ei_genome=default_ei_genome,
        bed_genome="hg19"
    threads:
        adaptive_proc(anno_dir_grch37, "tissue_specific_expression", max_proc)
    shell:
        """
        ./extrainsight_vs_insight_batch.R -b {input.analysis_bed_dir} \
        -o {output.analysis_out_dir} --max-sites {max_sites} -g {params.bed_genome} \
        --extrainsight-genome {params.ei_genome} --max-processes {threads} --exclude-chr {exclude_chr}
        """

rule tissue_group_exclusivity:
    input:
        analysis_bed_dir=os.path.join(anno_dir_grch37, "tissue_group_exclusivity")
    output:
        analysis_out_dir=directory(os.path.join(out_dir_grch38, "tissue_group_exclusivity"))
    params:
        ei_genome=default_ei_genome,
        bed_genome="hg19"
    threads:
        adaptive_proc(anno_dir_grch37, "tissue_group_exclusivity", max_proc)
    shell:
        """
        ./extrainsight_vs_insight_batch.R -b {input.analysis_bed_dir} \
        -o {output.analysis_out_dir} --max-sites {max_sites} -g {params.bed_genome} \
        --extrainsight-genome {params.ei_genome} --max-processes {threads} --exclude-chr {exclude_chr}
        """

#########################################
## Genome-wide analysis - INSIGHT - hg19
## NOTE: These are not included in the "all" rule
#########################################

genomewide_results="../../../results/genomewide"
rule genomewide_coding_insight:
    input:
        bed=os.path.join(anno_dir_grch37, "coding_genomewide", "coding.bed.gz")
    output:
        unzipped_bed=temp(os.path.join(anno_dir_grch37, "coding_genomewide", "coding.bed")),
        out_dir=directory(os.path.join(genomewide_results, "INSIGHT", "coding"))
    shell:
        """
        gunzip -c {input.bed} | sort-bed - >  {output.unzipped_bed}
        ../../INSIGHT2/INSIGHT2.py -b {output.unzipped_bed} -o {output.out_dir}
        """

rule genomewide_noncoding_insight:
    input:
        bed=os.path.join(anno_dir_grch37, "noncoding_genomewide", "noncoding.bed.gz")
    output:
        unzipped_bed=temp(os.path.join(anno_dir_grch37, "noncoding_genomewide", "noncoding.bed")),
        out_dir=directory(os.path.join(genomewide_results, "INSIGHT", "noncoding"))
    shell:
        """
        gunzip -c {input.bed} | sort-bed - >  {output.unzipped_bed}
        ../../INSIGHT2/INSIGHT2.py -b {output.unzipped_bed} -o {output.out_dir}
        """

#########################################
## Genome-wide analysis - ExtRaINSIGHT - hg38
## NOTE: These are not included in the "all" rule
#########################################

genomewide_results="../../../results/genomewide"
ei_mutation_model_dir="../../../results/grch38/gnomad_v3.0/mutation_model"
rule genomewide_coding_extrainsight:
    input:
        bed=os.path.join(anno_dir_grch38, "coding_genomewide", "coding.bed.gz"),
        mutation_rates=os.path.join(ei_mutation_model_dir, "final_mutation_rates.bed.gz"),
        mutation_coverage=os.path.join(ei_mutation_model_dir, "final_mutation_site_coverage.bed.gz"),
    output:
        unzipped_bed=temp(os.path.join(anno_dir_grch38, "coding_genomewide", "coding.bed")),
        out_dir=directory(os.path.join(genomewide_results, "ExtRaINSIGHT", "coding"))
    shell:
        """
        gunzip -c {input.bed} | sort-bed - | bedops -i - <(zcat {input.mutation_coverage}) >  {output.unzipped_bed}
        ../../extraINSIGHT/ExtRaINSIGHT.py -b {output.unzipped_bed} -r {input.mutation_rates} -o {output.out_dir}
        """

## Have to subsample for this one to avoid memory limitations
rule genomewide_noncoding_extrainsight:
    input:
        bed=os.path.join(anno_dir_grch38, "noncoding_genomewide", "noncoding.bed.gz"),
        mutation_rates=os.path.join(ei_mutation_model_dir, "final_mutation_rates.bed.gz"),
        mutation_coverage=os.path.join(ei_mutation_model_dir, "final_mutation_site_coverage.bed.gz"),
    output:
        unzipped_bed=temp(os.path.join(anno_dir_grch38, "noncoding_genomewide", "noncoding.bed")),
        out_dir=directory(os.path.join(genomewide_results, "ExtRaINSIGHT", "noncoding"))
    params:
        keep = 0.005,
        seed = 2395473
    shell:
        """
        gunzip -c {input.bed} | sort-bed - | bedops -i - <(zcat {input.mutation_coverage}) | bedops -w 1 - |\
        awk 'BEGIN{{srand({params.seed})}}{{ if (rand() <= {params.keep}) print $0 }}'>  {output.unzipped_bed}
        ../../extraINSIGHT/ExtRaINSIGHT.py -T -b {output.unzipped_bed} -r {input.mutation_rates} -o {output.out_dir}
        """