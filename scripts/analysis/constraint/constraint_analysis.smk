import os.path

## Important directories
anno_dir_grch38="../../../data/grch38/annotation_bed/"
anno_dir_grch37="../../../data/grch37/annotation_bed/"
out_dir_grch38="../../../results/grch38/gnomad_v3.0/constraint/"

## Parameters
max_sites=1000000 # max number of sites to subsample for analysis
exclude_chr="X Y" # chromosomes to exclude, space delimited
max_proc=5 # max num. of analysis that will be run by batch script simultaineously
default_ei_genome="hg38"

## Run all analysis
rule all:
    input:
        os.path.join(anno_dir_grch38, "gene_annotations"),
        os.path.join(out_dir_grch38, "top_5_expressed"),
        os.path.join(out_dir_grch38, "expression_pleiotropy"),
        os.path.join(out_dir_grch38, "reactome"),
        os.path.join(out_dir_grch38, "noncoding_rna"),
        os.path.join(out_dir_grch38, "pli"),
        os.path.join(out_dir_grch38, "ucne"),
        os.path.join(out_dir_grch38, "mirna"),
        os.path.join(out_dir_grch38, "three_utr_decomposition")
        
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
        max_proc
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
        max_proc
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
        max_proc
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
        max_proc
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
        max_proc
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
        max_proc + 1
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
        max_proc
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
        max_proc
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
        max_proc
    shell:
        """
        ./extrainsight_vs_insight_batch.R -b {input.analysis_bed_dir} \
        -o {output.analysis_out_dir} --max-sites {max_sites} -g {params.bed_genome} \
        --extrainsight-genome {params.ei_genome} --max-processes {threads} --exclude-chr {exclude_chr}
        """
