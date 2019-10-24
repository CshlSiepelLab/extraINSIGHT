## Choose gnomad version for both vcf and coverage as sometimes multiple vcf versions correspong to the same
## coverage version (e.g. 2.1.1 and 2.1)
gnomad_version_vcf=3.0
gnomad_version_coverage=3.0
genome_name=grch38

if ! gsutil_loc="$(type -p gsutil)" || [[ -z $gsutil_loc ]]; then
    echo "gsutil must be installed"
else
    ## Set-up the data directories
    gnomad_vcf_data_dir=/local/storage/data/extraINSIGHT/no-backup/data/${genome_name}/gnomad/"v"${gnomad_version_vcf}
    gnomad_coverage_data_dir=/local/storage/data/extraINSIGHT/no-backup/data/${genome_name}/gnomad/"v"${gnomad_version_coverage}

    ## Create data directory
    mkdir -p ${gnomad_vcf_data_dir}/vcf/genomes
    mkdir -p ${gnomad_coverage_data_dir}/coverage/genomes
    
    ## Download the VCF data using gsutil
    nohup gsutil cp gs://gnomad-public/release/${gnomad_version_vcf}/vcf/genomes/*.sites.*.vcf.bgz ${gnomad_vcf_data_dir}/vcf/genomes &> download_vcf_"v"${gnomad_version_vcf}.log &
    
    ## Download the VCF data using gsutil
    nohup gsutil cp gs://gnomad-public/release/${gnomad_version_coverage}/coverage/genomes/gnomad.genomes.*.coverage.summary.tsv.bgz ${gnomad_coverage_data_dir}/coverage/genomes &> download_coverage_"v"${gnomad_version_coverage}.log &
fi
