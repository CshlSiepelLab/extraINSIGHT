out_dir=../../data/grch37/pli_v2.1.1
mkdir -p ${out_dir}

wget -b -O ${out_dir}/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz
wget -b -O ${out_dir}/gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz
