#!/usr/local/bin/perl
use strict;
use warnings;

my ($bed_file, $mutation_file, $neutral_file, $global_para_file, $flanking_size, $min_coverage, $freq_cutoff, $min_n_mutation, $context_frequency_table, $covariate_file) = @ARGV;

# Iterate over each window in the chunked bed file
open IN, $bed_file or die "$!\n";
while (<IN>){
    chomp;
    # Store the co-ordinates as variables
    my ($chr, $bin_start, $bin_end) = split /\t/;

    # Add flanking regions to co-ordinates
    my $start = $bin_start - $flanking_size;
    $start = 0 if $start < 0;
    my $end = $bin_end + $flanking_size;

    # tabix queries assume a 1 based coordinate system (when not feeding in an input file with a *.bed or *.bed.gz suffix) 
    # and our other co-ordinates are 0-based so we need to correct for that here so that the tabix queries are for the 
    # correct co-ordinates. Outside the bed case tabix also assumes that the intervals are inclusive so there is no need
    # to adjust the end co-ordinates.
    $start += 1; 

    # Construct a commans that extracts the mutations in the window +/- flanking regions and filters based on coverage and
    # frequency. The remaining mutations are then merged with the genomic covariates ans subsetted again to contain only
    # mutations in neutral regions. This subcommand returns a file with 13 columns:
    ## 1)chrom
    ## 2)start
    ## 3)end
    ## 4)context-alternate_allele
    ## 5)mutation_flag (0 = no mutation, 1 = mutation at that particular site)
    ## 6)average sequencing coverage
    ## 7)log(mutation probability) based on context from mutation_frequency_table.txt
    ## 8)GC_content
    ## 9)is_cpg_island
    ## 10)number of bases that this site overlaps with the the covariate bed file by
    ## 11)chrom (of overlapping neutral region, or "." if not in neutral region) 
    ## 12)start (of overlapping neutral region, or "-1" if not in neutral region)
    ## 13)end (of overlapping neutral region, or "-1" if not in neutral region)
    
    my $cmd = qq[tabix $mutation_file $chr:$start-$end | perl filter_mutation.pl $min_coverage $freq_cutoff | perl add_context_mutation_rate.pl $context_frequency_table | bedtools intersect -a - -b $covariate_file -sorted -wo | cut -f1-7,11- | intersectBed -loj -a - -b $neutral_file -sorted];

    # An R command that passes the filtered data, global parameter file, and the minimum number of muations that must be present in the region to estimate a local scaling factor to the R-script
    # It then filters the output for th correct region and prints it
    my $R_cmd = $cmd . "| Rscript --no-restore local_mutation_regression.R $global_para_file $min_n_mutation /dev/fd/0" . qq[ | grep -v NULL | awk '{if(\$1 == "$chr" && \$2 >= $bin_start && \$3 <= $bin_end) print \$0}'];
    #print $R_cmd
    system $R_cmd;
    last;
}

close IN;
