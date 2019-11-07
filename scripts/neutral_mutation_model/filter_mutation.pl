#!/usr/local/bin/perl
use strict;
use warnings;

my $min_coverage = shift;
my $freq_cutoff = shift;

while (<>){
    chomp;
    my @items = split /\t/;
    # The alternate allele
    my $alt = $items[4];
    # The 7-mer mutation context
    my $context = $items[5];
    # The frequency of the alternative allele
    my $allele_freq = $items[6];
    # The mean sequencing depth of the site
    my $coverage = $items[7];
    # Default binary mutation flag value set to 0
    my $mutation_flag = 0;

    # If there a rare mutation in the site set the mutation flag to 1
    if ($allele_freq > 0 and $allele_freq <= $freq_cutoff){
        $mutation_flag = 1;
    }
    
    # Replace the allele frequency with the binary mutation flag
    $items[6] = $mutation_flag;
    
    # Print out the items array with the following columns
    # chrom start end context_mutation mutation_flag mean_coverage
    if ($coverage >= $min_coverage){
        # print (join("\t", @items), "\n");
        print (join("\t", @items[0 .. 2], "$items[5]-$items[4]", @items[6 .. $#items]), "\n");
    }
}
