#!/usr/local/bin/perl
use strict;
use warnings;

my $min_coverage = shift;
my $freq_cutoff = shift;

while (<>){
    chomp;
    my @items = split /\t/;
    my $alt = $items[4];
    my $context = $items[5];
    my $allele_freq = $items[6];
    my $coverage = $items[7];
    my $mutation_flag = 0;

    if ($allele_freq > 0 and $allele_freq <= $freq_cutoff){
        $mutation_flag = 1;
    }
    $items[6] = $mutation_flag;
    if ($coverage >= $min_coverage){
        # print (join("\t", @items), "\n");
        print (join("\t", @items[0 .. 2], "$items[5]-$items[4]", @items[6 .. $#items]), "\n");
    }
}
