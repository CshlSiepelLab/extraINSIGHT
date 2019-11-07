#!/usr/local/bin/perl
use strict;
use warnings;

# Create a hashtable to hold mutation frequency counts
my %table;
while (<>){
    chomp;
    my @item = split /\t/;
    # Get context-alt_allele string
    my $mut = $item[3];
    # Get frequency of mutation
    my $freq = $item[4];
    # Increment value of $table{context-alt_allele}{total} by 1. This is the total number of sites
    # for which that context-alt_allele combination is possible. Note: barewords used as hashkeys 
    # are treated as strings
    $table{$mut}{total} ++;
    # If there is no 'mut' key defined to count the actual occurance rate of a mutation, create it.
    if (!defined $table{$mut}{mut}){
        $table{$mut}{mut} = 0;
    }
    # If the mutation occurs increment the mutation occurance counter by 1
    if ($freq > 0){
        $table{$mut}{mut} ++;
    }
    # print "$mut\t$freq\n";
}

# Sort the output and output the fraction of potential sites at which a mutation was observed
for (sort keys %table){
    my $f = $table{$_}{mut} / $table{$_}{total};
    print "$_\t$f\n";
}
