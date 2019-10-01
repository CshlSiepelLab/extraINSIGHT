#!/usr/local/bin/perl
use strict;
use warnings;

# This script takes in a file of ranges of the form chr:X-Y that are 7-mers and generates all possible
# SNPs at the focal site as well as checking to make sure that the 7-mer is only A,C, G, or T.
my @alphabet = qw/A C G T/;
while (<>){
    chomp;
    my ($id, $seq) = split /\s+/;
    $seq = uc $seq;
    my ($chr, $start, $end)  = split /[\:\-]/, $id;
    # Get position of focal site from end of 7-mer
    my $pos = $end - 3;
    $chr =~ s/^chr//;
    # check valid context (7-mer ACGT) else skip rest of loop
    next if $seq !~ /^[AGCT]{7}$/;
    my $ref = substr $seq, 3, 1;
    # print("$chr\t$start\t$end\t$seq\t$ref\n");

    for my $allele (@alphabet){
        next if $allele eq $ref;
        print (join("-", $chr . ("!" x (3 - length $chr)), sprintf("%011d", $pos), $ref, $allele), "\t", $seq, "\n");
    }
}
