#!/usr/local/bin/perl
use strict;
use warnings;

my @alphabet = qw/A C G T/;
while (<>){
    chomp;
    my ($id, $seq) = split /\s+/;
    $seq = uc $seq;
    my ($chr, $start, $end)  = split /[\:\-]/, $id;
    my $pos = $end - 3;
    $chr =~ s/^chr//;
    # check valid context
    next if $seq !~ /^[AGCT]{7}$/;
    my $ref = substr $seq, 3, 1;
    # print "$chr\t$start\t$end\t$seq\t$ref\n";

    for my $allele (@alphabet){
        next if $allele eq $ref;
        print (join("-", $chr . ("!" x (3 - length $chr)), sprintf("%011d", $pos), $ref, $allele), "\t", $seq, "\n");
    }
}
