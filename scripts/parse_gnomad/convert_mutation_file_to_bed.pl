#!/usr/local/bin/perl
use strict;
use warnings;

while (<>){
    chomp;
    my ($id, @info) = split /\s/;
    my ($chr, $pos, $ref, $alt) = split /[\-]/, $id;
    $chr =~ s/!//g;
    $chr = "chr$chr";
    $pos =~ s/^0*//;
    print (join("\t", $chr, $pos - 1, $pos, $ref, $alt, @info), "\n");
}
