#!/usr/local/bin/perl
use strict;
use warnings;

my $context_rate_table = shift;

my %table;

open IN, $context_rate_table or die "$!\n";
while (<IN>){
    chomp;
    my ($mut, $rate) = split /\t/;
    $table{$mut} = $rate;
}
close IN;

while (<>){
    chomp;
    my @item = split /\t/;
    my $mut = $item[3];
    my $core_3mer = substr $mut, 2, 3;
    if (exists $table{$mut}){
        # filter out CpG
        if ($core_3mer !~ m/CG/){
            my $rate = $table{$mut};
            my $logit = log($rate / (1. - $rate));
            print (join("\t", @item, $logit), "\n");
        }
    }
    else{
        die "$mut\n";
    }
}
