#!/usr/local/bin/perl
use strict;
use warnings;

my %table;
while (<>){
    chomp;
    my @item = split /\t/;
    my $mut = $item[3];
    my $freq = $item[4];
    $table{$mut}{total} ++;
    if (!defined $table{$mut}{mut}){
        $table{$mut}{mut} = 0;
    }
    if ($freq > 0){
        $table{$mut}{mut} ++;
    }
    # print "$mut\t$freq\n";
}

for (sort keys %table){
    my $f = $table{$_}{mut} / $table{$_}{total};
    print "$_\t$f\n";
}
