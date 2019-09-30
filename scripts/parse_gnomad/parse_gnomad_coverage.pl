#!/usr/local/bin/perl
use strict;
use warnings;

# input summary coverage file
my $file = shift;
chomp($file);

# Read in and unzip coverage file, removing header
open IN, qq{gunzip -c $file | grep -v "#" | grep -v "pos" | } or die qq{Failed to open "$file" for reading: $!};
while (<IN>){
    chomp;
    my @items = split /\t/;
    my $chr = $items[0];
    my $pos = $items[1];
    # Get the mean coverage
    my $coverage = $items[2];
    # Return the position as 0-indexed position in bed format
    print (join ("\t", "chr".$chr, $pos - 1, $pos, $coverage), "\n");
    # last;
}
close IN;

