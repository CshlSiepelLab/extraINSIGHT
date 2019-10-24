#!/usr/local/bin/perl
use strict;
use warnings;

# input summary coverage file
my $file = shift;
chomp($file);
my $head;

open IN, qq{gunzip -c $file | grep -v "#" | head -1 | } or die qq{Failed to open "$file" for reading: $!};
while(<IN>){
    chomp;
    my @header = split /\t/;
    $head = $header[0];
}
close IN;
    
# Read in and unzip coverage file, removing header
open IN, qq{gunzip -c $file | grep -v "#" | grep -v pos | grep -v locus | } or die qq{Failed to open "$file" for reading: $!};
while (<IN>){
    chomp;
    my @items = split(/\t/,$_);
    my $chr;
    my $pos;
    if($head eq "chrom"){
	$chr = $items[0];
	$pos = $items[1];
    }
    elsif ($head eq "locus"){
	my @loci = split(':',$items[0]);
	$chr = $loci[0] =~ s/chr//;
	$pos = $loci[1];
    }
    # Get the mean coverage
    my $coverage = $items[2];    
    # Return the position as 0-indexed position in bed format
    print (join ("\t", "chr".$chr, $pos - 1, $pos, $coverage), "\n");
    # last;
}
close IN;
