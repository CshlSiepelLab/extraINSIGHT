#!/usr/local/bin/perl
use strict;
use warnings;

# input vcf
my $file = shift @ARGV; 
# "/gnomad.exomes.r2.1.sites.vcf.bgz";
# print "# $key $file\n";
# next;
# filter out coding only variants
# get chromosome
# print "$file\n";

# Stream in unzipped file, removing all lines with '#' at the begining
open IN, qq{gunzip -c $file | grep -v "#" | } or die "$!\n";
while (<IN>){
    # Remove endline whitespace
    chomp;
    # Split line on tabs
    my @items = split /\t/;
    # Retrieve elements of line
    my $chr = $items[0]; # Chromosome
    $chr =~ s/chr//g;
    my $pos = $items[1]; # Position (1-indexed)
    my $ref = $items[3]; # Reference allele
    my $all_alt = $items[4]; # Alternate allele
    # remove all sites that are not A,G,G, or T in reference or does not pass quality checks
    next if $ref !~ m/^[AGCT]$/ or $items[6] ne "PASS";

    # Get variant frequency
    my @all_freqs;
    # Extract the string of allele frequencies and split into array
    if ($items[7] =~ m/AF=([^;]+)/){
        @all_freqs = split /,/, $1;
    }
    else{
        die "Allele information not available: $_\n";
    }
    
    my @alts = split /,/, $all_alt;
    # If the number of frequencies does not match the number of alternative alleles, stop
    if (scalar @alts != scalar @all_freqs){
        die "Unmatched allele information: $_\n";
    }

    # Iterate over alternative alleles
    @alts = sort {$a cmp $b} @alts;
    for (my $i = 0; $i < @alts; $i ++){
        my $allele = $alts[$i];
        my $freq = $all_freqs[$i];
        # filter out indels and variants with 0 allele frequency
        if ($allele =~ m/^[AGCT]$/ and $freq > 0.){
            # print (join("\t", "chr" . $chr, $pos - 1, $pos, $ref, $allele, $freq), "\n");
	    # Create output where first column is fixed with key that can be used for later joins including position,
	    # ref-alt alleles, and the second column is the allele frequency
            print (join("-", $chr . ("!" x (3 - length $chr)), sprintf("%011d", $pos), $ref, $allele), "\t", $freq, "\n");
        }
    }
    # last;
}
close IN;
