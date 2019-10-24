#!/usr/local/bin/perl
use strict;
use warnings;

my $min_effective_size = shift;

while (<>){
    chomp;
    my ($chr, $start, $end, $seq) = split /[\:\-\s]+/;
    # convert to upper cases
    $seq = uc $seq;

    # count
    my $total_count = ($seq =~ tr/AGCT/AGCT/);
    my $GC = ($seq =~ tr/GC/GC/);

    # filter
    if ($total_count >= $min_effective_size){
        my $GC_content = $GC / $total_count;
        # print "$chr|$start|$end|$seq|$total_count|$GC|$GC_content\n";
        print "$chr\t$start\t$end\t$GC_content\n";
    }
}
