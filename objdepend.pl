#!/usr/bin/perl -w
use strict;

my $pn = 'objdepend';

sub aggregate {
    my ($deps, $agg, $have, $entry, $depth) = @_;
    $have->{$entry} = 1;
    push @$agg, $entry;
    foreach my $dep (@{$deps->{$entry}}) {
        unless ( exists $have->{$dep} ) {
            aggregate($deps, $agg, $have, $dep, $depth+1);
        }
    }
}

my %deps;

my @targs;
while (1) {
    my $x = shift @ARGV;
    defined $x or die "no -- found";
    last if $x eq "--";
    push @targs, $x;
}

foreach my $fn (@ARGV) {
    
    open(IN,"<".$fn) or die $pn.": can't open file: '".$fn."'\n";

    while (<IN>) {
        my @toks = split /\s+:\s+/, $_;
        next if scalar(@toks) != 2;
        
        my @a = split /\s+/, $toks[0];
        my @b = split /\s+/, $toks[1];
        foreach my $aa (@a) {
            next unless $aa;
            foreach my $bb (@b) {
                next unless $bb;
                push @{$deps{$aa}}, $bb;
            }
        }
    } 
}

foreach my $x (@targs) {
    my %have;
    my @agg;
    aggregate(\%deps, \@agg, \%have, $x.".o", 0);
    
    my @objects  = grep { /.o$/ } @agg;
    print $x, " : ", join(" ", @objects), "\n";
}
