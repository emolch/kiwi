#!/usr/bin/perl -w

use strict;
use Getopt::Long;

my $pn = "fdepend.pl";

my $genfile;
my $incdep;
my @ignores;
GetOptions( "generate-file|g" => \$genfile,
            "include-depfile|d" => \$incdep,
            "ignore|i=s" => \@ignores )
    or die "usage: fdepend.pl [-d] [-g] [-i ignore.mod] source.f90 ...\n";

my %ignores = map { $_ => 1 } @ignores;
            
foreach my $fn (@ARGV) {
    open(IN,"<".$fn) or die $pn.": can't open file: '".$fn."'\n";
    my %deps;
    my %prov;
    while (<IN>) {
        next if /module\s+procedure\s+[A-Za-z0-9_]+/;
        if (/^\s*use\s+([A-Za-z0-9_]+)/) {
            $deps{$1.".mod"} = 1;
        }
        if (/^\s*module\s+([A-Za-z0-9_]+)/) {
            $prov{$1.".mod"} = 1;
        }
    }
    my @provides =  sort keys %prov;
    my @deps = grep { ! exists $ignores{$_} && ! exists $prov{$_} } sort keys %deps;
    my $fno = $fn;
    $fno =~ s/\.f90$/.o/;
    
    my $fnd = $fn;
    $fnd =~ s/\.f90$/.d/;
    my $str = $fno.($incdep?" ".$fnd:"")." : ".
              $fn."\n";
    if (@provides) {
        $str .= join(" ",@provides) . " : ". $fno."\n";
    }
    if (@deps) {
        $str .= $fno." : ".join(" ",@deps)."\n";
    }
    if ($genfile) {
        open(OUT,">".$fnd) or die $pn.": can't write to file: '".$fnd."'\n";
        print OUT $str;
        close OUT;
    } else {
       print $str
    }
    close IN;
}
