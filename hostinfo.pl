#!/usr/bin/perl -w

use strict;
use Getopt::Long;
my ($printmachine,$printos);
GetOptions( "machine|m" => \$printmachine,
            "os|o" => \$printos )
    or die "usage: hostinfo.pl [-m|-o]\n";
    
if (!$printmachine && !$printos) {
    $printmachine = 1;
    $printos = 1;
}    

my $uname_a = `uname -a`;

my ($machine,$os);
$machine = 'i386' if $uname_a =~ /i?86/i;
$machine = 'ppc' if $uname_a =~ /powerpc/i;

$os = 'linux' if $uname_a =~ /linux/i;
$os = 'macosx' if $uname_a =~ /darwin/i;

defined $machine or die("hostinfo.pl: cannot determine type");
defined $os or die("hostinfo.pl: cannot determine os type");

my @output;
push @output, $machine if $printmachine;
push @output, $os if $printos;

print join(" ",@output),"\n";
