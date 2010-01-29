#!/usr/bin/perl -w

use strict;


my $atdoc = 0;
while (<>) {

    if (! $atdoc && /^\s*!! /) {
        $atdoc = 1;
    }

    if ($atdoc) {
        if (/^\s*!!? ?(.*)/) {
            print $1,"\n";
        } else {
            print "\n";
            $atdoc = 0;

        }
    }

        
}
