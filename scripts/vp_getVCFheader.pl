#!/usr/bin/perl
#
#
use strict;
use warnings;
use File::Basename;

# usage
my $thisScript = basename $0;
die "perl $thisScript in.vcf\n" if(not @ARGV);

# command arguments
my $file_vcf = shift @ARGV;

#
open IN, $file_vcf or die "";
while(<IN>){
	last if($_!~m/^##/);
	print $_;
}

close IN;

