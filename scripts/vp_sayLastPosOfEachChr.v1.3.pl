#!/usr/bin/perl
# author: biotang
# version 1.1
use strict;
use warnings;
use File::Basename;

# usage
my $thisScript = basename $0;
my $USAGE = "Usage:
  perl $thisScript in.vcf
  gzip -c -d in.vcf.gz | perl $thisScript -
";
die "$USAGE\n" if(not @ARGV);

# command
my $file_vcf = shift @ARGV;

# global variable
my $CK_CHR;
my $MK_POS;
my $COUNT;

#
open IN, $file_vcf or die "";
# process headers and title
while(<IN>){
	if($_=~m/^#CHROM/){
		last;
	}
}
# vcf lines
while(<IN>){
	my ($chr, $pos) = (split /\t/, $_)[0,1];
	if(not defined $CK_CHR){
		print "$chr\t$pos";
		$CK_CHR = $chr;
		$COUNT = 0;
	} elsif($CK_CHR ne $chr){
		print "\t$MK_POS\t$COUNT\n";
		$CK_CHR = $chr;
		$COUNT = 0;
		print "$chr\t$pos";		
	}
	$MK_POS = $pos;
	$COUNT++;
}
print "\t$MK_POS\t$COUNT\n";
close IN;

