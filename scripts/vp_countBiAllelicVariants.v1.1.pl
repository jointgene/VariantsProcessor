#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use 5.010001;

# usage
my $thisScript = basename $0;
die "Usage: perl $thisScript in.vcf\n" if(not @ARGV);

# command arguments
# input
my $file_vcf  = shift @ARGV;

# counter
my $counter = 0;
my $snp_count = 0;
my $indel_count = 0;
my $other_count = 0;

# read
open FILE, $file_vcf or die "";
while(<FILE>){
	# chomp ;
	# header and title
	if($_=~m/^#/){
		next;
	}
	# lines
	my ($ref,$alt) = (split /\t/, $_)[3,4];
	#
	if($alt =~ m/,/){
		$other_count++;
	} elsif(length $ref == 1 and length $alt == 1){
		$snp_count++;
	} elsif(length $ref != length $alt){
		$indel_count++;
	} else{
		$other_count++;
	}
	$counter++;
}
close FILE;

# summary
printf "Bi-allelic SNP\t$snp_count\t%.3f\n", $snp_count/$counter if($counter>0);
printf "Bi-allelic InDel\t$indel_count\t%.3f\n", $indel_count/$counter if($counter>0);
printf "Others\t$other_count\t%.3f\n", $other_count/$counter if($counter>0);
print  "Total\t$counter\n";

