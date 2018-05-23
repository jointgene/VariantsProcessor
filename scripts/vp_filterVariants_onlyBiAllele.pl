#!/usr/bin/perl
# pick proper marker
# biotang
# 2016.8.1
use strict;
use warnings;
use File::Basename;
use 5.016;

# usage
my $script = basename $0;
#die "Usage: perl $script conf.txt in.vcf out.vcf\n" if(not @ARGV);
die "Usage: perl $script in.vcf out.vcf\n" if(not @ARGV);

# command line
# input 
my $file_vcf = shift @ARGV;
# output
my $file_out = shift @ARGV;

####################################################
# main, read and write
####################################################
# read
open FILE, $file_vcf or die "";
open OUT, ">$file_out" or die "";
while(<FILE>){
	chomp ;
	# read header info
	if($_=~m/^##/){
		#print OUT "$_\n";
		next;
	}
	# read title
	if($_=~m/^#(.*)/){
		# write out
		print OUT "$_\n";
		next;
	}
	# read vcf lines
	my $alt = (split /\t/, $_)[4];
	
	# filtered if not bi-allelic
	next if($alt =~ m/,/);

	# output
	print OUT "$_\n";
}
close FILE;
close OUT;

