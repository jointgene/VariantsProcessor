#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;

# usage
my $thisScript = basename $0;
my $USAGE = "Usage:
  perl $thisScript chr_list.txt in.vcf out.vcf
  perl $thisScript chr_list.txt in.vcf - | bgzip -c >out.vcf.gz
  gzip -c -d in.vcf.gz | perl $thisScript chr_list.txt - out.vcf
  gzip -c -d in.vcf.gz | perl $thisScript chr_list.txt - - | bgzip -c >out.vcf.gz

tabix for bgzip file:
  tabix -p vcf out.vcf.gz

Format of chr_list.txt: 
  the first column must be chr names and no header
";
die "$USAGE\n" if(not @ARGV);

# command
my $file_chr = shift @ARGV;
my $file_vcf = shift @ARGV;
my $file_out = shift @ARGV;

#
open CHR, $file_chr or die "$file_chr cannot be open\n";
my @chr_list = ();
my %chr_list = ();
while(<CHR>){
	chomp ;
	my $chr = (split /\t/, $_)[0];
	push @chr_list, $chr;
	$chr_list{$chr}++;
}
close CHR;

#
open IN, $file_vcf or die "$file_vcf cannot be open\n";
open OUT, ">$file_out" or die "$file_out cannot be writen\n";
# process headers and title
while(<IN>){
	if($_=~m/^#CHROM/){
		print OUT $_;
		last;
	}
	if($_=~m/^##contig=<ID=(\S+?),/){
		print OUT $_ if($1 and exists $chr_list{$1});
		next;
	}
	print OUT $_;
}
# process vcf lines
while(<IN>){
	my $chr = (split /\t/, $_)[0];
	print OUT $_ if($chr and exists $chr_list{$chr});
}
close IN;
close OUT;

