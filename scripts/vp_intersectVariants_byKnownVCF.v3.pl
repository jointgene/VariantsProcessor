#!/usr/bin/perl
use strict;
use warnings;
use 5.010;
use File::Basename;

# usage
my $script = basename $0;
die "Usage: perl $script known.vcf input.vcf out.vcf\n" if(not @ARGV);

# command line
# input
my $file_vcf1 = shift @ARGV;
my $file_vcf2 = shift @ARGV;
# output
my $file_out  = shift @ARGV;

#
open IN1, $file_vcf1 or die "";
open IN2, $file_vcf2 or die "";
open OUT, ">$file_out"  or die "";

# vcf1
while(my $vcf1=<IN1>){
	next if($vcf1=~m/^##/);
	last if($vcf1=~m/^#/);
}
# vcf2
while(my $vcf2=<IN2>){
	if($vcf2=~m/^##/){
		print OUT $vcf2;
		next;
	}
	if($vcf2=~m/^#/){
		print OUT $vcf2;
		last;
	}
}

# vcf lines
my $lineA  = "";
my $lineB  = "";
my $ck_chr = "";
while($lineA=<IN1>){
	my $keyA = &get_key($lineA);
	#
	last if(eof(IN2));
	$lineB=<IN2> if(not $lineB);
	my $keyB = &get_key($lineB);
	
	#
	my $chrA = (split /:/, $keyA)[0];
	my $chrB = (split /:/, $keyB)[0];
	$ck_chr = $chrA if($ck_chr ne $chrA and $ck_chr ne $chrB and $chrA eq $chrB);
	
	#
	while($keyA ne $keyB){
		my $ord = &check_order($lineA, $lineB, $ck_chr);
		last if($ord==0);
		last if(eof(IN2));
		$lineB=<IN2>;
		$keyB = &get_key($lineB);
		$chrB = (split /:/, $keyB)[0];
		$ck_chr = $chrA if($ck_chr ne $chrA and $ck_chr ne $chrB and $chrA eq $chrB);
	}
	#
	if($keyA eq $keyB){
		print OUT $lineB;
		$ck_chr = (split /:/, $keyA)[0];
		#
		$lineA="";
		$lineB="";
	}
}
close IN1;
close IN2;
close OUT;

# functions
sub get_key{
	my $line = shift @_;
	my ($chr, $loc, $id, $ref, $alt) = (split /\t/, $line)[0,1,2,3,4];
	my $key = "$chr:$loc:$ref:$alt";
	return $key;
}

sub check_order{
	my $lineA = shift @_;
	my $lineB = shift @_;
	my $ck_chr = shift @_;
	#
	my ($chr1, $loc1) = (split /\t/, $lineA)[0,1];
	my ($chr2, $loc2) = (split /\t/, $lineB)[0,1];
	#
	return 0 if($chr1 eq $ck_chr and $chr2 ne $ck_chr);
	return 1 if($chr2 eq $ck_chr and $chr1 ne $ck_chr);
	return 0 if($loc1 < $loc2);
	return 1 if($loc1 > $loc2);
}


