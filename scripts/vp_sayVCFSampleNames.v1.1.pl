#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use 5.010;

# usage
my $script = basename $0;
die "Usage: perl $script in.vcf\n" if(not @ARGV);

# command
my $file_vcf = shift @ARGV;

# read vcf
open VCF, $file_vcf or die $!;
my $line = "";
while(<VCF>){
	last if($_!~m/^#/);
	$line = $_;
}
close VCF;

# sample names
chomp $line;
my ($chr,$loc,$id,$ref,$alt,$qual,$flt,$inf,$tag,@samples) = (split /\t/, $line);
my $i=1;
foreach my $name (@samples){
	print "$i:\t$name\n";
	$i++;
}


