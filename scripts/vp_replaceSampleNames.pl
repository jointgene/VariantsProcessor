#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;

# usage
my $thisScript = basename $0;
my $USAGE = "Usage:
  perl $thisScript name_list.txt in.vcf out.vcf
  perl $thisScript name_list.txt in.vcf - | bgzip -c >out.vcf.gz
  gzip -c -d in.vcf.gz | perl $thisScript name_list.txt - out.vcf
  gzip -c -d in.vcf.gz | perl $thisScript name_list.txt - - | bgzip -c >out.vcf.gz

tabix for bgzip file:
  tabix -p vcf out.vcf.gz

Format of name_list.txt: 
  the first column must be old-name of samples
  the second column must be new-name of samples
";
die "$USAGE\n" if(not @ARGV);

# command
my $file_lst = shift @ARGV;
my $file_vcf = shift @ARGV;
my $file_out = shift @ARGV;


# read the list file
open LIST, $file_lst or die "$file_lst cannot be open\n";
my %name_list = ();
while(<LIST>){
	chomp ;
	my ($old, $new) = (split /\t/, $_)[0,1];
	$name_list{$old} = $new;
}
close LIST;

#
open IN, $file_vcf or die "$file_vcf cannot be open\n";
open OUT, ">$file_out" or die "$file_out cannot be writen\n";
# process headers and title
while(<IN>){
	chomp ;
	if($_=~m/^#CHROM/){
		my ($chr, $loc, $id, $ref, $alt, $qual, $flt, $inf, $tag, @samples) = (split /\t/, $_);
		#
		my @new_rg = ();
		foreach my $rg (@samples){
			push @new_rg, $name_list{$rg};
		}
		my $str = join "\t", @new_rg;
		#
		print OUT "$chr\t$loc\t$id\t$ref\t$alt\t$qual\t$flt\t$inf\t$tag\t$str\n";
		last;
	}
	if($_=~m/^##/){
		print OUT "$_\n";
		next;
	}
}

# process vcf lines
while(<IN>){
	chomp ;
	print OUT "$_\n";
}
close IN;
close OUT;

