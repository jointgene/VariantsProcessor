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
  the first column must be the name of samples
";
die "$USAGE\n" if(not @ARGV);

# command
my $file_lst = shift @ARGV;
my $file_vcf = shift @ARGV;
my $file_out = shift @ARGV;


# read the list file
open LIST, $file_lst or die "$file_lst cannot be open\n";
my @NEWORDER = ();
while(<LIST>){
	chomp ;
	my $name = (split /\t/, $_)[0];
	push @NEWORDER, $name;
}
close LIST;

#
# global var
my @RG=();

#
open IN, $file_vcf or die "$file_vcf cannot be open\n";
open OUT, ">$file_out" or die "$file_out cannot be writen\n";
# process headers and title
while(<IN>){
	chomp ;
	# read header info
	if($_=~m/^##/){
		print OUT "$_\n";
		next;
	}
	# read title	
	if($_=~m/^#CHROM/){
		my ($chr, $loc, $id, $ref, $alt, $qual, $flt, $inf, $tag, @samples) = (split /\t/, $_);
		#
		@RG = @samples;
		my %samples = &assign_samples(\@samples);
		die "error sample name!!!\n" if( &check_samples(\%samples, \@NEWORDER) );

		# re-order
		my @new_rg = ();
		foreach my $rg (@NEWORDER){
			push @new_rg, $samples{$rg};
		}
		my $str = join "\t", @new_rg;
		#
		print OUT "$chr\t$loc\t$id\t$ref\t$alt\t$qual\t$flt\t$inf\t$tag\t$str\n";
		next;
	}
	# read vcf lines
	my ($chr, $loc, $id, $ref, $alt, $qual, $flt, $inf, $tag, @samples) = (split /\t/, $_);
	my %samples = &assign_samples(\@samples);
	
	# re-order
	my @new_rg = ();
	foreach my $rg (@NEWORDER){
		push @new_rg, $samples{$rg};
	}
	my $str = join "\t", @new_rg;
	#
	print OUT "$chr\t$loc\t$id\t$ref\t$alt\t$qual\t$flt\t$inf\t$tag\t$str\n";
}

close IN;
close OUT;

####################################################
# functions: read vcf
####################################################
sub assign_samples{
	my $samples = shift @_;
	#
	my %hash = ();
	for(my $i=0;$i<scalar @RG;$i++){
		$hash{$RG[$i]}=$$samples[$i];
	}
	return %hash;
}

sub check_samples{
	my $samples = shift @_;
	my $rg_arr  = shift @_;
	my $err = 0;
	foreach my $rg (@$rg_arr){
		if(not exists $$samples{$rg}){
			print STDERR "$rg not exists\n";
			$err++;
		}
	}
	return $err;	
}

