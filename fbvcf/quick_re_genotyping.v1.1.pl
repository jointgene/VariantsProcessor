#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use 5.010;

# usage
my $script = basename $0;
die "Usage: perl $script raw_GT.fb.vcf new_GT.fb.vcf\n" if(not @ARGV);

# command line
my $file_vcf = shift @ARGV;# full or slim version
my $file_out = shift @ARGV;

# global options
my $MIS_LIMIT = 1;
my $MAF_LIMIT = 0.2;

# open filehandle
open IN, $file_vcf or die "";
open OUT,">$file_out" or die "";
while(<IN>){
	chomp ;
	# header
	if($_=~m/^##/){
		#print OUT "$_\n";
		next;
	}
	# title
	if($_=~m/^#/){
		print OUT "$_\n";
		next;
	}
	# lines
	my ($chr,$loc,$id,$ref,$alt,$qual,$flt,$inf,$tag,@samples) = (split /\t/, $_);
	&quick_re_genotyping($tag, \@samples);
	my $str=join "\t", @samples;
	print OUT "$chr\t$loc\t$id\t$ref\t$alt\t$qual\t$flt\t$inf\t$tag\t$str\n";
}

# close filehandle
close IN;
close OUT;

# functions
sub quick_re_genotyping{
	my $tag = shift @_;
	my $arr = shift @_;
	#
	my @tag = split /:/, $tag;
	my @new = ();
	foreach my $str (@$arr){
		if($str eq "."){
			push @new, ".";
			next;
		}
		#
		my %hash = &read_sample_into_hash($str,\@tag);
		my %hcov = &read_cov_into_hash("$hash{'RO'},$hash{'AO'}");# fb-style
		my %rate = &cal_allele_freq(\%hcov, $hash{'DP'});
		#
		my $gt = &simple_genotyping(\%rate, \%hcov);
		$str=~s/\Q$hash{'GT'}\E/$gt/;
		push @new,$str;
	}
}

sub read_sample_into_hash{
	my $str = shift @_;
	my $tag = shift @_;

	#
	my @arr=split /:/, $str;
	my %hash=();
	for(my $i=0;$i<scalar @$tag;$i++){
		$hash{$$tag[$i]}=$arr[$i];
	}
	return %hash;	
}

sub read_cov_into_hash{
	my $cov = shift @_;
	#
	my @cov = ();
	push @cov, (split /,/, $cov);
	#
	my %hcov = ();
	for(my $i=0;$i<scalar @cov;$i++){
		$hcov{$i}=$cov[$i];
	}
	return %hcov;
}

sub cal_allele_freq{
	my $hcov  = shift @_;
	my $depth = shift @_;
	my %rate  = ();
	foreach my $al (keys %$hcov){
		$rate{$al} = 0 if(not $depth);
		$rate{$al} = $$hcov{$al}/$depth;
	}
	return %rate;
}

sub simple_genotyping{
	my $haf = shift @_;
	my $hcov= shift @_;
	
	#
	my ($max_af,$minor_af) = (sort {$b<=>$a} values %$haf)[0,1];
	my ($max_allele,$minor_allele)=(0,0);
	
	foreach my $allele (keys %$haf){
		if($max_af == $$haf{$allele}){
			$max_allele=$allele;
			last;
		}
	}
	foreach my $allele (keys %$haf){
		if($minor_af == $$haf{$allele} and $allele ne $max_allele){
			$minor_allele=$allele;
		}
	}
	#
	my $minor_cov = $$hcov{$minor_allele};
	
	#
	my $gt="";
	my $a1=$max_allele;
	my $a2=$minor_allele;
	# Condition:
	if($minor_cov<=$MIS_LIMIT and $minor_af<$MAF_LIMIT){
		$gt="$a1/$a1";
	}
	else{
		if($a1>$a2){
			my $tmp=$a1;
			$a1=$a2;
			$a2=$tmp;
		}
		$gt="$a1/$a2";;
	}
	return $gt;
}


