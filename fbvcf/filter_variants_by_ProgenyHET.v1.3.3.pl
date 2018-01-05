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
die "Usage: perl $script conf.txt in.vcf out.vcf\n" if(not @ARGV);

# command line
# input 
my $file_conf= shift @ARGV;
my $file_vcf = shift @ARGV;
# output
my $file_out = shift @ARGV;

# read configure file
my %CONF = &read_conf($file_conf);

# check and over
die "Configure file is error: Progeny is missing\n" if(not exists $CONF{'Progeny'});

# global set: for special experiment
my @PROGENY = split /\s*,\s*/, $CONF{'Progeny'};

# global options
# DP
#my $DP_HET_TOTAL = 15; #
#$DP_HET_TOTAL = $CONF{'HET.TOTAL.DP'} if(exists $CONF{'HET.TOTAL.DP'});
my $DP_HET_EACH = 7;
$DP_HET_EACH =  $CONF{'HET.EACH.DP'} if(exists $CONF{'HET.EACH.DP'});
# AF
my $AF_HET_LOWER = 0.2;
$AF_HET_LOWER = $CONF{'AF.LOWER'} if(exists $CONF{'AF.LOWER'});
my $AF_HET_UPPER = 0.8;
$AF_HET_UPPER = $CONF{'AF.UPPER'} if(exists $CONF{'AF.UPPER'});
# Mis
my $MIS_HET_SAMPLE = 1;
$MIS_HET_SAMPLE = $CONF{'HET.Mis'} if(exists $CONF{'HET.Mis'});
# Fit
my $num_progeny = scalar @PROGENY;
my $FIT_HET_SAMPLE = $num_progeny;
$FIT_HET_SAMPLE = $CONF{'HET.Fit'} if(exists $CONF{'HET.Fit'});

#
say "conditions:";
#say "DP_HET_TOTAL   = $DP_HET_TOTAL";
say "DP_HET_EACH    = $DP_HET_EACH";
say "AF_HET_LOWER   = $AF_HET_LOWER";
say "AF_HET_UPPER   = $AF_HET_UPPER";
say "MIS_HET_SAMPLE = $MIS_HET_SAMPLE";
say "FIT_HET_SAMPLE = $FIT_HET_SAMPLE";
say "";

#
say "Progeny samples:";
say join ", ", @PROGENY;

####################################################
# main, read and write
####################################################
# global var
my @RG=();

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
		my ($chr,$loc,$id,$ref,$alt,$qual,$flt,$inf,$tag,@samples) = (split /\t/, $_);
		@RG = @samples;
		say "samples = ". join ", ", @RG;
		# check
		my %samples=&assign_samples(\@samples);
		die "Configure file error!!!\n" if( &check_samples(\%samples, \@PROGENY) );
		# write out
		print OUT "$_\n";
		next;
	}
	# read vcf lines
	my ($chr,$loc,$id,$ref,$alt,$qual,$flt,$inf,$tag,@samples) = (split /\t/, $_);
	my %samples = &assign_samples(\@samples);
	my %sam_tag = &split_sample_tag($tag, \%samples, \@RG);
	my %sam_cov = &get_allele_counts(\%sam_tag, \@RG);
	my %sam_frq = &cal_allele_freq(\%sam_tag, \%sam_cov, \@RG);

	# picking markers by the following conditions
	# condition 1
	my ($allele, $depth, $freq) = &call_minor_allele(\%sam_cov, \@PROGENY);
	# check allele
	next if($allele eq "NA");
	# condition 2
	my $het = 0;
	my $mis = 0;
	foreach my $rg (@PROGENY){
		do {$mis++; next;} if(not exists $sam_tag{$rg});
		do {$mis++; next;} if($sam_tag{$rg}{'DP'} < $DP_HET_EACH);
		next if($sam_frq{$rg}{$allele} < $AF_HET_LOWER);
		next if($sam_frq{$rg}{$allele} > $AF_HET_UPPER);
		$het++;
	}
	next if($mis > $MIS_HET_SAMPLE);
	next if($het < $FIT_HET_SAMPLE);
	# output
	print OUT "$_\n";
}
close FILE;
close OUT;

####################################################
# functions: read configure file
####################################################
sub read_conf{
	# read configure file
	my $file_conf = shift @_;
	my %CONF = ();
	print "read $file_conf\n";
	open CONF, "<", $file_conf or die "";
	while(<CONF>){
		chomp ;
		next if($_=~m/^#/);
		next if($_=~m/^\s/);
		next if($_ eq "");
		$_=~s/\s*#.*$//;
		$_=~s/\s*;\s*$//;
		$_=~s/\s*([=\t])\s*/$1/;
		my ($id, $item) = (split /[=\t]/, $_);
		next if(not defined $id or not defined $item);
		$CONF{$id}=$item;
	}
	close CONF;
	return %CONF;
}

####################################################
# functions: 
####################################################
sub call_minor_allele{
	my $sam_cov = shift @_;
	my $rg_arr  = shift @_;
	#
	my %cov = ();
	%cov = &sum_cov_total($sam_cov, $rg_arr);
	#
	my $allele = "NA";
	my $depth  = 0;
	my $freq   = 0;
	my $total  = 0;
	if(%cov){
		my @allele = sort {$cov{$b} <=> $cov{$a}} keys %cov; 
		my $total  = 0;
		foreach my $a (@allele){
			$total += $cov{$a};
		}
		if(scalar @allele >= 2){
			$allele = $allele[1];
			$depth  = $cov{$allele};
			$freq   = $depth/$total if($total>0);
		}
	}
	return ($allele, $depth, $freq);
}

# 
sub count_homozygous_sample{
	my $sam_tag = shift @_;
	my $rg_arr  = shift @_;
	my $allele  = shift @_;
	#
	my $right = 0;
	my $wrong = 0;
	foreach my $rg (@$rg_arr){
		next if(not exists $$sam_tag{$rg} or not $$sam_tag{$rg}{'GT'});
		my ($mk, $mis) = &is_specific_homozygous($$sam_tag{$rg}{'GT'}, $allele); 
		$right++ if($mk==1);
		$wrong++ if($mk==0);
	}
	return ($right, $wrong);
}

#
sub is_specific_homozygous{
	my $str_gt = shift @_;
	my $allele = shift @_;
	#
	my $mk = 1;
	my @gt = split /\//, $str_gt;
	$mk=0 if($gt[0] == $gt[1] and $gt[0] != $allele);
	$mk=0 if($gt[0] != $gt[1]);
	return $mk;
}

sub sum_cov_total{
	my $sam_cov = shift @_;
	my $rg_arr  = shift @_;
	#
	my %cov=();
	foreach my $rg (@$rg_arr){
		next if(not exists $$sam_cov{$rg});
		#
		foreach my $i (keys %{$$sam_cov{$rg} } ){
			$cov{$i}+=$$sam_cov{$rg}{$i};
		}
	}
	return %cov;
}

sub mean_frq_total{
	my $sam_frq = shift @_;
	my $rg_arr  = shift @_;
	#
	my %sum=();
	my %num=();
	foreach my $rg (@$rg_arr){
		next if(not exists $$sam_frq{$rg});
		#
		foreach my $i (keys %{$$sam_frq{$rg} } ){
			$sum{$i}+=$$sam_frq{$rg}{$i};
			$num{$i}++;
		}
	}
	my %frq=();
	foreach my $i (keys %sum){
		$frq{$i}=$sum{$i}/$num{$i};
	}
	return %frq;
}

#
sub cal_mean_AF{
	my $sam_tag = shift @_;
	my $sam_frq = shift @_;
	my $rg_arr  = shift @_;
	my $allele  = shift @_;
	#
	my $rg_num = scalar @$rg_arr;
	my $rg_mis = 0;
	my $total  = 0;
	foreach my $rg (@$rg_arr){
		if(not exists $$sam_frq{$rg} or $$sam_tag{$rg}{'DP'}<$DP_HET_EACH){
			$rg_mis++;
			next;
		}
		$total+=$$sam_frq{$rg}{$allele};
	}
	my $num = $rg_num - $rg_mis;
	my $mean_AF = "NA";
	$mean_AF = $total/$num if($num>0);
	return ($mean_AF, $rg_mis);
}

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

sub split_sample_tag{
	my $tag     = shift @_;
	my $samples = shift @_;
	my $rg_arr  = shift @_;
	#
	my @tag = split /:/, $tag;
	#
	my %hash=();
	foreach my $rg (@$rg_arr){
		next if($$samples{$rg} eq ".");
		my @arr = split /:/, $$samples{$rg};
		for(my $i=0;$i<scalar @tag;$i++){
			$hash{$rg}{$tag[$i]}=$arr[$i];
		}
	}
	return %hash;
}

sub get_allele_counts{
	my $sam_tag = shift @_;
	my $rg_arr  = shift @_;
	#
	my %hash = ();
	foreach my $rg (@$rg_arr){
		next if(not exists $$sam_tag{$rg});
		my @cov = ();
		push @cov, $$sam_tag{$rg}{'RO'};
		push @cov, (split /,/, $$sam_tag{$rg}{'AO'});
		for(my $i=0;$i<scalar @cov;$i++){
			$hash{$rg}{$i}=$cov[$i];
		}
	}
	return %hash;
}

sub cal_allele_freq{
	my $sam_tag = shift @_;
	my $sam_cov = shift @_;
	my $rg_arr  = shift @_;
	#
	my %hash = ();
	foreach my $rg (@$rg_arr){
		next if(not exists $$sam_cov{$rg} or not exists $$sam_tag{$rg});
		next if($$sam_tag{$rg}{'DP'}==0);
		foreach my $i (keys %{$$sam_cov{$rg} } ){
			$hash{$rg}{$i} = $$sam_cov{$rg}{$i}/$$sam_tag{$rg}{'DP'};
		}
	}
	return %hash;
}

#
sub sum{
	my $arr = shift @_;
	my $sum = 0;
	foreach my $value (@$arr){
		$sum+=$value;
	}
	return $sum;
}


