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

# keys of configure file
# progeny       = sampleNames # required
# lower-dp      = 7    # depth of each sample, default is 7
# upper-dp      = 200 # not required
# lower-af      = 0.2  # swith 1
# upper-af      = 0.8  # swith 1
# mis-num       = 0    # default is 0, (0 mean not allow missing)
# het-num       = 1    # default is the number of all samples
# lower-avgaf   = 0.25 # swith 2
# upper-avgaf   = 0.75 # swith 2
# swith 1 or swith 2, at least one is required,

# check required parameters
die "Configure file is error: Progeny is missing\n" if(not exists $CONF{'progeny'});

# global set: for special experiment
my @PROGENY = split /\s*,\s*/, $CONF{'progeny'};

# global options
# DP
my $DP_HET_EACH = 7;
$DP_HET_EACH = $CONF{'lower-dp'} if(exists $CONF{'lower-dp'});
my $DP_HET_EACH_UPPER = $CONF{'upper-dp'} if(exists $CONF{'upper-dp'});

# Mis
my $MIS_HET_SAMPLE = 0;
$MIS_HET_SAMPLE = $CONF{'mis-num'} if(exists $CONF{'mis-num'});

# swith 1
# AF
my $AF_HET_LOWER = $CONF{'lower-af'} if(exists $CONF{'lower-af'});
my $AF_HET_UPPER = $CONF{'upper-af'} if(exists $CONF{'upper-af'});
# Fit
my $num_progeny = scalar @PROGENY;
my $FIT_HET_SAMPLE = $num_progeny;
$FIT_HET_SAMPLE = $CONF{'het-num'} if(exists $CONF{'het-num'});

# swith 2
# AF
my $AvgAF_HET_LOWER = $CONF{'lower-avgaf'} if(exists $CONF{'lower-avgaf'});
my $AvgAF_HET_UPPER = $CONF{'upper-avgaf'} if(exists $CONF{'upper-avgaf'});

# if both are not defined, swith 2 will be open
$AvgAF_HET_LOWER = 0.3 if(not defined $AF_HET_LOWER and not defined $AvgAF_HET_LOWER);
$AvgAF_HET_UPPER = 0.7 if(not defined $AF_HET_UPPER and not defined $AvgAF_HET_UPPER);

#
say "Progeny samples:";
say join ", ", @PROGENY;
say "";

#
say "conditions:";
say "DP_HET_EACH       = $DP_HET_EACH";
say "DP_HET_EACH_UPPER = $DP_HET_EACH_UPPER" if(defined $DP_HET_EACH_UPPER);
say "MIS_HET_SAMPLE    = $MIS_HET_SAMPLE";
say "AF_HET_LOWER      = $AF_HET_LOWER"    if(defined $AF_HET_LOWER);
say "AF_HET_UPPER      = $AF_HET_UPPER"    if(defined $AF_HET_UPPER);
say "FIT_HET_SAMPLE    = $FIT_HET_SAMPLE"  if(defined $AF_HET_LOWER and defined $AF_HET_UPPER);
say "AvgAF_HET_LOWER   = $AvgAF_HET_LOWER" if(defined $AvgAF_HET_LOWER);
say "AvgAF_HET_UPPER   = $AvgAF_HET_UPPER" if(defined $AvgAF_HET_UPPER);
say "";

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
	my $het   = 0;
	my $mis   = 0;
	my $sumAF = 0;
	my $noMis = 0;
	foreach my $rg (@PROGENY){
		do {$mis++; next;} if(not exists $sam_tag{$rg});
		do {$mis++; next;} if($sam_tag{$rg}{'DP'} < $DP_HET_EACH);
		do {$mis++; next;} if(defined $DP_HET_EACH_UPPER and $sam_tag{$rg}{'DP'} > $DP_HET_EACH_UPPER);
		$sumAF+= $sam_frq{$rg}{$allele};
		$noMis++;
		#
		next if(defined $AF_HET_LOWER and $sam_frq{$rg}{$allele} < $AF_HET_LOWER);
		next if(defined $AF_HET_UPPER and $sam_frq{$rg}{$allele} > $AF_HET_UPPER);
		$het++;
	}
	# about missing data
	next if($mis > $MIS_HET_SAMPLE);

	# swith 1
	next if(defined $AF_HET_LOWER and defined $AF_HET_UPPER and $het < $FIT_HET_SAMPLE);

	# swith 2
	next if($noMis == 0);	
	my $avgAF = $sumAF / $noMis;
	next if(defined $AvgAF_HET_LOWER and defined $AvgAF_HET_UPPER and $avgAF < $AvgAF_HET_LOWER);
	next if(defined $AvgAF_HET_LOWER and defined $AvgAF_HET_UPPER and $avgAF > $AvgAF_HET_UPPER);

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
		next if($$samples{$rg}=~m"^\.");
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
		# from GATK calling
		if(exists $$sam_tag{$rg}{'AD'}){ 
			push @cov, (split /,/, $$sam_tag{$rg}{'AD'});
		}
		# from freebayes calling
		elsif(exists $$sam_tag{$rg}{'RO'} and exists $$sam_tag{$rg}{'AO'}){
			push @cov, $$sam_tag{$rg}{'RO'};
			push @cov, (split /,/, $$sam_tag{$rg}{'AO'});
		}
		# save in a hash
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


