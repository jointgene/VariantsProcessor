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
die "Usage: perl $script conf.txt in.vcf out.vcf\n" if(not @ARGV);

# command line
# input 
my $file_conf= shift @ARGV;
my $file_vcf = shift @ARGV;
# output
my $file_out = shift @ARGV;

# read configure file
my %CONF = &read_conf($file_conf);
# example
#my %CONF = (
#	"major-allele-hom"  => 4,  #
#	"minor-allele-hom"  => 4,  #
#	"lower-dp"          => 7,  # each sample
#	"upper-dp"          => 100 # each sample
#);

# global set: for special experiment
my @POPULATION = (); # default: include all samples in this population
@POPULATION = split /\s*,\s*/, $CONF{'population'} if(exists $CONF{'population'});

# check required parameters
die "Configure Error: 'major-allele-hom' not defined\n" if(not exists $CONF{'major-allele-hom'});
die "Configure Error: 'minor-allele-hom' not defined\n" if(not exists $CONF{'minor-allele-hom'});

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
		my ($chr, $loc, $id, $ref, $alt, $qual, $flt, $inf, $tag, @samples) = (split /\t/, $_);
		@RG = @samples;
		say "samples = ". join ", ", @RG;
		#
		@POPULATION = @RG if(scalar @POPULATION == 0);
		# check
		my %samples=&assign_samples(\@samples);
		die "Configure file error!!!\n" if( &check_samples(\%samples, \@POPULATION) );
		# write out
		print OUT "$_\n";
		next;
	}
	# read vcf lines
	my ($chr, $loc, $id, $ref, $alt, $qual, $flt, $inf, $tag, @samples) = (split /\t/, $_);
	my %samples = &assign_samples(\@samples);
	my %sam_tag = &split_sample_tag($tag, \%samples, \@POPULATION) if($tag and @POPULATION);
	my %sam_cov = &get_allele_counts(\%sam_tag, \@POPULATION);
	my %sam_frq = &cal_allele_freq(\%sam_tag, \%sam_cov, \@POPULATION);

	# picking markers by the following conditions
	# condition 1
	my ($major_allele, $major_depth, $major_freq, $minor_allele, $minor_depth, $minor_freq) = &call_major_and_minor_allele(\%sam_cov, \@POPULATION);
	# check major allele and minor allele existed
	next if($major_allele eq "NA");	
	next if($minor_allele eq "NA");
	## minor allele frequency
	#next if($minor_freq < $CONF{'minor-allele-freq'});

	# condition 2
	# major-allele-hom and minor-allele-hom
	my ($right_major, $wrong_major) = &count_homozygous_sample_with_depth_filter(\%sam_tag, \@POPULATION, $major_allele) if(defined $major_allele);
	my ($right_minor, $wrong_minor) = &count_homozygous_sample_with_depth_filter(\%sam_tag, \@POPULATION, $minor_allele) if(defined $minor_allele);
	#
	next if($right_major < $CONF{'major-allele-hom'});
	next if($right_minor < $CONF{'minor-allele-hom'});

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
sub call_major_and_minor_allele{
	my $sam_cov = shift @_;
	my $rg_arr  = shift @_;
	#
	my %cov = ();
	%cov = &sum_cov_total($sam_cov, $rg_arr);
	#
	my $major_allele = "NA";
	my $major_depth  = 0;
	my $major_freq   = 0;
	my $minor_allele = "NA";
	my $minor_depth  = 0;
	my $minor_freq   = 0;
	my $total  = 0;
	if(%cov){
		my @allele = sort {$cov{$b} <=> $cov{$a}} keys %cov; 
		my $total  = 0;
		foreach my $a (@allele){
			$total += $cov{$a};
		}
		# major
		$major_allele = $allele[0];
		$major_depth  = $cov{$major_allele};
		$major_freq   = $major_depth/$total if($total>0);
		# minor
		if(scalar @allele >= 2){
			$minor_allele = $allele[1];
			$minor_depth  = $cov{$minor_allele};
			$minor_freq   = $minor_depth/$total if($total>0);
		}
	}
	return ($major_allele, $major_depth, $major_freq, $minor_allele, $minor_depth, $minor_freq);
}

# 
sub count_homozygous_sample_with_depth_filter{
	my $sam_tag = shift @_;
	my $rg_arr  = shift @_;
	my $allele  = shift @_;
	#
	my $right = 0;
	my $wrong = 0;
	foreach my $rg (@$rg_arr){
		next if(not exists $$sam_tag{$rg} or not $$sam_tag{$rg}{'GT'});
		next if(exists $CONF{'lower-dp'} and $$sam_tag{$rg}{'DP'} < $CONF{'lower-dp'});
		next if(exists $CONF{'upper-dp'} and $$sam_tag{$rg}{'DP'} > $CONF{'upper-dp'});
		#
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


