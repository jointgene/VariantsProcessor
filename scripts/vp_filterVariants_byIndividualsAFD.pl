#!/usr/bin/perl
# pick proper marker
# biotang
# 2019.5.4
# version: 1.1.3
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
#	"individuals-afd"   => 0.5,
#	"lower-dp"          => 3,
#	"upper-dp"          => 100,
#);

# global set: for special experiment
my @POPULATION = (); # default: include all samples in this population
@POPULATION = split /\s*,\s*/, $CONF{'population'} if(exists $CONF{'population'});
my @POPULATION2 = (); # default: include all samples in this population
@POPULATION2 = split /\s*,\s*/, $CONF{'population2'} if(exists $CONF{'population2'});

#
my $LOWER_DP;
my $UPPER_DP;
$LOWER_DP = $CONF{'lower-dp'} if(exists $CONF{'lower-dp'});
$UPPER_DP = $CONF{'upper-dp'} if(exists $CONF{'upper-dp'});

# check required parameters
die "Configure Error: 'individuals-afd' not defined\n" if(not exists $CONF{'individuals-afd'});
die "Configure Error: 'population' not defined\n" if(not exists $CONF{'population'});
die "Configure Error: 'population2' not defined\n" if(not exists $CONF{'population2'});

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
		print OUT "$_\n";
		next;
	}
	# read title
	if($_=~m/^#(.*)/){
		my ($chr, $loc, $id, $ref, $alt, $qual, $flt, $inf, $tag, @samples) = (split /\t/, $_);
		@RG = @samples;
		#say "samples = ". join ", ", @RG;
		# check
		my %samples = &assign_samples(\@samples);
		die "Sample error!!!\n" if( &check_samples(\%samples, \@POPULATION) );
		die "Sample error!!!\n" if( &check_samples(\%samples, \@POPULATION2) );
		# write out
		print OUT "$_\n";
		next;
	}
	# read vcf lines
	my ($chr, $loc, $id, $ref, $alt, $qual, $flt, $inf, $tag, @samples) = (split /\t/, $_);
	my %samples = &assign_samples(\@samples);
	my %sam_tag = &split_sample_tag($tag, \%samples, \@RG) if($tag and @RG);

	# picking markers by the following conditions
	my @alt = split /,/, $alt;
	my @allele = ($ref, @alt);
	my $num_allele = scalar @allele;
	my %allele_freq  = &calculate_individuals_af(\%sam_tag, \@POPULATION, $num_allele);
	my %allele_freq2 = &calculate_individuals_af(\%sam_tag, \@POPULATION2, $num_allele);
	#
	my $mk = 0;
	#my @frq1 = ();
	#my @frq2 = ();
	for(my $i=0;$i<$num_allele;$i++){
		$mk = 1 if( abs($allele_freq{$i} - $allele_freq2{$i}) > $CONF{'individuals-afd'});
		#push @frq1, $allele_freq{$i};
		#push @frq2, $allele_freq2{$i};
	}
	next if($mk == 0);

	# output
	print OUT "$_\n";

	#
	#my $str_frq1 = join ",", @frq1;
	#my $str_frq2 = join ",", @frq2;
	#print "$chr\t$loc\t$ref\t$alt\t".$str_frq1."\t".$str_frq2."\n";
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
	#print "read $file_conf\n";
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
sub calculate_individuals_af{
	my $sam_tag = shift @_;
	my $rg_arr  = shift @_;
	my $number  = shift @_;
	#
	my $individual_number = 0;
	my %allele_counts = ();
	#
	for(my $i=0;$i<$number;$i++){
		$allele_counts{$i}=0;
	}
	#
	foreach my $sample (@$rg_arr){
		next if(not exists $$sam_tag{$sample} or not exists $$sam_tag{$sample}{'GT'});
		next if(defined $LOWER_DP and $$sam_tag{$sample}{'DP'} < $LOWER_DP);
		next if(defined $UPPER_DP and $$sam_tag{$sample}{'DP'} > $UPPER_DP);
		#
		$individual_number++;
		#
		my @gt = split /[\/\|]/, $$sam_tag{$sample}{'GT'};
		foreach my $a (@gt){
			$allele_counts{$a} += 0.5;
		}
	}
	#
	my %allele_freq = ();
	foreach my $gt (sort {$a <=> $b} keys %allele_counts){
		my $counts = $allele_counts{$gt};
		my $freq = 0;
		$freq = $counts / $individual_number if($individual_number); 
		$allele_freq{$gt} = $freq;
	}
	return %allele_freq;
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

