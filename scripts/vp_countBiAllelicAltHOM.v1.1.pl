#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use 5.010001;

# usage
my $script = basename $0;
die "Usage: perl $script conf.txt in.vcf out.xls\n" if(not @ARGV);

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
#	"lower-dp"          => 3,
#	"upper-dp"          => 100
#);

# global set: for special experiment
my @POPULATION = (); # default: include all samples in this population
@POPULATION = split /\s*,\s*/, $CONF{'population'} if(exists $CONF{'population'});

#
my $LOWER_DP;
my $UPPER_DP;
$LOWER_DP = $CONF{'lower-dp'} if(exists $CONF{'lower-dp'});
$UPPER_DP = $CONF{'upper-dp'} if(exists $CONF{'upper-dp'});

####################################################
# main, read and write
####################################################
# global var
my @RG=();
my %ALTHOM = ();

# read
open FILE, $file_vcf or die "";
while(<FILE>){
	chomp ;
	# read header info
	if($_=~m/^##/){
		next;
	}
	# read title
	if($_=~m/^#(.*)/){
		my ($chr, $loc, $id, $ref, $alt, $qual, $flt, $inf, $tag, @samples) = (split /\t/, $_);
		@RG = @samples;
		#say "samples = ". join ", ", @RG;
		#
		@POPULATION = @RG if(scalar @POPULATION == 0);
		# check
		my %samples = &assign_samples(\@samples);
		die "Configure file error!!!\n" if( &check_samples(\%samples, \@POPULATION) );
		next;
	}
	# read vcf lines
	my ($chr, $loc, $id, $ref, $alt, $qual, $flt, $inf, $tag, @samples) = (split /\t/, $_);
	
	# if not Bi-Allelic, next line
	next if($alt =~ m/,/);

	#
	my %samples = &assign_samples(\@samples);
	my %sam_tag = &split_sample_tag($tag, \%samples, \@POPULATION);

	# picking markers by the following conditions
	foreach my $sample (@POPULATION){
		next if(not exists $sam_tag{$sample});
		next if(defined $LOWER_DP and $sam_tag{$sample}{'DP'} < $LOWER_DP);
		next if(defined $UPPER_DP and $sam_tag{$sample}{'DP'} > $UPPER_DP);
		$ALTHOM{$sample}++ if($sam_tag{$sample}{'GT'} eq "1/1");
	}
}
close FILE;

#
open OUT, ">$file_out" or die "";
foreach my $sample (@POPULATION){
	print OUT "$sample\t$ALTHOM{$sample}\n";
}
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

