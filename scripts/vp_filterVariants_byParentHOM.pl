#!/usr/bin/perl
# pick proper marker
# biotang
# 2016.8.1
# version: 2.1.5
use strict;
use warnings;
use 5.014;
use File::Basename;

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
# Parent1  = sampleNames
# Parent2  = sampleNames
# P1.Depth = 10
# P2.Depth = 10
# P1.Depth_upper = 100
# P2.Depth_upper = 100
# Total.Depth = 25
# Total.Depth_upper = 200
# P1.GT_WRONG = 0
# P2.GT_WRONG = 0
# P1.GT_RIGHT = 1
# P2.GT_RIGHT = 1

# check and over
#die "Configure file is error: Parent1 is missing\n" if(not exists $CONF{'Parent1'});
#die "Configure file is error: Parent2 is missing\n" if(not exists $CONF{'Parent2'});
die "Configure file is error: at least need one Parent\n" if(not exists $CONF{'parent1'} and not exists $CONF{'parent2'});

# global set: for special experiment
my @PARENT1 = ();
my @PARENT2 = ();
@PARENT1 = split /\s*,\s*/, $CONF{'parent1'} if(exists $CONF{'parent1'});
@PARENT2 = split /\s*,\s*/, $CONF{'parent2'} if(exists $CONF{'parent2'});

# global options: lower depth
my $DP_PAR1_TOTAL = 10; #
my $DP_PAR2_TOTAL = 10; #
$DP_PAR1_TOTAL = $CONF{'p1-lower-dp'} if(exists $CONF{'p1-lower-dp'});
$DP_PAR2_TOTAL = $CONF{'p2-lower-dp'} if(exists $CONF{'p2-lower-dp'});

# global options: upper depth
my $DP_PAR1_UPPER = $CONF{'p1-upper-dp'} if(exists $CONF{'p1-upper-dp'});
my $DP_PAR2_UPPER = $CONF{'p2-upper-dp'} if(exists $CONF{'p2-upper-dp'});

# global options: lower depth
my $TOTAL_DEPTH = $CONF{'total-lower-dp'} if(exists $CONF{'total-lower-dp'});

# global options: upper depth
my $TOTAL_DEPTH_UPPER = $CONF{'total-upper-dp'} if(exists $CONF{'total-upper-dp'});

# global options: limitation of genotype
my $P1GT_WRONG = 0;
my $P2GT_WRONG = 0;
my $n1 = scalar @PARENT1;
my $n2 = scalar @PARENT2;
my $P1GT_RIGHT = $n1;
my $P2GT_RIGHT = $n2;
$P1GT_RIGHT = $n1 - 1 if($n1 -1 >= 1);
$P2GT_RIGHT = $n2 - 1 if($n2 -1 >= 1);
$P1GT_WRONG =  $CONF{'p1-hom-wrong'} if(exists $CONF{'p1-hom-wrong'});
$P2GT_WRONG =  $CONF{'p2-hom-wrong'} if(exists $CONF{'p2-hom-wrong'});
$P1GT_RIGHT =  $CONF{'p1-hom-right'} if(exists $CONF{'p1-hom-right'});
$P2GT_RIGHT =  $CONF{'p2-hom-right'} if(exists $CONF{'p2-hom-right'});

#
#say "conditions:";
#say "DP_PAR1_TOTAL     = $DP_PAR1_TOTAL"       if(exists $CONF{'p1-lower-dp'});
#say "DP_PAR2_TOTAL     = $DP_PAR2_TOTAL"       if(exists $CONF{'p2-lower-dp'});
#say "DP_PAR1_UPPER     = $DP_PAR1_UPPER"       if(defined $DP_PAR1_UPPER);
#say "DP_PAR2_UPPER     = $DP_PAR2_UPPER"       if(defined $DP_PAR2_UPPER);
#say "TOTAL_DEPTH       = $TOTAL_DEPTH"         if(defined $TOTAL_DEPTH);
#say "TOTAL_DEPTH_UPPER = $TOTAL_DEPTH_UPPER"   if(defined $TOTAL_DEPTH_UPPER);

#say "P1.GT_WRONG = $P1GT_WRONG";
#say "P2.GT_WRONG = $P2GT_WRONG";
#say "P1.GT_RIGHT = $P1GT_RIGHT";
#say "P2.GT_RIGHT = $P2GT_RIGHT";

#say "Parent1 samples = ". join ", ", @PARENT1 if(exists $CONF{'parent1'});
#say "Parent2 samples = ". join ", ", @PARENT2 if(exists $CONF{'parent2'});
#say "";

####################################################
# main, read and write
####################################################
# global var
my @RG=();

#
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
		my ($chr,$loc,$id,$ref,$alt,$qual,$flt,$inf,$tag,@samples) = (split /\t/, $_);
		@RG = @samples;
		#say "samples = ". join ", ", @RG;
		# check
		my %samples = &assign_samples(\@samples);
		die "error sample name!!!\n\tParent1 missing.\n" if( &check_samples(\%samples, \@PARENT1) );
		die "error sample name!!!\n\tParent2 missing.\n" if( &check_samples(\%samples, \@PARENT2) );
		# write out
		print OUT "$_\n";
		next;
	}
	# read vcf lines
	my ($chr,$loc,$id,$ref,$alt,$qual,$flt,$inf,$tag,@samples) = (split /\t/, $_);
	my %pos_inf = &split_info($inf);
	my %samples = &assign_samples(\@samples);
	my %sam_tag = &split_sample_tag($tag, \%samples, \@RG);
	my %sam_cov = &get_allele_counts(\%sam_tag, \@RG);
	my %sam_frq = &cal_allele_freq(\%sam_tag, \%sam_cov, \@RG);

	# picking markers by the following conditions
	# condition 1
	# identify the main allele of parent
	my ($p1_allele, $p1_depth, $p1_freq) = &call_main_allele(\%sam_cov, \@PARENT1) if(scalar @PARENT1 > 0);
	my ($p2_allele, $p2_depth, $p2_freq) = &call_main_allele(\%sam_cov, \@PARENT2) if(scalar @PARENT2 > 0);

	# check allele
	next if(defined $p1_allele and $p1_allele eq "NA");
	next if(defined $p2_allele and $p2_allele eq "NA");
	next if(defined $p1_depth and $p1_depth < $DP_PAR1_TOTAL);
	next if(defined $p2_depth and $p2_depth < $DP_PAR2_TOTAL);
	next if(defined $p1_allele and defined $p2_allele and $p1_allele eq $p2_allele);
	next if(defined $DP_PAR1_UPPER and defined $p1_depth and $p1_depth > $DP_PAR1_UPPER);
	next if(defined $DP_PAR2_UPPER and defined $p2_depth and $p2_depth > $DP_PAR2_UPPER);
	next if(defined $TOTAL_DEPTH and defined $p1_depth and defined $p2_depth and $p1_depth + $p2_depth < $TOTAL_DEPTH);
	next if(defined $TOTAL_DEPTH_UPPER and defined $p1_depth and defined $p2_depth and $p1_depth + $p2_depth > $TOTAL_DEPTH_UPPER);

	# condition 2
	# check genotype (homozygous)
	my ($right1, $wrong1) = &count_homozygous_sample(\%sam_tag, \@PARENT1, $p1_allele) if(defined $p1_allele);
	my ($right2, $wrong2) = &count_homozygous_sample(\%sam_tag, \@PARENT2, $p2_allele) if(defined $p2_allele);
	# under strict limitation of samples number
	next if(defined $wrong1 and $wrong1 > $P1GT_WRONG);
	next if(defined $wrong2 and $wrong2 > $P2GT_WRONG);
	next if(defined $right1 and $right1 < $P1GT_RIGHT);
	next if(defined $right2 and $right2 < $P2GT_RIGHT);

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
	#print STDERR "read $file_conf\n";
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
		next  if(not $id or not $item);
		$CONF{$id}=$item;
	}
	close CONF;
	return %CONF;
}

####################################################
# functions: 
####################################################
sub call_main_allele{
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
	#
	my $allele = 0;
	my $depth  = 0;
	my $freq   = 0;
	my $total  = 0;
	if(not %cov){
		$allele = "NA";
		$depth  = 0;
		$freq   = 0;
	} else{
		foreach my $a (keys %cov){
			$total+=$cov{$a};
			if($depth<$cov{$a}){
				$allele=$a;
				$depth =$cov{$a};
			}
		}
		$freq = $depth/$total if($total>0);
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

sub split_info{
	my $inf =shift @_;
	
	my @info = split /;/,$inf;
	my %hash;
	for (@info){
		my ($key, $value) = split /=/, $_;
		$hash{$key} = $value;
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


