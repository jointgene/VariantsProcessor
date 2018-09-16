#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use 5.010;

# usage
my $script = basename $0;
die "Usage: perl $script list.tab in.vcf out.vcf\n" if(not @ARGV);

# command
my $file_tab = shift @ARGV;
my $file_vcf = shift @ARGV;
my $file_out = shift @ARGV;

# global var
my @RG=();
my %CHANGE = ();

# read sample names
open TAB, $file_tab or die $!;
while(<TAB>){
	chomp ;
	my ($name, $change) = split /\t/, $_;
	$CHANGE{$name} = $change;
}
close TAB;

# read vcf
open VCF, $file_vcf or die $!;
open OUT, ">$file_out" or die $!;
while(<VCF>){
	chomp ;
	if($_=~m/^##/){
		print OUT "$_\n";
	} elsif($_=~m/^#CHROM/){
		my ($chr,$loc,$id,$ref,$alt,$qual,$flt,$inf,$tag,@samples) = (split /\t/, $_);
		@RG = @samples;
		say "samples = ". join ", ", @RG;
		# check
		my %samples=&assign_samples(\@samples);
		my @name_array = keys %CHANGE;
		die "Sample Name error!!!\n" if( &check_samples(\%samples, \@name_array ) );
		#
		my @newsamples = ();
		foreach my $name (@samples){
			$name = $CHANGE{$name} if(exists $CHANGE{$name});
			push @newsamples, $name;
		}
		@RG = @newsamples;
		# print
		my $str = join "\t", @newsamples;
		print OUT "$chr\t$loc\t$id\t$ref\t$alt\t$qual\t$flt\t$inf\t$tag\t$str\n";

	} else{
		print OUT "$_\n";
	}
}
close VCF;

#
my $i=1;
foreach my $name (@RG){
	print "$i:\t$name\n";
	$i++;
}

####################################################
# functions: 
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


