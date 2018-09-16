#!/usr/bin/perl
use strict;
use warnings;
use 5.010;
use File::Basename;

# usage
my $script = basename $0;
die "Usage: perl $script raw.fb.vcf slim.fb.vcf\n" if(not @ARGV);

# command line
my $file_vcf = shift @ARGV;
my $file_out = shift @ARGV;

#
open FILE, $file_vcf or die "";
open OUT, ">$file_out" or die "";
while(<FILE>){
	chomp ;
	if($_=~m/^##/){
		#print OUT "$_\n";
		next;
	}
	if($_=~m/^#/){
		print OUT "$_\n";
		next;
	}
	#
	my ($chr,$loc,$id,$ref,$alt,$qual,$flt,$inf,$tag,@samples) = (split /\t/, $_);
	my $new_inf=&slim_info(\$inf);
	my $new_tag="GT:DP:AD";
	my $new_samples=&slim_tag(\$tag,\@samples);
	print OUT "$chr\t$loc\t$id\t$ref\t$alt\t$qual\t$flt\t$new_inf\t".
	"$new_tag\t".$new_samples."\n";
}
close FILE;
close OUT;

sub slim_info{
	my $inf = shift @_;
	# per-base coverage depth, combined across samples
	my $dp = $1 if($$inf=~m/\bDP=(\S+?)[;\t]/);
	# Allele count in genotype, for each ALT allele
	my $ac = $1 if($$inf=~m/\bAC=(\S+?)[;\t]/);
	# Allele frequency in genotype, for each ALT allele
	my $af = $1 if($$inf=~m/\bAF=(\S+?)[;\t]/);
	#
	my $new="DP=$dp;AC=$ac;AF=$af;";	
	return $new;
}

sub slim_tag{
	my $tag = shift @_;
	my $spl = shift @_;

	my @new = ();
	my @key = split /:/, $$tag;
	foreach my $str (@$spl){
		if($str =~ m"^\."){
			push @new, ".";
		}
		else{
			my @value = split /:/, $str;
			my %tag = ();
			for(my $i=0;$i<scalar @key;$i++){
				$tag{$key[$i]}=$value[$i];
			}
			push @new, "$tag{'GT'}:$tag{'DP'}:$tag{'AD'}";
		}
	}
	my $new_samples=join "\t",@new;
	return $new_samples;
}

