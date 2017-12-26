#!/usr/bin/perl
# author: biotangweiqi
# 2017.12.26
use strict;
use warnings;
use File::Basename;
use 5.010;

# usage
my $thisScript = basename $0;
die "Usage: perl $thisScript <in.bed> <in.vcf> <out.txt>\n" if(not @ARGV);

# command arguments
my $file_bed = shift @ARGV;
my $file_vcf = shift @ARGV;
my $file_out = shift @ARGV;

# global variables
my %bedlines  = ();
my %bed_start = ();
my %bed_end   = ();
my $num_bedline = 0;

# read bed into the hash
open BED, $file_bed or die "BED file $file_bed cannnot open\n";
while(<BED>){
	chomp ;
	my ($chr, $start, $end) = (split /\t/, $_)[0,1,2];
	push @{$bed_start{$chr} }, $start;
	push @{$bed_end{$chr} }, $end;
	$bedlines{"$chr:$start-$end"} = $_;
	$num_bedline++;
}
close BED;

# sort bed-lines for each chr in hash
foreach my $chr (keys %bed_start){
	my $size = scalar @{$bed_start{$chr} };
	my @order = sort {
		$bed_start{$chr}[$a] <=> $bed_start{$chr}[$b] or
		$bed_end{$chr}[$a] <=> $bed_end{$chr}[$b] } (0 .. $size-1);
	#
	my @start = ();
	my @end   = ();
	foreach my $i (@order){
		push @start, $bed_start{$chr}[$i];
		push @end,   $bed_end{$chr}[$i];
	}
	@{$bed_start{$chr} } = @start;
	@{$bed_end{$chr} }   = @end;
}

# global variable
my $VCF_TITLE = "";

#
my $ck_chr = "";
my $mk_idx = 0;
my %vcf_in_mem = ();

# read vcf and write out
open VCF, $file_vcf or die "VCF file $file_vcf cannot open\n";
open OUT, ">$file_out" or die "output file $file_out cannot write\n";
while(<VCF>){
	chomp ;
	next if($_=~m/^##/);
	do {$VCF_TITLE = $_; next} if($_=~m/^#CHROM/);
	
	# read vcf lines
	my ($chr, $pos, $ref) = (split /\t/, $_)[0,1,3];
	my $end = $pos + length($ref) - 1;
	
	# read next vcf line
	next if(not exists $bed_start{$chr});
	
	# new chr
	if($chr ne $ck_chr){
		if($ck_chr){ 
			# write out and clean
			foreach my $i (sort {$a <=> $b} keys %vcf_in_mem){
				&write_bed_vcf($bedlines{"$ck_chr:$bed_start{$ck_chr}[$i]-$bed_end{$ck_chr}[$i]"}, $vcf_in_mem{$i});
				delete $vcf_in_mem{$i};
			}
		}
		# initial
		$ck_chr = $chr;
		$mk_idx = 0;
	}

	# before, read next vcf line
	next if($end < $bed_start{$chr}[$mk_idx]);

	# within, process
	for(my $i=$mk_idx;$i<scalar @{$bed_start{$chr} };$i++){
		# before
		last if($end < $bed_start{$chr}[$i]);
		# within, process (record into mem)
		if($pos >= $bed_start{$chr}[$i] and $end <= $bed_end{$chr}[$i]){
			push @{$vcf_in_mem{$i} }, $_;
		}
	}

	# change the index of mark, and write out
	for(my $i=$mk_idx;$i<scalar @{$bed_start{$chr} };$i++){
		# before or within
		do {$mk_idx=$i; last} if($pos <= $bed_end{$chr}[$i]);
		# after, the bed segment is over time
		if($pos > $bed_end{$chr}[$i]){
			$mk_idx=$i;
			# write out and clean
			&write_bed_vcf($bedlines{"$chr:$bed_start{$chr}[$i]-$bed_end{$chr}[$i]"}, $vcf_in_mem{$i});
			delete $vcf_in_mem{$i};
		}
	}
}
close VCF;
close OUT;

# sub-functions
sub write_bed_vcf{
	my $bedline  = shift @_;
	my $vcflines = shift @_;
	#
	print OUT "BED\t$bedline\n";
	print OUT "VCF\t$VCF_TITLE\n";
	my $count = 0;
	foreach my $line (@$vcflines){
		print OUT "VCF\t$line\n";
		$count++;
	}
	print OUT "END\t$count\n";
}

