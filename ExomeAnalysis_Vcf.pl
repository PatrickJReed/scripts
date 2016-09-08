#!/usr/bin/perl -w
use strict;
use warnings;


# Globals:
my $ref_fa 		= "/glusterfs/SEQreference/hg19/hg19_Ordered.fa";
my $ref_index 	= "/glusterfs/SEQreference/hg19/hg19_Ordered.fa.fai";
my $samtools	= "/usr/local/tools/samtools-0.1.18/samtools";
my $bcftools 	= "/usr/local/tools/samtools-0.1.18/bcftools/bcftools";
my $vcfutils 	= "/usr/local/tools/samtools-0.1.18/bcftools/vcfutils.pl";
my $ID			= "K87_NEW_7";
my $bcf 		= $ID.".raw.bcf";
my $vcf			= $ID.".Exome_hg19.vcf";

	
	
	`$samtools mpileup -Dgu -f $ref_fa ./*.recalibrated_final.bam | $bcftools view -bvcg - > $bcf`;  

	`$bcftools view $bcf | $vcfutils varFilter -D1000000 > $vcf`;
	
	
