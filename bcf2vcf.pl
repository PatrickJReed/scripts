#!/usr/bin/perl -w
use strict;
use warnings;


# Globals:
my $samtools	= "./samtools-0.1.18";
my $bcftools 	= "./samtools-0.1.18/bcftools/bcftools";
my $vcfutils 	= "./samtools-0.1.18/bcftools/vcfutils.pl";
my $ID			= "SZ_All";
my $bcf 		= $ID.".raw.bcf";
my $vcf			= $ID.".Exome_hg19.vcf";

	
	`$bcftools view $bcf | $vcfutils varFilter -D1000000 > $vcf`;