#!/usr/bin/perl -w
use strict;
use warnings;

my $ID				= "2013-196_All";
my $samtools		= "/path_to/samtools-1.0/samtools";
my $gatk 			= "/path_to/GenomeAnalysisTK.jar";
my $ref_fa			= "/path_to/hg19_Ordered.fa";
my $dbsnp_vcf   	= "/path_to/dbsnp138_hg19_Ordered.vcf";
my $bam_FINAL		= $ID.".recalibrated_final.bam";
my $recal_file  	= $ID.".recal_data_test.grp";
my $recal1_log   	= $ID.".recal1.log";
my $recal2_log   	= $ID.".recal2.log";

	
########## Base Recalibration ##############
	`$samtools index $ID.srt.real.bam`;
		
	`java -Xmx8g -Djava.io.tmpdir=./tmp -jar $gatk -R $ref_fa -T BaseRecalibrator -log $recal1_log -dfrac 0.1 --knownSites $dbsnp_vcf -I $ID.srt.real.bam --out $recal_file`;
	
	
 	`java -Xmx8g -Djava.io.tmpdir=./tmp -jar $gatk -R $ref_fa -T PrintReads -I $ID.srt.real.bam -log $recal2_log --out $bam_FINAL -BQSR $recal_file`;





