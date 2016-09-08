#!/usr/bin/perl -w
use strict;
use warnings;

##Used GATK version 3.3-0 and got good results.

my $ID				= "2013-196_All";
my $samtools		= "/path_to/samtools-1.0/samtools";
my $picard			= "/path_to/picard-tools-1.118";
my $seqprep     	= "/path_to/SeqPrep-1.1/SeqPrep";
my $ref_fa			= "/path_to/hg19_Ordered.fa";
my $gatk 			= "/path_to/GenomeAnalysisTK.jar";
my $dbsnp_vcf   	= "/path_to/dbsnp_138_hg19_Ordered.vcf";
my $log1 			= $ID.".intervals.log.txt";
my $log2 			= $ID.".indelrealign.log.txt";
my $intervals_out 	= $ID.".intervals";	
		
	`samtools index $ID.srtrg.bam`;
			

########## Correct for Indels ##############


	`java -Xmx8g -Djava.io.tmpdir=./tmp -jar $gatk -R $ref_fa -T RealignerTargetCreator -I $ID.srtrg.bam -o $intervals_out -known $dbsnp_vcf -maxInterval 1000 -log $log1`;
		
	`java -Xmx8g -Djava.io.tmpdir=./tmp -jar $gatk -R $ref_fa -T IndelRealigner -LOD 2 -maxPosMove 400 -maxConsensuses 100 -greedy 1000 --maxReadsForRealignment 1000000 -targetIntervals $intervals_out -known $dbsnp_vcf -I $ID.srtrg.bam -o $ID.real.bam -log $log2`;
		
	`java -Xmx8g -Djava.io.tmpdir=./tmp -jar $picard/SortSam.jar I=$ID.real.bam O=$ID.srt.real.bam SO=coordinate`;
		
	`$samtools index $ID.srt.real.bam`;

	






