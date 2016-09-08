#!/usr/bin/perl -w
use strict;
use warnings;

# Very basic pipeline for expression analysis of RNAseq Data

my $tophat_DIR	= "/glusterfs/users/preed/tophat-2.0.9.Linux_x86_64";
my $cuff_DIR 	= "/glusterfs/users/preed/cufflinks-2.1.1.Linux_x86_64";
my $mm10_EBWT	= "/glusterfs/SEQreference/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome";
my $mm10_fa 	= "/glusterfs/SEQreference/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa";
my $mm10_GTF 	= "/glusterfs/SEQreference/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf";
my $Chrm_Dir	= "/glusterfs/SEQreference/Mus_musculus/UCSC/mm10/Sequence/Chromosomes/";
my $Ref_GTF		= "/glusterfs/SEQreference/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf";
my $MRG_GTF		= "/glusterfs/users/preed/KARA/ESC_EPI.combined.gtf";
my $ESC 		= "/glusterfs/users/preed/KARA/ESC.txt.gz"; 	
my $EPI			= "/glusterfs/users/preed/KARA/EpiSC.txt.gz";
my $EPI2ESC		= "/glusterfs/users/preed/KARA/EpiSC_ESC_P8.txt.gz";

#align all files against mm10



`$tophat_DIR/tophat2 -o ./ESC_mm10 -p 8 --report-secondary-alignments --coverage-search --segment-length 15 --keep-tmp $mm10_EBWT $ESC`; 


`$tophat_DIR/tophat2 -o ./EPI_mm10 -p 8 --report-secondary-alignments --coverage-search --segment-length 15 --keep-tmp $mm10_EBWT $EPI`;


`$tophat_DIR/tophat2 -o ./EPI2ESC_mm10 -p 8 --report-secondary-alignments --coverage-search --segment-length 15 --keep-tmp $mm10_EBWT $EPI2ESC`; 


`$cuff_DIR/cufflinks -o ./ESC_mm10 -p 8 -g $Ref_GTF -b $mm10_fa -u --no-update-check ./ESC_mm10/accepted_hits.bam`;


`$cuff_DIR/cufflinks -o ./EPI_mm10 -p 8 -g $Ref_GTF -b $mm10_fa -u --no-update-check ./EPI_mm10/accepted_hits.bam`;


`$cuff_DIR/cufflinks -o ./EPI2ESC_mm10 -p 8 -g $Ref_GTF -b $mm10_fa -u --no-update-check ./EPI2ESC_mm10/accepted_hits.bam`;


`$cuff_DIR/cuffcompare -o ESC_EPI -r $Ref_GTF -s $Chrm_Dir ./ESC_mm10/*.gtf ./EPI_mm10/*.gtf ./EPI2ESC_mm10/*.gtf`;


`$cuff_DIR/cuffdiff -L ESC,EPI,EPI2ESC -p 8 -b $mm10_fa -u $MRG_GTF ./ESC_mm10/accepted_hits.bam ./EPI_mm10/accepted_hits.bam ./EPI2ESC_mm10/accepted_hits.bam`;
