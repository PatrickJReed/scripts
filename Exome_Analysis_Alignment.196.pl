#!/usr/bin/perl -w
use strict;
use warnings;



my $bwa			= "/glusterfs/users/preed/bwa-0.7.10/bwa";
my $samtools	= "/glusterfs/users/preed/samtools-1.0/samtools";
my $picard		= "/glusterfs/users/preed/picard-tools-1.118";
my $seqprep     = "/glusterfs/users/preed/SeqPrep-1.1/SeqPrep";
my $ref_fa		= "/glusterfs/SEQreference/hg19/hg19_Ordered.fa";
my $maxMem		= 8000000000;
my $dir			= "/glusterfs/users/preed/HaloScreen/";
my $sample1_1 	= "2013-196_1_sequence.txt";
my $sample1_2 	= "2013-196_2_sequence.txt";

my $ID = "2013-196_All";
my @Fractions 	= (.9,.8,.7,.6,.5,.4,.3,.2,.1);
my @Replicate	= (1,2,3,4,5,6,7,8,9,10);

#prep PE reads, merge and clip
foreach my $fraction_n (@Fractions) {
		foreach my $replicate_n (@Replicate) {	


`$seqprep -6 -L 70 -A AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -B AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -f 2013-196-$fraction_n.$replicate_n._1_sequence.txt -r 2013-196-$fraction_n.$replicate_n._2_sequence.txt -1 $ID-$fraction_n.$replicate_n._1_sequence.txt -2 $ID-$fraction_n.$replicate_n._2_sequence.txt -3 $ID-$fraction_n.$replicate_n._1_sequence.x.txt -4 $ID-$fraction_n.$replicate_n._2_sequence.x.txt -s $ID-$fraction_n.$replicate_n.sequence.txt >& $ID-$fraction_n.$replicate_n.sp.log`;

	
	# Pair bwa alignments, convert to BAM, throw out unaligned reads, sort output by coords
		`$bwa mem -t 4 -M $ref_fa $ID-$fraction_n.$replicate_n._1_sequence.txt $ID-$fraction_n.$replicate_n._2_sequence.txt | 
		 $samtools view -bS - > $ID-$fraction_n.$replicate_n.pe.bam`;
		
		`$bwa mem -t 4 -M $ref_fa $ID-$fraction_n.$replicate_n._1_sequence.x.txt $ID-$fraction_n.$replicate_n._2_sequence.x.txt | 
		$samtools view -bS - > $ID-$fraction_n.$replicate_n.pe.x.bam`;
		
		`$bwa mem -t 4 -M $ref_fa $ID-$fraction_n.$replicate_n.sequence.txt | $samtools view -bS - > $ID-$fraction_n.$replicate_n.se.bam`;
		  
		  
		`java -Xmx8g -Djava.io.tmpdir=./tmp -jar $picard/SortSam.jar I=$ID-$fraction_n.$replicate_n.pe.bam O=$ID-$fraction_n.$replicate_n.pe.srt.bam SO=coordinate`;				
		
		`java -Xmx8g -Djava.io.tmpdir=./tmp -jar $picard/SortSam.jar I=$ID-$fraction_n.$replicate_n.pe.x.bam O=$ID-$fraction_n.$replicate_n.pe.x.srt.bam SO=coordinate`;				

		`java -Xmx8g -Djava.io.tmpdir=./tmp -jar $picard/SortSam.jar I=$ID-$fraction_n.$replicate_n.se.bam O=$ID-$fraction_n.$replicate_n.se.srt.bam SO=coordinate`;			



# Merge clipped and merged alignments
	
	`java -Xmx8g -Djava.io.tmpdir=./tmp -jar $picard/MergeSamFiles.jar I=$ID-$fraction_n.$replicate_n.pe.srt.bam I=$ID-$fraction_n.$replicate_n.pe.x.srt.bam I=$ID-$fraction_n.$replicate_n.se.srt.bam O=$ID-$fraction_n.$replicate_n.srt.bam SO=coordinate AS=true USE_THREADING=true`;
	
	#Add read group to bam header
	
	#`java -Xmx8g -Djava.io.tmpdir=./tmp -jar $picard/MarkDuplicates.jar I=$ID.srt.bam O=$ID.srt.mkdup.bam M=$ID.mkdups.out AS=true`;
	
	`java -Xmx8g -Djava.io.tmpdir=./tmp -jar $picard/AddOrReplaceReadGroups.jar I=$ID-$fraction_n.$replicate_n.srt.bam O=$ID-$fraction_n.$replicate_n.srtrg.bam RGID=All RGLB=Screen RGPL=Illumina RGPU=HISEQ2000 RGSM=$ID.$fraction_n.$replicate_n`;
			}
		}