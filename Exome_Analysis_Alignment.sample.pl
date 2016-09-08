#!/usr/bin/perl -w
use strict;
use warnings;

#reference genome should be indexed
#sequence files can be gziped (I think)
#look at seqprep docs for how to find adaptor_sequence to use or use defaults 

my $bwa			= "/path_to/bwa-0.7.10/bwa";
my $samtools	= "/path_to/samtools-1.0/samtools";
my $picard		= "/path_to/picard-tools-1.118";
my $seqprep     = "/path_to/SeqPrep-1.1/SeqPrep";
my $ref_fa		= "/path_to/hg19_Ordered.fa";
my $dir			= "/path_to/HaloScreen/";
my $sample1_1 	= "/path_to/2013-196_1_sequence.txt";
my $sample1_2 	= "/path_to/2013-196_2_sequence.txt";
my $ID = "2013-196_All";

`$seqprep -6 -L 70 -A adaptor_sequence -B adaptor_sequence -f 2013-196_1_sequence.txt -r 2013-196_2_sequence.txt -1 $ID_1_sequence.txt -2 $ID_2_sequence.txt -3 $ID_1_sequence.x.txt -4 $ID_2_sequence.x.txt -s $ID_sequence.txt >& $ID_sp.log`;

	
	# Pair bwa alignments, convert to BAM, throw out unaligned reads, sort output by coords
		`$bwa mem -t 8 -M $ref_fa $ID_1_sequence.txt $ID_2_sequence.txt | 
		 $samtools view -bS - > $ID.pe.bam`;
		
		`$bwa mem -t 8 -M $ref_fa $ID_1_sequence.x.txt $ID_2_sequence.x.txt | 
		$samtools view -bS - > $ID.pe.x.bam`;
		
		`$bwa mem -t 8 -M $ref_fa $ID.sequence.txt | $samtools view -bS - > $ID.se.bam`;  
		  
		`java -Xmx8g -Djava.io.tmpdir=./tmp -jar $picard/SortSam.jar I=$ID.pe.bam O=$ID.pe.srt.bam SO=coordinate`;				
		
		`java -Xmx8g -Djava.io.tmpdir=./tmp -jar $picard/SortSam.jar I=$ID.pe.x.bam O=$ID.pe.x.srt.bam SO=coordinate`;				

		`java -Xmx8g -Djava.io.tmpdir=./tmp -jar $picard/SortSam.jar I=$ID.se.bam O=$ID.se.srt.bam SO=coordinate`;			



# Merge clipped and overlapping alignments
	
	`java -Xmx8g -Djava.io.tmpdir=./tmp -jar $picard/MergeSamFiles.jar I=$ID.pe.srt.bam I=$ID.pe.x.srt.bam I=$ID.se.srt.bam O=$ID.srt.bam SO=coordinate AS=true USE_THREADING=true`;

#Mark Duplicates	
	
	`java -Xmx8g -Djava.io.tmpdir=./tmp -jar $picard/MarkDuplicates.jar I=$ID.srt.bam O=$ID.srt.mkdup.bam M=$ID.mkdups.out AS=true`;

	#Add read group to bam header
	
	`java -Xmx8g -Djava.io.tmpdir=./tmp -jar $picard/AddOrReplaceReadGroups.jar I=$ID.srt.mkdup.bam O=$ID.srtrg.bam RGID=All RGLB=Screen RGPL=Illumina RGPU=HISEQ2000 RGSM=$ID`;


