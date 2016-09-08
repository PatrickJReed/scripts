#!/usr/bin/perl -w
use strict;
use warnings;



my $lobSTR		= "/glusterfs/users/preed/lobSTR-bin-Linux-x86_64-2.0.5/bin/lobSTR";
my $allelotype	= "/glusterfs/users/preed/lobSTR-bin-Linux-x86_64-2.0.5/bin/allelotype";
my $lobSTR_idx	= "/glusterfs/users/preed/hg19/lobstr_v2.0.3_hg19_ref/lobSTR_";
my $lobSTR_str	= "/glusterfs/users/preed/hg19/lobstr_v2.0.3_hg19_strinfo.tab";
my $lobSTR_mdl	= "/glusterfs/users/preed/lobSTR-bin-Linux-x86_64-2.0.5/share/lobSTR/models";
my $samtools	= "/usr/local/tools/samtools-0.1.18/samtools";
my $picard		= "/glusterfs/users/preed/picard-tools-1.102";
my $gatk 		= "/glusterfs/users/preed/GenomeAnalysisTK-2.1-13/GenomeAnalysisTK.jar";
my $ref_fa		= "/glusterfs/SEQreference/hg19/hg19_Ordered.fa";
my $annovar		= "/glusterfs/users/preed/annovar/";
my $bedtools	= "/glusterfs/users/preed/BEDTools-Version-2.12.0/bin";
my $Halo_Exons  = "HaloScreen_Covered.bed";

my $sample1_1 	= "2013-196_130418_SN484_0205_BD1L1AACXX_7_1_sequence.txt.gz";
my $sample1_2 	= "2013-196_130418_SN484_0205_BD1L1AACXX_7_2_sequence.txt.gz";
my $sample2_1 	= "2013-196_131025_SN484_0229_AC2Y21ACXX_6_1_sequence.txt.gz";
my $sample2_2 	= "2013-196_131025_SN484_0229_AC2Y21ACXX_6_2_sequence.txt.gz";
my $sample3_1 	= "2013-196_131203_SN1070_0144_AH7LY5ADXX_1_1_sequence.txt.gz";
my $sample3_2 	= "2013-196_131203_SN1070_0144_AH7LY5ADXX_1_2_sequence.txt.gz";

my $ID = "2013-196_lobSTR";

#`$lobSTR --p1 $sample1_1,$sample2_1,$sample3_1 --p2 $sample1_2,$sample2_2,$sample3_2 -q --gzip -o $ID --index-prefix $lobSTR_idx --rg-sample $ID --rg-lib $ID -p 8`;

#`$samtools sort $ID.aligned.bam $ID.srt`;

#`$bedtools/intersectBed -abam $ID.srt.bam -b $Halo_Exons > $ID.Final.bam`;

#`$samtools index $ID.Final.bam`; 

#`$allelotype --command classify --no-rmdup --noweb --bam $ID.Final.bam --noise_model $lobSTR_mdl/illumina_v2.0.3 --out $ID --strinfo $lobSTR_str --index-prefix $lobSTR_idx`;

`$annovar/convert2annovar.pl -format vcf4old --includeinfo --withzyg $ID.vcf > $ID.avinput`;

`$annovar/table_annovar.pl $ID.avinput $annovar/humandb/ -buildver hg19 --remove --otherinfo --outfile $ID -protocol refGene,esp6500si_all,1000g2012apr_all,snp138NonFlagged,clinvar_20140211,ljb23_all,caddgt10,OMIM_Regions,tfbs,gerp++gt2,gwasCatalog,wgRna,targetScanS,pubsBlat,genomicSuperDups -operation g,f,f,f,f,f,f,r,r,f,r,r,r,r,r`;