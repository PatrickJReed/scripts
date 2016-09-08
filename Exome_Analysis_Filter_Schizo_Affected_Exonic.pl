#!/usr/bin/perl -w
use strict;
use warnings;
my $Input_File1		= "./Schizo_Affected_Final.recal_Both.vcf.hg19_multianno.txt";
my $Output_File		= "./Schizo_Affected_Final.recal_Both.vcf.hg19_Coding.txt";
my (@tmp1, @IDs, %DATA1);

open (LIST, $Input_File1) || die "File not found\n";     
     while (<LIST>) {
     chomp;
     my $Line =$_;
         @tmp1 = split(/\t/, $_);
         if ($tmp1[1] ne "Start"){
		 $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"Chr"} = $tmp1[0];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"Start"} = $tmp1[1];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"End"} = $tmp1[2];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"Ref"} = $tmp1[3];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"Alt"} = $tmp1[4];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"Func.refGene"} = $tmp1[5];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"Gene.refGene"} = $tmp1[6];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"ExonicFunc.refGene"} = $tmp1[7];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"AAChange.refGene"} = $tmp1[8];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"esp6500si_all"} = $tmp1[9];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"1000g2012apr_all"} = $tmp1[10];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"snp138NonFlagged"} = $tmp1[11];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"clinvar_20140211"} = $tmp1[12];         
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"ljb23_all"} = $tmp1[13];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"caddgt10"} = $tmp1[14];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"OMIM_Regions"} = $tmp1[15];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"tfbs"} = $tmp1[16];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"gerp++gt2"} = $tmp1[17];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"gwasCatalog"} = $tmp1[18];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"wgRna"} = $tmp1[19];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"targetScanS"} = $tmp1[20];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"pubsBlat"} = $tmp1[21];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"DGVMerged"} = $tmp1[22];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"OtherInfo"} = $tmp1[23];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"All"} = $Line;
         

      	} 
     if ($DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"Func.refGene"} ne "exonic") {delete $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"};}
	elsif ($DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"ExonicFunc.refGene"} ne "nonsynonymous SNV") {delete $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"};}
	elsif ($DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"1000g2012apr_all"}) 
		{if ($DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"1000g2012apr_all"} > .1) {delete $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"};}}
};
close(LIST);


open(MYOUTPUTFILE, ">$Output_File");
@IDs =keys(%DATA1);
print MYOUTPUTFILE "Chr\tStart\tEnd\tRef\tAlt\tFunc.refGene\tGene.refGene\tExonicFunc.refGene\tAAChange.refGene\tesp6500si_all\t1000g2012apr_all\tsnp138NonFlagged\tclinvar_20140211\tljb23_all\tcaddgt10\tOMIM_Regions\ttfbs\tgerp++gt2\tgwasCatalog\twgRna\ttargetScanS\tpubsBlat\tdgvMerged\tOtherData\n"; 
foreach (@IDs) { 
print MYOUTPUTFILE $DATA1{$_}{"All"}."\n";
  	};		
close(MYOUTPUTFILE);
