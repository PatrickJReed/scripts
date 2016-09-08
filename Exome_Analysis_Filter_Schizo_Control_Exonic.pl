#!/usr/bin/perl -w
use strict;
use warnings;
my $Input_File1		= "/Users/patrickreed/VCF_Files/2013-25_NEW.recal_both_haplotypes.AtaxiaExons.specific.vcf";
my $Input_File2		= "/Users/patrickreed/VCF_Files/2013-25_X.recal_both_haplotypes.AtaxiaExons.specific.vcf";
my $Input_File2		= "/Users/patrickreed/VCF_Files/2013-25_Overlap_QD.txt";
my (@tmp1, @tmp2, @IDs, %DATA);

open (LIST, $Input_File1) || die "File not found\n";
     while (<LIST>) {
     chomp;
         @tmp1 = split(/[;\t\/]/, $_);
		 $DATA{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$tmp1[4]"}{"X"}{"Quality"} = $tmp1[5];
		 
         foreach my $info (@tmp1) {
         @tmp2 = split('=',$info);
         if($tmp2[0] eq "DP") {
         $DATA{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$tmp1[4]"}{"X"}{"Depth"} = $tmp2[1];
    	}
	}	
};
close(LIST);

open (LIST, $Input_File2) || die "File not found\n";
     while (<LIST>) {
     chomp;
         @tmp1 = split(/[;\t\/]/, $_);
		 $DATA{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$tmp1[4]"}{"NEW"}{"Quality"} = $tmp1[5];
		 
         foreach my $info (@tmp1) {
         @tmp2 = split('=',$info);
         if($tmp2[0] eq "DP") {
         $DATA{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$tmp1[4]"}{"NEW"}{"Depth"} = $tmp2[1];
    	}
	}	
};
close(LIST);



open(MYOUTPUTFILE, ">$Output_File");
@IDs =keys(%DATA);
print MYOUTPUTFILE "Quality_X\tDepth_X\tQuality_NEW\tDepth_NEW\n";
foreach (@IDs) {
print MYOUTPUTFILE $DATA1{$_}{"X"}{"Quality"}."\t";
print MYOUTPUTFILE $DATA1{$_}{"X"}{"Depth"}."\t";
print MYOUTPUTFILE $DATA1{$_}{"NEW"}{"Quality"}."\t";
print MYOUTPUTFILE $DATA1{$_}{"NEW"}{"Depth"}."\n";
  	};
close(MYOUTPUTFILE);
