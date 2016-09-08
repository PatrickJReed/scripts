#!/usr/bin/perl -w
use warnings;

 my $GPL					= "./GPL96_Info.txt";
 my $Transcriptome 			= "./Tissue_Transcriptome.txt";
 my $Final_Transcriptome 	= "./Tissue_Expression_Final.txt";
 my ($S1, $S2, %EData, @Line);


#Load GPL Annotation Data
open (LIST, $GPL) || die "File not found\n";     
     while (<LIST>) {
         ($S1,$S2) = split($_);
         print "$S1";
         print"\n";
         #$EData{"$S1"}{"Affy_ID"} = $S1;
         #$EData{"$S1"}{"Gene_Symbol"} = $S2;  
     };
close(LIST);
exit;
#Load Expression Data
#open (LIST, $Transcriptome) || die "File not found\n";     
#     while (<LIST>) {
#      @Line = split(/\t/, $_);
      #print @Line."\n";
#         if ($EData{$Line[0]}{"Gene_Symbol"}) {
#         		$EData{$Line[0]}{"Entire_Line"} = "@Line";
#         				}
#         			};            
#close(LIST);

my @IDs =keys(%EData);
print "@IDs";
open(MYOUTPUTFILE, ">$Final_Transcriptome");
	
foreach (@IDs) { 
	print MYOUTPUTFILE $EData{$_}{"Affy_ID"}."\t";
	print MYOUTPUTFILE $EData{$_}{"Gene_Symbol"}."\n";
	#print MYOUTPUTFILE $EData{$_}{"Entire_Line"}."\n";

  }		

close(MYOUTPUTFILE);

 
