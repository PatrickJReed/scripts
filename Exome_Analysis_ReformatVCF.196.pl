#!/usr/bin/perl -w
use strict;
use warnings;
use Bio::Phenotype::OMIM::OMIMparser;
my $Input_File1		= "./2013-196_final.hg19_multianno.txt";
my $Input_File2		= "../pubsArticle.txt";
my $Input_File3		= "../genemap";
my $Input_File4		= "../omim.txt";
my $Output_File		= "./2013-196_final.Annotated.txt";
my (%DATA1, %DATA2, %OMIM_DATA, @tmp1, @tmp2, @tmp_omim, $HEAD, $article, @titles, @links, @Dumb, @IDs);
my $omim_parser = Bio::Phenotype::OMIM::OMIMparser->new( -genemap  => $Input_File3, -omimtext => $Input_File4 );
my %Genes; 
$Genes{"ABHD12"} = "True";
$Genes{"ADCK3"} = "True";
$Genes{"AFG3L2"} = "True";
$Genes{"ANO10"} = "True";
$Genes{"APTX"} = "True";
$Genes{"ATCAY"} = "True";
$Genes{"ATM"} = "True";
$Genes{"ATXN1"} = "True";
$Genes{"ATXN2"} = "True";
$Genes{"ATXN3"} = "True";
$Genes{"ATXN7"} = "True";
$Genes{"ATXN10"} = "True";
$Genes{"BEAN1"} = "True";
$Genes{"C10orf2"} = "True";
$Genes{"CA8"} = "True";
$Genes{"CACNA1A"} = "True";
$Genes{"CACNB4"} = "True";
$Genes{"COL18A1"} = "True";
$Genes{"CSTB"} = "True";
$Genes{"DNAJC19"} = "True";
$Genes{"EEF2"} = "True";
$Genes{"FGF14"} = "True";
$Genes{"FLVCR1"} = "True";
$Genes{"FXN"} = "True";
$Genes{"GOSR2"} = "True";
$Genes{"ITPR1"} = "True";
$Genes{"KCNA1"} = "True";
$Genes{"KCNC3"} = "True";
$Genes{"KCNJ10"} = "True";
$Genes{"KIAA0226"} = "True";
$Genes{"KLHL1"} = "True";
$Genes{"MRE11A"} = "True";
$Genes{"NOP56"} = "True";
$Genes{"PDYN"} = "True";
$Genes{"PEX10"} = "True";
$Genes{"POLG"} = "True";
$Genes{"PPP2R2B"} = "True";
$Genes{"PRICKLE1"} = "True";
$Genes{"PRKCG"} = "True";
$Genes{"SACS"} = "True";
$Genes{"SETX"} = "True";
$Genes{"SIL1"} = "True";
$Genes{"SLC1A3"} = "True";
$Genes{"SPTBN2"} = "True";
$Genes{"SYNE1"} = "True";
$Genes{"SYT14"} = "True";
$Genes{"TBP"} = "True";
$Genes{"TDP1"} = "True";
$Genes{"TGM6"} = "True";
$Genes{"TTBK2"} = "True";
$Genes{"TTPA"} = "True";
$Genes{"VLDLR"} = "True";
$Genes{"WDR81"} = "True"; 

open (LIST, $Input_File1) || die "File not found\n";     
     while (<LIST>) {
     chomp;
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
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"snp137"} = $tmp1[11];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"phastConsElements46way"} = $tmp1[12];         
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"LJB2_SIFT"} = $tmp1[13];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"LJB2_PolyPhen2_HDIV"} = $tmp1[14];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"LJB2_PolyPhen2_HVAR"} = $tmp1[16];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"LJB2_LRT"} = $tmp1[18];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"LJB2_MutationTaster"} = $tmp1[20];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"LJB_MutationAssessor"} = $tmp1[22];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"LJB2_FATHMM"} = $tmp1[24];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"LJB2_GERP++"} = $tmp1[25];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"LJB2_PhyloP"} = $tmp1[26];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"LJB2_SiPhy"} = $tmp1[27];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"avsift"} = $tmp1[28];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"OMIM"} = $tmp1[29];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"tfbs"} = $tmp1[30];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"gerp++gt2"} = $tmp1[31];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"gwasCatalog"} = $tmp1[32];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"evofold"} = $tmp1[33];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"wgRna"} = $tmp1[34];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"pubsBlat"} = $tmp1[35];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"miRNA_BindingSites"} = $tmp1[36];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"Info"} = $tmp1[45];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"Format"} = $tmp1[46];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"Data"} = $tmp1[47];
      	} 	 
     if ($DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"Func.refGene"} eq "intergenic") {delete $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"};}
	elsif ($DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"Func.refGene"} eq "ncRNA_intronic") {delete $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"};}
	elsif ($DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"Func.refGene"} eq "intronic") {delete $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"};}
	elsif ($DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"Func.refGene"} eq "upstream") {delete $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"};}
	elsif ($DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"Func.refGene"} eq "downstream") {delete $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"};}
	elsif ($DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"Func.refGene"} eq "UTR3") {delete $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"};}
	elsif ($DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"Func.refGene"} eq "UTR5") {delete $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"};}
	elsif ($DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"ExonicFunc.refGene"} eq "synonymous SNV") {delete $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"};}
	elsif ($DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"1000g2012apr_all"}) 
		{if ($DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"1000g2012apr_all"} > .1) {delete $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"};}}		                      
    elsif ($DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"Func.refGene"} eq "exonic") 
    	{unless ($Genes{$DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"Gene.refGene"}}) {delete $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"};}}  
     };
close(LIST);


#my $Covered=0;
#foreach (keys %DATA1) {
#my @tmp = split(':', $DATA1{$_}{"Data"});
#if ($tmp[2] < 8) { delete $DATA1{$_}};
#};

my $varcount = scalar keys %DATA1;
print $varcount."\n";


open (LIST, $Input_File2) || die "File not found\n";     
     while (<LIST>) {
     chomp;
         @tmp2 = split(/\t/, $_);
         $DATA2{$tmp2[0]}{"Article"} = "@tmp2";                    
     };
close(LIST);
print "2\n";
@IDs =keys(%DATA1);
foreach (@IDs) {
my @Pub_Articles=();
@Dumb = split('=', $DATA1{$_}{"pubsBlat"});  
	if ($Dumb[1]) {
	my @Articles = split(/\,/, $Dumb[1]);
		foreach my $article (@Articles) {
			push(@Pub_Articles,$DATA2{$article}{"Article"});		
			}
 		}
 	$DATA1{$_}{"Pub_Articles"} = "@Pub_Articles";
	};
print "3\n";
while ( my $omim_entry = $omim_parser->next_phenotype() ) {

    my $numb  = $omim_entry->MIM_number();                     # *FIELD* NO
    my $title = $omim_entry->title();                          # *FIELD* TI - first line
    my $desc  = $omim_entry->description();                    # *FIELD* TX
    my $cs    = $omim_entry->clinical_symptoms();              # *FIELD* CS
	
	#print Dumper(%$cs);
   	my @Classes = keys%{$cs};
   	my @Traits;
   		foreach (@Classes) {
   		my $Trait_Count = scalar(@{%{$cs}->{$_}}); 
   		my $count = 0;
   			while ($count < $Trait_Count) {
   				push(@Traits, @{%{$cs}->{$_}}->[$count]); 
   				$count++
   			}	
   		};		
    $OMIM_DATA{$numb}{"Title"} = $title;
    my @DESCs = split(/\n+/,$desc);
    $OMIM_DATA{$numb}{"Description"} = "@DESCs";
    $OMIM_DATA{$numb}{"ClinicalSymptoms"} = "@Traits";
	 
  };
print "4\n";
@IDs =keys(%DATA1);
foreach (@IDs) {
my @OMIM_Titles=();
my @OMIM_Descriptions=();
my @OMIM_Symptoms=();
my @Dumb2 = split('=', $DATA1{$_}{"OMIM"}); 
#print "@Dumb2\n"; 
	if ($Dumb2[1]) {
	my @MIMs = split(/\,/, $Dumb2[1]);
	my @OMIMs = ();
		foreach my $mim (@MIMs) {
		my @tmp3 = split(':', $mim);
		push(@OMIMs, $tmp3[1]);
		}
	#print "@OMIMs\n";
		foreach my $mim (@OMIMs) {
			push(@OMIM_Titles,$OMIM_DATA{$mim}{"Title"});		
			push(@OMIM_Descriptions,$OMIM_DATA{$mim}{"Description"});
			push(@OMIM_Symptoms,$OMIM_DATA{$mim}{"ClinicalSymptoms"});
			}
 		}
 	$DATA1{$_}{"OMIM_Titles"} = "@OMIM_Titles";
 	$DATA1{$_}{"OMIM_Descriptions"} = "@OMIM_Descriptions";
 	$DATA1{$_}{"OMIM_Symptoms"} = "@OMIM_Symptoms";	
	};



print "5\n";
@IDs =keys(%DATA1);
open(MYOUTPUTFILE, ">$Output_File");

print "6\n";
@IDs =keys(%DATA1);
print MYOUTPUTFILE "Chr\tStart\tEnd\tRef\tAlt\tFunc.refGene\tGene.refGene\tExonicFunc.refGene\tAAChange.refGene\tesp6500si_all\t1000g2012apr_all\tsnp137\tphastConsElements46way\tSIFT\tPolyPhen2_HDIV\tPolyPhen2_HVAR\tLRT\tMutationTaster\tMutationAssessor\tFATHMM\tGERP++\tPhyloP\tSiPhy\tavsift\tTFBS\tgerp++gt2\tgwasCatalog\tevofold\twgRna\tmiRNA_BindingSites\tPub_Articles\tOMIM_Titles\tOMIM_Descriptions\tOMIM_Symptoms\tInfo\tFormat\tData\n";
foreach (@IDs) { 
print MYOUTPUTFILE	$DATA1{$_}{"Chr"}."\t";
print MYOUTPUTFILE  $DATA1{$_}{"Start"}."\t";
print MYOUTPUTFILE  $DATA1{$_}{"End"}."\t";
print MYOUTPUTFILE  $DATA1{$_}{"Ref"}."\t";
print MYOUTPUTFILE  $DATA1{$_}{"Alt"}."\t";
print MYOUTPUTFILE  $DATA1{$_}{"Func.refGene"}."\t";
print MYOUTPUTFILE  $DATA1{$_}{"Gene.refGene"}."\t";
print MYOUTPUTFILE  $DATA1{$_}{"ExonicFunc.refGene"}."\t";
print MYOUTPUTFILE  $DATA1{$_}{"AAChange.refGene"}."\t";
print MYOUTPUTFILE  $DATA1{$_}{"esp6500si_all"}."\t";
print MYOUTPUTFILE  $DATA1{$_}{"1000g2012apr_all"}."\t";
print MYOUTPUTFILE  $DATA1{$_}{"snp137"}."\t";
print MYOUTPUTFILE  $DATA1{$_}{"phastConsElements46way"}."\t";         
print MYOUTPUTFILE  $DATA1{$_}{"LJB2_SIFT"}."\t";
print MYOUTPUTFILE  $DATA1{$_}{"LJB2_PolyPhen2_HDIV"}."\t";
print MYOUTPUTFILE  $DATA1{$_}{"LJB2_PolyPhen2_HVAR"}."\t";
print MYOUTPUTFILE  $DATA1{$_}{"LJB2_LRT"}."\t";
print MYOUTPUTFILE  $DATA1{$_}{"LJB2_MutationTaster"}."\t";
print MYOUTPUTFILE  $DATA1{$_}{"LJB_MutationAssessor"}."\t";
print MYOUTPUTFILE  $DATA1{$_}{"LJB2_FATHMM"}."\t";
print MYOUTPUTFILE  $DATA1{$_}{"LJB2_GERP++"}."\t";
print MYOUTPUTFILE  $DATA1{$_}{"LJB2_PhyloP"}."\t";
print MYOUTPUTFILE  $DATA1{$_}{"LJB2_SiPhy"}."\t";
print MYOUTPUTFILE  $DATA1{$_}{"avsift"}."\t";
print MYOUTPUTFILE  $DATA1{$_}{"tfbs"}."\t";
print MYOUTPUTFILE  $DATA1{$_}{"gerp++gt2"}."\t";
print MYOUTPUTFILE  $DATA1{$_}{"gwasCatalog"}."\t";
print MYOUTPUTFILE  $DATA1{$_}{"evofold"}."\t";
print MYOUTPUTFILE  $DATA1{$_}{"wgRna"}."\t";
print MYOUTPUTFILE  $DATA1{$_}{"miRNA_BindingSites"}."\t";
print MYOUTPUTFILE  $DATA1{$_}{"Pub_Articles"}."\t";
print MYOUTPUTFILE  $DATA1{$_}{"OMIM_Titles"}."\t";
print MYOUTPUTFILE  $DATA1{$_}{"OMIM_Descriptions"}."\t";
print MYOUTPUTFILE  $DATA1{$_}{"OMIM_Symptoms"}."\t";
print MYOUTPUTFILE  $DATA1{$_}{"Info"}."\t";
print MYOUTPUTFILE  $DATA1{$_}{"Format"}."\t";
print MYOUTPUTFILE  $DATA1{$_}{"Data"}."\n";
  	};		
close(MYOUTPUTFILE);

 
