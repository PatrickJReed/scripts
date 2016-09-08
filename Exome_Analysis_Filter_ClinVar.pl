#!/usr/bin/perl -w
use strict;
#use warnings;
my $Variants		= "./ExAC.r0.3.sites.ClinVar.txt";
my $MedGen_OMIM		= "./MedGen_HPO_OMIM_Mapping.txt";
my $rsID_OMIM		= "./OMIM_ClinVar_Variants3.txt";
my $MG_REL			= "./MG_REL.txt";
my $MG_NAMES		= "./MG_NAMES.txt";
my $MG_STY			= "./MG_STY.txt";

my $Output_File		= "./ClinVar_20141124.vcf.annovar.aa.hg19_multianno.Filtered.txt";
my (@tmp1, @tmp2, @IDs);
my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, $S1, $S2, $S3, $S4, $S5, $S6, $S7, %DATA1, @varID, $var, $PL, $GT, $DP, $GQ, $pPL, $pGT, $pDP, $pGQ, $aGT, $aPL, $aDP, $aGQ, $uGT, $uPL, $uDP, $uGQ, $dum, $depth, %Genotypes, $Ale1, $Ale2, $pAle1, $pAle2, $aAle1, $aAle2);


##Parse Variant VCF File 
open (LIST, $Variants) || die "File not found\n";     
     while (<LIST>) {
     chomp;
         @tmp1 = split(/\t/, $_);
         if ($tmp1[1] ne "Start"){
		 $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]\t$tmp1[4]"}{"Chr"} = $tmp1[0];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]\t$tmp1[4]"}{"Start"} = $tmp1[1];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]\t$tmp1[4]"}{"End"} = $tmp1[2];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]\t$tmp1[4]"}{"Ref"} = $tmp1[3];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]\t$tmp1[4]"}{"Alt"} = $tmp1[4];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]\t$tmp1[4]"}{"Func.refGene"} = $tmp1[5];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]\t$tmp1[4]"}{"Gene.refGene"} = $tmp1[6];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]\t$tmp1[4]"}{"ExonicFunc.refGene"} = $tmp1[7];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]\t$tmp1[4]"}{"AAChange.refGene"} = $tmp18];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]\t$tmp1[4]"}{"ESP6500si_ALL"} = $tmp1[9];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]\t$tmp1[4]"}{"1000g2014sep_all"} = $tmp1[10]; 
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]\t$tmp1[4]"}{"snp141"} = $tmp1[11];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]\t$tmp1[4]"}{"SIFT_score"} = $tmp1[12];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]\t$tmp1[4]"}{"GERP++_RS"} = $tmp1[13];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]\t$tmp1[4]"}{"gerp++gt2"} = $tmp1[14];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]\t$tmp1[4]"}{"dbSNP"} = $tmp1[15];         
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]\t$tmp1[4]"}{"Quality"} = $tmp1[16];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]\t$tmp1[4]"}{"VQSR"} = $tmp1[17];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]\t$tmp1[4]"}{"Data"} = $tmp1[18];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]\t$tmp1[4]"}{"CLINSIG"} = $tmp1[19];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]\t$tmp1[4]"}{"CLNDBN"} = $tmp1[20];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]\t$tmp1[4]"}{"CLNDSDB"} = $tmp1[21];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]\t$tmp1[4]"}{"CLNDSDBID"} = $tmp1[22];

      	};       
	
};
close(LIST);

##Parse dbSNP_OMIM relations
open (LIST, $rsID_OMIM) || die "File not found\n";     
     while (<LIST>) {
     chomp;
         @tmp1 = split(/\t/, $_);
         if ($tmp1[1] ne "Type"){
		 $DATA2{"$tmp1[6]"}{"AlleleID"} = $tmp1[0];
         $DATA2{"$tmp1[6]"}{"Type"} = $tmp1[1];
         $DATA2{"$tmp1[6]"}{"Name"} = $tmp1[2];
         $DATA2{"$tmp1[6]"}{"GeneID"} = $tmp1[3];
         $DATA2{"$tmp1[6]"}{"GeneSymbol"} = $tmp1[4];
         $DATA2{"$tmp1[6]"}{"ClinicalSignificance"} = $tmp1[5];
         $DATA2{"$tmp1[6]"}{"RS_dbSNP"} = $tmp1[6];
         $DATA2{"$tmp1[6]"}{"dbVar"} = $tmp1[7];
         $DATA2{"$tmp1[6]"}{"TestedInGTR"} = $tmp1[8];
         $DATA2{"$tmp1[6]"}{"Origin"} = $tmp1[9];
         $DATA2{"$tmp1[6]"}{"Assembly"} = $tmp1[10];
         $DATA2{"$tmp1[6]"}{"Chromosome"} = $tmp1[11];
         $DATA2{"$tmp1[6]"}{"Start"} = $tmp1[12];
         $DATA2{"$tmp1[6]"}{"Stop"} = $tmp1[13];
         $DATA2{"$tmp1[6]"}{"Cytogenetic"} = $tmp1[14];
         $DATA2{"$tmp1[6]"}{"ReviewStatus"} = $tmp1[15];
         $DATA2{"$tmp1[6]"}{"HGVS_c"} = $tmp1[16];
         $DATA2{"$tmp1[6]"}{"HGVS_p"} = $tmp1[17];
         $DATA2{"$tmp1[6]"}{"NumberSubmitters"} = $tmp1[18];
         $DATA2{"$tmp1[6]"}{"LastEvaluated"} = $tmp1[19];
         $DATA2{"$tmp1[6]"}{"OtherIDs"} = $tmp1[20];
         $DATA2{"$tmp1[6]"}{"VariantID1"} = $tmp1[21];
         $DATA2{"$tmp1[6]"}{"VariantID2"} = $tmp1[22];
         $DATA2{"$tmp1[6]"}{"VariantID3"} = $tmp1[23];
         $DATA2{"$tmp1[6]"}{"VariantID4"} = $tmp1[24];
         $DATA2{"$tmp1[6]"}{"VariantID5"} = $tmp1[25];
         $DATA2{"$tmp1[6]"}{"VariantID6"} = $tmp1[26];
         $DATA2{"$tmp1[6]"}{"VariantID7"} = $tmp1[27];
         $DATA2{"$tmp1[6]"}{"VariantID8"} = $tmp1[28];
         $DATA2{"$tmp1[6]"}{"MedGen1"} = $tmp1[29];
		 $DATA2{"$tmp1[6]"}{"MedGen2"} = $tmp1[30];
		 $DATA2{"$tmp1[6]"}{"MedGen3"} = $tmp1[31];
		 $DATA2{"$tmp1[6]"}{"MedGen4"} = $tmp1[32];
		 $DATA2{"$tmp1[6]"}{"MedGen5"} = $tmp1[33];
		 $DATA2{"$tmp1[6]"}{"MedGen6"} = $tmp1[34];
      	};       
	
};
close(LIST);

open (LIST, $MG_NAMES) || die "File not found\n";     
     while (<LIST>) {
     chomp;
         @tmp1 = split(/\|/, $_);
         if ($tmp1[1] ne "name"){
		 $NAMES{"$tmp1[0]"}{"CUI"} = $tmp1[0];
         $NAMES{"$tmp1[0]"}{"name"} = $tmp1[1];
         $NAMES{"$tmp1[0]"}{"source"} = $tmp1[2];
         $NAMES{"$tmp1[0]"}{"SUPPRESS"} = $tmp1[3];
      	};       
	
};
close(LIST);

my @VARs = keys(%DATA1);
foreach (@VARs) {
	if ($DATA2{$DDATA1{$_}{"snp141"}}) {
	my $RS = $DDATA1{$_}{"snp141"};
		 $DATA1{$_}{"AlleleID"} = $DATA2{$RS}{"AlleleID"};
         $DATA1{$_}{"Type"} = $DATA2{$RS}{"Type"};
         $DATA1{$_}{"Name"} = $DATA2{$RS}{"Name"};
         $DATA1{$_}{"GeneID"} = $DATA2{$RS}{"GeneID"};
         $DATA1{$_}{"GeneSymbol"} = $DATA2{$RS}{"GeneSymbol"};
         $DATA1{$_}{"ClinicalSignificance"} = $DATA2{$RS}{"ClinicalSignificance"};
         $DATA1{$_}{"RS_dbSNP"} = $DATA2{$RS}{"RS_dbSNP"};
         $DATA1{$_}{"dbVar"} = $DATA2{$RS}{"dbVar"};
         $DATA1{$_}{"TestedInGTR"} = $DATA2{$RS}{"TestedInGTR"};
         $DATA1{$_}{"Origin"} = $DATA2{$RS}{"Origin"};
         $DATA1{$_}{"Assembly"} = $DATA2{$RS}{"Assembly"};
         $DATA1{$_}{"Chromosome"} = $DATA2{$RS}{"Chromosome"};
         $DATA1{$_}{"Start"} = $DATA2{$RS}{"Start"};
         $DATA1{$_}{"Stop"} = $DATA2{$RS}{"Stop"};
         $DATA1{$_}{"Cytogenetic"} = $DATA2{$RS}{"Cytogenetic"};
         $DATA1{$_}{"ReviewStatus"} = $DATA2{$RS}{"ReviewStatus"};
         $DATA1{$_}{"HGVS_c"} = $DATA2{$RS}{"HGVS_c"};
         $DATA1{$_}{"HGVS_p"} = $DATA2{$RS}{"HGVS_p"};
         $DATA1{$_}{"NumberSubmitters"} = $DATA2{$RS}{"NumberSubmitters"};
         $DATA1{$_}{"LastEvaluated"} = $DATA2{$RS}{"LastEvaluated"};
         $DATA1{$_}{"OtherIDs"} = $DATA2{$RS}{"OtherIDs"};
         $DATA1{$_}{"VariantID1"} = $DATA2{$RS}{"VariantID1"};
         $DATA1{$_}{"VariantID2"} = $DATA2{$RS}{"VariantID2"};
         $DATA1{$_}{"VariantID3"} = $DATA2{$RS}{"VariantID3"};
         $DATA1{$_}{"VariantID4"} = $DATA2{$RS}{"VariantID4"};
         $DATA1{$_}{"VariantID5"} = $DATA2{$RS}{"VariantID5"};
         $DATA1{$_}{"VariantID6"} = $DATA2{$RS}{"VariantID6"};
         $DATA1{$_}{"VariantID7"} = $DATA2{$RS}{"VariantID7"};
         $DATA1{$_}{"VariantID8"} = $DATA2{$RS}{"VariantID8"};
         $DATA1{$_}{"MedGen1"} = $DATA2{$RS}{"MedGen1"};
		 $DATA1{$_}{"MedGen2"} = $DATA2{$RS}{"MedGen2"};
		 $DATA1{$_}{"MedGen3"} = $DATA2{$RS}{"MedGen3"};
		 $DATA1{$_}{"MedGen4"} = $DATA2{$RS}{"MedGen4"};
		 $DATA1{$_}{"MedGen5"} = $DATA2{$RS}{"MedGen5"};
		 $DATA1{$_}{"MedGen6"} = $DATA2{$RS}{"MedGen6"};
	};
};

##Parse MedGen_OMIM Relations
open (LIST, $MedGen_OMIM) || die "File not found\n";     
     while (<LIST>) {
     chomp;
         @tmp1 = split(/\|/, $_);
         if ($tmp1[1] ne "MIM_number"){
		 $DATA3{"$tmp1[0]"}{"$tmp1[0]\t$tmp[3]\t$tmp1[4]"}{"OMIM_CUI"} = $tmp1[0];
         $DATA3{"$tmp1[0]"}{"$tmp1[0]\t$tmp[3]\t$tmp1[4]"}{"MIM_number"} = $tmp1[1];
         $DATA3{"$tmp1[0]"}{"$tmp1[0]\t$tmp[3]\t$tmp1[4]"}{"OMIM_name"} = $tmp1[2];
         $DATA3{"$tmp1[0]"}{"$tmp1[0]\t$tmp[3]\t$tmp1[4]"}{"relationship"} = $tmp1[3];
         $DATA3{"$tmp1[0]"}{"$tmp1[0]\t$tmp[3]\t$tmp1[4]"}{"HPO_CUI"} = $tmp1[4];
         $DATA3{"$tmp1[0]"}{"$tmp1[0]\t$tmp[3]\t$tmp1[4]"}{"OMIM_Data"} = "$tmp1[1].:.$tmp1[2]";
      	};       
	
};
close(LIST);


my @VARs = keys(%DATA1);
foreach (@VARs) {
	if ($DATA3{$DATA1{$_}{"MedGen1"}}) {
	my $CUI_1 = $DATA1{$_}{"MedGen1"};
	my @HPO_CUIs = keys %{$DATA3{$CUI_1}};
	foreach my $CUI_2 (@HPO_CUIs) {
		$DATA1{$_}{"$CUI_2"} = $DATA3{$CUI_1}{$CUI_2}{"OMIM_Data"};
		};
	};
	elsif ($DATA3{$DATA1{$_}{"MedGen2"}}) {
	my $CUI_1 = $DATA1{$_}{"MedGen2"};
	my @HPO_CUIs = keys %{$DATA3{$CUI_1}};
	foreach my $CUI_2 (@HPO_CUIs) {
		$DATA1{$_}{"$CUI_2"} = $DATA3{$CUI_1}{$CUI_2}{"OMIM_Data"};
		};
	};
	elsif ($DATA3{$DATA1{$_}{"MedGen3"}})
	my $CUI_1 = $DATA1{$_}{"MedGen3"};
	my @HPO_CUIs = keys %{$DATA3{$CUI_1}};
	foreach my $CUI_2 (@HPO_CUIs) {
		$DATA1{$_}{"$CUI_2"} = $DATA3{$CUI_1}{$CUI_2}{"OMIM_Data"};
		};
	};
	elsif ($DATA3{$DATA1{$_}{"MedGen4"}})
	my $CUI_1 = $DATA1{$_}{"MedGen4"};
	my @HPO_CUIs = keys %{$DATA3{$CUI_1}};
	foreach my $CUI_2 (@HPO_CUIs) {
		$DATA1{$_}{"$CUI_2"} = $DATA3{$CUI_1}{$CUI_2}{"OMIM_Data"};
		};
	};
	elsif ($DATA3{$DATA1{$_}{"MedGen5"}})
	my $CUI_1 = $DATA1{$_}{"MedGen5"};
	my @HPO_CUIs = keys %{$DATA3{$CUI_1}};
	foreach my $CUI_2 (@HPO_CUIs) {
		$DATA1{$_}{"$CUI_2"} = $DATA3{$CUI_1}{$CUI_2}{"OMIM_Data"};
		};
	};
	elsif ($DATA3{$DATA1{$_}{"MedGen6"}})
	my $CUI_1 = $DATA1{$_}{"MedGen6"};
	my @HPO_CUIs = keys %{$DATA3{$CUI_1}};
	foreach my $CUI_2 (@HPO_CUIs) {
		$DATA1{$_}{"$CUI_2"} = $DATA3{$CUI_1}{$CUI_2}{"OMIM_Data"};
		};
	};


##Parse MedGen Ref Databases
open (LIST, $MG_REL) || die "File not found\n";     
     while (<LIST>) {
     chomp;
         @tmp1 = split(/\|/, $_);
         if ($tmp1[1] ne "AUI1"){
		 $DATA4{"$tmp1[0]"}{"$tmp1[0]\t$tmp1[6]\t$tmp1[4]"}{"CUI1"} = $tmp1[0];
         $DATA4{"$tmp1[0]"}{"$tmp1[0]\t$tmp1[6]\t$tmp1[4]"}{"AUI1"} = $tmp1[1];
         $DATA4{"$tmp1[0]"}{"$tmp1[0]\t$tmp1[6]\t$tmp1[4]"}{"STYPE1"} = $tmp1[2];
         $DATA1{"$tmp1[0]"}{"$tmp1[0]\t$tmp1[6]\t$tmp1[4]"}{"REL"} = $tmp1[3];
         $DATA1{"$tmp1[0]"}{"$tmp1[0]\t$tmp1[6]\t$tmp1[4]"}{"CUI2"} = $tmp1[4];
         $DATA1{"$tmp1[0]"}{"$tmp1[0]\t$tmp1[6]\t$tmp1[4]"}{"AUI2"} = $tmp1[5];
         $DATA1{"$tmp1[0]"}{"$tmp1[0]\t$tmp1[6]\t$tmp1[4]"}{"RELA"} = $tmp1[6];
         $DATA1{"$tmp1[0]"}{"$tmp1[0]\t$tmp1[6]\t$tmp1[4]"}{"RUI"} = $tmp1[7];
         $DATA1{"$tmp1[0]"}{"$tmp1[0]\t$tmp1[6]\t$tmp1[4]"}{"SAB"} = $tmp1[8];
         $DATA1{"$tmp1[0]"}{"$tmp1[0]\t$tmp1[6]\t$tmp1[4]"}{"SL"} = $tmp1[9];
         $DATA1{"$tmp1[0]"}{"$tmp1[0]\t$tmp1[6]\t$tmp1[4]"}{"SUPPRESS"} = $tmp1[10];		
      	};       
	
};
close(LIST);

my @VARs = keys(%DATA1);
foreach (@VARs) {
	if ($DATA4{$DATA1{$_}{"VariantID1"}}) {
	my $CUI_1 = $DATA1{$_}{"VariantID1"};
	my @RELAs = keys($DATA4{$CUI_1});
	foreach my $ID2 (@RELAs) {
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"CUI1"} = $DATA4{"$CUI_1"}{"$ID2"}{"CUI1"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"AUI1"} = $DATA4{"$CUI_1"}{"$ID2"}{"AUI1"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"STYPE1"} = $DATA4{"$CUI_1"}{"$ID2"}{"STYPE1"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"REL"} = $DATA4{"$CUI_1"}{"$ID2"}{"REL"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"CUI2"} = $DATA4{"$CUI_1"}{"$ID2"}{"CUI2"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"AUI2"} = $DATA4{"$CUI_1"}{"$ID2"}{"AUI2"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"RELA"} = $DATA4{"$CUI_1"}{"$ID2"}{"RELA"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"RUI"} = $DATA4{"$CUI_1"}{"$ID2"}{"RUI"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"SAB"} = $DATA4{"$CUI_1"}{"$ID2"}{"SAB"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"SL"} = $DATA4{"$CUI_1"}{"$ID2"}{"SL"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"SUPPRESS"} = $DATA4{"$CUI_1"}{"$ID2"}{"SUPPRESS"};
		};
	}:
	elsif($DATA4{$DATA1{$_}{"VariantID2"}}) {
	my $CUI_1 = $DATA1{$_}{"VariantID2"};
	my @RELAs = keys($DATA4{$CUI_1});
	foreach my $ID2 (@RELAs) {
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"CUI1"} = $DATA4{"$CUI_1"}{"$ID2"}{"CUI1"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"AUI1"} = $DATA4{"$CUI_1"}{"$ID2"}{"AUI1"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"STYPE1"} = $DATA4{"$CUI_1"}{"$ID2"}{"STYPE1"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"REL"} = $DATA4{"$CUI_1"}{"$ID2"}{"REL"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"CUI2"} = $DATA4{"$CUI_1"}{"$ID2"}{"CUI2"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"AUI2"} = $DATA4{"$CUI_1"}{"$ID2"}{"AUI2"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"RELA"} = $DATA4{"$CUI_1"}{"$ID2"}{"RELA"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"RUI"} = $DATA4{"$CUI_1"}{"$ID2"}{"RUI"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"SAB"} = $DATA4{"$CUI_1"}{"$ID2"}{"SAB"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"SL"} = $DATA4{"$CUI_1"}{"$ID2"}{"SL"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"SUPPRESS"} = $DATA4{"$CUI_1"}{"$ID2"}{"SUPPRESS"};
		};
	}:
	elsif($DATA4{$DATA1{$_}{"VariantID3"}}) {
	my $CUI_1 = $DATA1{$_}{"VariantID3"};
	my @RELAs = keys($DATA4{$CUI_1});
	foreach my $ID2 (@RELAs) {
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"CUI1"} = $DATA4{"$CUI_1"}{"$ID2"}{"CUI1"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"AUI1"} = $DATA4{"$CUI_1"}{"$ID2"}{"AUI1"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"STYPE1"} = $DATA4{"$CUI_1"}{"$ID2"}{"STYPE1"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"REL"} = $DATA4{"$CUI_1"}{"$ID2"}{"REL"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"CUI2"} = $DATA4{"$CUI_1"}{"$ID2"}{"CUI2"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"AUI2"} = $DATA4{"$CUI_1"}{"$ID2"}{"AUI2"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"RELA"} = $DATA4{"$CUI_1"}{"$ID2"}{"RELA"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"RUI"} = $DATA4{"$CUI_1"}{"$ID2"}{"RUI"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"SAB"} = $DATA4{"$CUI_1"}{"$ID2"}{"SAB"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"SL"} = $DATA4{"$CUI_1"}{"$ID2"}{"SL"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"SUPPRESS"} = $DATA4{"$CUI_1"}{"$ID2"}{"SUPPRESS"};
		};
	}:
	elsif($DATA4{$DATA1{$_}{"VariantID4"}}) {
	my $CUI_1 = $DATA1{$_}{"VariantID24"};
	my @RELAs = keys($DATA4{$CUI_1});
	foreach my $ID2 (@RELAs) {
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"CUI1"} = $DATA4{"$CUI_1"}{"$ID2"}{"CUI1"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"AUI1"} = $DATA4{"$CUI_1"}{"$ID2"}{"AUI1"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"STYPE1"} = $DATA4{"$CUI_1"}{"$ID2"}{"STYPE1"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"REL"} = $DATA4{"$CUI_1"}{"$ID2"}{"REL"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"CUI2"} = $DATA4{"$CUI_1"}{"$ID2"}{"CUI2"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"AUI2"} = $DATA4{"$CUI_1"}{"$ID2"}{"AUI2"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"RELA"} = $DATA4{"$CUI_1"}{"$ID2"}{"RELA"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"RUI"} = $DATA4{"$CUI_1"}{"$ID2"}{"RUI"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"SAB"} = $DATA4{"$CUI_1"}{"$ID2"}{"SAB"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"SL"} = $DATA4{"$CUI_1"}{"$ID2"}{"SL"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"SUPPRESS"} = $DATA4{"$CUI_1"}{"$ID2"}{"SUPPRESS"};
		};
	}:
	elsif($DATA4{$DATA1{$_}{"VariantID5"}}) {
	my $CUI_1 = $DATA1{$_}{"VariantID5"};
	my @RELAs = keys($DATA4{$CUI_1});
	foreach my $ID2 (@RELAs) {
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"CUI1"} = $DATA4{"$CUI_1"}{"$ID2"}{"CUI1"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"AUI1"} = $DATA4{"$CUI_1"}{"$ID2"}{"AUI1"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"STYPE1"} = $DATA4{"$CUI_1"}{"$ID2"}{"STYPE1"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"REL"} = $DATA4{"$CUI_1"}{"$ID2"}{"REL"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"CUI2"} = $DATA4{"$CUI_1"}{"$ID2"}{"CUI2"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"AUI2"} = $DATA4{"$CUI_1"}{"$ID2"}{"AUI2"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"RELA"} = $DATA4{"$CUI_1"}{"$ID2"}{"RELA"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"RUI"} = $DATA4{"$CUI_1"}{"$ID2"}{"RUI"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"SAB"} = $DATA4{"$CUI_1"}{"$ID2"}{"SAB"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"SL"} = $DATA4{"$CUI_1"}{"$ID2"}{"SL"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"SUPPRESS"} = $DATA4{"$CUI_1"}{"$ID2"}{"SUPPRESS"};
		};
	}:
	elsif($DATA4{$DATA1{$_}{"VariantID6"}}) {
	my $CUI_1 = $DATA1{$_}{"VariantID6"};
	my @RELAs = keys($DATA4{$CUI_1});
	foreach my $ID2 (@RELAs) {
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"CUI1"} = $DATA4{"$CUI_1"}{"$ID2"}{"CUI1"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"AUI1"} = $DATA4{"$CUI_1"}{"$ID2"}{"AUI1"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"STYPE1"} = $DATA4{"$CUI_1"}{"$ID2"}{"STYPE1"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"REL"} = $DATA4{"$CUI_1"}{"$ID2"}{"REL"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"CUI2"} = $DATA4{"$CUI_1"}{"$ID2"}{"CUI2"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"AUI2"} = $DATA4{"$CUI_1"}{"$ID2"}{"AUI2"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"RELA"} = $DATA4{"$CUI_1"}{"$ID2"}{"RELA"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"RUI"} = $DATA4{"$CUI_1"}{"$ID2"}{"RUI"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"SAB"} = $DATA4{"$CUI_1"}{"$ID2"}{"SAB"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"SL"} = $DATA4{"$CUI_1"}{"$ID2"}{"SL"};
		$DATA1{$_}{"$CUI_1"}{"$ID2"}{"SUPPRESS"} = $DATA4{"$CUI_1"}{"$ID2"}{"SUPPRESS"};
		};
	}:







open(MYOUTPUTFILE, ">$Output_File");
@IDs =keys(%DATA1);
print MYOUTPUTFILE 
foreach (@IDs) {
	if ($DATA1{$_}{"Chr"}) { 
print MYOUTPUTFILE $DATA1{$_}{"Chr"}."\t";
print MYOUTPUTFILE $DATA1{$_}{"Start"}."\t";
print MYOUTPUTFILE $DATA1{$_}{"End"}."\t"; 
print MYOUTPUTFILE $DATA1{$_}{"Ref"}."\t"; 
print MYOUTPUTFILE $DATA1{$_}{"Alt"}."\t"; 
print MYOUTPUTFILE $DATA1{$_}{"Func.refGene"}."\t"; 
print MYOUTPUTFILE $DATA1{$_}{"Gene.refGene"}."\t"; 
print MYOUTPUTFILE $DATA1{$_}{"GeneDetail.refGene"}."\t"; 
print MYOUTPUTFILE $DATA1{$_}{"ExonicFunc.refGene"}."\t"; 
print MYOUTPUTFILE $DATA1{$_}{"AAChange.refGene"}."\t"; 
print MYOUTPUTFILE $DATA1{$_}{"exac02"}."\t";
print MYOUTPUTFILE $DATA1{$_}{"CG69"}."\t";
print MYOUTPUTFILE $DATA1{$_}{"ESP6500si_ALL"}."\t";  
print MYOUTPUTFILE $DATA1{$_}{"1000g2014sep_all"}."\t"; 
print MYOUTPUTFILE $DATA1{$_}{"snp141"}."\t"; 
print MYOUTPUTFILE $DATA1{$_}{"snp138NonFlagged"}."\t"; 
print MYOUTPUTFILE $DATA1{$_}{"clinvar_20140211"}."\t"; 
print MYOUTPUTFILE $DATA1{$_}{"OMIM_Regions"}."\t";
print MYOUTPUTFILE $DATA1{$_}{"pubsBlat"}."\t";
print MYOUTPUTFILE $DATA1{$_}{"gerp++gt2"}."\t";
print MYOUTPUTFILE $DATA1{$_}{"gwasCatalog"}."\n"; 
  		}
  	};		
close(MYOUTPUTFILE);

##Parse Variant VCF File 

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
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"GeneDetail.refGene"} = $tmp1[7];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"ExonicFunc.refGene"} = $tmp1[8];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"AAChange.refGene"} = $tmp1[9];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"exac02"} = $tmp1[10];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"CG69"} = $tmp1[11];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"ESP6500si_ALL"} = $tmp1[12];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"1000g2014sep_all"} = $tmp1[13]; 
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"snp141"} = $tmp1[14];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"snp138NonFlagged"} = $tmp1[15];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"clinvar_20140211"} = $tmp1[16];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"OMIM_Regions"} = $tmp1[17];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"pubsBlat"} = $tmp1[18];         
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"gerp++gt2"} = $tmp1[19];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"gwasCatalog"} = $tmp1[20];

      	};       
	
};
close(LIST);

 
